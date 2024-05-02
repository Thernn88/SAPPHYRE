from collections import defaultdict
from multiprocessing import Pool
from os import listdir, makedirs, path
from pathlib import Path
from shutil import move
import copy
from msgspec import Struct, json
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    constrained_distance,
)
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta


def min_aa_check(sequence: list, minimum: int) -> bool:
    """
    Checks if the given sequence has at least the minimum amount of data characters.
    '-' and 'X' are not counted as data. All other characters are counted as data.
    Given sequence must be ascii characters.
    Returns the answer as a boolean.
    If enough data is found, returns True. Otherwise, returns False.
    """
    valid_count = len(sequence) - sequence.count("-") - sequence.count("X")
    return valid_count >= minimum


def first_last_X(string: str) -> tuple:
    """
    Finds the first and last occurrence of 'X' in the given string.
    If found, returns a tuple containing the slice indices [first, last).
    First in inclusive, last is exclusive.
    If no 'X' is found in the string, returns (None, None).
    """
    first = None
    last = None
    for i, bp in enumerate(string):
        if bp == "X":
            first = i
            break
    for i in range(len(string) - 1, -1, -1):
        if string[i] == "X":
            last = i
            break
    if first is None or last is None:
        return None, None
    return first, last + 1

def check_bad_regions(
    consensus: str, limit: float, initial_window=16, offset=0, _rig=False
) -> list:
    """
    Checks a consensus sequence for regions containing too many 'X' characters.

    Iterates over the consensus using a sliding window.
    If the required number of 'X' are found, expands the end of the window until
    the ratio of X falls below the given threshold, then reduces the window to the
    indices of the first and last 'X' to make the final result.

    If the required number of 'X' are not found, the current window is dumped and the
    exclusive ending index of the current window becomes the inclusive starting index of the
    next window.

    Consensus is the region of the consensus sequence to be searched.

    Limit is a float between 0.0 and 1.0 that represents the minimum ratio of
    X to all characters to flag a region.

    Initial_window is an int that sets the starting length of the window.

    Offset is the index of the first character to be searched. By passing this to
    the function and adding it to the results, we can search an arbitrary slice of the
    consensus sequence and return indices that are usable without modification.

    Rig is a debugging flag which can be ignored by a regular user.
    """
    if not _rig:
        first, last = 0, len(consensus)
    else:
        first, last = find_index_pair(consensus, "X")
    if first is None or last is None:
        return []
    left, right = first, first + initial_window
    num_x = consensus[left:right].count("X")
    begin = None
    output = {}

    while right <= last:
        ratio = num_x / (right - left)
        if ratio > limit:  # if ratio is above limit, this is a bad region
            if begin is None:  # if begin is None, this is the start of a region
                begin = left
            if right == len(consensus):
                # This break stops an index error. The region will still be
                # caught by the last check in this function.
                break
            # increment the X count, then expand the window
            # without this we have to recount all the characters in the window
            num_x += consensus[right] == "X"
            right += 1
            continue
        # If the following block executes, it means the ratio of X is below the threshold.
        if begin is not None:  # end of region, so make output
            a, b = first_last_X(consensus[begin : right + 1])
            a, b = a + begin + offset, b + begin + offset
            if b not in output:
                output[b] = (a, b)
        # since the current window is below the threshold, discard it and make a new one
        left = right
        right = left + initial_window
        begin = None
        num_x = consensus[left:right].count("X")
    # this check catches bad regions at the end of the string
    if begin is not None:
        a, b = first_last_X(consensus[begin:last])
        a, b = a + begin + offset, b + begin + offset
        if a is not None and b is not None:
            if b not in output:
                output[b] = (a, b)
    return [*output.values()]

def check_covered_bad_regions(consensus, min_ambiguous, ambig_char='X', max_distance=18):
    x_regions = []
    x_indices = []
    current_group = []

    start, stop = find_index_pair(consensus, "X")

    for i, base in enumerate(consensus[start:stop], start):
        if base == ambig_char:
            x_indices.append(i)

    for num in x_indices:
        if not current_group:
            current_group.append(num)
        elif num - current_group[-1] <= max_distance:
            current_group.append(num)
        else:
            if len(current_group) >= min_ambiguous:
                return (current_group[0], current_group[-1] + 1)
            current_group = [num]

    if current_group:
        if len(current_group) >= min_ambiguous:
            return (current_group[0], current_group[-1] + 1)


    return

def bundle_seqs_and_dupes(sequences: list, prepare_dupe_counts, reporter_dupe_counts):
    output = []
    for header, seq in sequences:
        node = header.split("|")[3]
        dupes = prepare_dupe_counts.get(node, 1) + sum(
            prepare_dupe_counts.get(node, 1)
            for node in reporter_dupe_counts.get(node, [])
        )
        output.append((seq, dupes))
    return output


def make_duped_consensus(
    raw_sequences: list, prepare_dupes: dict, reporter_dupes: dict, threshold: float
) -> str:
    seqs = [(header, seq) for header, seq in raw_sequences if header[-1] != "."]
    bundled_seqs = bundle_seqs_and_dupes(seqs, prepare_dupes, reporter_dupes)
    return dumb_consensus_dupe(bundled_seqs, threshold, 0)


class NODE(Struct):
    header: str
    sequence: str
    nt_sequence: str
    start: int
    end: int
    children: list

    def extend(self, node_2, overlap_coord):
        """
        Merges two nodes together, extending the current node to include the other node.
        Merges only occur if the overlapping kmer is a perfect match.

        Args:
        ----
            node_2 (NODE): The node to be merged into the current node
            overlap_coord (int): The index of the first matching character in the overlapping kmer
        Returns:
        -------
            None
        """
        # If node_2 is contained inside self
        if node_2.start >= self.start and node_2.end <= self.end:
            self.sequence = (
                self.sequence[:overlap_coord // 3]
                + node_2.sequence[overlap_coord // 3 : node_2.end]
                + self.sequence[node_2.end :]
            )

            self.nt_sequence = (
                self.nt_sequence[:overlap_coord // 3]
                + node_2.nt_sequence[overlap_coord // 3 : node_2.end * 3]
                + self.nt_sequence[node_2.end * 3 :]
            )

        # If node_2 contains self
        elif self.start >= node_2.start and self.end <= node_2.end:
            self.sequence = (
                node_2.sequence[:overlap_coord // 3]
                + self.sequence[overlap_coord // 3 : self.end]
                + node_2.sequence[self.end :]
            )

            self.nt_sequence = (
                node_2.nt_sequence[:overlap_coord // 3]
                + self.nt_sequence[overlap_coord // 3 : self.end * 3]
                + node_2.nt_sequence[self.end * 3 :]
            )

            self.start = node_2.start
            self.end = node_2.end
        # If node_2 is to the right of self
        elif node_2.start >= self.start:
            self.sequence = (
                self.sequence[:overlap_coord // 3] + node_2.sequence[overlap_coord // 3:]
            )

            self.nt_sequence = (
                self.nt_sequence[:overlap_coord // 3] + node_2.nt_sequence[overlap_coord // 3:]
            )

            self.end = node_2.end
        # If node_2 is to the left of self
        else:
            self.sequence = (
                node_2.sequence[:overlap_coord // 3] + self.sequence[overlap_coord // 3:]
            )

            self.nt_sequence = (
                node_2.nt_sequence[:overlap_coord // 3] + self.nt_sequence[overlap_coord // 3:]
            )

            self.start = node_2.start

        # Save node_2 and the children of node_2 to self
        self.children.append(node_2.header)
        self.children.extend(node_2.children)

    def clone(self):
        return NODE(self.header, self.sequence, self.nt_sequence, self.start, self.end, self.children)
    
    def contig_header(self):
        """
        Generates a header containg the contig node and all children nodes for debug

        Returns:
        -------
            str: The contig header
        """
        contig_node = self.header.split("|")[3]
        children_nodes = "|".join([i.split("|")[3] for i in self.children])
        return f"CONTIG_{contig_node}|{children_nodes}"


def simple_assembly(nodes, min_merge_overlap_percent):
    merged = set()
    for i, node in enumerate(nodes):
        if i in merged:
            continue
        merge_occured = True
        while merge_occured:
            merge_occured = False
            for j, node_b in enumerate(nodes):
                if j in merged:
                    continue
                if i == j:
                    continue

                overlap_coords=  get_overlap(node.start*3, node.end*3, node_b.start*3, node_b.end*3, 1)
                if overlap_coords:
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    overlap_percent = overlap_amount / (node.end - node.start)
                    if overlap_percent < min_merge_overlap_percent:
                        continue


                    kmer_a = node.nt_sequence[overlap_coords[0]:overlap_coords[1]]
                    kmer_b = node_b.nt_sequence[overlap_coords[0]:overlap_coords[1]]

                    if not is_same_kmer(kmer_a, kmer_b):
                        continue

                    merged.add(j)

                    overlap_coord = overlap_coords[0]
                    merge_occured = True
                    node.extend(node_b, overlap_coord)

                    nodes[j] = None
    
    nodes = [node for node in nodes if node is not None]

    return nodes
        

def identity(nodes_in_region, best_index, allowed_mismatches):
    indices = [i for i, node in enumerate(nodes_in_region) if i != best_index]
    best_node = nodes_in_region[best_index]
    #best_length = len(best_node.sequence) - best_node.sequence.count("-")
    extend_region = set()
    for i in indices:
        node = nodes_in_region[i]
        overlap_coords = get_overlap(best_node.start, best_node.end, node.start, node.end, 1)
        if overlap_coords is None:
            continue
        kmer_node = node.sequence[overlap_coords[0]:overlap_coords[1]]
        kmer_best = best_node.sequence[overlap_coords[0]:overlap_coords[1]]
        distance = constrained_distance(kmer_node, kmer_best)
        # node_length = len(node.sequence) - node.sequence.count("-")

        # if node_length < best_length * (1 - length_percentage):
        #     continue

        if distance < allowed_mismatches:
            extend_region.add(i)

    return extend_region


def indices_that_resolve(nodes, sequences_out_of_region, merge_percent):
    """
    This function is used to find the longest sequence in the list of sequences_out_of_region
    that can be merged into the list of nodes. If a sequence is found that can be merged,
    it is removed from the list of sequences_out_of_region and added to the list of nodes.
    """
    indices = []
    for i, node in enumerate(nodes):
        this_seq_out = copy.deepcopy(sequences_out_of_region)
        has_unambig_overlap = False
        for node_b in this_seq_out:
            overlap_coords = get_overlap(node.start, node.end, node_b.start, node_b.end, 1)
            
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]
                overlap_percent = overlap_amount / (node.end - node.start)
                
                if overlap_percent > merge_percent:
                    has_unambig_overlap = True
                    break

        if has_unambig_overlap:
            indices.append(i)

    return indices
        

def del_cols(sequence, columns, nt=False):
    if nt:
        seq = [sequence[i: i+3] for i in range(0, len(sequence), 3)]
        for i in columns:
            seq[i] = "---"
        return "".join(seq)
    seq = list(sequence)
    for i in columns:
        seq[i] = "-"
    return "".join(seq)


def get_coverage(aa_seqs, ref_avg_len):
    data_cols = 0
    for i in range(len(aa_seqs[0])):
        for seq in aa_seqs:
            if seq[i] != "-":
                data_cols += 1
                break

    return data_cols / ref_avg_len


def cluster(cluster_raw, cluster_threshold):
    clusters = []

    current_cluster = []
    for child_index, child_full in cluster_raw:
        if not current_cluster:
            current_cluster.append(child_full)
            current_index = child_index
        else:
            if child_index - current_index <= cluster_threshold:
                current_cluster.append(child_full)
                current_index = child_index
            else:
                clusters.append(current_cluster)
                current_cluster = [child_full]
                current_index = child_index
    
    if current_cluster:
        clusters.append(current_cluster)

    return clusters


def calculate_split(node_a: str, node_b: str, overlapping_coords: tuple, ref_consensus: dict) -> int:
    """Iterates over each position in the overlap range of sequence A and sequence B and
    creates a frankenstein sequence of sequence A + Sequence B joined at each
    position in the overlap.

    Final split position = pos in overlap with the highest score.

    Score is determined by the amount of characters that are the same between each
    position in the frankenstein sequence and the comparison sequence.
    """
    overlap_start, overlap_end = overlapping_coords

    highest_score = -1
    highest_scoring_pos = overlap_start

    for i in range(overlap_start, overlap_end + 1):
        prev_kmer = node_a.sequence[overlap_start:i]
        next_kmer = node_b.sequence[i:overlap_end]

        window = prev_kmer + next_kmer

        this_window_score = 0
        for x, let in enumerate(window, overlap_start):
            if let in ref_consensus[x]:
                this_window_score += ref_consensus[x].count(let)
            else:
                this_window_score -= 1
            if let == "-":
                this_window_score -= 1
                
        if this_window_score >= highest_score:
            highest_score = this_window_score
            highest_scoring_pos = i

    return highest_scoring_pos


def do_trim(aa_nodes, get_id, cluster_sets, x_positions, ref_consensus, kicked_headers, prepare_dupes, reporter_dupes, excise_trim_consensus):
    for cluster_set in cluster_sets:
        
        sub_aa_nodes = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or get_id(node.header) in cluster_set)]

        aa_sequences = [node.sequence for node in sub_aa_nodes]
        if aa_sequences:
            if prepare_dupes and reporter_dupes:
                current_raw_aa = [(node.header, node.sequence) for node in sub_aa_nodes]
                consensus_seq = make_duped_consensus(
                    current_raw_aa, prepare_dupes, reporter_dupes, excise_trim_consensus
                )
            else:
                consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)

            cstart, cend = find_index_pair(consensus_seq, "X")
            cstart, cend = find_index_pair(consensus_seq, "X")

            for i, maj_bp in enumerate(consensus_seq[cstart:cend], cstart):
                if maj_bp != "X":
                    continue
            for i, maj_bp in enumerate(consensus_seq[cstart:cend], cstart):
                if maj_bp != "X":
                    continue

                in_region = []
                out_of_region = []
                for x, node in enumerate(sub_aa_nodes):
                    # within 3 bp
                    within_left = i >= node.start and i <= node.start + 3
                    within_right = i <= node.end and i >= node.end - 3

                    if within_left or within_right:
                        in_region.append((x, node.sequence[i], within_right))
                    elif i >= node.start and i <= node.end:
                        out_of_region.append((x, node.sequence[i], within_right))
                    if within_left or within_right:
                        in_region.append((x, node.sequence[i], within_right))
                    elif i >= node.start and i <= node.end:
                        out_of_region.append((x, node.sequence[i], within_right))

                if not out_of_region and not in_region:
                    continue
                if not out_of_region and not in_region:
                    continue

                if not out_of_region and in_region:
                    for node_index, bp, on_end in in_region:
                        if bp in ref_consensus[i]:
                            continue
                        
                        if on_end:
                            for x in range(i, sub_aa_nodes[node_index].end):
                                x_positions[sub_aa_nodes[node_index].header].add(x)
                        else:
                            for x in range(sub_aa_nodes[node_index].start, i + 1):
                                x_positions[sub_aa_nodes[node_index].header].add(x)

                if out_of_region and in_region:
                    for node_index, bp, on_end in in_region:
                        if on_end:
                            for x in range(i, sub_aa_nodes[node_index].end):
                                x_positions[sub_aa_nodes[node_index].header].add(x)
                        else:
                            for x in range(sub_aa_nodes[node_index].start, i + 1):
                                x_positions[sub_aa_nodes[node_index].header].add(x)


            #refresh aa
            if x_positions:
                for node in sub_aa_nodes:
                    node.sequence = del_cols(node.sequence, x_positions[node.header])
                    node.start, node.end = find_index_pair(node.sequence, "-")

            if prepare_dupes and reporter_dupes:
                current_raw_aa = [(node.header, node.sequence) for node in sub_aa_nodes if node.header not in kicked_headers]
                consensus_seq = make_duped_consensus(
                    current_raw_aa, prepare_dupes, reporter_dupes, excise_trim_consensus
                )
            else:
                aa_sequences = [x.sequence for x in sub_aa_nodes if x.header not in kicked_headers]
                consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)

            for node in sub_aa_nodes:
                i = None
                for poss_i in range(node.start, node.start + 3):
                    if node.sequence[poss_i] != consensus_seq[poss_i]:
                        i = poss_i

                if not i is None:
                    for x in range(node.start , i + 1):
                        x_positions[node.header].add(x)

                i = None
                for poss_i in range(node.end -1, node.end - 4, -1):
                    if node.sequence[poss_i] != consensus_seq[poss_i]:
                        i = poss_i

                if not i is None:
                    for x in range(i, node.end):
                        x_positions[node.header].add(x)


def do_cluster(ids, ref_coords, max_distance=100):
    clusters = []
    ids.sort(key = lambda x: x[0])

    req_seq_coverage = 0.5

    current_cluster = []
    for i, (child_index, seq_coords) in enumerate(ids):

        coverage = len(seq_coords.intersection(ref_coords)) / len(ref_coords)


        if not current_cluster:
            current_cluster.append((child_index, coverage, i))
            current_index = child_index
        else:
            if child_index - current_index <= max_distance:
                current_cluster.append((child_index, coverage, i))
                current_index = child_index
            else:
                if len(current_cluster) >= 2:
                    cluster_data_cols = set()
                    for _, _, index in current_cluster:
                        cluster_data_cols.update(ids[index][1])
                        
                    cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

                    clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage))
                elif len(current_cluster) == 1:
                    if current_cluster[0][1] > req_seq_coverage:
                        cluster_coverage = current_cluster[0][1]
                        clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage))
                        
                current_cluster = [(child_index, coverage, i)]
                current_index = child_index

    if current_cluster:
        if len(current_cluster) >= 2:
            cluster_data_cols = set()
            for _, _, index in current_cluster:
                cluster_data_cols.update(ids[index][1])
                
            cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

            clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage))
        elif len(current_cluster) == 1:
            if current_cluster[0][1] > req_seq_coverage:
                cluster_coverage = current_cluster[0][1]
                clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage))
                
    return clusters



def log_excised_consensus(
    gene: str,
    is_assembly_or_genome: bool,
    is_genome: bool,
    input_path: Path,
    output_path: Path,
    compress_intermediates: bool,
    excise_overlap_merge,
    excise_overlap_ambig,
    excise_region_overlap,
    excise_consensus,
    excise_maximum_depth,
    excise_minimum_ambig,
    allowed_distance,
    excise_rescue_match,
    prepare_dupes: dict,
    reporter_dupes: dict,
    excise_trim_consensus,
    true_cluster_threshold = 24,
):
    """
    By default, this does non-dupe consensus. If you pass "dupes=True", it will call
    the dupe-weighted version of consensus. If dupes is enable, the program needs dupe counts
    in the standard rocksDB location.

    Consensus_threshold is a value between 0.0 and 1.0 which represents a ratio.
    At each location, a bp is chosen as the consensus bp if the ratio of occurrences/total_characters_at_location
    exceeds the consensus_threshold. Raising this value makes the selection more strict. If you want this to match
    the default in outlier.py, set this to 0.65

    After the consensus sequence is made, the tail is checked for excessive X placeholder characters.
    Excise_threshold is a value between 0.0 and 1.0, which represents the maximum allowable ratio of X/total_characters.
    The tail is checked in blocks of 16 characters. The block checks are cumulative, so the percentage at each block is
    affected directly by all the preceding blocks. For example, if there are 8 X in the first block, and 10 X in the
    second block, the percentage at the second block is calculated as (10+8)/(16+16) = 18/32 = 0.5625 .
    Lowering this value makes the check more strict. Small changes in the excise_threshold can result in
    disproportionately large increases in truncation. My best results happen around 0.40 so far.

    The first field in the log file is the gene name.
    The second field in the log file is the cut-index in regular slice notation. It starts at the index of the
    first removed character, inclusive. The end is the length of the string, which is exclusive.
    """
    log_output = []
    this_rescues = []

    aa_in = input_path.joinpath("aa", gene)
    aa_out = output_path.joinpath("aa", gene)

    nt_in = input_path.joinpath("nt", gene.replace(".aa.", ".nt."))
    nt_out = output_path.joinpath("nt", gene.replace(".aa.", ".nt."))

    x_positions = defaultdict(set)

    bp_count = lambda x: len(x) - x.count("-")

    raw_aa = list(parseFasta(str(aa_in)))

    ref_lens = []
    aa_nodes = []
    reference_cluster_data = set()
    ref_consensus = defaultdict(list)
    for header, seq in raw_aa:
        if header.endswith('.'):
            start, end = find_index_pair(seq, "-")
            for i, bp in enumerate(seq[start:end], start):
                if bp != "-":
                    reference_cluster_data.add(i)
                ref_consensus[i].append(bp)
                
            ref_lens.append(bp_count(seq))
            continue

        aa_nodes.append(NODE(header, seq, None, *find_index_pair(seq, "-"), []))

    ref_avg_len = sum(ref_lens) / len(ref_lens)
    kicked_headers = set()

    cluster_sets = [None]
    get_id = lambda x: int(x.split("|")[3].split("_")[1])
    if is_genome:
        ids = []
        for node in aa_nodes:
            if node.header not in kicked_headers:
                this_id = get_id(node.header)
                start, end = find_index_pair(node.sequence, "-")
                data_cols = {i for i, let in enumerate(node.sequence[start:end], start) if let != "-"}
                ids.append((this_id, data_cols))
    
        clusters = do_cluster(ids, reference_cluster_data)
        if clusters:
            cluster_sets = [set(range(a, b+1)) for a, b, _ in clusters]

    if not is_genome:
        do_trim(aa_nodes, get_id, cluster_sets, x_positions, ref_consensus, kicked_headers, prepare_dupes, reporter_dupes, excise_trim_consensus)

    aa_sequence = {}
    for node in aa_nodes:

        node.sequence = del_cols(node.sequence, x_positions[node.header])
        node.start, node.end = find_index_pair(node.sequence, "-")

        aa_sequence[node.header] = node.sequence
        node_kmer = node.sequence[node.start:node.end]
        data_len = bp_count(node_kmer)
        if data_len < 15:
            kicked_headers.add(node.header)

    raw_sequences = {header: del_cols(seq, x_positions[header], True) for header, seq in parseFasta(str(nt_in))}
    for node in aa_nodes:
        if node.header in kicked_headers:
            continue
        nt_seq = raw_sequences[node.header]
        node.nt_sequence = nt_seq
    
    true_cluster_raw = []

    for header in raw_sequences:
        true_cluster_raw.append((int(header.split("|")[3].split("_")[1]), header))

    true_cluster_raw.sort(key = lambda x: x[0])
    before_true_clusters = cluster(true_cluster_raw, true_cluster_threshold)


    # Search for slices of the consensus seq with a high ratio of 'X' to total characters
    region = True
    had_region = False
    recursion_max = 5 * len(cluster_sets)
    last_region = {}
    while region:
        region = None

        for cluster_i, cluster_set in enumerate(cluster_sets):

            aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or get_id(node.header) in cluster_set)]

            sequences = [node.nt_sequence for node in aa_subset]
            if not sequences:
                break

            if prepare_dupes and reporter_dupes:
                nt_sequences = [(node.header, node.nt_sequence) for node in aa_subset]
                consensus_seq = make_duped_consensus(
                    nt_sequences, prepare_dupes, reporter_dupes, excise_consensus
                )
            else:
                consensus_seq = dumb_consensus(sequences, excise_consensus, 0)

            consensus_seq = convert_consensus(sequences, consensus_seq)
            region = check_covered_bad_regions(consensus_seq, excise_minimum_ambig)
            if region == last_region.get(cluster_i, -1):
                recursion_max -= 1
                if recursion_max <= 0:
                    region = None
                    continue
            last_region[cluster_i] = region

            if region:
                had_region = True
                sequences_in_region = []
                sequences_out_of_region = []
                a, b = region         

                for i, node in enumerate(aa_subset):

                    if node.header in kicked_headers:
                        continue

                    overlap_coords = get_overlap(a, b, node.start * 3, node.end * 3, 1)
                    if overlap_coords:
                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        overlap_percent = overlap_amount / (node.end - node.start)
                        if overlap_percent >= excise_region_overlap: # Adjustable percent
                            sequences_in_region.append(aa_subset[i])
                        else:
                            sequences_out_of_region.append(aa_subset[i])

                if len(sequences_in_region) > excise_maximum_depth:
                    continue

                cluster_resolve_failed = True
                nodes_in_region = None
                if is_genome:
                    tagged_in_region = [(int(node.header.split("|")[3].split("_")[1]), node) for node in sequences_in_region]
                    tagged_in_region.sort(key=lambda x: x[0])
                    clusters = cluster(tagged_in_region, true_cluster_threshold)

                    log_output.append(f">{gene}_ambig_{a}:{b}\n{consensus_seq}")

                    for clust in clusters:
                        if len(clust) <= 1:
                            continue

                        clust.sort(key = lambda node: node.start)
                        for i, node in enumerate(clust):
                            if i == 0:
                                continue
                            
                            prev_node = clust[i-1]
                            overlapping_coords = get_overlap(node.start, node.end, prev_node.start, prev_node.end, 1)
                            if overlapping_coords:
                                kmer = node.sequence[overlapping_coords[0]:overlapping_coords[1]]
                                prev_kmer = prev_node.sequence[overlapping_coords[0]:overlapping_coords[1]]

                                if is_same_kmer(kmer, prev_kmer):
                                    continue

                                splice_index = calculate_split(prev_node, node, overlapping_coords, ref_consensus)
                                prev_positions = set()
                                node_positions = set()

                                for x in range(node.start, splice_index):
                                    node_positions.add(x)

                                for x in range(splice_index, prev_node.end):
                                    prev_positions.add(x)

                                log_output.append(f">{prev_node.header} vs {node.header}")
                                log_output.append(f"Split at {splice_index}")

                                prev_node.sequence = del_cols(prev_node.sequence, prev_positions)
                                node.sequence = del_cols(node.sequence, node_positions)

                                cluster_resolve_failed = False

                                either_kicked = False
                                if len(prev_node.sequence) - prev_node.sequence.count("-") < 15:
                                    kicked_headers.add(prev_node.header)
                                    either_kicked = True
                                    

                                if len(node.sequence) - node.sequence.count("-") < 15:
                                    kicked_headers.add(node.header)
                                    either_kicked = True

                                if either_kicked:
                                    break

                                prev_node.nt_sequence = del_cols(prev_node.nt_sequence, prev_positions, True)
                                node.nt_sequence = del_cols(node.nt_sequence, node_positions, True)

                                node.start, node.end = find_index_pair(node.sequence, "-")
                                prev_node.start, prev_node.end = find_index_pair(prev_node.sequence, "-")

                                x_positions[node.header].update(node_positions)
                                x_positions[prev_node.header].update(prev_positions)
                
                if cluster_resolve_failed:
                    
                    sequences_in_region = copy.deepcopy(sequences_in_region)
                    nodes_in_region = simple_assembly(sequences_in_region, excise_overlap_ambig)

                    node_indices = indices_that_resolve(nodes_in_region, sequences_out_of_region, excise_overlap_merge)


                    keep_indices = set()
                    if node_indices:
                        best_index = max(node_indices, key=lambda x: len(nodes_in_region[x].sequence) - nodes_in_region[x].sequence.count("-"))
                        keep_indices.update([best_index])

                        similar_indices = identity(nodes_in_region, best_index, allowed_distance)
                        keep_indices.update(similar_indices)

                    for i, node in enumerate(nodes_in_region):
                        if i in keep_indices:
                            continue

                        kicked_headers.add(node.header)
                        kicked_headers.update(node.children)

                if sequences_in_region:
                    log_output.append(f">{gene}_ambig_{a}:{b}\n{consensus_seq}")
                    if nodes_in_region:
                        log_output.extend([f">{node.contig_header()}_{'kept' if i in keep_indices else 'kicked'}\n{node.nt_sequence}" for i, node in enumerate(nodes_in_region)])
                    else:
                        log_output.extend([f">{node.header}_{'kept' if node.header not in kicked_headers else 'kicked'}\n{node.nt_sequence}" for node in sequences_in_region])
                    log_output.append("\n")
        if recursion_max <= 0:
            break

    if is_genome:
        do_trim(aa_nodes, get_id, cluster_sets, x_positions, ref_consensus, kicked_headers, prepare_dupes, reporter_dupes, excise_trim_consensus)

    if had_region:
        after_data = []
        for node in aa_nodes:
            if node.header in kicked_headers:
                continue
            after_data.append((node.header.split("|")[3].split("_")[1], node.header))

        after_data.sort(key = lambda x: x[0])
        after_true_clusters = []

        current_cluster = []
        for child_index, child_full in true_cluster_raw:
            if not current_cluster:
                current_cluster.append(child_full)
                current_index = child_index
            else:
                if child_index - current_index <= true_cluster_threshold:
                    current_cluster.append(child_full)
                    current_index = child_index
                else:
                    after_true_clusters.append(current_cluster)
                    current_cluster = [child_full]
                    current_index = child_index
        
        if is_assembly_or_genome:
            if current_cluster:
                after_true_clusters.append(current_cluster)

            after_true_cluster = set(max(after_true_clusters, key=lambda x: len(x)))
            matches = []
            for i, before_cluster in enumerate(before_true_clusters):
                matches.append((len(after_true_cluster.intersection(set(before_cluster))) / len(before_cluster), i))
            matches.sort(key=lambda x: x[0], reverse= True)
            
            best_match = max(matches, key=lambda x: x[0])[1]
            this_rescues = [f"Rescued in gene: {gene}"]
            for header in before_true_clusters[best_match]:
                if header in kicked_headers:

                    this_sequence = aa_sequence[header]
                    start, end = find_index_pair(this_sequence, "-")
                    matching_char = 0
                    for i, let in enumerate(this_sequence[start: end], start):
                        if let in ref_consensus[i]:
                            matching_char += 1

                    if matching_char / (end - start) < excise_rescue_match:
                        continue

                    kicked_headers.remove(header)
                    this_rescues.append(header)

    aa_output = [(header, del_cols(seq, x_positions[header])) for header, seq in raw_aa if header not in kicked_headers]

    aa_has_candidate = False
    for header, _ in aa_output:
        if not header.endswith("."):
            aa_has_candidate = True
            break

    if aa_has_candidate:
        gene_coverage = get_coverage([seq for header, seq in aa_output if header[-1] != "."], ref_avg_len)
        req_coverage = 0.4 if is_assembly_or_genome else 0.01
        if gene_coverage < req_coverage:
            log_output.append(f">{gene}_kicked_coverage_{gene_coverage}_of_{req_coverage}\n{consensus_seq}")
            return log_output, False, False, gene, False, len(aa_nodes), this_rescues
        
        writeFasta(aa_out, aa_output, compress_intermediates)
        nt_output = [(header, del_cols(seq, x_positions[header], True)) for header, seq in raw_sequences.items() if header not in kicked_headers]
        writeFasta(nt_out, nt_output, compress_intermediates)

        return log_output, had_region, False, False, gene, len(kicked_headers), this_rescues
    return log_output, had_region, gene, False, None, len(kicked_headers), this_rescues


def load_dupes(rocks_db_path: Path):
    rocksdb_db = RocksDB(str(rocks_db_path))
    prepare_dupe_counts = json.decode(
        rocksdb_db.get("getall:gene_dupes"), type=dict[str, dict[str, int]]
    )
    reporter_dupe_counts = json.decode(
        rocksdb_db.get("getall:reporter_dupes"), type=dict[str, dict[str, list]]
    )
    return prepare_dupe_counts, reporter_dupe_counts


### USED BY __main__
def do_move(from_, to_):
    move(from_, to_)


def move_flagged(to_move, processes):
    if processes <= 1:
        for from_, to_ in to_move:
            move(from_, to_)
    else:
        with Pool(processes) as pool:
            pool.starmap(do_move, to_move)


###


def main(args, sub_dir, is_genome, is_assembly_or_genome):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if not (0 < args.excise_overlap_merge < 1.0):
        if 0 < args.excise_overlap_merge <= 100:
            args.excise_overlap_merge = args.excise_overlap_merge / 100
        else:
            raise ValueError(
                "Cannot convert excise_overlap_merge to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.excise_overlap_ambig < 1.0):
        if 0 < args.excise_overlap_ambig <= 100:
            args.excise_overlap_ambig = args.excise_overlap_ambig / 100
        else:
            raise ValueError(
                "Cannot convert excise_overlap_ambig to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.excise_consensus < 1.0):
        if 0 < args.excise_consensus <= 100:
            args.excise_consensus = args.excise_consensus / 100
        else:
            raise ValueError(
                "Cannot convert excise_consensus to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.excise_region_overlap < 1.0):
        if 0 < args.excise_region_overlap <= 100:
            args.excise_region_overlap = args.excise_region_overlap / 100
        else:
            raise ValueError(
                "Cannot convert excise_region_overlap to a percent. Use a decimal or a whole number between 0 and 100"
            )

    folder = args.INPUT
    input_folder = Path(folder, "outlier", sub_dir)
    if not input_folder.exists():
        input_folder = Path(folder, sub_dir)

    output_folder = Path(folder, "outlier", "excise")

    output_aa_folder = output_folder.joinpath("aa")
    output_nt_folder = output_folder.joinpath("nt")

    if not path.exists(output_aa_folder):
        makedirs(str(output_aa_folder), exist_ok=True)
        makedirs(str(output_nt_folder), exist_ok=True)

    aa_input = input_folder.joinpath("aa")

    rocksdb_path = Path(folder, "rocksdb", "sequences", "nt")
    if args.no_dupes:
        prepare_dupes, reporter_dupes = {}, {}
    else:
        if not rocksdb_path.exists():
            err = f"cannot find rocksdb for {folder}"
            raise FileNotFoundError(err)
        prepare_dupes, reporter_dupes = load_dupes(rocksdb_path)

    compress = not args.uncompress_intermediates or args.compress

    genes = [fasta for fasta in listdir(aa_input) if ".fa" in fasta]
    log_path = Path(output_folder, "excise_regions.txt")
    gene_log_path = Path(output_folder, "excise_genes.txt")
    if args.processes > 1:
        arguments = [
            (
                gene,
                is_assembly_or_genome,
                is_genome,
                input_folder,
                output_folder,
                compress,
                args.excise_overlap_merge,
                args.excise_overlap_ambig,
                args.excise_region_overlap,
                args.excise_consensus,
                args.excise_maximum_depth,
                args.excise_minimum_ambig,
                args.excise_allowed_distance,
                args.excise_rescue_match,
                prepare_dupes.get(gene.split(".")[0], {}),
                reporter_dupes.get(gene.split(".")[0], {}),
                args.excise_trim_consensus,
            )
            for gene in genes
        ]
        with Pool(args.processes) as pool:
            results = pool.starmap(log_excised_consensus, arguments)
    else:
        results = []
        for gene in genes:
            results.append(
                log_excised_consensus(
                    gene,
                    is_assembly_or_genome,
                    is_genome,
                    input_folder,
                    output_folder,
                    compress,
                    args.excise_overlap_merge,
                    args.excise_overlap_ambig,
                    args.excise_region_overlap,
                    args.excise_consensus,
                    args.excise_maximum_depth,
                    args.excise_minimum_ambig,
                    args.excise_allowed_distance,
                    args.excise_rescue_match,
                    prepare_dupes.get(gene.split(".")[0], {}),
                    reporter_dupes.get(gene.split(".")[0], {}),
                    args.excise_trim_consensus,
                )
            )

    kicked_sequences = 0
    log_output = []
    loci_containing_bad_regions = 0
    kicked_no_resolve = []
    kicked_coverage = []
    has_resolution = []
    no_ambig = []
    rescues = []

    for glog, ghas_ambig, ghas_no_resolution, gcoverage_kick, g_has_resolution, gkicked_seq, grescues in results:
        log_output.extend(glog)
        rescues.extend(grescues)
        if ghas_ambig:
            loci_containing_bad_regions += 1
        kicked_sequences += gkicked_seq
        if ghas_no_resolution:
            kicked_no_resolve.append(ghas_no_resolution)
        if gcoverage_kick:
            kicked_coverage.append(gcoverage_kick)
        if g_has_resolution:
            if ghas_ambig:
                has_resolution.append(g_has_resolution)
            else:
                no_ambig.append(g_has_resolution)



    printv(
        f"{input_folder}: {loci_containing_bad_regions} ambiguous loci found. {len(kicked_coverage)} genes kicked due to low coverage. Kicked {kicked_sequences} sequences total.",
        args.verbose,
    )

    if args.debug:
        with open(log_path, "w") as f:
            f.write("Gene,Cut-Indices,Consensus Sequence\n")
            f.write("\n".join(log_output))
        with open(gene_log_path, "w") as f:
            f.write(f"{len(kicked_no_resolve)} gene(s) kicked due to no seqs with resolution:\n")
            f.write("\n".join(kicked_no_resolve))
            f.write("\n\n")
            f.write(f"{len(kicked_coverage)} gene(s) kicked due to low coverage:\n")
            f.write("\n".join(kicked_coverage))
            f.write("\n\n")
            f.write(f"{len(has_resolution)} gene(s) with resolution and ambig:\n")
            f.write("\n".join(has_resolution))
            f.write("\n\n")
            f.write(f"{len(no_ambig)} gene(s) resolve with no ambig:\n")
            f.write("\n".join(no_ambig))
            f.write("\n\n")
        with open(Path(output_folder, "rescues.txt"), "w") as f:
            f.write("\n".join(rescues))


    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)
    if len(genes) == 0:
        printv("WARNING: No genes in output.", args.verbose, 0)
        return True, True
    return True, loci_containing_bad_regions / len(genes) >= args.majority_excise


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
