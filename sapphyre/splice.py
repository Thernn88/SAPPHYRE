from collections import defaultdict
from itertools import combinations, product
from math import ceil, floor
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
    bio_revcomp,
    join_with_exclusions,
    join_triplets_with_exclusions,
)
from .directional_cluster import cluster_ids, within_distance, node_to_ids, quick_rec
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

def check_covered_bad_regions(nodes, consensus, min_ambiguous, ambig_char='X', max_distance=30):
    x_indices = []
    current_group = []
    regions = []

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
                regions.append((current_group[0], current_group[-1] + 1))
            current_group = [num]

    if current_group:
        if len(current_group) >= min_ambiguous:
            regions.append((current_group[0], current_group[-1] + 1))

    # MERGE
    merge_occured = True
    while merge_occured:
        merge_occured = False
        for (i, region_a), (j, region_b) in combinations(enumerate(regions), 2):
            #gap_between = region_a[1], region_b[0]
            aa_gap_region = floor(region_a[1] / 3), ceil(region_b[0] / 3)
            seqs_in_gap = [(node.start, node.end) for node in nodes if get_overlap(node.start, node.end, aa_gap_region[0], aa_gap_region[1], 1)]
            if seqs_in_gap:
                all_overlap = True
                cr_start, cr_end = seqs_in_gap[0]
                for start, end in seqs_in_gap[1:]:
                    if get_overlap(cr_start, cr_end, start, end, 1) is None:       
                        all_overlap = False
                        break
                    else:
                        cr_start, cr_end = min(cr_start, start), max(cr_end, end)
                    
                if all_overlap:
                    regions[i] = (region_a[0], region_b[1])
                    regions[j] = None
                    merge_occured = True

        regions = [i for i in regions if i]
        
    if regions:
        return regions

    return None, None


def make_duped_consensus(
    raw_sequences: list, threshold: float
) -> str:
    bundled_seqs = [(seq, int(header.split("|")[5])) for header, seq in raw_sequences if header[-1] != "."]
    return dumb_consensus_dupe(bundled_seqs, threshold, 0)


class NODE(Struct):
    header: str
    frame: int
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
                self.nt_sequence[:overlap_coord]
                + node_2.nt_sequence[overlap_coord : node_2.end * 3]
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
                node_2.nt_sequence[:overlap_coord]
                + self.nt_sequence[overlap_coord : self.end * 3]
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
                self.nt_sequence[:overlap_coord] + node_2.nt_sequence[overlap_coord:]
            )

            self.end = node_2.end
        # If node_2 is to the left of self
        else:
            self.sequence = (
                node_2.sequence[:overlap_coord // 3] + self.sequence[overlap_coord // 3:]
            )

            self.nt_sequence = (
                node_2.nt_sequence[:overlap_coord] + self.nt_sequence[overlap_coord:]
            )

            self.start = node_2.start

        if (self.start*3, self.end*3) != find_index_pair(self.nt_sequence, "-"):
            print("ERROR")
            print(self.header)
            print(self.nt_sequence)
            print(self.start, self.end,"/",self.start*3, self.end*3)
            print(find_index_pair(self.nt_sequence, "-"))
            print("ERROR")

        # Save node_2 and the children of node_2 to self
        self.children.append(node_2.header)
        self.children.extend(node_2.children)

    def clone(self):
        return NODE(self.header, self.frame, self.sequence, self.nt_sequence, self.start, self.end, self.children)
    
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
    nodes.sort(key = lambda x: int(x.header.split("|")[5]), reverse=True)
    merged = set()
    node_merge_occured = True
    while node_merge_occured:
        node_merge_occured = False
        for (i, node), (j, node_b) in combinations(enumerate(nodes), 2):
            if i in merged:
                continue
            if j in merged:
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
                node_merge_occured = True
                node.extend(node_b, overlap_coord)

                nodes[j] = None

            # if node.contig_header() == "CONTIG_NODE_111539475|NODE_142545853|NODE_20325325|NODE_17140544|NODE_138951782|NODE_12108535" or node.contig_header() == "CONTIG_NODE_2277847|NODE_11185482|NODE_10065300|NODE_10609136|NODE_13367829|NODE_112796552|NODE_122624644|NODE_118141251|NODE_67349755":
            #     print(i)
    
    nodes = [node for node in nodes if node is not None]

    return nodes
        

def identity(reads_in_region, best_contig, allowed_mismatches):
    # best_contig = contigs_in_region[best_index]
    #best_length = len(best_node.sequence) - best_node.sequence.count("-")
    extend_region = set()
    for i, read in enumerate(reads_in_region):
        
        overlap_coords = get_overlap(best_contig.start * 3, best_contig.end * 3, read.start * 3, read.end * 3, 1)
        if overlap_coords is None:
            extend_region.add(i)
            continue
        kmer_node = read.nt_sequence[overlap_coords[0]:overlap_coords[1]]
        kmer_best = best_contig.nt_sequence[overlap_coords[0]:overlap_coords[1]]
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
        return join_triplets_with_exclusions(sequence, set(), columns)

    return join_with_exclusions(sequence, columns)


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


def do_trim(aa_nodes, cluster_sets, x_positions, ref_consensus, kicked_headers, no_dupes, excise_trim_consensus):
    for cluster_set in cluster_sets:
        
        sub_aa_nodes = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
        sub_x_positions = defaultdict(set)
        aa_sequences = [node.sequence for node in sub_aa_nodes]
        if aa_sequences:
            if no_dupes:
                consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)
            else:
                current_raw_aa = [(node.header, node.sequence) for node in sub_aa_nodes]
                consensus_seq = make_duped_consensus(
                    current_raw_aa, excise_trim_consensus
                )

            cstart, cend = find_index_pair(consensus_seq, "X")

            for i, maj_bp in enumerate(consensus_seq[cstart:cend], cstart):
                if maj_bp != "X":
                    continue

                in_region = []
                out_of_region = False
                for x, node in enumerate(sub_aa_nodes):
                    # within 3 bp
                    if node.start <= i <= node.start + 3:
                        in_region.append((x, node.sequence[i], False))
                    elif node.end - 3 <= i <= node.end:
                        in_region.append((x, node.sequence[i], True))
                    elif i >= node.start and i <= node.end:
                        out_of_region = True
                        
                if not out_of_region and not in_region:
                    continue

                if not out_of_region and in_region:
                    for node_index, bp, on_end in in_region:
                        if bp in ref_consensus[i]:
                            continue
                        
                        if on_end:
                            for x in range(i, sub_aa_nodes[node_index].end):
                                sub_x_positions[sub_aa_nodes[node_index].header].add(x * 3)
                        else:
                            for x in range(sub_aa_nodes[node_index].start, i + 1):
                                sub_x_positions[sub_aa_nodes[node_index].header].add(x * 3)

                if out_of_region and in_region:
                    for node_index, bp, on_end in in_region:
                        if on_end:
                            for x in range(i, sub_aa_nodes[node_index].end):
                                sub_x_positions[sub_aa_nodes[node_index].header].add(x * 3)
                        else:
                            for x in range(sub_aa_nodes[node_index].start, i + 1):
                                sub_x_positions[sub_aa_nodes[node_index].header].add(x * 3)


            #refresh aa
            if sub_x_positions:
                for node in sub_aa_nodes:
                    node.sequence = del_cols(node.sequence, sub_x_positions[node.header])
                    node.start, node.end = find_index_pair(node.sequence, "-")
                    x_positions[node.header].update(sub_x_positions[node.header])

            if no_dupes:
                aa_sequences = [x.sequence for x in sub_aa_nodes if x.header not in kicked_headers]
                consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)
            else:
                current_raw_aa = [(node.header, node.sequence) for node in sub_aa_nodes if node.header not in kicked_headers]
                consensus_seq = make_duped_consensus(
                    current_raw_aa, excise_trim_consensus
                )
  

            for node in sub_aa_nodes:
                i = None
                for poss_i in range(node.start, node.start + 3):
                    if node.sequence[poss_i] != consensus_seq[poss_i]:
                        i = poss_i

                if not i is None:
                    for x in range(node.start , i + 1):
                        x_positions[node.header].add(x * 3)

                i = None
                for poss_i in range(node.end -1, node.end - 4, -1):
                    if node.sequence[poss_i] != consensus_seq[poss_i]:
                        i = poss_i

                if not i is None:
                    for x in range(i, node.end):
                        x_positions[node.header].add(x * 3)


def insert_gaps(input_string, positions, offset):
    input_string = list(input_string)
    
    for coord in positions:
        input_string.insert(offset + coord, "-")

    return "".join(input_string)


def get_combo_results(gt_positions, ag_positions, prev_node, node, FRANKENSTEIN_PENALTY, DELETION_PENALTY, GC_PENALTY, DNA_CODONS, ref_gaps, minimum_bp_for_splice = 15):
    this_results = []
                
    for (gt_size, act_gt_index, gt_index, this_prev_extensions, is_gc), (ag_size, act_ag_index_rev, ag_index_rev, this_node_extensions) in product(gt_positions, ag_positions):
        prev_deletions = []
        node_deletions = []

        prev_nt_seq = list(prev_node.nt_sequence)
        node_nt_seq = list(node.nt_sequence)
        
        for i, let in this_prev_extensions.items():
            prev_nt_seq[i] = let
        for i, let in this_node_extensions.items():
            node_nt_seq[i] = let
        
        for i in range(act_gt_index, act_gt_index + gt_size):
            prev_nt_seq[i] = "-"
            # prev_nt_seq[i] = prev_nt_seq[i].lower()
            
        for i in range(act_ag_index_rev, act_ag_index_rev + ag_size):
            node_nt_seq[i] = "-"
            # node_nt_seq[i] = node_nt_seq[i].lower()

        if prev_node.end * 3 > act_gt_index:
            for x in range(act_gt_index + gt_size, prev_node.end * 3):
                if x in this_prev_extensions:
                    continue
                prev_nt_seq[x] = "-"
                prev_deletions.append(x)
        
        if node.start * 3 < act_ag_index_rev:
            for x in range(node.start * 3, act_ag_index_rev):
                if x in this_node_extensions:
                    continue
                node_nt_seq[x] = "-"
                node_deletions.append(x)
                
        
        prev_nt_seq = [x for x in prev_nt_seq if x != ""]
        node_nt_seq = [x for x in node_nt_seq if x != ""]
                
        node_nt_start, node_nt_end = find_index_pair("".join(node_nt_seq), "-")
        length = node_nt_end - node_nt_start
        right_end_codon = node_nt_start - (3 - (length % 3))
        
        prev_nt_start, prev_nt_end = find_index_pair("".join(prev_nt_seq), "-")
        length = prev_nt_end - prev_nt_start
        left_last_codon = prev_nt_start + length - (length % 3)
        
        prev_gap_insertions = []
        node_gap_insertions = []
        
        if ag_size > 2 and left_last_codon in range(act_ag_index_rev, act_ag_index_rev + ag_size):
            for i in range(act_ag_index_rev + 1, act_ag_index_rev + ag_size - 1):
                if prev_nt_seq[i] != "-":
                    prev_gap_insertions.append(i)
                    prev_nt_seq.insert(i, "-")
                    prev_nt_seq.pop(-1)

        if gt_size > 2 and right_end_codon in range(act_gt_index, act_gt_index + gt_size):
            for i in range(act_gt_index + gt_size - 1, act_gt_index + 1, -1):
                if node_nt_seq[i - 1] != "-":
                    node_gap_insertions.append(i)
                    node_nt_seq.insert(i, "-")
                    node_nt_seq.pop(0)
                
        node_nt_start, node_nt_end = find_index_pair("".join(node_nt_seq), "-")
        length = node_nt_end - node_nt_start
        right_end_codon = node_nt_start - (3 - (length % 3))
        
        prev_nt_start, prev_nt_end = find_index_pair("".join(prev_nt_seq), "-")
        length = prev_nt_end - prev_nt_start
        left_last_codon = prev_nt_start + length - (length % 3)
        
        this_score = 0
        
        if is_gc:
            this_score += GC_PENALTY
        
        right_codon = node_nt_seq[right_end_codon: right_end_codon + 3]
        left_codon = prev_nt_seq[left_last_codon: left_last_codon + 3]

        if right_end_codon == left_last_codon:
            orphan_codon = []
            for i in range(3):
                if left_codon[i] == "-":
                    orphan_codon.append(right_codon[i])
                else:
                    orphan_codon.append(left_codon[i])
                    
            if "-" in orphan_codon:
                continue

            joined = "".join(orphan_codon)
            if joined not in DNA_CODONS:
                this_score += FRANKENSTEIN_PENALTY
            elif DNA_CODONS[joined] == "*":
                continue
        
        else:
            orphan_codon = None
            right_incomplete = "-" in right_codon and right_codon.count("-") != 3
            right_has_ref_gap = right_end_codon // 3 in ref_gaps
                    
            left_incomplete = "-" in left_codon and left_codon.count("-") != 3
            left_has_ref_gap = left_last_codon // 3 in ref_gaps
            
            deletion_possible = 3 - left_codon.count("-") == right_codon.count("-")
            
            # Has incomplete and either right or left is in the ref gap but not both
            if (right_incomplete or left_incomplete) and (((right_has_ref_gap or left_has_ref_gap) and (right_has_ref_gap != left_has_ref_gap)) or deletion_possible):
                if right_has_ref_gap:
                    while right_end_codon // 3 in ref_gaps:
                        for i in range(right_end_codon + 3,right_end_codon, -1):
                            node_gap_insertions.append((i))
                            node_nt_seq.insert(i, "-")
                            node_nt_seq.pop(0)
                        right_end_codon -= 3
                elif left_has_ref_gap:
                    # Extend left by ref gap and see if it meets the right
                    while left_last_codon // 3 in ref_gaps:
                        for i in range(left_last_codon, left_last_codon + 3):
                            prev_gap_insertions.append((i))
                            prev_nt_seq.insert(i, "-")
                            prev_nt_seq.pop(-1)
                        left_last_codon += 3
                elif deletion_possible:
                    if (right_end_codon - 3 - left_last_codon) % 3 == 0:
                        for i in range(0, right_end_codon - left_last_codon):
                            if prev_nt_seq[left_last_codon + i] != "-":
                                prev_gap_insertions.append((left_last_codon + i))
                                prev_nt_seq.insert(left_last_codon + i, "-")
                                prev_nt_seq.pop(-1)
                                this_score += DELETION_PENALTY
                            
                node_nt_start, node_nt_end = find_index_pair("".join(node_nt_seq), "-")
                length = node_nt_end - node_nt_start
                right_end_codon = node_nt_start - (3 - (length % 3))
                
                prev_nt_start, prev_nt_end = find_index_pair("".join(prev_nt_seq), "-")
                length = prev_nt_end - prev_nt_start
                left_last_codon = prev_nt_start + length - (length % 3)
                
                right_codon = node_nt_seq[right_end_codon: right_end_codon + 3]
                left_codon = prev_nt_seq[left_last_codon: left_last_codon + 3]
                
                if right_end_codon == left_last_codon:
                    orphan_codon = []
                    for i in range(3):
                        if left_codon[i] == "-":
                            orphan_codon.append(right_codon[i])
                        else:
                            orphan_codon.append(left_codon[i])
                    if "-" in orphan_codon:
                        continue
                    joined = "".join(orphan_codon)
                    if joined not in DNA_CODONS:
                        this_score += FRANKENSTEIN_PENALTY
                    elif DNA_CODONS[joined] == "*":
                        continue
                else:
                    continue
                
            else:
                continue
                

        distance = node_nt_start - prev_nt_end
        if not (0 <= distance <= 2):
            if prev_nt_end > node_nt_start:
                this_range = range(node_nt_start, prev_nt_end)
            else:
                this_range = range(prev_nt_end, node_nt_start)
            for i in this_range:
                if node_nt_seq[i] == "-" and prev_nt_seq[i] == "-" and i // 3 in ref_gaps:
                    continue
                this_score -= 1
                
        if len(prev_nt_seq) - prev_nt_seq.count("-") < minimum_bp_for_splice:
            continue
        
        if len(node_nt_seq) - node_nt_seq.count("-") < minimum_bp_for_splice:
            continue

        this_results.append((is_gc, this_score, gt_size, act_gt_index, gt_index, ag_size, act_ag_index_rev, ag_index_rev, prev_deletions, node_deletions, prev_gap_insertions, node_gap_insertions, this_prev_extensions, this_node_extensions, left_last_codon, right_end_codon, orphan_codon))

    return this_results

def find_gt_ag(prev_node, node, prev_start_index, prev_end_index, node_start_index, node_end_index, DNA_CODONS, prev_og, node_og, prev_nt_seq, node_nt_seq, ref_gaps):
    act_gt_index = None
    gt_index = None
    act_ag_index_rev = None
    ag_index_rev = None
    prev_extensions = {}
    node_extensions = {}
    gt_positions = []
    ag_positions = []
    
    EXTEND_WINDOW = 45
    
    overlapping_region = get_overlap(prev_node.start * 3, prev_node.end * 3, node.start * 3, node.end * 3, 1)
    overlap_amount = 0
    if overlapping_region:
        overlap_amount = overlapping_region[1] - overlapping_region[0]
        
    x = -1
    act_offset = (prev_node.end * 3) - overlap_amount - EXTEND_WINDOW
    prev_start_scan = prev_end_index - overlap_amount - EXTEND_WINDOW
    if prev_start_scan < prev_start_index:
        act_offset = (prev_node.start * 3) + 3
        prev_start_scan = prev_start_index + 3
    
    for i in range(prev_start_scan, len(prev_og)):
        x += 1
        prev_act_coord = act_offset + x
            
        # Get last codon
        
        if prev_act_coord <= 0:
            continue
        
        if x > 0 and i >= 3 and i < len(prev_og) - 3 and x % 3 == 0 and i > prev_end_index:
            last_codon = prev_og[i - 3: i]
            
            if last_codon not in DNA_CODONS:
                #Bandaid TODO
                continue
                
            
            if DNA_CODONS[last_codon] == "*":
                break

        if i + 1 >= len(prev_og) or prev_act_coord >= len(prev_nt_seq) - 1:
            break
        
        if prev_nt_seq[prev_act_coord + 1] == "-":
            prev_extensions[prev_act_coord + 1] = prev_og[i + 1]

        scan_index = 0
        if prev_og[i] == "G":
            while True:
                scan_index += 1
                if i + scan_index >= len(prev_og):
                    break
                
                if prev_og[i + scan_index] == "T" or prev_og[i + scan_index] == "C":
                    act_gt_index = prev_act_coord
                    gt_index = i + 1
                    gt_size = scan_index + 1
                    gt_positions.append((gt_size, act_gt_index, gt_index, copy.deepcopy(prev_extensions), prev_og[i + scan_index] == "C"))

                if prev_og[i + scan_index] != "-":
                    break

                if (prev_act_coord + scan_index)//3 not in ref_gaps:
                    break            
    # Iterate in reverse from the start of the kmer to the start of the original sequence
    x = -1
    
    offset_distance = 0
    if node.start < prev_node.start:
        offset_distance = (prev_node.start - node.start) * 3
    
    start_scan = node_start_index + overlap_amount + EXTEND_WINDOW + offset_distance
    act_offset = (node.start * 3) + overlap_amount + EXTEND_WINDOW + offset_distance
    if start_scan > node_end_index:
        start_scan = node_end_index
        act_offset = (node.end * 3)
        
    for i in range(start_scan - 1, -1, -1):
        x += 1
        node_act_coord = act_offset - 1 - x
        # Get last codon
        
        if node_act_coord >= len(node_nt_seq) - 1 or i >= len(node_og) - 1:
            continue
        
        if node_act_coord < 0 or node_act_coord < (prev_node.start * 3):
            break
        
        if x > 0 and i < len(node_og) - 3 and x % 3 == 0 and i < node_start_index:
            last_codon = node_og[i + 1: i + 4]
            if last_codon not in DNA_CODONS:
                #Bandaid TODO
                continue
            
            
            if DNA_CODONS[last_codon] == "*":
                break
        
        # Check if the next nucleotide (i + 1) is "A" and the current is "G"
        scan_index = 0
        
        if node_nt_seq[node_act_coord] == "-":
            node_extensions[node_act_coord] = node_og[i]
        
        if node_og[i] == "A":
            while True:
                scan_index += 1

                if i + scan_index >= len(node_og):
                    break
                
                if node_og[i + scan_index] == "G":
                    act_ag_index_rev = node_act_coord
                    ag_index_rev = i
                    ag_size = scan_index + 1
                    ag_positions.append((ag_size, act_ag_index_rev, ag_index_rev, copy.deepcopy(node_extensions)))
                if node_og[i + scan_index] != "-":
                    break

                if (node_act_coord + scan_index)//3 not in ref_gaps:
                    break

    return gt_positions, ag_positions


def splice_combo(add_results,
                 print_extra,
                 formed_seqs,
                 this_result,
                 prev_node,
                 node,
                 prev_og,
                 node_og,
                 DNA_CODONS,
                 scan_log,
                 combo_log,
                 prev_start_index,
                 node_start_index,
                 kmer,
                 kmer_internal_gaps,
                 prev_internal_gaps,
                 gff_coords
                 ):
    # if "223372" not in node.header:
    #     return
    
    prev_nt_seq = list(prev_node.nt_sequence)
    node_nt_seq = list(node.nt_sequence)
    
    prev_seq = list(prev_node.sequence)
    node_seq = list(node.sequence)
    
    is_gc, this_score, gt_size, act_gt_index, gt_index, ag_size, act_ag_index_rev, ag_index_rev, prev_deletions, node_deletions, final_prev_insertions, final_node_insertions, final_prev_extensions, final_node_extensions, left_last_codon, right_end_codon, orphan_codon = this_result
    for i, let in final_prev_extensions.items():
        prev_nt_seq[i] = let
    for i, let in final_node_extensions.items():
        node_nt_seq[i] = let
    for i in range(act_gt_index, act_gt_index + gt_size):
        if add_results:
            prev_seq[i//3] = "-"
        prev_nt_seq[i] = "-"
    for i in range(act_ag_index_rev, act_ag_index_rev + ag_size):
        if add_results:
            node_seq[i//3] = "-"
        node_nt_seq[i] = "-"
    for x in prev_deletions:
        if x in final_prev_extensions:
            continue
        if add_results:
            prev_seq[x//3] = "-"
        prev_nt_seq[x] = "-"
    for x in node_deletions:
        if x in final_node_extensions:
            continue
        if add_results:
            node_seq[x//3] = "-"
        node_nt_seq[x] = "-"
        
    prev_insertions = 0
    node_insertions = 0
        
    for i in final_prev_insertions:
        # if add_results:
            # if i % 3 == 0:
            #     gap_insertions_aa[prev_node.header].append((i//3, -1))
        prev_nt_seq.insert(i, "-")
        prev_nt_seq.pop(-1)
        prev_insertions += 1
    
    for i in final_node_insertions:
        # if add_results:
            # if i % 3 == 0:
            #     gap_insertions_aa[node.header].append((i//3, 0))
        node_nt_seq.insert(i, "-")
        node_nt_seq.pop(0)
        node_insertions += 1
        
    node_hit = node_og[ag_index_rev: (node_start_index + len(kmer) + len(kmer_internal_gaps))]
    prev_hit = prev_og[prev_start_index: gt_index + gt_size - 1]
        
    if orphan_codon:
        joined = "".join(orphan_codon)
        if joined not in DNA_CODONS:
            return None, None
        
        orphan_aa = DNA_CODONS[joined]
        if add_results:
            prev_seq[left_last_codon//3] = orphan_aa
            node_seq[right_end_codon//3] = orphan_aa

        prev_og = list(prev_og)
        node_og = list(node_og)

        for i, x in enumerate(range(left_last_codon, left_last_codon + 3)):
            prev_nt_seq[x] = orphan_codon[i]
            if add_results:
                prev_og[prev_start_index + x - (prev_node.start * 3) - prev_insertions] = orphan_codon[i]

        for i, x in enumerate(range(right_end_codon, right_end_codon + 3)):
            node_nt_seq[x] = orphan_codon[i]
            if add_results:
                node_og[node_start_index + x - (node.start * 3) - node_insertions] = orphan_codon[i]


        prev_og = "".join(prev_og)
        node_og = "".join(node_og)
        
        formed_seqs[prev_node.header] = prev_og.replace("-", "")
        formed_seqs[node.header] = node_og.replace("-", "")
        
    node_nt_seq = "".join(node_nt_seq)
    prev_nt_seq = "".join(prev_nt_seq)

    node_start, node_end = find_index_pair(node_nt_seq, "-")
    prev_start, prev_end = find_index_pair(prev_nt_seq, "-")
    
    
    if node_start % 3 != 0:
        node_start -= node_start % 3
    
    if prev_end % 3 != 0:
        prev_end += 3 - (prev_end % 3)
    
    if add_results:
        for i in range((prev_node.end * 3) - 3, prev_end, 3):
            codon = prev_nt_seq[i:i+3]
            if codon in DNA_CODONS:
                prev_seq[i//3] = DNA_CODONS[codon]
            
        for i in range(node_start, node.start * 3, 3):
            codon = node_nt_seq[i:i+3]
            if codon in DNA_CODONS:
                node_seq[i//3] = DNA_CODONS[codon]

    # Extend one codon to the left and right for debug to show end of sequence
    if prev_start >= 3:
        prev_start -= 3
    if node_end + 3 <= len(node_nt_seq):
        node_end += 3
        
    node_region = node.nt_sequence[prev_start: node_end]
    node_region_start, _ = find_index_pair(node_region, "-")
    
    final_prev_start, final_prev_end = find_index_pair(prev_nt_seq, "-")
    final_node_start, final_node_end = find_index_pair(node_nt_seq, "-")
    
    smallest_change = min(abs(final_prev_start - (prev_node.start*3)) + abs(final_prev_end - (prev_node.end*3)), abs(final_node_start - (node.start*3)) + abs(final_node_end - (node.end*3)))
    if print_extra: 
        scan_log.append(f">{prev_node.header}_orf")

        scan_log.append(prev_og[prev_start_index - 3 :][:node_end])  
        scan_log.append(f">{node.header}_orf")  
        scan_log.append(node_og[node_start_index - node_region_start :][:node_end])  
        scan_log.append(f">{prev_node.header}_excise_output")
        scan_log.append(prev_node.nt_sequence[prev_start: node_end])
        scan_log.append(f">{node.header}_excise_output")
        scan_log.append(node_region)
        scan_log.append("")
                
    scan_log.append(f">{prev_node.header}_spliced_s{this_score}_c{smallest_change}")
    scan_log.append(prev_nt_seq[prev_start: node_end])
    scan_log.append(f">{node.header}_spliced_s{this_score}_c{smallest_change}")
    scan_log.append(node_nt_seq[prev_start: node_end])
    scan_log.append("")    
    
    
    # lowercase
    prev_hit = prev_hit[:-gt_size] + prev_hit[-gt_size:].lower()
    node_hit = node_hit[:ag_size].lower() + node_hit[ag_size:]

    scan_log.append(f">{prev_node.header}_orf_scan")
    scan_log.append((("-" * (prev_node.start * 3)) + prev_hit)[prev_start: node_end])
    scan_log.append(f">{node.header}_orf_scan")
    scan_log.append((("-" * ((node.end * 3) - len(node_hit))) + node_hit)[prev_start: node_end] )
    scan_log.append("") 
    
    if add_results:
        combo_log.append(f"{prev_node.header} vs {node.header}: {this_score} " + ("GC-AG" if is_gc else "GT-AG"))
        
        prev_kmer = prev_hit.replace("-", "")
        gff_coord_prev = (
            prev_start_index + 1,
            prev_start_index + len(prev_kmer) - 2
        )

        gff_coord_node = (
            ag_index_rev + 3,
            node_start_index + len(kmer)
        )
        
        if prev_node.header in gff_coords:
            existing_gff = gff_coords[prev_node.header]
            gff_coord_prev = (existing_gff[0], gff_coord_prev[1])
            gff_coords[prev_node.header] = gff_coord_prev
        else:
            gff_coords[prev_node.header] = gff_coord_prev
            
        if node.header in gff_coords:
            existing_gff = gff_coords[node.header]
            gff_coord_node = (gff_coord_node[0], existing_gff[1])
            gff_coords[node.header] = gff_coord_node
        else:
            gff_coords[node.header] = gff_coord_node
        
        if prev_node.frame < 0:
            og_length = len(prev_og) - len(prev_internal_gaps) + 1
            gff_coord_prev = (
                og_length - gff_coord_prev[1],
                og_length - gff_coord_prev[0]
            )
            
        if node.frame < 0:
            og_length = len(node_og) - len(kmer_internal_gaps) + 1
            gff_coord_node = (
                og_length - gff_coord_node[1],
                og_length - gff_coord_node[0]
            )
            
        prev_node.nt_sequence = prev_nt_seq
        node.nt_sequence = node_nt_seq
        
        node.sequence = "".join(node_seq)
        prev_node.sequence = "".join(prev_seq)
        
        prev_node.start, prev_node.end = final_prev_start//3, final_prev_end//3
        node.start, node.end = final_node_start//3, final_node_end//3
        
        return (gff_coord_prev, gff_coord_node), smallest_change
    return None, smallest_change


def edge_check(aa_nodes, cluster_sets, kicked_headers, log_output, min_overlap_amount = 0.85, min_matching_percent = 0.95):
    for cluster_set in cluster_sets:
        aa_subset = [node for node in aa_nodes if (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
        if aa_subset[0].frame < 0:
            aa_subset.sort(key=lambda x: x.start, reverse=True)
        else:
            aa_subset.sort(key=lambda x: x.start)
            
        if len(aa_subset) < 2:
            continue
        
        left_most = aa_subset[0]
        left_next = aa_subset[1]
       
        left_overlap = get_overlap(left_most.start, left_most.end, left_next.start, left_next.end, 1)
        if left_overlap:
            left_overlap_amt = left_overlap[1] - left_overlap[0]
            left_overlap_percent = left_overlap_amt / (left_most.end - left_most.start)
            if left_overlap_percent >= min_overlap_amount:
                left_most_kmer = left_most.sequence[left_overlap[0]: left_overlap[1]]
                left_next_kmer = left_next.sequence[left_overlap[0]: left_overlap[1]]
                
                distance = constrained_distance(left_most_kmer, left_next_kmer)
                matching_percent = distance / len(left_most_kmer)
                
                
                if matching_percent <= min_matching_percent:
                    if left_next.start == left_most.start and len(left_next.sequence) < len(left_most.sequence):
                        kicked_headers.add(left_next.header)
                        log_output.append(f"Kicked {left_next.header} as it greatly overlapped on the edge")
                    else:
                        kicked_headers.add(left_most.header)
                        log_output.append(f"Kicked {left_most.header} as it greatly overlapped on the edge")
                    
        right_most = aa_subset[-1]
        right_next = aa_subset[-2]
        
        right_overlap = get_overlap(right_most.start, right_most.end, right_next.start, right_next.end, 1)
        if right_overlap:
            right_overlap_amt = right_overlap[1] - right_overlap[0]
            right_overlap_percent = right_overlap_amt / (right_most.end - right_most.start)
            if right_overlap_percent >= min_overlap_amount:
                right_most_kmer = right_most.sequence[right_overlap[0]: right_overlap[1]]
                right_next_kmer = right_next.sequence[right_overlap[0]: right_overlap[1]]
                
                distance = constrained_distance(right_most_kmer, right_next_kmer)
                matching_percent = distance / len(right_most_kmer)
                
                if matching_percent <= min_matching_percent:
                    if right_next.end == right_most.end and len(right_next.sequence) < len(right_most.sequence):
                        kicked_headers.add(right_next.header)
                        log_output.append(f"Kicked {right_next.header} as it greatly overlapped on the edge")
                    else:
                        kicked_headers.add(right_most.header)
                        log_output.append(f"Kicked {right_most.header} as it greatly overlapped on the edge")
                    
            
        


def log_excised_consensus(
    verbose: int,
    gene: str,
    input_path: Path,
    output_path: Path,
    compress_intermediates: bool,
    excise_consensus,
    excise_minimum_ambig,
    excise_rescue_match,
    no_dupes,
    head_to_seq,
    original_coords,
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
    printv(f"Processing {gene}", verbose, 2)
    log_output = []
    scan_log = []
    combo_log = []
    multi_log = []
    this_rescues = []
    rescue_jank_log = set()

    aa_in = input_path.joinpath("aa", gene)
    aa_out = output_path.joinpath("aa", gene)

    nt_in = input_path.joinpath("nt", gene.replace(".aa.", ".nt."))
    nt_out = output_path.joinpath("nt", gene.replace(".aa.", ".nt."))

    x_positions = defaultdict(set)

    bp_count = lambda x: len(x) - x.count("-")

    ref_lens = []
    aa_nodes = []
    reference_cluster_data = set()
    ref_consensus = defaultdict(list)
    ref_gaps = set()
    for header, seq in parseFasta(str(aa_in)):
        if header.endswith('.'):
            start, end = find_index_pair(seq, "-")
            for i, bp in enumerate(seq[start:end], start):
                ref_consensus[i].append(bp)
                
            ref_lens.append(bp_count(seq))
            continue

        frame = int(header.split("|")[4])
        aa_nodes.append(NODE(header, frame, seq, None, *find_index_pair(seq, "-"), []))

    ref_gap_percent = 0.7
    gap_set = {"-"}
    for i, lets in ref_consensus.items():
        if gap_set != set(lets):
            reference_cluster_data.add(i)
            
        if lets.count("-") / len(lets) > ref_gap_percent:
            ref_gaps.add(i)

    kicked_headers = set()

    cluster_sets = [None]
    get_id = lambda header: header.split("|")[3].replace("NODE_","")
    ids = []
    for node in aa_nodes:
        if node.header not in kicked_headers:
            start, end = find_index_pair(node.sequence, "-")
            ids.append(quick_rec(node.header.split("|")[3], node.frame, node.sequence, start, end))

    max_gap_size = round(len(aa_nodes[0].sequence) * 0.3) # Half MSA length

    clusters, _ = cluster_ids(ids, 100, max_gap_size, reference_cluster_data) #TODO: Make distance an arg

    if clusters:
        cluster_sets = [set(range(a, b+1)) for a, b, _ in clusters]
        
    edge_check(aa_nodes, cluster_sets, kicked_headers, log_output)

    head_to_node = {}
    for node in aa_nodes:

        node.sequence = del_cols(node.sequence, x_positions[node.header])
        node.start, node.end = find_index_pair(node.sequence, "-")

        head_to_node[node.header] = node
        node_kmer = node.sequence[node.start:node.end]
        data_bp = bp_count(node_kmer) * 3
        if data_bp < 15:
            log_output.append(f"Kicking {node.header} due to < 15 bp after trimming")
            kicked_headers.add(node.header)
            
    aa_nodes = [node for node in aa_nodes if node.header not in kicked_headers]

    raw_sequences = {header: del_cols(seq, x_positions[header], True) for header, seq in parseFasta(str(nt_in)) if header not in kicked_headers}
    for node in aa_nodes:
        nt_seq = raw_sequences[node.header]
        node.nt_sequence = nt_seq
    
    true_cluster_raw = []

    for header in raw_sequences:
        true_cluster_raw.append((int(header.split("|")[3].split("&&")[0].split("_")[1]), header))

    true_cluster_raw.sort(key = lambda x: x[0])
    before_true_clusters = cluster(true_cluster_raw, true_cluster_threshold)

    # Search for slices of the consensus seq with a high ratio of 'X' to total characters
    has_region = True
    had_region = False
    recursion_max = 5
    last_region = {}
    consensus_seq = ""
    while has_region:
        has_region = False
        
        for cluster_i, cluster_set in enumerate(cluster_sets):

            aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]

            sequences = [node.nt_sequence for node in aa_subset]
            if not sequences:
                break

            if no_dupes:
                consensus_seq = dumb_consensus(sequences, excise_consensus, 0)
            else:
                nt_sequences = [(node.header, node.nt_sequence) for node in aa_subset]
                consensus_seq = make_duped_consensus(
                    nt_sequences, excise_consensus
                )    

            consensus_seq = convert_consensus(sequences, consensus_seq)
            regions = check_covered_bad_regions(aa_subset, consensus_seq, excise_minimum_ambig)
            if regions[0] is None:
                continue
            for ri, (region_start, region_end) in enumerate(regions):
                if ri == 0:
                    if region_start == last_region.get(cluster_i, -1):
                        recursion_max -= 1
                        if recursion_max <= 0:
                            region_start = None
                            continue
                    last_region[cluster_i] = region_start

                if region_start:
                    has_region = True
                    had_region = True
                    break
                
        if recursion_max <= 0:
            break
        
    # do_trim(aa_nodes, cluster_sets, x_positions, ref_consensus, kicked_headers, no_dupes, excise_trim_consensus)
    # for node in aa_nodes:
    #     if node.header in kicked_headers:
    #         continue
        
    #     node.sequence = del_cols(node.sequence, x_positions[node.header])
    #     node.nt_sequence = del_cols(node.nt_sequence, x_positions[node.header], True)
    #     node.start, node.end = find_index_pair(node.sequence, "-")
          
    move_dict = defaultdict(list)
    debug_out = []
    merged = []

    # Detect ref gap regoins
    integers = list(ref_gaps)
    if integers:
        integers.sort()
        start = integers[0]
        end = integers[0]

        for num in integers[1:]:
            if num == end + 1:
                end = num
            else:
                merged.append((start, end + 1))
                start = num
                end = num

        merged.append((start, end + 1))
        
    for cluster_i, cluster_set in enumerate(cluster_sets):
        aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
        
        for start, end in merged:
            gap_size = end - start
            if gap_size == 1:
                continue
            
            for node in aa_nodes:
                overlap = get_overlap(start, end, node.start, node.end, 1)

                if overlap is None:
                    continue
                
                amount = overlap[1] - overlap[0]

                gap_region = node.sequence[start: end]
                gaps_in_region = gap_region.count("-")
                if gaps_in_region == len(gap_region) or gaps_in_region == 0:
                    continue
                
                gap_region_on_seq = node.sequence[overlap[0]: overlap[1]]
                internal_gap_offset = find_index_pair(gap_region_on_seq, "-")
                data_in_gap = gap_region_on_seq[internal_gap_offset[0]:internal_gap_offset[1]]
                # right side or left side of the gap
                
                this_seq = list(node.sequence)
                this_nt_seq = list(node.nt_sequence)

                if overlap[0] == start and internal_gap_offset[0] == 0: # Left side
                    match = True
                    for i, let in enumerate(data_in_gap, overlap[0]+gap_size):
                        if node.sequence[i] != "-" and let != node.sequence[i]:
                            match = False
                            break
                        if not let in ref_consensus[i]:
                            match = False
                            break

                    if match:
                        for new_i, original_i in enumerate(range(overlap[0], overlap[0]+len(data_in_gap)), overlap[0] + gap_size):
                            this_seq[new_i] = this_seq[original_i]
                            this_seq[original_i] = "-"
                            
                            this_nt_seq[(new_i*3):(new_i*3)+3] = this_nt_seq[(original_i*3):(original_i*3)+3]
                            this_nt_seq[(original_i*3):(original_i*3)+3] = ["-","-","-"]
                            
                            move_dict[node.header].append((original_i, new_i))          
                else: # Right side
                    match = True
                    for i, let in enumerate(data_in_gap, start-len(data_in_gap)):
                        
                        if node.sequence[i] != "-" and let != node.sequence[i]:
                            match = False
                            break
                        
                        if not let in ref_consensus[i]:
                            match = False
                            break
                    
                    if match:
                        for new_i, original_i in enumerate(range(overlap[0] + internal_gap_offset[0], overlap[0] + internal_gap_offset[1]), start - len(data_in_gap)):
                            this_seq[new_i] = this_seq[original_i]
                            this_seq[original_i] = "-"
                            
                            
                            this_nt_seq[(new_i*3):(new_i*3)+3] = this_nt_seq[(original_i*3):(original_i*3)+3]
                            this_nt_seq[(original_i*3):(original_i*3)+3] = ["-","-","-"]
                            
                            move_dict[node.header].append((original_i, new_i))

                node.sequence = "".join(this_seq)
                node.nt_sequence = "".join(this_nt_seq)
                node.start, node.end = find_index_pair(node.sequence, "-")

    if had_region:
        after_data = []
        for node in aa_nodes:
            if node.header in kicked_headers:
                continue
            after_data.append((node.header.split("|")[3].split("&&")[0].split("_")[1], node.header))

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
        
        if current_cluster:
            after_true_clusters.append(current_cluster)

        after_true_cluster = set(max(after_true_clusters, key=len))
        matches = []
        for i, before_cluster in enumerate(before_true_clusters):
            matches.append((len(after_true_cluster.intersection(set(before_cluster))) / len(before_cluster), i))
        matches.sort(key=lambda x: x[0], reverse= True)
        
        best_match = max(matches, key=lambda x: x[0])[1]
        this_rescues = [f"Rescued in gene: {gene}"]
        for header in before_true_clusters[best_match]:
            if header in kicked_headers:

                node = head_to_node[header]
                
                if not node.nt_sequence:
                    continue
                
                if len(node.nt_sequence) - node.nt_sequence.count("-") < 15:
                    continue
                
                this_sequence = node.sequence
                start, end = find_index_pair(this_sequence, "-")
                matching_char = 0
                for i, let in enumerate(this_sequence[start: end], start):
                    if let in ref_consensus[i]:
                        matching_char += 1

                if matching_char / (end - start) < excise_rescue_match:
                    continue

                kicked_headers.remove(header)
                this_rescues.append(header)
                    
    aa_nodes = [node for node in aa_nodes if node.header not in kicked_headers]
    ends = {}
    gff_out = defaultdict(dict)
    gff_coords = {}
    formed_seqs = {}
    has_exisiting_result = defaultdict(dict)

    DNA_CODONS = {
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "TGT": "C",
        "TGC": "C",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "TTT": "F",
        "TTC": "F",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
        "CAT": "H",
        "CAC": "H",
        "ATA": "I",
        "ATT": "I",
        "ATC": "I",
        "AAA": "K",
        "AAG": "K",
        "TTA": "L",
        "TTG": "L",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "ATG": "M",
        "AAT": "N",
        "AAC": "N",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGA": "R",
        "AGG": "R",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "AGT": "S",
        "AGC": "S",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "TGG": "W",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGA": "*",
        "---": "-",
    }
    gene_name = gene.split(".")[0]
    FRANKENSTEIN_PENALTY = -20
    GC_PENALTY = -20
    SIMILARITY_SKIP = 0.95
    DELETION_PENALTY = 0
    int_first_id = lambda x: int(x.split("_")[0])
    for cluster_i, cluster_set in enumerate(cluster_sets):
        aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
        aa_subset.sort(key = lambda x: x.start)
        for prev_node, node in combinations(aa_subset, 2):
            overlapping_coords = get_overlap(node.start, node.end, prev_node.start, prev_node.end, -40)
            if overlapping_coords:
                amount = overlapping_coords[1] - overlapping_coords[0]
                this_flipped = ""
                if amount > 1:
                    early_start = min(prev_node.start, node.start) * 3
                    late_end = max(prev_node.end, node.end) * 3
                    
                    distance = constrained_distance(prev_node.nt_sequence[early_start : late_end], node.nt_sequence[early_start : late_end])
                    percent_matching = 1 - (distance / (late_end - early_start))
                    if percent_matching >= SIMILARITY_SKIP:
                        continue
                    
                    if amount == (node.end - node.start) or amount == (prev_node.end - prev_node.start):
                        #Containment
                        
                        node_bp_indices = [i for i,let in enumerate(node.nt_sequence) if let != "-"]
                        node_avg_bp_index = sum(node_bp_indices) / len(node_bp_indices)
                        
                        prev_bp_indices = [i for i,let in enumerate(prev_node.nt_sequence) if let != "-"]
                        prev_avg_bp_index = sum(prev_bp_indices) / len(prev_bp_indices)
                        
                        if node_avg_bp_index < prev_avg_bp_index:
                            prev_node, node = node, prev_node
                            this_flipped = " (flipped)"
                            
                scan_log.append("")
                scan_log.append("")
                scan_log.append(f"Comparing {prev_node.header} vs {node.header}" + this_flipped)

                prev_kmer = prev_node.nt_sequence[prev_node.start * 3 : prev_node.end * 3]
                prev_internal_gaps = [i for i, let in enumerate(prev_kmer) if let == "-"]
                prev_kmer = prev_kmer.replace("-","")
                
                kmer = node.nt_sequence[node.start * 3 : node.end * 3]
                kmer_internal_gaps = [i for i, let in enumerate(kmer) if let == "-"]
                kmer = kmer.replace("-","")
                
                if prev_node.header not in formed_seqs:
                    children = list(map(int_first_id, prev_node.header.split("|")[3].replace("NODE_", "").split("&&")))
                    prev_og = head_to_seq[children[0]]
                    for i, child in enumerate(children):
                        if i == 0:
                            continue
                        prev_og += head_to_seq[child][250:]

                    if prev_node.frame < 0:
                        prev_og = bio_revcomp(prev_og)
                
                    formed_seqs[prev_node.header] = prev_og
                else:
                    prev_og = formed_seqs[prev_node.header]
                    
                if node.header not in formed_seqs:
                    children = list(map(int_first_id, node.header.split("|")[3].replace("NODE_", "").split("&&")))
                    node_og = head_to_seq[children[0]]
                    for i, child in enumerate(children):
                        if i == 0:
                            continue
                        node_og += head_to_seq[child][250:]

                    if node.frame < 0:
                        node_og = bio_revcomp(node_og)
                    formed_seqs[node.header] = node_og
                else:
                    node_og = formed_seqs[node.header]
                
                prev_start_index = prev_og.find(prev_kmer)
                if prev_start_index == -1:
                    continue

                prev_og = insert_gaps(prev_og, prev_internal_gaps, prev_start_index)
                prev_end_index = prev_start_index + len(prev_kmer) + len(prev_internal_gaps)

                node_start_index = node_og.find(kmer)
                if node_start_index == -1:
                    continue

                
                node_og = insert_gaps(node_og, kmer_internal_gaps, node_start_index)
                node_end_index = node_start_index + len(kmer) + len(kmer_internal_gaps)

                prev_nt_seq = list(prev_node.nt_sequence)
                node_nt_seq = list(node.nt_sequence)
                
                gt_positions, ag_positions = find_gt_ag(
                    prev_node, 
                    node, 
                    prev_start_index, 
                    prev_end_index, 
                    node_start_index, 
                    node_end_index, 
                    DNA_CODONS, 
                    prev_og, 
                    node_og, 
                    prev_nt_seq, 
                    node_nt_seq, 
                    ref_gaps
                    )

                # if "2A|AglaOr12CTE|splice_fix|NODE_343534&&343535|-3|1" not in prev_node.header:
                #     continue

                this_results = get_combo_results(gt_positions, ag_positions, prev_node, node, FRANKENSTEIN_PENALTY, GC_PENALTY, DELETION_PENALTY, DNA_CODONS, ref_gaps)
                splice_found = False
                this_best_splice = None
                if this_results:
                    this_best_score = max(i[1] for i in this_results)
                    highest_results = [i for i in this_results if i[1] == this_best_score]
                    if len(highest_results) > 1:
                        best_change = None
                        for i, result in enumerate(highest_results):
                            _, smallest_coords_change = splice_combo(False, i == 0, formed_seqs, result, prev_node, node, prev_og, node_og, DNA_CODONS, multi_log, [], prev_start_index, node_start_index, kmer, kmer_internal_gaps, prev_internal_gaps, gff_coords)
                            # print(smallest_coords_change)
                            if smallest_coords_change is not None and (best_change is None or smallest_coords_change < best_change):
                                best_change = smallest_coords_change
                                this_best_splice = result
                    else:
                        this_best_splice = highest_results[0]
                    
                    if this_best_splice is not None:
                        if "Right" in has_exisiting_result[prev_node.header]:
                            if has_exisiting_result[prev_node.header]["Right"] >= this_best_splice[1]:
                                continue
                        if "Left" in has_exisiting_result[node.header]:
                            if has_exisiting_result[node.header]["Left"] >= this_best_splice[1]:
                                continue
                        
                        gff, _ = splice_combo(True, True, formed_seqs, this_best_splice, prev_node, node, prev_og, node_og, DNA_CODONS, scan_log, combo_log, prev_start_index, node_start_index, kmer, kmer_internal_gaps, prev_internal_gaps, gff_coords)
                        if gff:
                            has_exisiting_result[prev_node.header]["Right"] = this_best_splice[1]
                            has_exisiting_result[node.header]["Left"] = this_best_splice[1]
                            
                            splice_found = True
                            prev_gff, node_gff = gff
                            
                            prev_id = get_id(prev_node.header)
                            tup = original_coords.get(prev_id.split("&&")[0].split("_")[0], None)
                            if tup:
                                parent, chomp_start, _, input_len, _ = tup
                                
                                prev_start = prev_gff[0] + chomp_start
                                prev_end = prev_gff[1] + chomp_start
                                
                                if parent not in ends:
                                    ends[parent] = input_len
                                    
                                strand = "+" if prev_node.frame > 0 else "-"
                                gff_out[parent][prev_id] = ((prev_start), f"{parent}\tSapphyre\texon\t{prev_start}\t{prev_end}\t.\t{strand}\t.\tDescription={prev_id};Name={gene_name};ID={gene_name};Note={prev_node.frame};")
                                
                            node_id = get_id(node.header)
                            tup = original_coords.get(node_id.split("&&")[0].split("_")[0], None)
                            if tup:
                                parent, chomp_start, _, input_len, _ = tup
                                
                                node_start = node_gff[0] + chomp_start
                                node_end = node_gff[1] + chomp_start
                                
                                if parent not in ends:
                                    ends[parent] = input_len
                                    
                                strand = "+" if node.frame > 0 else "-"
                                gff_out[parent][node_id] = ((node_start), f"{parent}\tSapphyre\texon\t{node_start}\t{node_end}\t.\t{strand}\t.\tDescription={node_id};Name={gene_name};ID={gene_name};Note={node.frame};")
                if True and not splice_found:    
                    scan_log.append(f">{prev_node.header}_orf")
                    # print(prev_start_index, node_end_index)
                    # input()
                    node_start, node_end = find_index_pair("".join(node_nt_seq), "-")
                    prev_start, prev_end = find_index_pair("".join(prev_nt_seq), "-")
                    prev_start -= 3
                    add_gaps = 0
                    if prev_start < 0:
                        prev_start = 0
                        add_gaps = 3
                    node_end += 3
                    node_region = node.nt_sequence[prev_start: node_end]
                    node_region_start, _ = find_index_pair(node_region, "-")
                    
                    
                    this_prev = prev_og[prev_start_index - 3 :][:node_end]
                    
                    
                    if node_start_index > prev_start_index:
                        this_node = node_og[prev_start_index - 3 :][:node_end]
                        gaps_needed = node_region_start - (node_start_index - prev_start_index) - 3
                    else:
                        gaps_needed = 0
                        this_node = node_og[node_start_index - node_region_start :][:node_end]
                    
                    scan_log.append(this_prev)  
                    scan_log.append(f">{node.header}_orf")  
                    scan_log.append(('-' * gaps_needed) + this_node)  
                    scan_log.append(f">{prev_node.header}_excise_output")
                    # scan_log.append("".join(prev_nt_seq))
                    # scan_log.append(prev_node.nt_sequence)
                    scan_log.append(("-" * add_gaps) + prev_node.nt_sequence[prev_start: node_end])
                    scan_log.append(f">{node.header}_excise_output")
                    scan_log.append(node_region)
                    scan_log.append("")
                    scan_log.append("No result found.")
                    scan_log.append("")
                    
                    
    aa_output = [(header, seq) for header, seq in parseFasta(str(aa_in)) if header.endswith(".")]
    for node in aa_nodes:
        if node.header in kicked_headers:
            continue
        header = node.header
        seq = node.sequence
        
        this_id = get_id(header)
        tup = original_coords.get(this_id.split("&&")[0].split("_")[0], None)
        if tup:
            parent, _, _, input_len, _ = tup
            if parent not in ends:
                ends[parent] = input_len
                
            if this_id not in gff_out[parent]:
                gff_out[parent][this_id] = None
                
        aa_output.append((header, seq))
    
    internal_headers = {}
    if merged:
        ref_gaps_over_ten = [group for group in merged if group[1] - group[0] >= 10]
        for header, seq in aa_output:
            if header.endswith("."):
                continue
        
            start, end = find_index_pair(seq, "-")
            for region in ref_gaps_over_ten:
                this_overlap = get_overlap(start, end, region[0], region[1], 1)
                if this_overlap:
                    amt_non_gap = 0
                    kmer = seq[start:end]
                    for i, let in enumerate(kmer, start):
                        if let != "-":
                            amt_non_gap += 1
                            
                    if amt_non_gap / (end - start) >= 0.9:
                        internal_headers[header] = region[0], region[1]
                        break
    rescue_jank_log = [f"{gene} - {header}" for header in rescue_jank_log]
    aa_has_candidate = False
    for header, _ in aa_output:
        if not header.endswith("."):
            aa_has_candidate = True
            break

    if aa_has_candidate:
        gene_coverage = 1
        req_coverage = 0.4
        if gene_coverage < req_coverage:
            log_output.append(f">{gene}_kicked_coverage_{gene_coverage}_of_{req_coverage}\n")
            return log_output, False, False, gene, False, len(aa_nodes), this_rescues, scan_log, multi_log, ends, gff_out, debug_out, rescue_jank_log
        
        writeFasta(aa_out, aa_output, compress_intermediates)
        EXTEND_WINDOW = 45
        final_nt_out = []
        
        nt_output = [(node.header, node.nt_sequence) for node in aa_nodes if node.header not in kicked_headers]
            
        for header, seq in nt_output:
            if header in internal_headers:
                region_start, region_end = internal_headers[header]
                gt_positions = []
                ag_positions = []
                for i in range(max(0, (region_start * 3) - EXTEND_WINDOW), min(len(seq), (region_end * 3) + EXTEND_WINDOW), 2):
                    if seq[i:i+2] == "GT":
                        gt_positions.append(i)
                    
                    if seq[i:i+2] == "AG":
                        ag_positions.append(i)
                
                for gt_i, ag_i in product(gt_positions, ag_positions):
                    if gt_i > ag_i:
                        continue
                    
                    gt_ag_kmer = seq[gt_i:ag_i+2]
                    difference = len(gt_ag_kmer) - ((region_end - region_start) * 3)
                    debug_out.append(">"+header+f" - {region_start}:{region_end} - {difference}")
                    debug_out.append(gt_ag_kmer)
            if (len(seq) - seq.count("-")) % 3 != 0:
                print(header)
            final_nt_out.append((header, seq))
        writeFasta(nt_out, final_nt_out, compress_intermediates)

        return log_output, had_region, False, False, gene, len(kicked_headers), this_rescues, scan_log, combo_log, multi_log, ends, gff_out, debug_out, rescue_jank_log
    return log_output, had_region, gene, False, None, len(kicked_headers), this_rescues, scan_log, combo_log, multi_log, ends, gff_out, debug_out, rescue_jank_log

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


def get_args(args, genes, head_to_seq, input_folder, output_folder, compress, no_dupes, original_coords):
    get_id = lambda x: int(x.split("_")[0]) 
    for gene in genes:
        this_headers = []
        for header, _ in parseFasta(str(Path(input_folder, "aa", gene))):
            if header.endswith("."):
                continue
            this_headers.extend(
                map(get_id, header.split("|")[3].replace("NODE_", "").split("&&"))
            )

            this_seqs = {i: head_to_seq[i] for i in set(this_headers)}
            this_original_coords = {str(i): original_coords[str(i)] for i in set(this_headers)}
    
        yield (
            args.verbose,
            gene,
            input_folder,
            output_folder,
            compress,
            args.excise_consensus,
            args.excise_minimum_ambig,
            args.excise_rescue_match,
            no_dupes,
            this_seqs,
            this_original_coords,
        )
    
def get_head_to_seq(nt_db):
    """Get a dictionary of headers to sequences.

    Args:
    ----
        nt_db (RocksDB): The NT rocksdb database.
        recipe (list[str]): A list of each batch index.
    Returns:
    -------
        dict[str, str]: A dictionary of headers to sequences.
    """
    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
        
    head_to_seq = {}
    for i in recipe:
        lines = nt_db.get_bytes(f"ntbatch:{i}").decode().splitlines()
        head_to_seq.update(
            {
                int(lines[i][1:]): lines[i + 1]
                for i in range(0, len(lines), 2)
                if lines[i] != ""
            },
        )

    return head_to_seq


def main(args, sub_dir):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if args.excise_consensus > 1.0:
        if 0 < args.excise_consensus <= 100:
            args.excise_consensus = args.excise_consensus / 100
        else:
            raise ValueError(
                "Cannot convert excise_consensus to a percent. Use a decimal or a whole number between 0 and 100"
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

    rocksdb_path = path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(rocksdb_path)

    head_to_seq = {}
    original_coords = {}
    raw_data = nt_db.get("getall:original_coords")
    
    if raw_data:
        original_coords = json.decode(raw_data, type = dict[str, tuple[str, int, int, int, int]])
    head_to_seq = get_head_to_seq(nt_db)
    del nt_db

    genes = [fasta for fasta in listdir(aa_input) if ".fa" in fasta]
    log_path = Path(output_folder, "excise_regions.txt")
    scan_log_path = Path(output_folder, "gt_ag_scan.txt")
    combo_log_path = Path(output_folder, "combo_log.txt")
    multi_log_path = Path(output_folder, "multi_log.txt")
    internal_path = Path(output_folder, "internal.txt")
    gene_log_path = Path(output_folder, "excise_genes.txt")
    new_rescue_path = Path(output_folder, "new_rescues.txt")
    coords_path = Path(folder, "coords")
    arguments = get_args(args, genes, head_to_seq, input_folder, output_folder, args.compress, args.no_dupes, original_coords)
    if args.processes > 1:
        with Pool(args.processes) as pool:
            results = pool.starmap(log_excised_consensus, arguments)
    else:
        results = []
        for arg in arguments:
            results.append(
                log_excised_consensus(
                    *arg
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
    gt_ag_scan_log = []
    multi_scan_log = []
    combo_scan_log = []
    this_debug_out = []
    jank_log = []
    gt_hits = 0
    gc_hits = 0

    parent_gff_output = defaultdict(dict)
    end_bp = {}

    for glog, ghas_ambig, ghas_no_resolution, gcoverage_kick, g_has_resolution, gkicked_seq, grescues, slog, clog, dlog, input_lengths, gff_result, debug_lines, jlog in results:
        for parent, node_values in gff_result.items():
            for id, value in node_values.items():
                parent_gff_output[parent][id] = value
                
        end_bp.update(input_lengths)
        
        jank_log.extend(jlog)
        log_output.extend(glog)
        gt_ag_scan_log.extend(slog)
        this_debug_out.extend(debug_lines)
        multi_scan_log.extend(dlog)
        for line in clog:
            if line:
                if "GT-AG" in line:
                    gt_hits += 1
                if "GC-AG" in line:
                    gc_hits += 1
                combo_scan_log.append(line)
                
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

    reporter_coords_path = coords_path.joinpath("coords.gff")
    if reporter_coords_path.exists():
        with open(reporter_coords_path, "r") as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                line = line.strip().split("\t")
                #AGOUTI_SCAF_51|6429119BP|CTG001940_1,CTG001110_1,CTG004120_1	Sapphyre	exon	4815540	4815717	.	-	.	Name=136854;Parent=1.aa.fa;Note=-2;
                parent = line[0]
                start = int(line[3])
                
                fields = line[-1].split(";")
                for field in fields:
                    if field.startswith("Description="):
                        id = field.split("=")[1]
                        break
                
                if not parent in parent_gff_output:
                    continue
                
                if not id in parent_gff_output[parent]:
                    continue
                
                if parent_gff_output[parent][id] is None:
                    parent_gff_output[parent][id] = (start, "\t".join(line))
    else:
        printv("No reporter coords found. Unable to fill in the blank.", args.verbose, 0)
    
    gff_output= []
    for parent, rows in parent_gff_output.items():
            
        rows = [i for i in rows.values() if i]
        end = end_bp[parent]
        gff_output.append(f"##sequence-region\t{parent}\t{1}\t{end}")
        rows.sort(key = lambda x: (x[0]))
        gff_output.extend(i[1] for i in rows)
    
    if gff_output:
        with open(path.join(coords_path, "splice.gff"), "w") as fp:
            fp.write("\n".join(gff_output))


    printv(
        f"{input_folder}: {loci_containing_bad_regions} ambiguous loci found. {len(kicked_coverage)} genes kicked due to low coverage. Kicked {kicked_sequences} sequences total.",
        args.verbose,
    )

    if args.debug:
        combo_scan_log.insert(0, f"GT-AG Hits: {gt_hits}\nGC-AG Hits: {gc_hits}\n")
        with open(log_path, "w") as f:
            f.write("Gene,Cut-Indices,Consensus Sequence\n")
            f.write("\n".join(log_output))
        with open(scan_log_path, "w") as f:
            f.write("\n".join(gt_ag_scan_log))
        with open(combo_log_path, "w") as f:
            f.write("\n".join(combo_scan_log))
        with open(multi_log_path, "w") as f:
            f.write("\n".join(multi_scan_log))
        with open(internal_path, "w") as f:
            f.write("\n".join(this_debug_out))
        with open(new_rescue_path, "w") as f:
            f.write("\n".join(jank_log))
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
        return True
    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )