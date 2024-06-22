from collections import defaultdict
from functools import cached_property
from itertools import combinations, product
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


    return None, None

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


def do_trim(aa_nodes, cluster_sets, x_positions, ref_consensus, kicked_headers, prepare_dupes, reporter_dupes, excise_trim_consensus):
    for cluster_set in cluster_sets:
        
        sub_aa_nodes = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]

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
                    within_left = node.start <= i <= node.start + 3
                    within_right = node.end - 3 <= i <= node.end

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

### Clustering code

def determine_direction(start, end, current_start, current_end, current_direction, debug = False):
    this_direction = None
    if start == current_start and end == current_end:
        this_direction =  current_direction
    else:
        if start == current_start or current_direction == "reverse":
            if end >= current_end:
                this_direction = "forward"
            else:
                this_direction = "reverse"
        else:
            if start >= current_start:
                this_direction = "forward"
            else:
                this_direction = "reverse"

    if current_direction == "bi" or this_direction == "bi" or this_direction == current_direction:
        return this_direction
    return None


def node_to_ids(node):
    if "NODE_" in node:
        node = node.replace("NODE_","")
    return list(map(lambda x: int(x.split("_")[0]), node.split("&&")))

class cluster_rec:
    def __init__(self, node, start, end, seq_data_coords, strand) -> None:
        self.node = node
        self.start = start
        self.end = end
        self.seq_data_coords = seq_data_coords
        self.strand = strand
    
    @cached_property
    def get_ids(self): 
        return node_to_ids(self.node)
    
    def get_first_id(self):
        return min(self.get_ids)
    
   
def get_min_distance(a: set, b: set) -> int:
    """Get the minimum distance between two sets."""
    min_distance = float("inf")
    for a_obj in a:
        for b_obj in b:
            distance = abs(a_obj - b_obj)
            if distance < min_distance:
                min_distance = distance
    return min_distance
    
    
def within_distance(a: set, b: set, distance: int) -> bool:
    """Check if any object from set a is within distance of set b."""
    for a_obj in a:
        for b_obj in b:
            if abs(a_obj - b_obj) <= distance:
                return True
    return False


def finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage):

    if len(current_cluster) >= 2:
        cluster_data_cols = set()
        for rec in current_cluster:
            cluster_data_cols.update(rec.seq_data_coords)
            
        cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)
        clusters.append((min(current_indices), max(current_indices), cluster_coverage))
    elif len(current_cluster) == 1:
        cluster_rec = current_cluster[0]
        cluster_coverage = len(cluster_rec.seq_data_coords.intersection(ref_coords)) / len(ref_coords)
        
        if cluster_coverage > req_seq_coverage:
            clusters.append((min(current_indices), max(current_indices), cluster_coverage))
        else:
            kicks.add(cluster_rec.node)


def cluster_ids(ids, max_id_distance, max_gap, ref_coords):
    clusters = []
    kicks = set()
    ids.sort(key=lambda x: x.get_first_id())
    req_seq_coverage = 0.5
    current_cluster = None
    debug = False
    
    for x, rec in enumerate(ids):
        
        if current_cluster is None:
            current_cluster = [rec]
            current_indices = set(rec.get_ids)
            current_direction = "bi"
        else:   
            passed = False
            cluster_distance = get_min_distance(current_indices, rec.get_ids)
            passed_id_distance = None
            passed_distance = None
            if rec.strand == current_cluster[-1].strand:
                if cluster_distance <= max_id_distance:
                    current_rec = current_cluster[-1]
                    cluster_overlap = get_overlap(rec.start, rec.end, current_rec.start, current_rec.end, -max_gap)
                    
                    if cluster_overlap:
                        cluster_pos_distance = min(
                            abs(rec.start - current_rec.start), 
                            abs(rec.start - current_rec.end), 
                            abs(rec.end - current_rec.start), 
                            abs(rec.end - current_rec.end)
                        )

                        this_direction = determine_direction(rec.start, rec.end, current_rec.start, current_rec.end, current_direction)
                        if this_direction:
                            passed_id_distance = cluster_distance
                            passed_distance = cluster_pos_distance
                            passed_direction = this_direction
                            passed = True
                
            if passed:
                # not last id
                if x != len(ids) - 1:
                    next_rec = ids[x + 1]
                    if within_distance(rec.get_ids, next_rec.get_ids, max_id_distance):
                        next_overlap = get_overlap(rec.start, rec.end, next_rec.start, next_rec.end, -max_gap)
                        if next_overlap:
                            next_distance = get_min_distance(rec.get_ids, next_rec.get_ids)
                                
                            next_amount = min(
                                abs(rec.start - next_rec.start), 
                                abs(rec.start - next_rec.end), 
                                abs(rec.end - next_rec.start), 
                                abs(rec.end - next_rec.end)
                            )
                            
                            next_direction = determine_direction(next_rec.start, next_rec.end, rec.start, rec.end, current_direction)


                            if next_direction and passed_direction and next_direction != passed_direction and next_distance < passed_id_distance and next_amount < passed_distance:
                                passed = False
                    
            if passed:
                current_cluster.append(rec)
                current_indices.update(rec.get_ids)
                if passed_direction != "bi":
                    current_direction = passed_direction
            else:
                if cluster_distance == 0:
                    new_cluster = [rec]
                    new_indices = set(rec.get_ids)
                    kick_occured = True
                    while kick_occured:
                        kick_occured = False
                        for rec_in in current_cluster:
                            if get_min_distance(new_indices, rec_in.get_ids) == 0:
                                current_cluster.remove(rec_in)
                                current_indices.difference_update(rec_in.get_ids)
                                new_cluster.append(rec_in)
                                new_indices.update(rec_in.get_ids)
                                kick_occured = True
                            
                    if current_cluster:
                        finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage)
                    
                    current_cluster = new_cluster
                    current_indices = new_indices
                    current_direction = determine_direction(rec.start, rec.end, current_rec.start, current_rec.end, "bi")
                else:
                    finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage)
                    current_cluster = [rec]
                    current_indices = set(rec.get_ids)
                    current_direction = "bi"
    
    if current_cluster:
        finalize_cluster(current_cluster, current_indices, ref_coords, clusters, kicks, req_seq_coverage)
                
    return clusters, kicks

###


def insert_gaps(input_string, positions, offset):
    input_string = list(input_string)
    
    for coord in positions:
        input_string.insert(offset + coord, "-")

    return "".join(input_string)


def get_combo_results(gt_positions, ag_positions, prev_node, node, FRANKENSTEIN_PENALTY, INSERTION_PENALTY, DNA_CODONS, ref_gaps):
    this_results = []
                
    for (gt_size, act_gt_index, gt_index, this_prev_extensions), (ag_size, act_ag_index_rev, ag_index_rev, this_node_extensions) in product(gt_positions, ag_positions):
        prev_deletions = []
        node_deletions = []

        prev_nt_seq = list(prev_node.nt_sequence)
        node_seq = list(node.nt_sequence)
        
        for i, let in this_prev_extensions.items():
            prev_nt_seq[i] = let
        for i, let in this_node_extensions.items():
            node_seq[i] = let
        
        for i in range(act_gt_index, act_gt_index + gt_size):
            prev_nt_seq[i] = "-"
            # prev_nt_seq[i] = prev_nt_seq[i].lower()
            
        for i in range(act_ag_index_rev, act_ag_index_rev + ag_size):
            node_seq[i] = "-"
            # node_seq[i] = node_seq[i].lower()

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
                node_seq[x] = "-"
                node_deletions.append(x)
                
        
        prev_nt_seq = [x for x in prev_nt_seq if x != ""]
        node_seq = [x for x in node_seq if x != ""]
                
        node_nt_start, node_nt_end = find_index_pair("".join(node_seq), "-")
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
            print(node.header)
            for i in range(act_gt_index + 1, act_gt_index + gt_size - 1):
                if node_seq[i] != "-":
                    node_gap_insertions.append(i)
                    node_seq.insert(i, "-")
                    node_seq.pop(0)
                
        node_nt_start, node_nt_end = find_index_pair("".join(node_seq), "-")
        length = node_nt_end - node_nt_start
        right_end_codon = node_nt_start - (3 - (length % 3))
        
        prev_nt_start, prev_nt_end = find_index_pair("".join(prev_nt_seq), "-")
        length = prev_nt_end - prev_nt_start
        left_last_codon = prev_nt_start + length - (length % 3)
        
        this_score = 0
        right_codon = node_seq[right_end_codon: right_end_codon + 3]
        left_codon = prev_nt_seq[left_last_codon: left_last_codon + 3]

        if right_end_codon == left_last_codon:
            orphan_codon = []
            for i in range(3):
                if left_codon[i] == "-":
                    orphan_codon.append(right_codon[i])
                else:
                    orphan_codon.append(left_codon[i])

            joined = "".join(orphan_codon)
            if joined not in DNA_CODONS:
                this_score += FRANKENSTEIN_PENALTY
            elif DNA_CODONS[joined] == "*":
                continue
        
        else:
            if "-" in right_codon and right_codon.count("-") != 3:
                continue
            if "-" in left_codon and left_codon.count("-") != 3:
                continue 
                
            orphan_codon = None
            
        distance = node_nt_start - prev_nt_end
        if not (0 <= distance <= 2):
            if prev_nt_end > node_nt_start:
                this_range = range(node_nt_start, prev_nt_end)
            else:
                this_range = range(prev_nt_end, node_nt_start)
            for i in this_range:
                if node_seq[i] == "-" and prev_nt_seq[i] == "-" and i // 3 in ref_gaps:
                    continue
                this_score -= 1

        this_results.append((this_score, gt_size, act_gt_index, gt_index, ag_size, act_ag_index_rev, ag_index_rev, prev_deletions, node_deletions, prev_gap_insertions, node_gap_insertions, this_prev_extensions, this_node_extensions, left_last_codon, right_end_codon, orphan_codon))

    return this_results

def find_gt_ag(prev_node, node, prev_start_index, prev_end_index, node_start_index, node_end_index, DNA_CODONS, prev_og, node_og, prev_nt_seq, node_seq, ref_gaps):
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
        
    x = 0
    act_offset = (prev_node.end * 3) - overlap_amount - EXTEND_WINDOW
    prev_start_scan = prev_end_index - overlap_amount - EXTEND_WINDOW
    if prev_start_scan < prev_start_index:
        act_offset = (prev_node.start * 3) + 3
        prev_start_scan = prev_start_index + 3
    
    for i in range(prev_start_scan, len(prev_og)):
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

        scan_index = 0
        if prev_og[i] == "G":
            while True:
                scan_index += 1
                if i + scan_index >= len(prev_og):
                    break
                
                if prev_og[i + scan_index] == "T":
                    act_gt_index = prev_act_coord
                    gt_index = i + 1
                    gt_size = scan_index + 1
                    gt_positions.append((gt_size, act_gt_index, gt_index, copy.deepcopy(prev_extensions)))

                if prev_og[i + scan_index] != "-":
                    break

                if prev_act_coord + scan_index not in ref_gaps:
                    break            

        if prev_nt_seq[prev_act_coord] == "-":
            prev_extensions[prev_act_coord] = prev_og[i]
        
        x += 1
        
    # Iterate in reverse from the start of the kmer to the start of the original sequence
    x = 0
    
    start_scan = node_start_index + overlap_amount + EXTEND_WINDOW
    act_offset = (node.start * 3) + overlap_amount + EXTEND_WINDOW
    if start_scan > node_end_index:
        start_scan = node_end_index
        act_offset = (node.end * 3)
        
    for i in range(start_scan - 1, -1, -1):
        node_act_coord = act_offset - 1 - x
        # Get last codon
        
        if node_act_coord >= len(node_seq) - 1 or i >= len(node_og) - 1:
            continue
        
        if node_act_coord < 0 or i < 0:
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

        if node_seq[node_act_coord - 1] == "-":
            node_extensions[node_act_coord - 1] = node_og[i - 1]
            
        x += 1
            
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
                 replacements,
                 replacements_aa,
                 extensions,
                 extensions_aa,
                 gap_insertions_aa,
                 gap_insertions_nt,
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
    node_seq = list(node.nt_sequence)
    
    this_score, gt_size, act_gt_index, gt_index, ag_size, act_ag_index_rev, ag_index_rev, prev_deletions, node_deletions, final_prev_insertions, final_node_insertions, final_prev_extensions, final_node_extensions, left_last_codon, right_end_codon, orphan_codon = this_result
    for i, let in final_prev_extensions.items():
        prev_nt_seq[i] = let
    for i, let in final_node_extensions.items():
        node_seq[i] = let
    for i in range(act_gt_index, act_gt_index + gt_size):
        if add_results:
            replacements[prev_node.header][i] = "-"
            replacements_aa[prev_node.header][i//3] = "-"
        prev_nt_seq[i] = "-"
    for i in range(act_ag_index_rev, act_ag_index_rev + ag_size):
        if add_results:
            replacements[node.header][i] = "-"
            replacements_aa[node.header][i//3] = "-"
        node_seq[i] = "-"
    for x in prev_deletions:
        if x in final_prev_extensions:
            continue
        if add_results:
            replacements_aa[prev_node.header][x//3] = "-"#
            replacements[prev_node.header][x] = "-"
        prev_nt_seq[x] = "-"
    for x in node_deletions:
        if x in final_node_extensions:
            continue
        if add_results:
            replacements_aa[node.header][x//3] = "-"
            replacements[node.header][x] = "-"
        node_seq[x] = "-"
        
    for i in final_prev_insertions:
        if add_results:
            # if i % 3 == 0:
            #     gap_insertions_aa[prev_node.header].append((i//3, -1))
            gap_insertions_nt[prev_node.header].append((i, -1))
        prev_nt_seq.insert(i, "-")
        prev_nt_seq.pop(-1)
    
    for i in final_node_insertions:
        if add_results:
            # if i % 3 == 0:
            #     gap_insertions_aa[node.header].append((i//3, 0))
            gap_insertions_nt[node.header].append((i, 0))
        node_seq.insert(i, "-")
        node_seq.pop(0)

    if orphan_codon:
        joined = "".join(orphan_codon)
        if joined not in DNA_CODONS:
            return None, None
        
        orphan_aa = DNA_CODONS[joined]
        if add_results:
            replacements_aa[prev_node.header][left_last_codon//3] = orphan_aa
            replacements_aa[node.header][right_end_codon//3] = orphan_aa

        prev_og = list(prev_og)
        node_og = list(node_og)

        for i, x in enumerate(range(left_last_codon, left_last_codon + 3)):
            prev_nt_seq[x] = orphan_codon[i]
            if add_results:
                prev_og[prev_start_index + x - (prev_node.start * 3)] = orphan_codon[i]
                replacements[prev_node.header][x] = orphan_codon[i]

        for i, x in enumerate(range(right_end_codon, right_end_codon + 3)):
            node_seq[x] = orphan_codon[i]
            if add_results:
                node_og[node_start_index + x - (node.start * 3)] = orphan_codon[i]
                replacements[node.header][x] = orphan_codon[i]
                
        prev_og = "".join(prev_og)
        node_og = "".join(node_og)
        
        formed_seqs[prev_node.header] = prev_og.replace("-", "")
        formed_seqs[node.header] = node_og.replace("-", "")
        
    node_seq = "".join(node_seq)
    prev_nt_seq = "".join(prev_nt_seq)
        
    node_start, node_end = find_index_pair(node_seq, "-")
    prev_start, prev_end = find_index_pair(prev_nt_seq, "-")
    
    
    if node_start % 3 != 0:
        node_start -= node_start % 3
    
    if prev_end % 3 != 0:
        prev_end += 3 - (prev_end % 3)
    
    if add_results:
        for i in range((prev_node.end * 3) - 3, prev_end, 3):
            codon = prev_nt_seq[i:i+3]
            if codon in DNA_CODONS:
                extensions_aa[prev_node.header][i//3] = DNA_CODONS[codon]
            
        for i in range(node_start, node.start * 3, 3):
            codon = node_seq[i:i+3]
            if codon in DNA_CODONS:
                extensions_aa[node.header][i//3] = DNA_CODONS[codon]

        extensions[node.header].update(final_node_extensions)
        extensions[prev_node.header].update(final_prev_extensions)
        
    # Extend one codon to the left and right for debug to show end of sequence
    if prev_start >= 3:
        prev_start -= 3
    if node_end + 3 <= len(node_seq):
        node_end += 3
        
    node_region = node.nt_sequence[prev_start: node_end]
    node_region_start, _ = find_index_pair(node_region, "-")
    
    final_prev_start, final_prev_end = find_index_pair(prev_nt_seq, "-")
    final_node_start, final_node_end = find_index_pair(node_seq, "-")
    
    smallest_change = min(abs(final_prev_start - (prev_node.start*3)) + abs(final_prev_end - (prev_node.end*3)), abs(final_node_start - (node.start*3)) + abs(final_node_end - (node.end*3)))
    
    if print_extra: 
        scan_log.append("")
        scan_log.append("")   
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
    scan_log.append(node_seq[prev_start: node_end])
    scan_log.append("")    
    
    node_hit = node_og[ag_index_rev: (node_start_index + len(kmer) + len(kmer_internal_gaps))]
    prev_hit = prev_og[prev_start_index: gt_index + 1]
    # lowercase
    prev_hit = prev_hit[:-gt_size] + prev_hit[-gt_size:].lower()
    node_hit = node_hit[:ag_size].lower() + node_hit[ag_size:]

    scan_log.append(f">{prev_node.header}_orf_scan")
    scan_log.append((("-" * (prev_node.start * 3)) + prev_hit)[prev_start: node_end])
    scan_log.append(f">{node.header}_orf_scan")
    scan_log.append((("-" * ((node.end * 3) - len(node_hit))) + node_hit)[prev_start: node_end] )
    scan_log.append("") 
    
    if add_results:
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
        node.nt_sequence = node_seq
        
        prev_node.start, prev_node.end = final_prev_start//3, final_prev_end//3
        node.start, node.end = final_node_start//3, final_node_end//3
        
        return (gff_coord_prev, gff_coord_node), smallest_change
    return None, smallest_change

def log_excised_consensus(
    verbose: int,
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
    head_to_seq,
    original_coords,
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
    printv(f"Processing {gene}", verbose, 2)
    log_output = []
    scan_log = []
    multi_log = []
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
    ref_gaps = set()
    for header, seq in raw_aa:
        if header.endswith('.'):
            start, end = find_index_pair(seq, "-")
            for i, bp in enumerate(seq[start:end], start):
                if bp != "-":
                    reference_cluster_data.add(i)
                ref_consensus[i].append(bp)
                
            ref_lens.append(bp_count(seq))
            continue

        frame = int(header.split("|")[4])
        aa_nodes.append(NODE(header, frame, seq, None, *find_index_pair(seq, "-"), []))

    ref_gap_percent = 0.75
    for i, lets in ref_consensus.items():
        if lets.count("-") / len(lets) > ref_gap_percent:
            ref_gaps.add(i)

    ref_avg_len = sum(ref_lens) / len(ref_lens)
    kicked_headers = set()
    replacements = defaultdict(dict)
    replacements_aa = defaultdict(dict)

    cluster_sets = [None]
    get_id = lambda header: header.split("|")[3].replace("NODE_","")
    get_parent_id = lambda header: int(header.split("|")[3].split("&&")[0].split("_")[1])
    if is_genome:
        ids = []
        for node in aa_nodes:
            if node.header not in kicked_headers:
                start, end = find_index_pair(node.sequence, "-")
                data_cols = {i for i, let in enumerate(node.sequence[start:end], start) if let != "-"}
                ids.append(cluster_rec(node.header.split("|")[3], start, end, data_cols, "-" if node.frame < 0 else "+"))
    
        max_gap_size = round(len(aa_nodes[0].sequence) * 0.5) # Half MSA length
    
        clusters, kicks = cluster_ids(ids, 100, max_gap_size, reference_cluster_data) #TODO: Make distance an arg
        if kicks:
            for node in aa_nodes:
                this_parent_id = int(get_id(node.header).split("&&")[0].split("_")[0])
                if this_parent_id in kicks:
                    kicked_headers.add(node.header)
                    log_output.append(f"Kicking {node.header} due to low coverage")
        
        if clusters:
            cluster_sets = [set(range(a, b+1)) for a, b, _ in clusters]
            

    if not is_genome:
        do_trim(aa_nodes, cluster_sets, x_positions, ref_consensus, kicked_headers, prepare_dupes, reporter_dupes, excise_trim_consensus)

    aa_sequence = {}
    for node in aa_nodes:

        node.sequence = del_cols(node.sequence, x_positions[node.header])
        node.start, node.end = find_index_pair(node.sequence, "-")

        aa_sequence[node.header] = node.sequence
        node_kmer = node.sequence[node.start:node.end]
        data_len = bp_count(node_kmer)
        if data_len < 15:
            log_output.append(f"Kicking {node.header} due to < 15 bp after trimming")
            kicked_headers.add(node.header)

    raw_sequences = {header: del_cols(seq, x_positions[header], True) for header, seq in parseFasta(str(nt_in))}
    for node in aa_nodes:
        if node.header in kicked_headers:
            continue
        nt_seq = raw_sequences[node.header]
        node.nt_sequence = nt_seq
    
    true_cluster_raw = []

    for header in raw_sequences:
        true_cluster_raw.append((int(header.split("|")[3].split("&&")[0].split("_")[1]), header))

    true_cluster_raw.sort(key = lambda x: x[0])
    before_true_clusters = cluster(true_cluster_raw, true_cluster_threshold)

    # with open(gene+"_debug.txt", "w") as fp:
    #     for cluster_i, cluster_set in enumerate(cluster_sets):
    #         if 1482133 in cluster_set:
    #             for node in aa_nodes:
    #                 if "NODE_1482133&&1482134" in node.header:
                        
    #                     print(node.header in kicked_headers)
    #                     print(cluster_set)
    #                     print(node_to_ids(node.header.split("|")[3]))
    #                     print(within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))
                        
    #         fp.write(f"Cluster {cluster_i}\n")
    #         aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
    #         for node in aa_subset:
    #             if "NODE_1482133&&1482134" in node.header:
    #                 print("wot")
    #             fp.write(f">{node.header}\n{node.sequence}\n")
        

    # Search for slices of the consensus seq with a high ratio of 'X' to total characters
    has_region = True
    had_region = False
    recursion_max = 5 * len(cluster_sets)
    last_region = {}
    consensus_seq = ""
    while has_region:
        has_region = False

        for cluster_i, cluster_set in enumerate(cluster_sets):

            aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]

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
            region_start, region_end = check_covered_bad_regions(consensus_seq, excise_minimum_ambig)
            if region_start == last_region.get(cluster_i, -1):
                recursion_max -= 1
                if recursion_max <= 0:
                    region_start = None
                    continue
            last_region[cluster_i] = region_start

            if region_start:
                has_region = True
                had_region = True
                sequences_in_region = []
                sequences_out_of_region = []      

                for i, node in enumerate(aa_subset):

                    if node.header in kicked_headers:
                        continue

                    overlap_coords = get_overlap(region_start, region_end, node.start * 3, node.end * 3, 1)
                    if overlap_coords:
                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        overlap_percent = overlap_amount / (node.end - node.start)
                        if overlap_percent >= excise_region_overlap: # Adjustable percent
                            sequences_in_region.append(aa_subset[i])
                        else:
                            sequences_out_of_region.append(aa_subset[i])

                if len(sequences_in_region) > excise_maximum_depth:
                    continue
                
                nodes_in_region = None
                if is_genome:
                    continue
                    tagged_in_region = [(int(node.header.split("|")[3].split("&&")[0].split("_")[1]), node) for node in sequences_in_region]
                    tagged_in_region.sort(key=lambda x: x[0])
                    clusters = cluster(tagged_in_region, true_cluster_threshold)

                    log_output.append(f">{gene}_ambig_{region_start}:{region_end}\n{consensus_seq}")

                    for clust in clusters:
                        if len(clust) <= 1:
                            continue

                        clust.sort(key = lambda node: node.start)
                        for (_, prev_node), (i, node) in combinations(enumerate(clust), 2):
                            overlapping_coords = get_overlap(node.start, node.end, prev_node.start, prev_node.end, 1)
                            if overlapping_coords:
                                kmer = node.sequence[overlapping_coords[0]:overlapping_coords[1]]
                                prev_kmer = prev_node.sequence[overlapping_coords[0]:overlapping_coords[1]]

                                if is_same_kmer(kmer, prev_kmer):
                                    continue

                                splice_index = calculate_split(prev_node, node, overlapping_coords, ref_consensus)
                                prev_positions = set()
                                node_positions = set()
                                
                                consecutive_match = True
                                for x in range(splice_index - 1, node.start - 1, -1):
                                    if consecutive_match and node.sequence[x] == prev_node.sequence[x]:
                                        continue
                                    
                                    consecutive_match = False
                                    node_positions.add(x)

                                for x in range(splice_index, prev_node.end):
                                    prev_positions.add(x)

                                log_output.append(f">{prev_node.header} vs {node.header}")
                                log_output.append(f"Split at {splice_index}")

                                prev_node.sequence = del_cols(prev_node.sequence, prev_positions)
                                node.sequence = del_cols(node.sequence, node_positions)

                                if len(prev_node.sequence) - prev_node.sequence.count("-") < 15:
                                    log_output.append(f"Kicking {prev_node.header} due to < 15 bp after splice")
                                    kicked_headers.add(prev_node.header)
                                    

                                if len(node.sequence) - node.sequence.count("-") < 15:
                                    log_output.append(f"Kicking {node.header} due to < 15 bp after splice")
                                    kicked_headers.add(node.header)

                                # if either_kicked:
                                #     break

                                prev_node.nt_sequence = del_cols(prev_node.nt_sequence, prev_positions, True)
                                node.nt_sequence = del_cols(node.nt_sequence, node_positions, True)

                                node.start, node.end = find_index_pair(node.sequence, "-")
                                prev_node.start, prev_node.end = find_index_pair(prev_node.sequence, "-")

                                x_positions[node.header].update(node_positions)
                                x_positions[prev_node.header].update(prev_positions)
                
                else:
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
                    log_output.append(f">{gene}_ambig_{region_start}:{region_end}\n{consensus_seq}")
                    if nodes_in_region:
                        log_output.extend([f">{node.contig_header()}_{'kept' if i in keep_indices else 'kicked'}\n{node.nt_sequence}" for i, node in enumerate(nodes_in_region)])
                    else:
                        log_output.extend([f">{node.header}_{'kept' if node.header not in kicked_headers else 'kicked'}\n{node.nt_sequence}" for node in sequences_in_region])
                    log_output.append("\n")
                
        if recursion_max <= 0:
            break
        
    # if is_genome:
    #     do_trim(aa_nodes, cluster_sets, x_positions, ref_consensus, kicked_headers, prepare_dupes, reporter_dupes, excise_trim_consensus)
    #     for node in aa_nodes:
    #         if node.header in kicked_headers:
    #             continue
            
    #         node.sequence = del_cols(node.sequence, x_positions[node.header])
    #         node.nt_sequence = del_cols(node.nt_sequence, x_positions[node.header], True)
    #         node.start, node.end = find_index_pair(node.sequence, "-")
                          
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
        
        if is_assembly_or_genome:
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
    
    FRANKENSTEIN_PENALTY = -20
    INSERTION_PENALTY = -1
    SIMILARITY_SKIP = 0.95
    int_first_id = lambda x: int(x.split("_")[0])
    extensions = defaultdict(dict)
    extensions_aa = defaultdict(dict)
    gap_insertions_aa = defaultdict(list)
    gap_insertions_nt = defaultdict(list)
    ends = {}
    gff_out = defaultdict(dict)
    gff_coords = {}
    formed_seqs = {}
    for cluster_i, cluster_set in enumerate(cluster_sets):
        aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
        aa_subset.sort(key = lambda x: x.start)
        og_starts = {}
        for prev_node, node in combinations(aa_subset, 2):
            overlapping_coords = get_overlap(node.start, node.end, prev_node.start, prev_node.end, -10)
            if overlapping_coords:
                amount = overlapping_coords[1] - overlapping_coords[0]
                if amount > 1:
                    early_start = min(prev_node.start, node.start) * 3
                    late_end = max(prev_node.end, node.end) * 3
                    
                    distance = constrained_distance(prev_node.nt_sequence[early_start : late_end], node.nt_sequence[early_start : late_end])
                    percent_matching = 1 - (distance / (late_end - early_start))
                    if percent_matching >= SIMILARITY_SKIP:
                        continue
                
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
                    # print(prev_node.header)
                    # print("FAIL")
                    continue

                prev_og = insert_gaps(prev_og, prev_internal_gaps, prev_start_index)
                prev_end_index = prev_start_index + len(prev_kmer) + len(prev_internal_gaps)

                node_start_index = node_og.find(kmer)
                if node_start_index == -1:
                    continue

                
                node_og = insert_gaps(node_og, kmer_internal_gaps, node_start_index)
                node_end_index = node_start_index + len(kmer) + len(kmer_internal_gaps)

                prev_nt_seq = list(prev_node.nt_sequence)
                node_seq = list(node.nt_sequence)
                
                gt_positions, ag_positions = find_gt_ag(prev_node, node, prev_start_index, prev_end_index, node_start_index, node_end_index, DNA_CODONS, prev_og, node_og, prev_nt_seq, node_seq, ref_gaps)

                # if "2A|AglaOr12CTE|splice_fix|NODE_343534&&343535|-3|1" not in prev_node.header:
                #     continue

                this_results = get_combo_results(gt_positions, ag_positions, prev_node, node, FRANKENSTEIN_PENALTY, INSERTION_PENALTY, DNA_CODONS, ref_gaps)
                splice_found = False
                if this_results:
                    this_best_score = max(i[0] for i in this_results)
                    highest_results = [i for i in this_results if i[0] == this_best_score]
                    if len(highest_results) > 1:
                        best_change = None
                        for i, result in enumerate(highest_results):
                            _, smallest_coords_change = splice_combo(False, i == 0, formed_seqs, result, prev_node, node, prev_og, node_og, DNA_CODONS, multi_log, replacements, replacements_aa, extensions, extensions_aa, gap_insertions_aa, gap_insertions_nt, prev_start_index, node_start_index, kmer, kmer_internal_gaps, prev_internal_gaps, gff_coords)
                            if best_change is None or smallest_coords_change < best_change:
                                best_change = smallest_coords_change
                                this_best_splice = result
                    else:
                        this_best_splice = highest_results[0]
                    
                    gff, _ = splice_combo(True, True, formed_seqs, this_best_splice, prev_node, node, prev_og, node_og, DNA_CODONS, scan_log, replacements, replacements_aa, extensions, extensions_aa, gap_insertions_aa, gap_insertions_nt, prev_start_index, node_start_index, kmer, kmer_internal_gaps, prev_internal_gaps, gff_coords)
                    if gff:
                        splice_found = True
                        prev_gff, node_gff = gff
                        
                        prev_id = get_id(prev_node.header)
                        tup = original_coords.get(prev_id.split("&&")[0].split("_")[0], None)
                        if tup:
                            parent, chomp_start, chomp_end, input_len, chomp_len = tup
                            
                            prev_start = prev_gff[0] + chomp_start
                            prev_end = prev_gff[1] + chomp_start
                            
                            if parent not in ends:
                                ends[parent] = input_len
                                
                            strand = "+" if prev_node.frame > 0 else "-"
                            gff_out[parent][prev_id] = ((prev_start), f"{parent}\tSapphyre\texon\t{prev_start}\t{prev_end}\t.\t{strand}\t.\tID={prev_id};Parent={gene};Note={prev_node.frame};")
                            
                        node_id = get_id(node.header)
                        tup = original_coords.get(node_id.split("&&")[0].split("_")[0], None)
                        if tup:
                            parent, chomp_start, chomp_end, input_len, chomp_len = tup
                            
                            node_start = node_gff[0] + chomp_start
                            node_end = node_gff[1] + chomp_start
                            
                            if parent not in ends:
                                ends[parent] = input_len
                                
                            strand = "+" if node.frame > 0 else "-"
                            gff_out[parent][node_id] = ((node_start), f"{parent}\tSapphyre\texon\t{node_start}\t{node_end}\t.\t{strand}\t.\tID={node_id};Parent={gene};Note={node.frame};")
                if True and not splice_found:
                    scan_log.append("")
                    scan_log.append("")    
                    scan_log.append(f">{prev_node.header}_orf")
                    # print(prev_start_index, node_end_index)
                    # input()
                    node_start, node_end = find_index_pair("".join(node_seq), "-")
                    prev_start, prev_end = find_index_pair("".join(prev_nt_seq), "-")
                    prev_start -= 3
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
                    scan_log.append(prev_node.nt_sequence[prev_start: node_end])
                    scan_log.append(f">{node.header}_excise_output")
                    scan_log.append(node_region)
                    scan_log.append("")
                    scan_log.append("No result found.")
                    scan_log.append("")
                    
                    
    aa_raw_output = [(header, del_cols(seq, x_positions[header])) for header, seq in raw_aa if header not in kicked_headers]
    aa_output = []
    for header, seq in aa_raw_output:
        this_id = get_id(header)
        tup = original_coords.get(this_id.split("&&")[0].split("_")[0], None)
        if tup:
            parent, _, _, input_len, _ = tup
            if parent not in ends:
                ends[parent] = input_len
                
            if this_id not in gff_out[parent]:
                gff_out[parent][this_id] = None
        
        if header in replacements_aa or header in extensions_aa or header in gap_insertions_aa:
            seq = list(seq)
            if header in extensions_aa:
                for i, bp in extensions_aa[header].items():
                    seq[i] = bp
            # if header in gap_insertions_aa:
            #     for insert_i, pop_i in gap_insertions_aa[header]:
            #         seq.insert(insert_i, "-")
            #         seq.pop(pop_i)
            if header in replacements_aa:
                for i, bp in replacements_aa[header].items():
                    seq[i] = bp
                    
            seq = "".join(seq)
        aa_output.append((header, seq))
        
    integers = list(ref_gaps)
    debug_out = []
    merged = []
    if integers:
        integers.sort()
        start = integers[0]
        end = integers[0]

        for num in integers[1:]:
            if num == end + 1:
                end = num
            else:
                if end - start >= 10:
                    merged.append((start, end))
                start = num
                end = num

        if end - start >= 10:
            merged.append((start, end))
    
    internal_headers = {}
    if merged:
        for header, seq in aa_output:
            if header.endswith("."):
                continue
        
            start, end = find_index_pair(seq, "-")
            for region in merged:
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

    aa_has_candidate = False
    for header, _ in aa_output:
        if not header.endswith("."):
            aa_has_candidate = True
            break

    if aa_has_candidate:
        gene_coverage = 1#get_coverage([seq for header, seq in aa_output if header[-1] != "."], ref_avg_len)
        req_coverage = 0.4 if is_assembly_or_genome else 0.01
        if gene_coverage < req_coverage:
            log_output.append(f">{gene}_kicked_coverage_{gene_coverage}_of_{req_coverage}\n{consensus_seq}")
            return log_output, False, False, gene, False, len(aa_nodes), this_rescues, scan_log, multi_log, ends, gff_out, debug_out
        
        writeFasta(aa_out, aa_output, compress_intermediates)
        nt_output = [(header, del_cols(seq, x_positions[header], True)) for header, seq in raw_sequences.items() if header not in kicked_headers]
        EXTEND_WINDOW = 45
        final_nt_out = []
            
        for header, seq in nt_output:
            if header in replacements or header in extensions or header in gap_insertions_nt:
                seq = list(seq)
                if header in extensions: 
                    for i, bp in extensions[header].items():
                        seq[i] = bp
                if header in gap_insertions_nt:
                    gap_inserts = sorted(list(gap_insertions_nt[header]))
                    for insert_i, pop_i in gap_inserts:
                        seq.insert(insert_i, "-")
                        seq.pop(pop_i)
                if header in replacements:
                    for i, bp in replacements[header].items():
                        seq[i] = bp
                seq = "".join(seq)
                
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
            
            final_nt_out.append((header, seq))
        writeFasta(nt_out, final_nt_out, compress_intermediates)

        return log_output, had_region, False, False, gene, len(kicked_headers), this_rescues, scan_log, multi_log, ends, gff_out, debug_out
    return log_output, had_region, gene, False, None, len(kicked_headers), this_rescues, scan_log, multi_log, ends, gff_out, debug_out


def load_dupes(rocksdb_db):
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


def get_args(args, genes, head_to_seq, is_assembly_or_genome, is_genome, input_folder, output_folder, compress, prepare_dupes, reporter_dupes, original_coords):

    get_id = lambda x: int(x.split("_")[0]) 
    for gene in genes:
        this_prepare_dupes = prepare_dupes.get(gene.split(".")[0], {})
        this_reporter_dupes = reporter_dupes.get(gene.split(".")[0], {})
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
            this_prepare_dupes,
            this_reporter_dupes,
            this_seqs,
            this_original_coords,
            args.excise_trim_consensus,
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


def main(args, sub_dir, is_genome, is_assembly_or_genome):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if args.excise_overlap_merge > 1.0:
        if 0 < args.excise_overlap_merge <= 100:
            args.excise_overlap_merge = args.excise_overlap_merge / 100
        else:
            raise ValueError(
                "Cannot convert excise_overlap_merge to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if args.excise_overlap_ambig > 1.0:
        if 0 < args.excise_overlap_ambig <= 100:
            args.excise_overlap_ambig = args.excise_overlap_ambig / 100
        else:
            raise ValueError(
                "Cannot convert excise_overlap_ambig to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if args.excise_consensus > 1.0:
        if 0 < args.excise_consensus <= 100:
            args.excise_consensus = args.excise_consensus / 100
        else:
            raise ValueError(
                "Cannot convert excise_consensus to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if args.excise_region_overlap > 1.0:
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

    rocksdb_path = path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(rocksdb_path)
    
    head_to_seq = get_head_to_seq(nt_db)
    
    original_coords = {}
    if is_genome:
        raw_data = nt_db.get("getall:original_coords")
        
        if raw_data:
            original_coords = json.decode(raw_data, type = dict[str, tuple[str, int, int, int, int]])
    
    if args.no_dupes:
        prepare_dupes, reporter_dupes = {}, {}
    else:
        if not path.exists(rocksdb_path):
            err = f"cannot find rocksdb for {folder}"
            raise FileNotFoundError(err)
        prepare_dupes, reporter_dupes = load_dupes(nt_db)
        
    del nt_db

    compress = not args.uncompress_intermediates or args.compress

    genes = [fasta for fasta in listdir(aa_input) if ".fa" in fasta]
    log_path = Path(output_folder, "excise_regions.txt")
    scan_log_path = Path(output_folder, "gt_ag_scan.txt")
    multi_log_path = Path(output_folder, "multi_log.txt")
    internal_path = Path(output_folder, "internal.txt")
    gene_log_path = Path(output_folder, "excise_genes.txt")
    coords_path = Path(folder, "coords")
    arguments = get_args(args, genes, head_to_seq, is_assembly_or_genome, is_genome, input_folder, output_folder, compress, prepare_dupes, reporter_dupes, original_coords)
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
    this_debug_out = []

    parent_gff_output = defaultdict(dict)
    end_bp = {}

    for glog, ghas_ambig, ghas_no_resolution, gcoverage_kick, g_has_resolution, gkicked_seq, grescues, slog, dlog, input_lengths, gff_result, debug_lines in results:
        for parent, node_values in gff_result.items():
            for id, value in node_values.items():
                parent_gff_output[parent][id] = value
                
        end_bp.update(input_lengths)
        
        log_output.extend(glog)
        gt_ag_scan_log.extend(slog)
        this_debug_out.extend(debug_lines)
        multi_scan_log.extend(dlog)
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

    if is_genome:
        
        reporter_coords_path = coords_path.joinpath("coords.gff")
        if reporter_coords_path.exists():
            with open(reporter_coords_path, "r") as fp:
                for line in fp:
                    if line.startswith("#"):
                        continue
                    line = line.strip().split("\t")
                    #AGOUTI_SCAF_51|6429119BP|CTG001940_1,CTG001110_1,CTG004120_1	Sapphyre	exon	4815540	4815717	.	-	.	ID=136854;Parent=1.aa.fa;Note=-2;
                    parent = line[0]
                    id = line[-1].split(";")[0].split("=")[1]
                    start = int(line[3])
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
        with open(log_path, "w") as f:
            f.write("Gene,Cut-Indices,Consensus Sequence\n")
            f.write("\n".join(log_output))
        with open(scan_log_path, "w") as f:
            f.write("\n".join(gt_ag_scan_log))
        with open(multi_log_path, "w") as f:
            f.write("\n".join(multi_scan_log))
        with open(internal_path, "w") as f:
            f.write("\n".join(this_debug_out))
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
