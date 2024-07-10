"""Outlier Check."""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from multiprocessing.pool import Pool
from os import mkdir, path, remove
from pathlib import Path
from sys import exit

from msgspec import Struct
from numpy import float16, isnan, nanpercentile,nanmedian
from sapphyre_tools import (
    blosum62_candidate_to_reference,
    blosum62_distance,
    delete_empty_columns,
    find_index_pair,
    get_overlap
)
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, write2Line2Fasta


class Record(Struct):
    id: str
    raw: str
    start: int
    end: int
    sequence: str = None
    upper_bound: float16 = None
    iqr: float16 = None
    mean_distance: float16 = None
    grade: str = None

    def __str__(self) -> str:
        return self.sequence

    def get_result(self):
        return (
            self.id
            + ", "
            + str(self.mean_distance)
            + ", "
            + str(self.upper_bound)
            + ", "
            + self.grade
            + ", "
            + str(self.iqr)
            + "\n"
        )


def nan_check(iterable) -> bool:
    """Checks elements in iterable for numeric values.
    Returns True if numeric value is found. Otherwise returns False.
    """
    return any(not isnan(element) for element in iterable)


def folder_check(taxa_path: Path, debug: bool) -> None:
    """Create subfolders 'aa' and 'nt' to given path."""
    aa_folder = Path(taxa_path, "aa")
    nt_folder = Path(taxa_path, "nt")

    aa_folder.mkdir(parents=True, exist_ok=True)
    nt_folder.mkdir(parents=True, exist_ok=True)

    if debug:
        logs_folder = Path(taxa_path, "logs")
        logs_folder.mkdir(parents=True, exist_ok=True)


def split_sequences(gene_path: str) -> tuple:
    """Reads over a fasta record in the given list and returns a tuple of two smaller lists.
    The first returned list is the reference sequences found, the second returned list
    is the candidate sequences found.
    """
    references = []
    candidates = []
    end_of_references = False
    ref_check = set()
    try:
        for header, sequence in parseFasta(gene_path):
            header = ">" + header
            start, end = find_index_pair(sequence, '-')
            if end_of_references is False:
                # The reference header identifier is present in the header
                if header[-1] == ".":
                    if header[-9] == ":":
                        ref_check.add(header[:-9])

                    references.append(Record(header, sequence, start, end))
                else:
                    end_of_references = True

            if end_of_references is True:
                candidates.append(Record(header, sequence, start, end))
    except ValueError as e:
        print(f"Error in file: {path}, Problem with {header},\n{e}")
        exit(1)
    except TypeError as e:
        print(f"Wrong IO type: {path},\n{e}")
        exit(1)
    return references, candidates, ref_check


def find_index_groups(candidates: list) -> dict:
    """Iterate over a list of candidate fastas as lines of text and finds their start
    and stop indices. Makes a tuple out of the pairs, then uses the
    tuple as a key in two dictionaries. One dictionary stores lists of
    candidates with identical indices, and the other dictionary stores
    the ref set after constraining to those indices.
    """

    def lst():
        return []

    candidate_dict = defaultdict(lst)
    for candidate in candidates:
        start, stop = find_index_pair(candidate.raw, "-")
        candidate.sequence = candidate.raw[start:stop]
        # candidate.id = candidate.id + f"$${start}$${stop}"
        candidate_dict[(start, stop)].append(candidate)
    return candidate_dict

def candidate_pairwise_calls(candidate: Record, refs: list) -> list:
    """Calls calc._pairwise on a candidate and each ref sequence in the list.
    Returns the distances as list. Used to avoid recalculating the ref distances
    for any given candidate index.
    """
    result = []
    for ref in refs:
        result.append(blosum62_candidate_to_reference(candidate.sequence, ref.sequence))
    result.append(0.0)
    return result


def has_minimum_data(
    seq: str, cand_rejected_indices: set, min_data: float, start_offset: int, gap="-"
):
    data_chars = []
    for raw_i, character in enumerate(seq):
        i = raw_i + start_offset
        if i not in cand_rejected_indices:
            data_chars.append(character != gap)

    if len(data_chars) == 0:
        return False

    return sum(data_chars) / len(data_chars) >= min_data


def is_same_variant(header1, header2) -> bool:
    return header1[0:-8] == header2[0:-8]


def calc_ref_distance(ref1, ref2, start, stop) -> float:
    if ref1.start > start:
        start = ref1.start
    if ref2.start > start:
        start = ref2.start

    if ref1.end < stop:
        stop = ref1.end
    if ref2.end < stop:
        stop = ref2.end

    return blosum62_distance(ref1.raw[start:stop], ref2.raw[start:stop])


def compare_means(
    references: list,
    regulars: list,
    candidates_dict: dict,
    threshold: float,
    refs_in_file: int,
    rejected_indices: set,
    index_group_min_bp: int,
    ref_gap_percent: float,
    ref_min_percent: float,
    top_refs: set,
) -> tuple:
    """For each candidate record, finds the index of the first non-gap bp and makes
    matching cuts in the reference sequences. Afterwards finds the mean of the trimmed
    data.
    """
    passing = []
    failing = []
    for index_pair, candidates_at_index in candidates_dict.items():
        # constrain refs here to avoid deepcopy
        for ref in references:
            ref.sequence = ref.raw[index_pair[0] : index_pair[1]]
        # get first candidate in group and check the amount of bp left after
        # column cull
        first_candidate = str(candidates_at_index[0])
        bp_count = 0
        for raw_i, _ in enumerate(first_candidate):
            i = raw_i + index_pair[0]
            if i in rejected_indices:
                continue
            if first_candidate[raw_i] != "-":
                bp_count += 1
            if bp_count >= index_group_min_bp:
                break
        # this line kicks the index group if it fails the bp requirement
        if bp_count < index_group_min_bp:
            for candidate in candidates_at_index:
                mean_distance = "No refs"
                candidate.mean_distance = mean_distance
                candidate.upper_bound = "N/A"
                candidate.grade = "Fail"
                candidate.iqr = "min_candidate_bp"
                failing.append(candidate)
            continue

        # first we have to calculate the reference distances to make the ref mean
        ref_alignments = [
            ref
            for ref in references
            if has_minimum_data(
                ref.sequence, rejected_indices, ref_gap_percent, index_pair[0]
            )
        ]

        ref_distances = []

        # find number of unique ref variants remaining after bp kick
        if references[0].id[-9] == ":":
            found_set = set()
            for ref in references:
                found_set.add(ref.id[:-9])
            found = len(found_set)
        else:
            found = len(references)
        # hadrcode a lower bound in case there are too few refs
        refs_needed = max(2, refs_in_file * ref_min_percent)
        if found < refs_needed:
            has_ref_distances = False
        else:
            ref_distances = [
                # blosum62_distance(ref1.sequence, ref2.sequence)
                calc_ref_distance(ref1, ref2, index_pair[0], index_pair[1])
                for ref1, ref2 in combinations(ref_alignments, 2)
                if not is_same_variant(ref1.id, ref2.id)
            ]
            has_ref_distances = nan_check(ref_distances)
        if has_ref_distances:
            # First quartile (Q1)
            try:
                Q1 = nanpercentile(ref_distances, 25, method="midpoint")
            except IndexError:
                Q1 = 0.0
                print(
                    f'Q1 Runtime Error caused by references in {ref_alignments[0].id.split("|")[0]}',
                )
            except RuntimeError:
                Q1 = 0.0
                print(
                    f'Index Error caused by references in {ref_alignments[0].id.split("|")[0]}',
                )
            # Third quartile (Q3)
            try:
                Q3 = nanpercentile(ref_distances, 75, method="midpoint")
            except IndexError:
                Q3 = 0.0
                print(
                    f'Index Error caused by references in {ref_alignments[0].id.split("|")[0]}',
                )
            except RuntimeError:
                Q3 = 0.0
                print(
                    f'Runtime Error caused by references in {ref_alignments[0].id.split("|")[0]}',
                )
            # Interquartile range (IQR)
            IQR = Q3 - Q1
            margin = 0.02
            #if IQR <= .2: margin = .025
            #if IQR <= .1: margin = .05
            IQR = max(IQR, 0.05)
            upper_bound = Q3 + (threshold * IQR) + margin
        else:  # if no ref_distances, reject
            upper_bound = "N/A"
            IQR = "N/A"
        for candidate in candidates_at_index:
            mean_distance = "No refs"
            candidate.grade = "Ref Fail"
            cupper_bound = upper_bound * 2 if candidate.id.split("|")[1] in top_refs else upper_bound

            if has_ref_distances:
                candidate_distances = candidate_pairwise_calls(
                    candidate,
                    ref_alignments,
                )
                candidate.mean_distance = nanmedian(candidate_distances)
                candidate.iqr = IQR

                candidate.upper_bound = cupper_bound
                if candidate.mean_distance <= cupper_bound:
                    candidate.grade = "Pass"
                    passing.append(candidate)
                else:
                    failing.append(candidate)
            else:
                candidate.mean_distance = mean_distance
                candidate.iqr = IQR
                candidate.upper_bound = cupper_bound
                failing.append(candidate)

    return regulars, passing, failing


def align_col_removal(raw_fed_sequences: list, positions_to_keep: list) -> list:
    """Iterates over each sequence and deletes columns
    that were removed in the empty column removal.
    """
    raw_sequences = [
        i.replace("\n", "") for i in raw_fed_sequences if i.replace("\n", "") != ""
    ]

    result = []

    for i in range(0, len(raw_sequences), 2):
        result.append(raw_sequences[i])

        sequence = raw_sequences[i + 1]

        sequence = [sequence[i * 3 : (i * 3) + 3] for i in positions_to_keep]

        result.append("".join(sequence))

    return result


def aa_to_nt(index: int):
    return [index * 3, index * 3 + 1, index * 3 + 2]


def remove_excluded_sequences(lines: list, excluded: set) -> list:
    """Given a list of fasta lines and a set of headers to be excluded from output,
    returns a list of all valid headers and sequences. Use before the delete_column
    call in the nt portion.
    """
    output = []
    for i in range(0, len(lines), 2):
        if lines[i].strip() not in excluded:
            output.append(lines[i])
            output.append(lines[i + 1])
    return output


def original_order_sort(original: list, candidate_records: list) -> list:
    """Accepts a list of headers and a list of Records. Returns a list of
    Records in the order of the original list.
    """
    candidates = {cand.id: cand for cand in candidate_records}
    output = [candidates.get(header, False) for header in original]
    return [x for x in output if x]



def find_ref_len(ref: str) -> int:
    start, stop = find_index_pair(ref, '-')
    return stop - start


def cull_reference_outliers(reference_records: list) -> list:
    """
    Removes reference sequences which have an unusually large mean
    blosum distance. Finds the constrained blosum distance between
    each reference pull any reference with a mean 1.5x higher than
    the group mean. Returns the remaining references in a list.
    """
    distances_by_index = defaultdict(list)
    all_distances = []
    indices = {i:find_index_pair(reference_records[i].raw, '-') for i in range(len(reference_records))}
    # generate reference distances
    for i, ref1 in enumerate(reference_records[:-1]):
        start1, stop1 = indices[i]
        for j, ref2 in enumerate(reference_records[i+1:],i+1):
            start2, stop2 = indices[j]
            start = max(start1, start2)
            stop = min(stop1, stop2)
            # avoid the occasional rust-nan result
            if start >= stop:
                distances_by_index[i].append(1)
                distances_by_index[j].append(1)
                continue
            dist = blosum62_distance(ref1.raw[start:stop], ref2.raw[start:stop])
            distances_by_index[i].append(dist)
            distances_by_index[j].append(dist)
            all_distances.append(dist)

    total_mean = sum(all_distances) / len(all_distances)
    ALLOWABLE_COEFFICENT = 2
    allowable = max(total_mean * ALLOWABLE_COEFFICENT, 0.3)

    # if a record's mean is too high, cull it
    filtered = []
    for index, distances in distances_by_index.items():
        mean = sum(distances) / len(distances)
        if mean > allowable:
            distances_by_index[index] = None
            filtered.append( (reference_records[index], mean) )
        # else:
        #     distances_by_index[index] = mean
    # get all remaining records
    output = [reference_records[i] for i in range(len(reference_records)) if distances_by_index[i] is not None]
    return output, filtered, total_mean


def report_overlaps(records, true_cluster_headers):
    THRESHOLD = 0
    true, not_true = [], []
    reported = []

    for rec in records:
        if rec.id in true_cluster_headers:
            true.append((rec, *find_index_pair(rec.sequence, "-")))
        else:
            not_true.append((rec, *find_index_pair(rec.sequence, "-")))

    for node_2, start_2, end_2 in not_true:
        for node, start, end in true:
            overlap_coords = get_overlap(
                start_2, end_2, start, end, 1
            )
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]
                percent = overlap_amount / min((end_2 - start_2), (end - start))

                if percent > THRESHOLD:
                    reported.append(f"{node.id},{node_2.id},{percent}")
                    break

    return reported


def do_cluster(ids, ref_coords, max_distance=100):
    clusters = []
    ids.sort(key = lambda x: x[0])

    req_seq_coverage = 0.5

    current_cluster = []
    for i, (child_index, seq_coords, passed) in enumerate(ids):

        coverage = len(seq_coords.intersection(ref_coords)) / len(ref_coords)


        if not current_cluster:
            current_cluster.append((child_index, coverage, i, passed))
            current_index = child_index
        else:
            if child_index - current_index <= max_distance:
                current_cluster.append((child_index, coverage, i, passed))
                current_index = child_index
            else:
                amt_passed = sum([passed for _, _, _, passed in current_cluster])
                if len(current_cluster) >= 2:
                    cluster_data_cols = set()
                    for _, _, index, _ in current_cluster:
                        cluster_data_cols.update(ids[index][1])
                        
                    cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

                    clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage, amt_passed / len(current_cluster)))
                elif len(current_cluster) == 1:
                    if current_cluster[0][1] > req_seq_coverage:
                        cluster_coverage = current_cluster[0][1]
                        clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage, amt_passed / len(current_cluster)))
                        
                current_cluster = [(child_index, coverage, i, passed)]
                current_index = child_index

    if current_cluster:
        amt_passed = sum([passed for _, _, _, passed in current_cluster])
        if len(current_cluster) >= 2:
            cluster_data_cols = set()
            for _, _, index, _ in current_cluster:
                cluster_data_cols.update(ids[index][1])
                
            cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

            clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage, amt_passed / len(current_cluster)))
        elif len(current_cluster) == 1:
            if current_cluster[0][1] > req_seq_coverage:
                cluster_coverage = current_cluster[0][1]
                clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage, amt_passed / len(current_cluster)))
                
    return clusters

# def do_cluster(ids, ref_coords, id_chomp_distance=100, max_distance=120):
#     clusters = []
#     ids.sort(key = lambda x: x[0])
#     grouped_ids = defaultdict(list)
#     for i, (child_index, seq_coords, start, end, child_passed) in enumerate(ids):
#         id = int(child_index.split("_")[0])
#         grouped_ids[id].append((i, child_index, seq_coords, start, end, child_passed))
        
#     ids_ascending = sorted(grouped_ids.keys())
    

#     req_seq_coverage = 0.5

#     current_cluster = []
    
#     for id in ids_ascending:
#         seq_list = grouped_ids[id]
#         if not current_cluster:
#             current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i, cpassed) for i, _, seq_coords, _, _, cpassed in seq_list]
#             current_index = id
#             current_seqs = seq_list
#             current_direction = "bi"
#         else:
#             passed = False
#             passing_direction = None
            
#             if id - current_index <= id_chomp_distance:
#                 for i, child_index, seq_coords, start, end, _ in seq_list:
#                     for _, _, _, current_start, current_end, _ in current_seqs:
#                         this_direction = None
#                         if start == current_start and end == current_end:
#                             this_direction = "bi"
#                         else:
#                             if start == current_start:
#                                 if end >= current_end:
#                                     this_direction = "forward"
#                                 else:
#                                     this_direction = "reverse"
#                             else:
#                                 if start >= current_start:
#                                     this_direction = "forward"
#                                 else:
#                                     this_direction = "reverse"
                     
#                         if current_direction == "bi" or this_direction == "bi" or this_direction == current_direction:
#                             # distance = get_overlap(start, end, current_start, current_end, -max_distance)
#                             # if distance is not None:
#                             #     distance = abs(distance[1] - distance[0])
#                             # if distance is not None and distance < max_distance:
#                             passed = True
#                             passing_direction = this_direction
#                             break
#                     if passed:
#                         break
            
            

                
#             if passed:
#                 current_cluster.extend([(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i, cpassed) for i, _, seq_coords, _, _, cpassed in seq_list])
#                 current_index = id
#                 current_seqs = seq_list
#                 if passing_direction != "bi":
#                     current_direction = passing_direction
#             else:
#                 amt_passed = sum([passed for _, _, _, passed in current_cluster])
#                 if len(current_cluster) >= 2:
#                     cluster_data_cols = set()
#                     for _, _, index, _ in current_cluster:
#                         cluster_data_cols.update(ids[index][1])
                        
#                     cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

#                     clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage, amt_passed / len(current_cluster)))
#                 elif len(current_cluster) == 1:
#                     if current_cluster[0][1] > req_seq_coverage:
#                         cluster_coverage = current_cluster[0][1]
#                         clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage, amt_passed / len(current_cluster)))
                        
#                 current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i, cpassed) for i, _, seq_coords, _, _, cpassed in seq_list]
#                 current_index = id
#                 current_seqs = seq_list
#                 current_direction = "bi"
    
#     if current_cluster:
#         amt_passed = sum([passed for _, _, _, passed in current_cluster])
#         if len(current_cluster) >= 2:
#             cluster_data_cols = set()
#             for _, _, index, _ in current_cluster:
#                 cluster_data_cols.update(ids[index][1])
                
#             cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

#             clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage, amt_passed / len(current_cluster)))
#         elif len(current_cluster) == 1:
#             if current_cluster[0][1] > req_seq_coverage:
#                 cluster_coverage = current_cluster[0][1]
#                 clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage, amt_passed / len(current_cluster)))
                
#     return clusters


def main_process(
    args_input,
    nt_input,
    args_output,
    args_threshold,
    args_references,
    nt_output_path: str,
    debug: bool,
    verbose: int,
    compress: bool,
    true_cluster_threshold: int,
    col_cull_percent: float,
    index_group_min_bp: int,
    ref_gap_percent: float,
    ref_min_percent: int,
    is_genome: bool,
    top_refs: set,
    passing_rescue_percent: float,
    rescue_consensus_percent: float,
):
    keep_refs = not args_references

    file_input = args_input
    filename = path.basename(file_input)

    printv(f"Doing: {filename}", verbose, 2)

    threshold = args_threshold
    aa_output = path.join(args_output, "aa")
    aa_output = path.join(aa_output, filename.rstrip(".gz"))
    reference_records, candidate_records, ref_check = split_sequences(file_input)
    
    original_order = [candidate.id for candidate in candidate_records]

    candidates_dict = find_index_groups(candidate_records)
    # if not assembly:
    #     candidates_dict = find_index_groups(candidate_records)
    # else:
    #     candidates_dict = find_asm_index_groups(candidate_records)

    # calculate indices that have valid data columns
    rejected_indices = set()
    # ref_seqs = reference_sequences[1::2]
    ref_seqs = [ref.raw for ref in reference_records]
    for i in range(len(ref_seqs[0])):
        percent_of_non_dash = len([ref[i] for ref in ref_seqs if ref[i] != "-"]) / len(
            ref_seqs,
        )
        if percent_of_non_dash <= col_cull_percent:
            rejected_indices.add(i)
    # find number of unique reference variants in file, use for refs_in_file
    if ref_check:
        refs_in_file = len(ref_check)
    else:
        refs_in_file = len(ref_seqs)
    regulars = []
    if keep_refs:
        for ref in reference_records:
            regulars.append(ref.id)
            regulars.append(ref.raw)

    raw_regulars, passing, failing = compare_means(
        reference_records,
        regulars,
        candidates_dict,
        threshold,
        refs_in_file,
        rejected_indices,
        index_group_min_bp,
        ref_gap_percent,
        ref_min_percent,
        top_refs,
    )
    logs = []
    # if passing:
    #     if assembly:
    #         # save any failed subseqs if the original seq had a passing segment
    #         any_passed = get_passing_headers(passing)
    #         to_save, failing = save_partial_fails(failing, any_passed)
    #         passing.extend(to_save)
    #         passing, header_to_indices = remake_introns(passing)

    ids = []
    cluster_out = None
    saved_fails = set()
    save_data_cols = {}
    if is_genome: # Rescue
        passed = {candidate.id for candidate in passing}
        ref_coords = set()
        # get_id = lambda header: header.split("|")[3].replace("NODE_","")
        get_id = lambda header: int(header.split("|")[3].split("&&")[0].split("_")[1])
        for ref in reference_records:
            start, end = find_index_pair(ref.raw, "-")
            for i, let in enumerate(ref.raw[start:end], start):
                if let != "-":
                    ref_coords.add(i)
        candidate_flex_consensus = defaultdict(set)
        for candidate in candidate_records:
            data_cols = set()  
            for i, let in enumerate(candidate.raw[candidate.start:candidate.end], candidate.start):
                if let != "-":
                    if candidate.id in passed:
                        candidate_flex_consensus[i].add(let)
                    data_cols.add(i)
                    
            save_data_cols[candidate.id] = data_cols
            ids.append((get_id(candidate.id), data_cols, candidate.id in passed))

        clusters = []
        cluster_sets = set()
        if ids:
            clusters = do_cluster(ids, ref_coords, true_cluster_threshold)
            cluster_sets = [set(range(start, end+1)) for start, end, _, perc_passed in clusters if perc_passed > passing_rescue_percent]
            if cluster_sets:
                flattened_set = set.union(*cluster_sets)

        if cluster_sets:
            for i, candidate in enumerate(failing):
                if get_id(candidate.id) in flattened_set:
                    
                    matches = 0
                    for x, let in enumerate(candidate.raw[candidate.start:candidate.end], candidate.start):
                        if not candidate_flex_consensus[x] or let in candidate_flex_consensus[x]:
                            matches += 1
                            
                    if matches / (candidate.end - candidate.start) < rescue_consensus_percent:
                        continue
                    
                    candidate.grade = "Saved By Cluster / " + candidate.grade
                    passing.append(candidate)
                    saved_fails.add(i)

    passing = original_order_sort(original_order, passing)
    for candidate in passing:
        raw_regulars.extend([candidate.id, candidate.raw])

    # logging
    if debug:
        logs = [candidate.get_result() for candidate in failing]

    if saved_fails:
        failing = [candidate for i, candidate in enumerate(failing) if i not in saved_fails]

    #update cluster
    if is_genome:
        ids = []
        for candidate in passing:
            data_cols = save_data_cols[candidate.id]
            # data_cols = {i for i, let in enumerate(candidate.raw[candidate.start:candidate.end], candidate.start) if let != "-"}
            ids.append((get_id(candidate.id), data_cols, candidate.id in passed))
            
        clusters = do_cluster(ids, ref_coords, true_cluster_threshold)

        cluster_string = ", ".join([f"{cluster[0]}-{cluster[1]} {(cluster[2]*100):.2f}%" for cluster in clusters])         
        gene = filename.split(".")[0]
        cluster_out = (gene, f"{gene},{len(ids)},{len(clusters)},{cluster_string}")

    # regulars, allowed_columns = delete_empty_columns(raw_regulars, verbose)
    regulars, allowed_columns = delete_empty_columns(raw_regulars)
    # if assembly:
    #     to_be_excluded = make_asm_exclusions(passing, failing)
    # else:

    to_be_excluded = {candidate.id for candidate in failing}

    if passing:  # If candidate added to fasta
        write2Line2Fasta(aa_output, regulars, compress)

        nt_file = filename.replace(".aa.", ".nt.")
        nt_input_path = path.join(nt_input, nt_file)
        if not path.exists(nt_output_path):
            mkdir(nt_output_path)
        nt_output_path = path.join(nt_output_path, nt_file.rstrip(".gz"))
        lines = []
        for header, sequence in parseFasta(nt_input_path):
            lines.append(">" + header)
            lines.append(sequence)

        non_empty_lines = remove_excluded_sequences(lines, to_be_excluded)
        # if assembly:
        #     non_empty_lines = align_intron_removal(non_empty_lines, header_to_indices)
        non_empty_lines = align_col_removal(non_empty_lines, allowed_columns)

        write2Line2Fasta(nt_output_path, non_empty_lines, compress)

    
    
    return logs, cluster_out


def do_folder(folder, args, is_genome, gene_source):
    ALLOWED_EXTENSIONS = {".fa", ".fas", ".fasta", ".gz", ".fq", ".fastq"}
    reference_kick_path = Path(folder, "reference_kicks.log")
    if path.exists(reference_kick_path):
        remove(reference_kick_path)

    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    wanted_aa_path = Path(folder, "outlier", gene_source, "aa")
    if not args.map and wanted_aa_path.exists():
        aa_input = wanted_aa_path
        nt_input = Path(folder, "outlier", gene_source, "nt")
    elif gene_source == "trimmed" or gene_source == "motif":
        aa_input = Path(folder, gene_source, "aa")
        nt_input = Path(folder, gene_source, "nt")
        if not aa_input.exists():
            if gene_source == "motif":
                aa_input = Path(folder, "trimmed", "aa")
                nt_input = Path(folder, "trimmed", "nt")
            else:
                aa_input = Path(folder, "align")
                nt_input = Path(folder, "nt_aligned")
    rocks_db_path = Path(folder, "rocksdb", "sequences", "nt")
    if not rocks_db_path.exists():
        err = f"cannot find dupe databases for {folder}"
        raise FileNotFoundError(err)
    rocksdb_db = RocksDB(str(rocks_db_path))
    top_refs = set(rocksdb_db.get("getall:valid_refs").split(",")[:1])

    # debugging lines, uncomment to manually set a flag
    # is_assembly = False
    # is_genome = False
    # top_refs = set()
    del rocksdb_db

    if not aa_input.exists():  # exit early
        printv(
            f"WARNING: Can't find aa folder for taxa {folder}: '{wanted_aa_path}'. Aborting",
            args.verbose,
            0,
        )
        return True

    file_inputs = [
        gene
        for gene in aa_input.iterdir() 
        if ".aa" in gene.suffixes and gene.suffix in ALLOWED_EXTENSIONS
    ]
    output_path = Path(folder, "outlier", "blosum")
    nt_output_path = path.join(output_path, "nt")
    folder_check(output_path, args.debug)

    compress = not args.uncompress_intermediates or args.compress

    #  convert threshold to percent
    args.threshold = args.threshold / 100
    if not 0 < args.col_cull_percent < 1.0:
        if 0 < args.col_cull_percent <= 100:
            args.col_cull_percent = args.col_cull_percent / 100
        else:
            raise ValueError(
                "Cannot convert column cull percent to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not 0 < args.ref_gap_percent < 1.0:
        if 0 < args.ref_gap_percent <= 100:
            args.ref_gap_percent = args.ref_gap_percent / 100
        else:
            raise ValueError(
                "Cannot convert ref gap percent to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not 0 < args.ref_min_percent < 1.0:
        if 0 < args.ref_min_percent <= 100:
            args.ref_min_percent = args.ref_min_percent / 100
        else:
            raise ValueError(
                "Cannot convert ref min percent to a percent. Use a decimal or a whole number between 0 and 100"
            )

    file_inputs.sort(key=lambda x: x.stat().st_size, reverse=True)
    if args.processes > 1:
        arguments = []
        for gene in file_inputs:
            arguments.append(
                (
                    gene,
                    nt_input,
                    output_path,
                    args.threshold,
                    args.no_references,
                    nt_output_path,
                    args.debug,
                    args.verbose,
                    compress,
                    args.true_cluster_threshold,
                    args.col_cull_percent,
                    args.index_group_min_bp,
                    args.ref_gap_percent,
                    args.ref_min_percent,
                    # args.internal_consensus_threshold,
                    # args.internal_kick_threshold,
                    is_genome,
                    top_refs,
                    args.rescue_passing_cluster,
                    args.rescue_consensus_percent,
                ),
            )

        with Pool(args.processes) as pool:
            process_data = pool.starmap(main_process, arguments, chunksize=1)
    else:
        process_data = []
        for gene in file_inputs:
            process_data.append(
                main_process(
                    gene,
                    nt_input,
                    output_path,
                    args.threshold,
                    args.no_references,
                    nt_output_path,
                    args.debug,
                    args.verbose,
                    compress,
                    args.true_cluster_threshold,
                    args.col_cull_percent,
                    args.index_group_min_bp,
                    args.ref_gap_percent,
                    args.ref_min_percent,
                    # args.internal_consensus_threshold,
                    # args.internal_kick_threshold,
                    is_genome,
                    top_refs,
                    args.rescue_passing_cluster,
                    args.rescue_consensus_percent,
                ),
            )
    if args.debug:
        log_folder_path = path.join(output_path, "logs")
        global_csv_path = path.join(log_folder_path, "outliers_global.csv")
        with open(global_csv_path, "w", encoding="UTF-8") as global_csv:
            global_csv.write("Header,Candidate_Median,Upper_Bound, Grade, IQR\n")
            cluster_data = []
            for log_data, cluster_line in process_data:
                cluster_data.append(cluster_line)
                for line in log_data:
                    if line.split(",")[-2] != "Pass":
                        if line[-1] != "\n":
                            line = f"{line}\n"
                        global_csv.write(line)
        
        if any(cluster_data):
            cluster_data = [x for x in cluster_data if x]
            cluster_data.sort(key=lambda x: x[0])
            cluster_data = [x[1] for x in cluster_data]
                
            with open(path.join(folder, "blosum_clusters.csv"), "w") as f:
                f.write("Gene,Seq count,Cluster count,Cluster ranges\n")
                f.write("\n".join(cluster_data))

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose)
    return True


def main(args, is_genome = False, gene_source = "trimmed"):
    success = False
    if isinstance(args.INPUT, list):
        success = all(do_folder(Path(folder), args, is_genome, gene_source)[0] for folder in args.INPUT)
    elif isinstance(args.INPUT, str):
        success = do_folder(Path(args.INPUT), args, is_genome, gene_source)
    return success


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise RuntimeError(
        MSG,
    )