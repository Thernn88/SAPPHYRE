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
    ref_seq_len: int,
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
            if IQR <= .2: margin = .05
            if IQR <= .1: margin = .1
            upper_bound = Q3 + (threshold * IQR) + margin
        else:  # if no ref_distances, this is an orthograph, so reject
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


def grab_index_cluster(true_cluster_threshold, true_cluster_raw):
    before_true_clusters = []

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
                before_true_clusters.append(current_cluster)
                current_cluster = [child_full]
                current_index = child_index
        
    if current_cluster:
        before_true_clusters.append(current_cluster)
    return before_true_clusters


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
    assembly: bool,
    top_refs: set,
    ref_kick_path,
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

    true_cluster_raw = [(int(header.split("|")[3].split("_")[1]), header[1:]) for header in original_order]
    true_cluster_raw.sort(key=lambda x: x[0])
    before_true_clusters = grab_index_cluster(true_cluster_threshold, true_cluster_raw)

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
    ref_seq_len = sum(find_ref_len(ref) for ref in ref_seqs) / len(ref_seqs)
    # ref_seq_len = len(ref_seqs[0]) - len(rejected_indices)
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
        ref_seq_len,
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
    passing = original_order_sort(original_order, passing)

    after_data = []
    for candidate in passing:
        after_data.append((int(candidate.id.split("|")[3].split("_")[1]), candidate.id[1:]))
        raw_regulars.extend([candidate.id, candidate.raw])

    reported = []
    if after_data:
        after_data.sort(key=lambda x: x[0])
        after_true_clusters = grab_index_cluster(true_cluster_threshold, after_data)
        after_true_cluster = set(max(after_true_clusters, key=lambda x: len(x)))
        matches = []
        for i, cluster in enumerate(before_true_clusters):
            matches.append((len(after_true_cluster.intersection(set(cluster))) / len( cluster), i))
        matches.sort(key=lambda x: x[0], reverse= True)
        
        best_match = max(matches, key=lambda x: x[0])[1]

        reported = report_overlaps(passing, before_true_clusters[best_match])
    
    # regulars, allowed_columns = delete_empty_columns(raw_regulars, verbose)
    regulars, allowed_columns = delete_empty_columns(raw_regulars)
    # if assembly:
    #     to_be_excluded = make_asm_exclusions(passing, failing)
    # else:
    to_be_excluded = {candidate.id for candidate in failing}
    ids = []
    ref_lens = []
    if passing:  # If candidate added to fasta
        for ref in reference_records:
            ref_lens.append(len(ref.raw) - ref.raw.count("-"))
        for candidate in passing:
            ids.append((int(candidate.id.split("|")[3].split("_")[1]), len(candidate.sequence) - candidate.sequence.count("-")))
        write2Line2Fasta(aa_output, regulars, compress)

    # logging
    if debug:
        logs = [candidate.get_result() for candidate in failing]

    # if valid candidates found, do nt output
    if passing:
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

    clusters = []
    if ids:
        ids.sort(key = lambda x: x[0])

        req_seq_coverage = 0.5

        average_ref_len = sum(ref_lens) / len(ref_lens)

        current_cluster = []
        for child_index, seq_len in ids:
            if not current_cluster:
                current_cluster.append((child_index, seq_len / average_ref_len))
                current_index = child_index
            else:
                if child_index - current_index <= 24:
                    current_cluster.append((child_index, seq_len / average_ref_len))
                    current_index = child_index
                else:
                    if len(current_cluster) >= 2:
                        avg_cluster_coverage = sum([x[1] for x in current_cluster]) / len(current_cluster)
                        clusters.append((current_cluster[0][0], current_cluster[-1][0], avg_cluster_coverage))
                    elif len(current_cluster) == 1:
                        if current_cluster[0][1] > req_seq_coverage:
                            avg_cluster_coverage = sum([x[1] for x in current_cluster]) / len(current_cluster)
                            clusters.append((current_cluster[0][0], current_cluster[0][0], avg_cluster_coverage))
                            
                    current_cluster = [(child_index, seq_len / average_ref_len)]
                    current_index = child_index

        if current_cluster:
            if len(current_cluster) >= 2:
                avg_cluster_coverage = sum([x[1] for x in current_cluster]) / len(current_cluster)
                clusters.append((current_cluster[0][0], current_cluster[-1][0], avg_cluster_coverage))
            elif len(current_cluster) == 1:
                if current_cluster[0][1] > req_seq_coverage:
                    avg_cluster_coverage = sum([x[1] for x in current_cluster]) / len(current_cluster)
                    clusters.append((current_cluster[0][0], current_cluster[0][0], avg_cluster_coverage))
                    
        clusters.sort(key=lambda x: x[2], reverse=True)

        cluster_string = ", ".join([f"{cluster[0]}-{cluster[1]} {(cluster[2]*100):.2f}%" for cluster in clusters])         
        gene = filename.split(".")[0]
        cluster_out = (clusters[0][2] if clusters else 0, f"{gene},{len(ids)},{len(clusters)},{cluster_string}")
    
    return logs, reported, cluster_out


def do_folder(folder, args):
    ALLOWED_EXTENSIONS = {".fa", ".fas", ".fasta", ".gz", ".fq", ".fastq"}
    reference_kick_path = Path(folder, "reference_kicks.log")
    if path.exists(reference_kick_path):
        remove(reference_kick_path)

    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    wanted_aa_path = Path(folder, "trimmed", "aa")
    if not args.map and wanted_aa_path.exists():
        aa_input = wanted_aa_path
        nt_input = Path(folder, "trimmed", "nt")
    else:
        aa_input = Path(folder, "align")
        nt_input = Path(folder, "nt_aligned")
    rocks_db_path = Path(folder, "rocksdb", "sequences", "nt")
    if not rocks_db_path.exists():
        err = f"cannot find dupe databases for {folder}"
        raise FileNotFoundError(err)
    rocksdb_db = RocksDB(str(rocks_db_path))
    is_assembly = rocksdb_db.get("get:isassembly")
    is_assembly = is_assembly == "True"
    is_genome = rocksdb_db.get("get:isgenome")
    is_genome = is_genome == "True"
    top_refs = set(rocksdb_db.get("getall:valid_refs").split(",")[:1])

    # debugging lines, uncomment to manually set a flag
    # is_assembly = False
    # is_genome = False
    del rocksdb_db

    if not aa_input.exists():  # exit early
        printv(
            f"WARNING: Can't find aa folder for taxa {folder}: '{wanted_aa_path}'. Aborting",
            args.verbose,
            0,
        )
        return True, is_assembly, is_genome

    file_inputs = [
        gene
        for gene in aa_input.iterdir()
        if ".aa" in gene.suffixes and gene.suffix in ALLOWED_EXTENSIONS
    ]
    output_path = Path(folder, "outlier", "blosum")
    nt_output_path = path.join(output_path, "nt")
    folder_check(output_path, args.debug)

    compress = not args.uncompress_intermediates or args.compress

    if not (0 < args.threshold < 1.0):
        if 0 < args.threshold <= 100:
            args.threshold = args.threshold / 100
        else:
            raise ValueError(
                "Cannot convert threshold to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.col_cull_percent < 1.0):
        if 0 < args.col_cull_percent <= 100:
            args.col_cull_percent = args.col_cull_percent / 100
        else:
            raise ValueError(
                "Cannot convert column cull percent to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.ref_gap_percent < 1.0):
        if 0 < args.ref_gap_percent <= 100:
            args.ref_gap_percent = args.ref_gap_percent / 100
        else:
            raise ValueError(
                "Cannot convert ref gap percent to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.ref_min_percent < 1.0):
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
                    is_genome or is_assembly,
                    top_refs,
                    reference_kick_path,
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
                    is_genome or is_assembly,
                    top_refs,
                    reference_kick_path,
                ),
            )
    if args.debug:
        log_folder_path = path.join(output_path, "logs")
        global_csv_path = path.join(log_folder_path, "outliers_global.csv")
        with open(global_csv_path, "w", encoding="UTF-8") as global_csv:
            global_csv.write("Gene,Header,Mean_Dist,Ref_Mean,IQR\n")
            reported_nodes = ["True Node,Node,Overlap Percent"]
            cluster_data = []
            for log_data, reported_log, cluster_line in process_data:
                reported_nodes.extend(reported_log)
                cluster_data.append(cluster_line)
                for line in log_data:
                    if line.split(",")[-2] != "Pass":
                        if line[-1] != "\n":
                            line = f"{line}\n"
                        global_csv.write(line)
        
        cluster_data.sort(key=lambda x: x[0], reverse=True)
        cluster_data = [x[1] for x in cluster_data]
            
        with open(path.join(folder, "blosum_clusters.csv"), "w") as f:
            f.write("Gene,Seq count,Cluster count,Cluster ranges\n")
            f.write("\n".join(cluster_data))

        with open(path.join(log_folder_path, "outliers_reported.csv"), "w", encoding="UTF-8") as reported_csv:
            reported_csv.write("\n".join(reported_nodes))

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose)
    return True, is_assembly, is_genome


def main(args):
    is_assembly = None
    is_genome = None
    success = False
    if isinstance(args.INPUT, list):
        success = all([do_folder(Path(folder), args)[0] for folder in args.INPUT])
    elif isinstance(args.INPUT, str):
        success, is_assembly, is_genome = do_folder(Path(args.INPUT), args)
    return success, is_assembly, is_genome


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise RuntimeError(
        MSG,
    )
