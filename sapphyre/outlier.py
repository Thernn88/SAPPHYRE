"""Outlier Check."""
from __future__ import annotations

import os
import sys
from copy import deepcopy
from itertools import combinations
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import rmtree

import numpy as np
import phymmr_tools as bd
from msgspec import Struct
from phymmr_tools import constrained_distance, dumb_consensus

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, write2Line2Fasta

ALLOWED_EXTENSIONS = (".fa", ".fas", ".fasta", ".fa", ".gz", ".fq", ".fastq")


class Record(Struct):
    id: str
    sequence: str
    raw: str
    upper_bound: np.float16 = None
    iqr: np.float16 = None
    mean_distance: np.float16 = None
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
    return any(not np.isnan(element) for element in iterable)


def fast_pop(array: list, i: int):
    temp = array[i]
    array[i] = array[-1]
    array[-1] = temp
    return array.pop()


def folder_check(path: Path, debug: bool) -> None:
    """Create subfolders 'aa' and 'nt' to given path."""
    aa_folder = Path(path, "aa")
    nt_folder = Path(path, "nt")

    rmtree(aa_folder, ignore_errors=True)
    rmtree(nt_folder, ignore_errors=True)

    aa_folder.mkdir(parents=True, exist_ok=True)
    nt_folder.mkdir(parents=True, exist_ok=True)

    if debug:
        logs_folder = Path(path, "logs")
        logs_folder.mkdir(parents=True, exist_ok=True)


def get_headers(lines: list) -> list:
    """Returns a list of every other line in the provided argument. Used to get
    header names from a list of sequences.
    """
    result = []
    for i in range(0, len(lines), 2):
        result.append(lines[i])
    return result


def split_sequences_ex(path: str, excluded: set) -> tuple:
    """Reads over a fasta record in the given list and returns a tuple of two smaller lists.
    The first returned list is the reference sequences found, the second returned list
    is the candidate sequences found.
    """
    references = []
    candidates = []
    end_of_references = False
    ref_check = set()
    try:
        for header, sequence in parseFasta(path):
            header = ">" + header

            if end_of_references is False:
                # The reference header identifier is present in the header
                if header[-1] == ".":
                    if header.split("|")[1].lower() in excluded:
                        continue
                    if header[-9] == ":":
                        ref_check.add(header[:-9])

                    references.append(header)
                    references.append(sequence)
                else:
                    end_of_references = True

            if end_of_references is True:
                candidates.append(header)
                candidates.append(sequence)
    except ValueError:
        print(f"Error in file: {path}, error reading {header}")
        sys.exit(1)
    except TypeError:
        print(f"Wrong IO type: {path}")
        sys.exit(1)
    return references, candidates, ref_check


def split_sequences(path: str) -> tuple:
    """Reads over a fasta record in the given list and returns a tuple of two smaller lists.
    The first returned list is the reference sequences found, the second returned list
    is the candidate sequences found.
    """
    references = []
    candidates = []
    end_of_references = False
    ref_check = set()
    try:
        for header, sequence in parseFasta(path):
            header = ">" + header

            if end_of_references is False:
                # The reference header identifier is present in the header
                if header[-1] == ".":
                    if header[-9] == ":":
                        ref_check.add(header[:-9])

                    references.append(header)
                    references.append(sequence)
                else:
                    end_of_references = True

            if end_of_references is True:
                candidates.append(header)
                candidates.append(sequence)
    except ValueError:
        print(f"Error in file: {path}, Problem with {header}")
        sys.exit(1)
    except TypeError:
        print(f"Wrong IO type: {path}")
        sys.exit(1)
    return references, candidates, ref_check


def make_indices(sequence: str, gap_character="-") -> tuple:
    """Finds the index of the first and last non-gap bp in a sequence.
    Returns the start value and the end values + 1 as a tuple.
    """
    start = None
    end = None
    for i, character in enumerate(sequence):
        if character != gap_character:
            start = i
            break
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != gap_character:
            end = i + 1
            break
    if start is None or end is None:
        raise ValueError
    return start, end


def constrain_data_lines(lines: list, start: int, end: int) -> tuple:
    """Given a start and end value, iterates over the list of sequences and
    trims the non-header lines to given values. No return, mutates the original data.
    """
    full = []
    heads = []
    for i in range(0, len(lines), 2):
        newline = lines[i + 1][start:end]
        # Do work if the string contains a non-gap character.
        if any(character != "-" for character in newline):
            full.append(lines[i])
            full.append(newline)
            heads.append(lines[i])
    return (full, heads)


def convert_to_record_objects(lines: list) -> list:
    """Given a list of stings from a fasta file, returns a list of Sequence objects
    from the biopython module. This allows us to make a MultipleSequenceAlignment
    object later.
    """
    return [Record(lines[i], lines[i + 1], None) for i in range(0, len(lines), 2)]


def find_index_groups(references: list, candidates: list) -> tuple:
    """Iterate over a list of candidate fastas as lines of text and finds their start
    and stop indices. Makes a tuple out of the pairs, then uses the
    tuple as a key in two dictionaries. One dictionary stores lists of
    candidates with identical indices, and the other dictionary stores
    the ref set after constraining to those indices.
    """
    candidate_dict = {}
    for i in range(0, len(candidates), 2):
        sequence = candidates[i + 1]
        raw_seq = sequence
        index_tuple = make_indices(sequence)
        start, stop = index_tuple
        lines = [candidates[i], candidates[i + 1]]
        lines, _ = constrain_data_lines(lines, start, stop)
        cand_seq = Record(lines[0], lines[1], raw_seq)
        made_already = candidate_dict.get(index_tuple, False)
        if not made_already:
            seq_list = []
            seq_list.append(cand_seq)
            candidate_dict[index_tuple] = seq_list
        else:
            made_already.append(cand_seq)
            candidate_dict[index_tuple] = made_already
    # after processing candidates, make appropriate ref sets
    reference_dict = {}
    raw_ref_dict = {}
    for key in candidate_dict:
        start, stop = key
        ref_lines = deepcopy(references)
        raw_ref_dict[key] = ref_lines
        ref_lines, _ = constrain_data_lines(ref_lines, start, stop)
        reference_dict[key] = ref_lines
    return reference_dict, candidate_dict


def make_ref_mean(matrix: list, ignore_zeros=False) -> float:
    """Iterates over a distance matrix and calculates the mean value of all found
    distances. Returns the value as a float. If ignore_zeros is enabled, ignores
    any distance value of zero.
    """
    isum = 0
    zeros_found = 0
    total_number = 0
    for row in matrix:
        for column in row:
            isum += column
            total_number += 1
            if column == 0:
                zeros_found += 1
    if ignore_zeros:
        total_number -= zeros_found
    return isum / total_number


def candidate_pairwise_calls(candidate: Record, refs: list) -> list:
    """Calls calc._pairwise on a candidate and each ref sequence in the list.
    Returns the distances as list. Used to avoid recalculating the ref distances
    for any given candidate index.
    """
    result = []
    for ref in refs:
        result.append(bd.blosum62_candidate_to_reference(str(candidate), str(ref)))
    result.append(0.0)
    return result


def has_minimum_data(
    seq: str, cand_rejected_indices: set, min_data: float, start_offset: int, gap="-",
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
    if header1[-9] == ":" and header2[-9] == ":" and header1[:-9] == header2[:-9]:
        return True
    return False


def compare_means(
    references: list,
    ref_dict: dict,
    candidates_dict: dict,
    threshold: float,
    keep_refs: bool,
    refs_in_file: int,
    rejected_indices: set,
    index_group_min_bp: int,
    ref_gap_percent: float,
    ref_min_percent: float,
) -> tuple:
    """For each candidate record, finds the index of the first non-gap bp and makes
    matching cuts in the reference sequences. Afterwards finds the mean of the trimmed
    data.
    """
    regulars = []
    passing = []
    failing = []
    if keep_refs:
        for line in references:
            regulars.append(line)
    for index_pair, current_refs in ref_dict.items():
        # get candidates in this index group
        candidates_at_index = candidates_dict[index_pair]

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
        ref_records = convert_to_record_objects(current_refs)
        ref_alignments = [
            seq
            for seq in ref_records
            if has_minimum_data(
                seq.sequence, rejected_indices, ref_gap_percent, index_pair[0],
            )
        ]

        ref_distances = []

        # find number of unique ref variants remaining after bp kick
        if ref_records[0].id[-9] == ":":
            found_set = set()
            for ref in ref_records:
                found_set.add(ref.id[:-9])
            found = len(found_set)
        else:
            found = len(ref_records)

        if found / refs_in_file < ref_min_percent:
            has_ref_distances = False
        else:
            for seq1, seq2 in combinations(ref_alignments, 2):
                if is_same_variant(seq1.id, seq2.id):
                    continue
                ref1 = str(seq1)
                ref2 = str(seq2)
                ref_distances.append(bd.blosum62_distance(ref1, ref2))
            has_ref_distances = nan_check(ref_distances)
        if has_ref_distances:
            # First quartile (Q1)
            try:
                Q1 = np.nanpercentile(ref_distances, 25, method="midpoint")
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
                Q3 = np.nanpercentile(ref_distances, 75, method="midpoint")
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
            upper_bound = Q3 + (threshold * IQR) + 0.02
        else:  # if no ref_distances, this is an orthograph, so reject
            upper_bound = "N/A"
            IQR = "N/A"
        for candidate in candidates_at_index:
            mean_distance = "No refs"
            candidate.grade = "Ref Fail"
            if has_ref_distances:
                candidate_distances = candidate_pairwise_calls(
                    candidate, ref_alignments,
                )
                candidate.mean_distance = np.nanmean(candidate_distances)
                candidate.iqr = IQR
                candidate.upper_bound = upper_bound
                if candidate.mean_distance <= upper_bound:
                    candidate.grade = "Pass"
                    passing.append(candidate)
                else:
                    failing.append(candidate)
            else:
                candidate.mean_distance = mean_distance
                candidate.iqr = IQR
                candidate.upper_bound = upper_bound
                failing.append(candidate)

    return regulars, passing, failing


def delete_empty_columns(raw_fed_sequences: list, verbose: bool) -> tuple[list, list]:
    """Iterates over each sequence and deletes columns
    that consist of 100% dashes.
    """
    result = []
    sequences = []
    raw_sequences = [
        i.replace("\n", "") for i in raw_fed_sequences if i.replace("\n", "") != ""
    ]

    for i in range(0, len(raw_sequences), 2):
        sequences.append(raw_sequences[i + 1])

    positions_to_keep = []
    if sequences:
        for i in range(len(sequences[0])):
            for sequence in sequences:
                if sequence[i] != "-":
                    positions_to_keep.append(i)
                    break

        for i in range(0, len(raw_sequences), 2):
            try:
                sequence = [raw_sequences[i + 1][x] for x in positions_to_keep]
                result.append(raw_sequences[i])
            except IndexError:
                printv(
                    f"WARNING: Sequence length is not the same as other sequences: {raw_sequences[i]}",
                    verbose,
                    0,
                )
                continue
            sequence = "".join(sequence)

            result.append(sequence)

    return result, positions_to_keep


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


def make_exclusion_set(path: str) -> set:
    """Reads a file at a given path and returns a set containing
    each line. Used to make a taxa exclusion list.
    """
    excluded = set()
    if not path:
        return excluded
    with open(path) as f:
        for line in f:
            excluded.add(line.rstrip())
    return excluded


def delete_excluded(lines: list, excluded: set) -> list:
    """Given a list of fasta lines and a set of headers to be excluded from output,
    returns a list of all valid headers and sequences.
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
    col_cull_percent: float,
    index_group_min_bp: int,
    ref_gap_percent: float,
    ref_min_percent: int,
    internal_consensus_threshold: float,
    internal_kick_threshold: int,
):
    keep_refs = not args_references

    file_input = args_input
    filename = os.path.basename(file_input)

    printv(f"Doing: {filename}", verbose, 2)

    threshold = args_threshold / 100
    aa_output = os.path.join(args_output, "aa")
    aa_output = os.path.join(aa_output, filename.rstrip(".gz"))

    reference_sequences, candidate_sequences, ref_check = split_sequences(file_input)
    original_order = list(candidate_sequences[0::2])
    ref_dict, candidates_dict = find_index_groups(
        reference_sequences, candidate_sequences,
    )

    # calculate indices that have valid data columns
    rejected_indices = set()
    ref_seqs = reference_sequences[1::2]
    for i in range(len(ref_seqs[0])):
        percent_of_non_dash = len([ref[i] for ref in ref_seqs if ref[i] != "-"]) / len(
            ref_seqs,
        )
        if percent_of_non_dash <= col_cull_percent:
            rejected_indices.add(i)

    # find number of unique reference variants in file, use for refs_in_file
    refs_in_file = len(ref_check) if ref_check else len(ref_seqs)
    raw_regulars, passing, failing = compare_means(
        reference_sequences,
        ref_dict,
        candidates_dict,
        threshold,
        keep_refs,
        refs_in_file,
        rejected_indices,
        index_group_min_bp,
        ref_gap_percent,
        ref_min_percent,
    )
    logs = []
    if passing:
        consensus = dumb_consensus(
            [cand.raw for cand in passing], internal_consensus_threshold,
        )
        for i, candidate in enumerate(passing):
            distance = constrained_distance(consensus, candidate.raw)
            if distance >= internal_kick_threshold:
                candidate.grade = "Internal Fail"
                failing.append(candidate)
                passing[i] = None
    passing = [x for x in passing if x is not None]
    passing = original_order_sort(original_order, passing)
    for candidate in passing:
        raw_regulars.extend([candidate.id, candidate.raw])
    regulars, allowed_columns = delete_empty_columns(raw_regulars, verbose)

    to_be_excluded = {candidate.id for candidate in failing}
    if passing:  # If candidate added to fasta
        write2Line2Fasta(aa_output, regulars, compress)

    # logging
    if debug:
        logs = [candidate.get_result() for candidate in failing]

    # if valid candidates found, do nt output
    if passing:
        nt_file = filename.replace(".aa.", ".nt.")
        nt_input_path = os.path.join(nt_input, nt_file)
        if not os.path.exists(nt_output_path):
            os.mkdir(nt_output_path)
        nt_output_path = os.path.join(nt_output_path, nt_file.rstrip(".gz"))
        lines = []
        for header, sequence in parseFasta(nt_input_path):
            lines.append(">" + header)
            lines.append(sequence)

        non_empty_lines = remove_excluded_sequences(lines, to_be_excluded)
        non_empty_lines = align_col_removal(non_empty_lines, allowed_columns)

        write2Line2Fasta(nt_output_path, non_empty_lines, compress)
    return logs


def do_folder(folder, args):
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    wanted_aa_path = Path(folder, "trimmed", "aa")
    if wanted_aa_path.exists():
        aa_input = wanted_aa_path
        nt_input = Path(folder, "trimmed", "nt")
    else:
        aa_input = Path(folder, "align")
        nt_input = Path(folder, "nt_aligned")

    printv(f"Processing: {os.path.basename(folder)}", args.verbose, 0)

    if not aa_input.exists():  # exit early
        printv(
            f"WARNING: Can't find aa folder for taxa {folder}: '{wanted_aa_path}'. Aborting",
            args.verbose,
            0,
        )
        return

    file_inputs = [
        gene
        for gene in aa_input.iterdir()
        if ".aa" in gene.suffixes and gene.suffix in ALLOWED_EXTENSIONS
    ]
    output_path = Path(folder, args.output)
    nt_output_path = os.path.join(output_path, "nt")
    folder_check(output_path, args.debug)

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
                    args.compress,
                    args.col_cull_percent,
                    args.index_group_min_bp,
                    args.ref_gap_percent,
                    args.ref_min_percent,
                    args.internal_consensus_threshold,
                    args.internal_kick_threshold,
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
                    args.compress,
                    args.col_cull_percent,
                    args.index_group_min_bp,
                    args.ref_gap_percent,
                    args.ref_min_percent,
                    args.internal_consensus_threshold,
                    args.internal_kick_threshold,
                ),
            )
    if args.debug:
        log_folder_path = os.path.join(output_path, "logs")
        global_csv_path = os.path.join(log_folder_path, "outliers_global.csv")
        global_csv = open(global_csv_path, "w", encoding="UTF-8")
        global_csv.write("Gene,Header,Mean_Dist,Ref_Mean,IQR\n")

    for log_data in process_data:
        if args.debug:
            for line in log_data:
                if line.split(",")[-2] != "Pass":
                    if line[-1] != "\n":
                        line = f"{line}\n"
                    global_csv.write(line)
    if args.debug:
        global_csv.close()

    printv(f"Done! Took {time_keeper.differential():.2f}s", args.verbose)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for folder in args.INPUT:
        do_folder(Path(folder), args)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise Exception(
        msg,
    )
