"""
Outlier Check
"""
from __future__ import annotations

import os
from copy import deepcopy
from itertools import combinations
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import rmtree
from statistics import mean
import numpy as np
import phymmr_tools as bd

from .utils import printv, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode
ALLOWED_EXTENSIONS = (".fa", ".fas", ".fasta", ".fa", ".gz", ".fq", ".fastq")

class Record:
    def __init__(self, head, seq, raw_seq=None):
        self.id = head
        self.sequence = seq
        if raw_seq is None:
            self.raw = seq
        else:
            self.raw = raw_seq

    def __hash__(self):
        return hash(self.id + self.sequence)

    def __str__(self):
        return self.sequence


def nan_check(iterable) -> bool:
    """
    Checks elements in iterable for numeric values.
    Returns True if numeric value is found. Otherwise returns False.
    """
    for element in iterable:
        if not np.isnan(element):
            return True
    return False


def taxa_sort(records: list) -> list:
    """
    Iterates over a list of candidates and creates Records. Sorts by
    taxa, then makes a fasta list. Returns the list.
    """
    records = []
    for header, sequence in records:
        specimen = Record(header, sequence)
        records.append(specimen)
    records.sort(key=lambda x: (x.id.split("|")[2], x.id.split("|")[3]))
    output = []
    for record in records:
        output.append(record.id)
        output.append(record.sequence)
    return output


def original_sort(headers, records) -> list:
    """
    Returns candidate sequences to their original order.
    """
    output = []
    record_dict = {}
    for header, sequence in records:
        record_dict[header] = sequence
    for header in headers:
        sequence = record_dict.get(header, False)
        if sequence:
            output.append((header,sequence))
    return output


def folder_check(path: Path, debug: bool) -> None:
    """Create subfolders 'aa' and 'nt' to given path."""
    aa_folder = Path(path, "aa")
    nt_folder = Path(path, "nt")

    rmtree(aa_folder, ignore_errors = True)
    rmtree(nt_folder, ignore_errors = True)

    aa_folder.mkdir(parents=True, exist_ok=True)
    nt_folder.mkdir(parents=True, exist_ok=True)

    if debug:
        logs_folder = Path(path, "logs")
        logs_folder.mkdir(parents=True, exist_ok=True)

def split_sequences(path: str, excluded: set) -> tuple:
    """
    Reads over a fasta record in the given list and returns a tuple of two smaller lists.
    The first returned list is the reference sequences found, the second returned list
    is the candidate sequences found.
    """
    bad_names = {"bombyx_mori", "danaus_plexippus"}
    references = []
    candidates = []

    end_of_references = False
    try:
        for header, sequence in parseFasta(path):
            if end_of_references is False:
                # The reference header identifier is present in the header
                if header[-1] == ".":
                    if header.split("|")[1].lower() in bad_names:
                        excluded.add(header)

                    references.append((header,sequence))
                else:
                    end_of_references = True

            if end_of_references is True:
                candidates.append((header,sequence))
    except ValueError:
        print(f'Error in file: {path}')
        raise ValueError("Found sequences of different length")
    except TypeError:
        print(f'Wrong IO type: {path}')
        raise TypeError
    return references, candidates


def make_indices(sequence: str, gap_character="-") -> tuple:
    """
    Finds the index of the first and last non-gap bp in a sequence.
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
        raise ValueError()
    return start, end


def constrain_data_lines(records: list, start: int, end: int) -> tuple:
    """
    Given a start and end value, iterates over the list of sequences and
    trims the non-header lines to given values. No return, mutates the original data.
    """
    full = []
    heads = []
    for header, sequence in records:
        newline = sequence[start:end]
        # Do work if the string contains a non-gap character.
        if any(character != "-" for character in newline):
            full.append((header, sequence))
            heads.append(header)
    return (full, heads)


def convert_to_record_objects(records: list) -> list:
    """
    Given a list of stings from a fasta file, returns a list of Sequence objects
    from the biopython module. This allows us to make a MultipleSequenceAlignment
    object later.
    """
    return [Record(header, sequence) for header, sequence in records]


def find_index_groups(references: list, candidates: list) -> tuple:
    """
    Iterate over a list of candidate fastas as lines of text and finds their start
    and stop indices. Makes a tuple out of the pairs, then uses the
    tuple as a key in two dictionaries. One dictionary stores lists of
    candidates with identical indices, and the other dictionary stores
    the ref set after constraining to those indices.
    """
    candidate_dict = {}
    for header, sequence in candidates:
        raw_seq = sequence
        index_tuple = make_indices(sequence)
        start, stop = index_tuple
        lines = [(header, sequence)]
        lines, _ = constrain_data_lines(lines, start, stop)
        cand_seq = Record(lines[0][0], lines[0][1], raw_seq)
        made_already = candidate_dict.get(index_tuple, False)
        if not made_already:
            seq_set = set()
            seq_set.add(cand_seq)
            candidate_dict[index_tuple] = seq_set
        else:
            made_already.add(cand_seq)
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
    """
    Iterates over a distance matrix and calculates the mean value of all found
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
    """
    Calls calc._pairwise on a candidate and each ref sequence in the list.
    Returns the distances as list. Used to avoid recalculating the ref distances
    for any given candidate index.
    """
    result = []
    for ref in refs:
        result.append(bd.blosum62_distance(str(candidate), str(ref)))
    result.append(0.0)
    return result


def compare_means(
    references: list,
    candidates: list,
    threshold: float,
    excluded_headers: set,
    keep_refs: bool,
    sort: str,
) -> tuple:
    """
    For each candidate record, finds the index of the first non-gap bp and makes
    matching cuts in the reference sequences. Afterwards finds the mean of the trimmed
    data.
    """
    regulars = []
    outliers = []
    if keep_refs:
        regulars.extend(references)
    ref_dict, candidates_dict = find_index_groups(references, candidates)
    to_add_later = []
    for index_pair, current_refs in ref_dict.items():
        # start, stop = index_pair
        candidates_at_index = candidates_dict[index_pair]
        # first we have to calculate the reference distances to make the ref mean
        ref_alignments = [
            seq
            for seq in convert_to_record_objects(current_refs)
            if seq.id not in excluded_headers
        ]
        ref_distances = []

        for seq1, seq2 in combinations(ref_alignments, 2):
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
                print(f'Q1 Runtime Error caused by references in {ref_alignments[0].id.split("|")[0]}')
            except RuntimeError:
                Q1 = 0.0
                print(f'Index Error caused by references in {ref_alignments[0].id.split("|")[0]}')
            # Third quartile (Q3)
            try:
                Q3 = np.nanpercentile(ref_distances, 75, method="midpoint")
            except IndexError:
                Q3 = 0.0
                print(f'Index Error caused by references in {ref_alignments[0].id.split("|")[0]}')
            except RuntimeError:
                Q3 = 0.0
                print(f'Runtime Error caused by references in {ref_alignments[0].id.split("|")[0]}')
            # Interquartile range (IQR)
            IQR = Q3 - Q1
            upper_bound = Q3 + (threshold * IQR) + 0.02
        else:  # if no ref_distances, this is an orthograph, so reject
            upper_bound = "N/A"
            IQR = "N/A"
        intermediate_list = []
        for candidate in candidates_at_index:
            mean_distance = "No refs"
            header = candidate.id
            raw_sequence = candidate.raw
            grade = "Fail"
            if has_ref_distances:
                candidate_distances = candidate_pairwise_calls(candidate, ref_alignments)
                mean_distance = np.nanmean(candidate_distances)

                if mean_distance <= upper_bound:
                    if sort == "original":
                        to_add_later.append((header,raw_sequence))
                    elif sort == "cluster":
                        intermediate_list.append((header,raw_sequence))
                    grade = "Pass"
            outliers.append((header, mean_distance, upper_bound, grade, IQR))
        if sort == "cluster":
            intermediate_list = taxa_sort(intermediate_list)
            to_add_later.extend(intermediate_list)

    return regulars, to_add_later, outliers


def delete_empty_columns(raw_fed_sequences: list, verbose: bool) -> tuple[list, list]:
    """
    Iterates over each sequence and deletes columns
    that consist of 100% dashes.
    """
    result = []
    sequences = [i[1] for i in raw_fed_sequences]

    positions_to_keep = []
    if sequences:
        for i in range(len(sequences[-1])):
            for sequence in sequences:
                if sequence[i] != "-":
                    positions_to_keep.append(i)
                    break
                
        for header, sequence in raw_fed_sequences:
            try:
                sequence = "".join([sequence[x] for x in positions_to_keep])
            except IndexError:
                printv(f"WARNING: Sequence length is not the same as other sequences: {header}", verbose, 0)
                continue

            result.append((header, sequence))

    return result, positions_to_keep

def align_col_removal(raw_fed_sequences: list, positions_to_keep: list) -> list:
    """
    Iterates over each sequence and deletes columns
    that were removed in the empty column removal.
    """
    result = []

    for header, sequence in raw_fed_sequences:
        sequence = "".join([sequence[i*3:(i*3)+3] for i in positions_to_keep])
        
        result.append((header, sequence))
    
    return result

def remove_excluded_sequences(records: list, excluded: set) -> list:
    """
    Given a list of fasta records and a set of headers to be excluded from output,
    returns a list of all valid headers and sequences. Use before the delete_column
    call in the nt portion.
    """
    output = []
    for header, sequence in records:
        if header not in excluded:
            output.append((header, sequence))
    return output


def main_process(
    args_input,
    nt_input,
    args_output,
    args_threshold,
    args_references,
    sort: str,
    nt_output_path: str,
    debug: bool,
    verbose: int,
    compress: bool,
):
    keep_refs = not args_references

    file_input = args_input
    filename = os.path.basename(file_input)

    printv(f"Doing: {filename}", verbose, 2)

    name = filename.split(".")[0]
    threshold = args_threshold / 100
    aa_output = os.path.join(args_output, "aa")
    aa_output = os.path.join(aa_output, filename.rstrip(".gz"))

    to_be_excluded = set()
    outliers_csv_path = os.path.join(args_output, "logs", "outliers_" + name + ".csv")
    reference_sequences, candidate_sequences = split_sequences(
        file_input, to_be_excluded
    )

    candidate_headers = [
        header for header, _ in candidate_sequences
    ]
    raw_regulars, to_add, outliers = compare_means(
        reference_sequences,
        candidate_sequences,
        threshold,
        to_be_excluded,
        keep_refs,
        sort,
    )
    if sort == "original":
        to_add = original_sort(candidate_headers, to_add)

    for line in to_add:
        raw_regulars.append(line)

    regulars, allowed_columns = delete_empty_columns(raw_regulars, verbose)

    if to_add:  # If candidate added to fasta
        writeFasta(aa_output, regulars, compress)

    to_be_excluded = set()
    if debug:
        with open(outliers_csv_path, "w", encoding="UTF-8") as outliers_csv:
            for outlier in outliers:
                header, distance, ref_dist, grade, iqr = outlier
                if grade == "Fail":
                    to_be_excluded.add(header)
                    header = header[1:]
                result = [header, str(distance), str(ref_dist), str(iqr), grade]
                outliers_csv.write(",".join(result) + "\n")
    else:
        for outlier in outliers:
            header, distance, ref_dist, grade, iqr = outlier
            if grade == "Fail":
                to_be_excluded.add(header)
    if to_add:
        nt_file = filename.replace(".aa.", ".nt.")
        nt_input_path = os.path.join(nt_input, nt_file)
        if not os.path.exists(nt_input_path):
            nt_input_path = os.path.join(nt_input, nt_file)
        if not os.path.exists(nt_output_path):
            os.mkdir(nt_output_path)
        nt_output_path = os.path.join(nt_output_path, nt_file.rstrip(".gz"))

        lines = parseFasta(nt_input_path)

        non_empty_lines = remove_excluded_sequences(
            lines, to_be_excluded
        )
        non_empty_lines = align_col_removal(non_empty_lines, allowed_columns)

        writeFasta(nt_output_path, non_empty_lines, compress)


def run_command(arg_tuple: tuple) -> None:
    main_process(*arg_tuple)


def do_folder(folder, args):
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    wanted_aa_path = Path(folder, "trimmed", "aa")
    if wanted_aa_path.exists():
        aa_input = wanted_aa_path
        nt_input = Path(folder, "trimmed", "nt")
    else:
        aa_input = Path(folder, "mafft")
        nt_input = Path(folder, "nt_aligned")

    printv(f"Processing: {os.path.basename(folder)}", args.verbose, 0)

    if not aa_input.exists():  # exit early
        printv(f"WARNING: Can't find aa folder for taxa {folder}: '{wanted_aa_path}'. Aborting", args.verbose, 0)
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
                    args.sort,
                    nt_output_path,
                    args.debug,
                    args.verbose,
                    args.compress
                )
            )

        with Pool(args.processes) as pool:
            pool.map(run_command, arguments, chunksize=1)
    else:
        for gene in file_inputs:
            main_process(
                gene,
                nt_input,
                output_path,
                args.threshold,
                args.no_references,
                args.sort,
                nt_output_path,
                args.debug,
                args.verbose,
                args.compress
            )
    if args.debug:
        log_folder_path = os.path.join(output_path, "logs")
        global_csv_path = os.path.join(log_folder_path, "outliers_global.csv")

        logs = [
            x
            for x in os.listdir(log_folder_path)
            if "outliers_" in x and "global" not in x
        ]
        with open(global_csv_path, "w", encoding="UTF-8") as global_csv:
            global_csv.write("Gene,Header,Mean_Dist,Ref_Mean,IQR\n")
            for log in logs:
                log_file_path = os.path.join(log_folder_path, log)
                with open(log_file_path, encoding="UTF-8") as log_f:
                    for line in log_f:
                        if line.strip().split(",")[-1] == "Fail":
                            if line[-1] != "\n":
                                line = f"{line}\n"
                            global_csv.write(line)

    printv(f"Done! Took {time_keeper.differential():.2f}s", args.verbose)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for folder in args.INPUT:
        do_folder(Path(folder), args)
    if len(args.INPUT) > 1 or not args.verbose :
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr OutlierCheck"
    )
