"""
Outlier Check

PyLint 8.99/10
"""

import argparse
import os
from copy import deepcopy
from statistics import mean
from multiprocessing.pool import Pool
from itertools import combinations
from time import time
import numpy as np
import phymmr_tools as bd

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


def taxa_sort(lines: list) -> list:
    """
    Iterates over a list of candidates and creates Records. Sorts by
    taxa, then makes a fasta list. Returns the list.
    """
    records = []
    for i in range(0, len(lines), 2):
        specimen = Record(lines[i], lines[i + 1])
        records.append(specimen)
    records.sort(key=lambda x: (x.id.split("|")[2], x.id.split("|")[3]))
    output = []
    for record in records:
        output.append(record.id)
        output.append(record.sequence)
    return output


def original_sort(headers, lines) -> list:
    """
    Returns candidate sequences to their original order.
    """
    output = []
    record_dict = {}
    for i in range(0, len(lines), 2):
        record_dict[lines[i]] = lines[i + 1]
    for header in headers:
        sequence = record_dict.get(header, False)
        if sequence:
            output.append(header)
            output.append(sequence)
    return output


def folder_check(path: str) -> None:

    aa_folder = os.path.join(path, "aa")
    nt_folder = os.path.join(path, "nt")
    logs_folder = os.path.join(path, "logs")
    os.makedirs(aa_folder, exist_ok=True)
    os.makedirs(nt_folder, exist_ok=True)
    os.makedirs(logs_folder, exist_ok=True)

def get_headers(lines: list) -> list:
    """
    Returns a list of every other line in the provided argument. Used to get
    header names from a list of sequences.
    """
    result = []
    for i in range(0, len(lines), 2):
        result.append(lines[i])
    return result


def split_sequences(lines: list, excluded: set) -> tuple:
    """
    Reads over a fasta record in the given list and returns a tuple of two smaller lists.
    The first returned list is the reference sequences found, the second returned list
    is the candidate sequences found.
    """
    bad_names = {"bombyx_mori", "danaus_plexippus"}
    references = []
    candidates = []

    end_of_references = False
    for i in range(0, len(lines), 2):
        header = lines[i].strip()
        sequence = lines[i + 1].strip()

        if end_of_references is False:
            # The reference header identifier is present in the header
            if header[-1] == '.':
                if header.split("|")[1].lower() in bad_names:
                    excluded.add(header)

                references.append(header)
                references.append(sequence)
            else:
                end_of_references = True

        if end_of_references is True:
            candidates.append(header)
            candidates.append(sequence)

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


def constrain_data_lines(lines: list, start: int, end: int) -> tuple:
    """
    Given a start and end value, iterates over the list of sequences and
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
    """
    Given a list of stings from a fasta file, returns a list of Sequence objects
    from the biopython module. This allows us to make a MultipleSequenceAlignment
    object later.
    """
    return [Record(lines[i], lines[i + 1]) for i in range(0, len(lines), 2)]


def find_index_groups(references: list, candidates: list) -> tuple:
    """
    Iterate over a list of candidate fastas as lines of text and finds their start
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
    sum = 0
    zeros_found = 0
    total_number = 0
    for row in matrix:
        for column in row:
            sum += column
            total_number += 1
            if column == 0:
                zeros_found += 1
    if ignore_zeros:
        total_number -= zeros_found
    mean = sum / total_number
    return mean


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
        for line in references:
            regulars.append(line)
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
        # First quartile (Q1)
        try:
            Q1 = np.percentile(ref_distances, 25, method="midpoint")
        except IndexError:
            Q1 = 0.0
        # Third quartile (Q3)
        try:
            Q3 = np.percentile(ref_distances, 75, method="midpoint")
        except IndexError:
            Q3 = 0.0
        # Interquartile range (IQR)
        IQR = Q3 - Q1
        upper_bound = Q3 + (threshold * IQR) + 0.02
        intermediate_list = []
        for candidate in candidates_at_index:
            candidate_distances = candidate_pairwise_calls(candidate, ref_alignments)
            mean_distance = mean(candidate_distances)
            header = candidate.id
            raw_sequence = candidate.raw
            grade = "Fail"
            if mean_distance <= upper_bound:
                if sort == "original":
                    to_add_later.append(header)
                    to_add_later.append(raw_sequence)
                elif sort == "cluster":
                    intermediate_list.append(header)
                    intermediate_list.append(raw_sequence)
                grade = "Pass"
            outliers.append((header, mean_distance, upper_bound, grade, IQR))
        if sort == "cluster":
            intermediate_list = taxa_sort(intermediate_list)
            to_add_later.extend(intermediate_list)
    return regulars, to_add_later, outliers


def make_nt_folder(path: str) -> str:  # FIXME: not used anywhere (but in a comment)
    head, tail = os.path.split(path)
    possible = os.listdir(head)
    result = None
    if tail == "mafft" and "nt_aligned" in possible:
        result = os.path.join(head, "nt_aligned")
    elif tail == "aa" and "nt" in possible:
        result = os.path.join(head, "nt")
    if result is None:
        raise ValueError("no valid nt folder found in input")
    return result


def make_nt_out_folder(output: str, path: str) -> str:  # FIXME: not used anywhere
    _, tail = os.path.split(path)
    return os.path.join(output, tail)


def deinterleave(fasta_lines: list) -> list:
    result = []
    this_out = []
    for line in fasta_lines:
        if line[0] == ">":
            if this_out:
                result.append("".join(this_out))
            result.append(line)
            this_out = []
        else:
            this_out.append(line.strip())
    if this_out:
        result.append("".join(this_out))
    return result


def delete_empty_columns(raw_fed_sequences: list) -> list:
    """
    Iterates over each sequence and deletes columns
    that consist of 100% dashes.
    """
    result = []
    sequences = []
    raw_sequences = [i.replace("\n", "") for i in raw_fed_sequences if i.replace("\n", "") != ""]

    for i in range(0, len(raw_sequences), 2):
        sequences.append(raw_sequences[i + 1])
    min_length = float("inf")
    for seq in sequences:
        min_length = min(len(seq), min_length)

    positions_to_keep = []
    if sequences:
        for i in range(len(sequences[0])):
            for sequence in sequences:
                if sequence[i] != "-":
                    positions_to_keep.append(i)
                    break
        for i in range(0, len(raw_sequences), 2):
            result.append(raw_sequences[i] + "\n")
            try:
                sequence = [raw_sequences[i + 1][x] for x in positions_to_keep]
            except IndexError:
                print(sequence)
            sequence.append("\n")
            sequence = "".join(sequence)

            result.append(sequence)

    return result


def remove_excluded_sequences(lines: list, excluded: set) -> list:
    """
    Given a list of fasta lines and a set of headers to be excluded from output,
    returns a list of all valid headers and sequences. Use before the delete_column
    call in the nt portion.
    """
    output = []
    for i in range(0, len(lines), 2):
        if lines[i].strip() not in excluded:
            output.append(lines[i])
            output.append(lines[i + 1])
    return output


def main_process(
    args_input,
    nt_input,
    args_output,
    args_threshold,
    args_references,
    sort: str,
    nt_output_path: str,
):

    keep_refs = not args_references

    file_input = args_input
    filename = os.path.basename(file_input)
    name = filename.split(".")[0]
    threshold = args_threshold / 100
    aa_output = os.path.join(args_output, "aa")
    aa_output = os.path.join(aa_output, filename)

    outliers_csv = os.path.join(args_output, "logs", "outliers_" + name + ".csv")
    with open(outliers_csv, "w+", encoding = "UTF-8") as outliers_csv:
        lines = []
        with open(file_input, encoding = "UTF-8") as fasta_in:
            lines = fasta_in.readlines()
            lines = deinterleave(lines)

        to_be_excluded = set()
        reference_sequences, candidate_sequences = split_sequences(lines, to_be_excluded)
        candidate_headers = [header for header in candidate_sequences if header[0] == ">"]
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

        regulars = delete_empty_columns(raw_regulars)

        if to_add:  # If candidate added to fasta
            with open(aa_output, "w+", encoding = "UTF-8") as aa_output:
                aa_output.writelines(regulars)

            to_be_excluded = set()
            for outlier in outliers:
                header, distance, ref_dist, grade, iqr = outlier
                if grade == "Fail":
                    to_be_excluded.add(header)

                header = header[1:]
                result = [header, str(distance), str(ref_dist), str(iqr), grade]
                outliers_csv.write(",".join(result) + "\n")

            nt_file = filename.replace(".aa.", ".nt.")
            nt_input_path = os.path.join(nt_input, nt_file)
            if not os.path.exists(nt_output_path):
                os.mkdir(nt_output_path)
            nt_output_path = os.path.join(nt_output_path, nt_file)

            with open(nt_output_path, "w+", encoding = "UTF-8") as nt_output_handle:
                with open(nt_input_path, encoding = "UTF-8") as nt_input_handle:
                    lines = nt_input_handle.readlines()
                    de_lines = deinterleave(lines)
                    non_empty_lines = remove_excluded_sequences(de_lines, to_be_excluded)
                    non_empty_lines = delete_empty_columns(non_empty_lines)

                for i in range(0, len(non_empty_lines), 2):
                    nt_output_handle.write(non_empty_lines[i])
                    nt_output_handle.write(non_empty_lines[i + 1])


def run_command(arg_tuple: tuple) -> None:
    input, nt_input, output, threshold, references_args, sort, nt = arg_tuple
    main_process(input, nt_input, output, threshold, references_args, sort, nt)


if __name__ == "__main__":
    start = time()
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", default="Taxa", help="Path to taxa")
    parser.add_argument("-o", "--output", default="outlier", help="Output folder")
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=0,
        help="Number of threads used to call processes.",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=50,
        help="Greater than reference mean to be counted as an outlier. Default is 2x.",
    )
    parser.add_argument(
        "--no-references",
        action="store_true",
        help="Disable output of reference sequences",
    )
    parser.add_argument(
        "-s",
        "--sort",
        choices=["cluster", "original"],
        default="original",
        help="Sort candidate output by cluster and taxa, or preserver original order.",
    )
    args = parser.parse_args()
    allowed_extensions = {"fa", "fas", "fasta"}

    for taxa in os.listdir(args.input):
        print(f"Doing taxa {taxa}")
        taxa_path = os.path.join(args.input, taxa)

        wanted_aa_path = os.path.join(taxa_path, "trimmed", "aa")
        if os.path.exists(wanted_aa_path):
            aa_input = wanted_aa_path
            nt_input = os.path.join(taxa_path, "trimmed", "nt")
        else:
            aa_input = os.path.join(taxa_path, "mafft")
            nt_input = os.path.join(taxa_path, "nt_aligned")

        if os.path.exists(aa_input):
            file_inputs = [
                os.path.join(aa_input, gene)
                for gene in os.listdir(os.path.join(aa_input))
                if ".aa" in gene and gene.split(".")[-1] in allowed_extensions
            ]
            output_path = os.path.join(taxa_path, args.output)
            nt_output_path = os.path.join(output_path, "nt")
            folder_check(output_path)
            # nt_folder = args.nt_input
            # if not nt_folder:
            #    nt_folder = make_nt_folder(args.aa_input)

            if args.processes:
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
                        )
                    )

                with Pool(args.processes) as pool:
                    pool.map(run_command, arguments, chunksize=1)
            else:
                for gene in file_inputs:
                    print(gene)
                    main_process(
                        gene,
                        nt_input,
                        output_path,
                        args.threshold,
                        args.no_references,
                        args.sort,
                        nt_output_path,
                    )

            log_folder_path = os.path.join(output_path, "logs")
            global_csv_path = os.path.join(log_folder_path, "outliers_global.csv")

            logs = [
                x
                for x in os.listdir(log_folder_path)
                if "outliers_" in x and "global" not in x
            ]
            with open(global_csv_path, "w", encoding = "UTF-8") as global_csv:
                global_csv.write("Gene,Header,Mean_Dist,Ref_Mean,IQR\n")
                for log in logs:
                    log_file_path = os.path.join(log_folder_path, log)
                    with open(log_file_path, encoding = "UTF-8") as log_f:
                        for line in log_f:
                            if line.strip().split(",")[-1] == "Fail":
                                global_csv.write(line)
                                if line[-1] != "\n":
                                    global_csv.write("\n")
            time_taken = time()
            time_taken = round(time_taken - start)

            print(f"Finished in {time_taken} seconds")

        else:
            print(f"Can't find aa folder for taxa {taxa}")
