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
import sys
import numpy as np
import phymmr_tools as bd

from Bio.Align.AlignInfo import SummaryInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from phymmr_tools import constrained_distance

from .utils import parseFasta, printv, write2Line2Fasta
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
    """
    Returns a list of every other line in the provided argument. Used to get
    header names from a list of sequences.
    """
    result = []
    for i in range(0, len(lines), 2):
        result.append(lines[i])
    return result


def split_sequences_ex(path: str, excluded: set) -> tuple:
    """
    Reads over a fasta record in the given list and returns a tuple of two smaller lists.
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
                    # fields = header.split('|')
                    # if header.split("|")[1].lower() in excluded: continue
                    if header.split("|")[1].lower() in excluded: continue
                    if header[-9] == ':':
                        ref_check.add(header[-9])

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
    """
    Reads over a fasta record in the given list and returns a tuple of two smaller lists.
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
                    # fields = header.split('|')
                    # if header.split("|")[1].lower() in excluded: continue
                    if header[-9] == ':':
                        ref_check.add(header[-9])

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


# function now skips any index which was previously rejected
# do the rejected indices make it to the distance calc? do I have to care at all?
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


def compare_means(
    references: list,
    ref_dict: dict,
    candidates_dict: dict,
    threshold: float,
    excluded_headers: set,
    keep_refs: bool,
    sort: str,
    refs_in_file: int,
    rejected_indices: set,
    index_group_min_bp: int,
    ref_gap_percent: float,
    ref_min_percent: float,
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
    to_add_later = []
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
                outliers.append(
                    (candidate.id, mean_distance, "N/A", "Fail", "min_cand_bp")
                )
            continue

        # first we have to calculate the reference distances to make the ref mean
        ref_records = convert_to_record_objects(current_refs)
        ref_alignments = [
            seq
            for seq in ref_records
            if seq.id not in excluded_headers
            and has_minimum_data(
                seq.sequence, rejected_indices, ref_gap_percent, index_pair[0]
            )
        ]

        ref_distances = []

        # find number of unique ref variants remaining after bp kick
        if ref_records[0].id[-9] == ':':
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
                    f'Q1 Runtime Error caused by references in {ref_alignments[0].id.split("|")[0]}'
                )
            except RuntimeError:
                Q1 = 0.0
                print(
                    f'Index Error caused by references in {ref_alignments[0].id.split("|")[0]}'
                )
            # Third quartile (Q3)
            try:
                Q3 = np.nanpercentile(ref_distances, 75, method="midpoint")
            except IndexError:
                Q3 = 0.0
                print(
                    f'Index Error caused by references in {ref_alignments[0].id.split("|")[0]}'
                )
            except RuntimeError:
                Q3 = 0.0
                print(
                    f'Runtime Error caused by references in {ref_alignments[0].id.split("|")[0]}'
                )
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
                candidate_distances = candidate_pairwise_calls(
                    candidate, ref_alignments
                )
                mean_distance = np.nanmean(candidate_distances)

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


def delete_empty_columns(raw_fed_sequences: list, verbose: bool) -> tuple[list, list]:
    """
    Iterates over each sequence and deletes columns
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
    """
    Iterates over each sequence and deletes columns
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


def make_exclusion_set(path: str) -> set:
    """
    Reads a file at a given path and returns a set containing
    each line. Used to make a taxa exclusion list.
    """
    excluded = set()
    if not path: return excluded
    with open(path) as f:
        for line in f:
            excluded.add(line.rstrip())
    return excluded


def alignment_from_2line(lines: list) -> list:
    """
    Given a list of header and sequences in 2-line fasta format,
    return list of BioPython SeqRecord objects.
    """
    result = []

    for i in range(0, len(lines), 2):
        header = lines[i]
        sequence = lines[i+1]
        result.append(
            SeqRecord(Seq(sequence.strip()),
                   id=header.strip()[1:])
        )
    
    return result

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
    col_cull_percent: float,
    index_group_min_bp: int,
    ref_gap_percent: float,
    ref_min_percent: int,
    internal_consensus_threshold: float,
    internal_kick_threshold: int,
    exclusion_file
):
    keep_refs = not args_references

    file_input = args_input
    filename = os.path.basename(file_input)

    printv(f"Doing: {filename}", verbose, 2)

    threshold = args_threshold / 100
    aa_output = os.path.join(args_output, "aa")
    aa_output = os.path.join(aa_output, filename.rstrip(".gz"))

    to_be_excluded = make_exclusion_set(exclusion_file)

    if to_be_excluded:
        reference_sequences, candidate_sequences, ref_check = split_sequences_ex(
        file_input, to_be_excluded
        )
    else:
        reference_sequences, candidate_sequences, ref_check = split_sequences(
            file_input
        )
    # refs_in_file = len(reference_sequences) / 2

    ref_dict, candidates_dict = find_index_groups(
        reference_sequences, candidate_sequences
    )

    # calculate indices that have valid data columns
    rejected_indices = set()
    ref_seqs = reference_sequences[1::2]
    for i in range(len(ref_seqs[0])):
        percent_of_non_dash = len([ref[i] for ref in ref_seqs if ref[i] != "-"]) / len(
            ref_seqs
        )
        if percent_of_non_dash <= col_cull_percent:
            rejected_indices.add(i)

    # find number of unique reference variants in file, use for refs_in_file
    if ref_check:
        refs_in_file = len(ref_check)
    else:
        refs_in_file = len(ref_seqs)
    candidate_headers = [header for header in candidate_sequences if header[0] == ">"]

    raw_regulars, to_add, outliers = compare_means(
        reference_sequences,
        ref_dict,
        candidates_dict,
        threshold,
        to_be_excluded,
        keep_refs,
        sort,
        refs_in_file,
        rejected_indices,
        index_group_min_bp,
        ref_gap_percent,
        ref_min_percent,
    )
    if sort == "original":
        to_add = original_sort(candidate_headers, to_add)
    to_be_excluded = set()
    logs = []
    internal_pass = []

    msr = alignment_from_2line(to_add)
    msa = MultipleSeqAlignment(msr)
    summary = SummaryInfo(msa)
    consensus = summary.dumb_consensus(threshold=internal_consensus_threshold)
    for seq in msa:
        header = seq.id
        distance = constrained_distance(consensus._data, seq.seq._data)
        if distance >= internal_kick_threshold:
            to_be_excluded.add(header)
            if debug:
                logs.append(f"{header},{distance},,,Internal Fail")
        else:
            internal_pass.append(f">{seq.id}\n")
            internal_pass.append(f"{seq.seq}\n")

    if internal_pass:
        for line in internal_pass:
            raw_regulars.append(line)

        regulars, allowed_columns = delete_empty_columns(raw_regulars, verbose)

        if to_add:  # If candidate added to fasta
            write2Line2Fasta(aa_output, regulars, compress)

        if debug:
            for outlier in outliers:
                header, distance, ref_dist, grade, iqr = outlier
                if grade == "Fail":
                    to_be_excluded.add(header)
                    header = header[1:]
                result = [header, str(distance), str(ref_dist), str(iqr), grade]
                logs.append(",".join(result) + "\n")
        else:
            for outlier in outliers:
                header, distance, ref_dist, grade, iqr = outlier
                if grade == "Fail":
                    to_be_excluded.add(header)
        if to_add:
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
        aa_input = Path(folder, "mafft")
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
                    args.sort,
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
                    args.exclude
                )
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
                    args.sort,
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
                    args.exclude
                )
            )
    if args.debug:
        log_folder_path = os.path.join(output_path, "logs")
        global_csv_path = os.path.join(log_folder_path, "outliers_global.csv")
        global_csv = open(global_csv_path, "w", encoding="UTF-8")
        global_csv.write("Gene,Header,Mean_Dist,Ref_Mean,IQR\n")

    for log_data in process_data:
        if args.debug:
            for line in log_data:
                if line.strip().split(",")[-1] == "Fail" or line.strip().split(",")[-1] == "Internal Fail":
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
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr OutlierCheck"
    )
