"""Outlier Check."""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from multiprocessing.pool import Pool
from os import mkdir, path
from pathlib import Path
from sys import exit

from msgspec import Struct
from numpy import float16, isnan, nanmean, nanpercentile
from phymmr_tools import (
    asm_index_split,
    blosum62_candidate_to_reference,
    blosum62_distance,
    delete_empty_columns,
    find_index_pair,
)
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, write2Line2Fasta


class Record(Struct):
    id: str
    raw: str
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

            if end_of_references is False:
                # The reference header identifier is present in the header
                if header[-1] == ".":
                    if header[-9] == ":":
                        ref_check.add(header[:-9])

                    references.append(Record(header, sequence))
                else:
                    end_of_references = True

            if end_of_references is True:
                candidates.append(Record(header, sequence))
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
        candidate_dict[(start, stop)].append(candidate)
    return candidate_dict


def find_asm_index_groups(candidates: list) -> dict:
    def lst():
        return []

    candidate_dict = defaultdict(lst)
    for candidate in candidates:
        indices = asm_index_split(candidate.raw)
        for index_pair in indices:
            start, stop = index_pair
            header = candidate.id + f"$${start}$${stop}"
            intron_candidate = Record(header, candidate.raw)
            intron_candidate.sequence = intron_candidate.raw[start:stop]
            candidate_dict[(start, stop)].append(intron_candidate)

    return candidate_dict


def remake_introns(passing: list) -> tuple:
    def lst():
        return []

    result = {}
    header_to_intron_records = defaultdict(lst)
    # first pass to find passing indices
    for record in passing:
        header, start, stop = record.id.split("$$")
        start = int(start)
        stop = int(stop)
        header_to_intron_records[header].append((start, stop))
    # second pass to reconstruct records with only valid indices
    for record in passing:
        header, start, stop = record.id.split("$$")
        if header in result:
            continue
        indices_list = header_to_intron_records[header]
        # make a single set of all sites to use
        valid_indices = set()
        for index_pair in indices_list:
            valid_indices.update(set(range(index_pair[0], index_pair[1])))
        seq = list(record.raw)
        for i, _ in enumerate(seq):
            if i not in valid_indices:
                seq[i] = "-"
        record.raw = "".join(seq)
        record.id = header
        result[header] = record
    return list(result.values()), header_to_intron_records


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
                blosum62_distance(ref1.sequence, ref2.sequence)
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
            upper_bound = Q3 + (threshold * IQR) + 0.02
        else:  # if no ref_distances, this is an orthograph, so reject
            upper_bound = "N/A"
            IQR = "N/A"
        for candidate in candidates_at_index:
            mean_distance = "No refs"
            candidate.grade = "Ref Fail"
            if has_ref_distances:
                candidate_distances = candidate_pairwise_calls(
                    candidate,
                    ref_alignments,
                )
                candidate.mean_distance = nanmean(candidate_distances)
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


def align_intron_removal(lines: list, header_to_indices: dict) -> list:
    """Replaces failing introns with hyphens.

    Iterates over the nt lines. For each two lines, pulls the header and
    checks the aa indices, translates those to nt indices, and then gaps
    the appropriate sites. Returns the altered list of nt lines. This should
    only be called when the assembly flag is found in the nt database.
    """

    def set_shim():
        return set()

    nt_indices = defaultdict(set_shim)
    # convert aa index pairs to nt indices
    for header, index_list in header_to_indices.items():
        for start, stop in index_list:
            index_set = sum([aa_to_nt(index) for index in range(start, stop)], [])
            nt_indices[header].update(index_set)
    # gap introns at any index not present in nt_indices
    for i in range(0, len(lines), 2):
        header = lines[i].strip()
        if header not in nt_indices:
            continue
        indices = nt_indices[header]
        sequence = list(lines[i + 1].strip())

        for j, _ in enumerate(sequence):
            if j not in indices:
                sequence[j] = "-"
        lines[i + 1] = "".join(sequence)
    return lines


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


def make_asm_exclusions(passing: list, failing: list) -> set:
    """Makes to_be_excluded set for assembly files.

    Finds failing - passing. These names need to be
    excluded from the nt files. This function is required because
    the intron location has been appended to the record's id, so
    it's possible to have a single name in both passing and failing.
    """
    failing_set = {record.id.split("$$")[0] for record in failing}
    passing_set = {record.id.split("$$")[0] for record in passing}

    return failing_set.difference(passing_set)


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
    assembly: bool,
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
    if not assembly:
        candidates_dict = find_index_groups(candidate_records)
    else:
        candidates_dict = find_asm_index_groups(candidate_records)
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
    # reference_records = [ref for ref in reference_records if has_data(ref.raw)]
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
    )
    logs = []
    if passing:
        if assembly:
            passing, header_to_indices = remake_introns(passing)
    passing = original_order_sort(original_order, passing)
    for candidate in passing:
        raw_regulars.extend([candidate.id, candidate.raw])
    # regulars, allowed_columns = delete_empty_columns(raw_regulars, verbose)
    regulars, allowed_columns = delete_empty_columns(raw_regulars)
    if assembly:
        to_be_excluded = make_asm_exclusions(passing, failing)
    else:
        to_be_excluded = {candidate.id for candidate in failing}

    if passing:  # If candidate added to fasta
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
        if assembly:
            non_empty_lines = align_intron_removal(non_empty_lines, header_to_indices)
        non_empty_lines = align_col_removal(non_empty_lines, allowed_columns)

        write2Line2Fasta(nt_output_path, non_empty_lines, compress)
    return logs


def do_folder(folder, args):
    ALLOWED_EXTENSIONS = {".fa", ".fas", ".fasta", ".fa", ".gz", ".fq", ".fastq"}

    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    wanted_aa_path = Path(folder, "trimmed", "aa")
    if wanted_aa_path.exists():
        aa_input = wanted_aa_path
        nt_input = Path(folder, "trimmed", "nt")
    else:
        aa_input = Path(folder, "align")
        nt_input = Path(folder, "nt_aligned")
    rocks_db_path = Path(folder, "rocksdb", "sequences", "nt")
    if rocks_db_path.exists():
        rocksdb_db = RocksDB(str(rocks_db_path))
        assembly = rocksdb_db.get("get:isassembly")
        assembly = assembly == "True"
        del rocksdb_db
    else:
        err = f"cannot find dupe databases for {folder}"
        raise FileNotFoundError(err)

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
                    args.col_cull_percent,
                    args.index_group_min_bp,
                    args.ref_gap_percent,
                    args.ref_min_percent,
                    # args.internal_consensus_threshold,
                    # args.internal_kick_threshold,
                    assembly,
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
                    args.col_cull_percent,
                    args.index_group_min_bp,
                    args.ref_gap_percent,
                    args.ref_min_percent,
                    # args.internal_consensus_threshold,
                    # args.internal_kick_threshold,
                    assembly,
                ),
            )
    if args.debug:
        log_folder_path = path.join(output_path, "logs")
        global_csv_path = path.join(log_folder_path, "outliers_global.csv")
        with open(global_csv_path, "w", encoding="UTF-8") as global_csv:
            global_csv.write("Gene,Header,Mean_Dist,Ref_Mean,IQR\n")

            for log_data in process_data:
                for line in log_data:
                    if line.split(",")[-2] != "Pass":
                        if line[-1] != "\n":
                            line = f"{line}\n"
                        global_csv.write(line)

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose)
    return True, assembly


def main(args):
    is_assembly = None
    success = False
    if isinstance(args.INPUT, list):
        success = all([do_folder(Path(folder), args)[0] for folder in args.INPUT])
    elif isinstance(args.INPUT, str):
        success, is_assembly = do_folder(Path(args.INPUT), args)
    return success, is_assembly


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise RuntimeError(
        MSG,
    )
