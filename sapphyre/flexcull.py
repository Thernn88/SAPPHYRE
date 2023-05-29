"""FlexCull Description Goes Here.

PyLint 9.81/10
"""
from __future__ import annotations

import os
from collections import Counter, namedtuple
from itertools import chain
from multiprocessing.pool import Pool
from phymmr_tools import (
    join_by_tripled_index,
    join_with_exclusions,
    join_triplets_with_exclusions,
)
from shutil import rmtree
import wrap_rocks
import numpy as np

import blosum as bl

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

MainArgs = namedtuple(
    "MainArgs",
    [
        "verbose",
        "processes",
        "debug",
        "INPUT",
        "output",
        "amino_acid",
        "nucleotide",
        "matches",
        "base_pair",
        "compress",
        "gap_threshold",
        "mismatches",
        "column_cull",
        "blosum_strictness",
    ],
)

FlexcullArgs = namedtuple(
    "FlexcullArgs",
    [
        "aa_input",
        "nt_input",
        "output",
        "amt_matches",
        "aa_file",
        "debug",
        "bp",
        "verbosity",
        "compress",
        "gap_threshold",
        "mismatches",
        "column_cull_percent",
        "filtered_mat",
        "is_assembly",
    ],
)

def align_col_removal(raw_fed_sequences: list, positions_to_keep: list) -> list:
    """Iterates over each sequence and deletes columns
    that were removed in the empty column removal.

    Args:
        raw_fed_sequences (list): List of tuples containing header and sequence
        positions_to_keep (list): List of positions to keep
    Returns:
        list: List of tuples containing header and sequence with empty columns removed
    """
    result = []
    for raw_sequence in raw_fed_sequences:
        sequence = raw_sequence[1]
        # sequence = ''.join(sequence[i*3:(i*3)+3] for i in positions_to_keep)
        sequence = join_by_tripled_index(sequence, positions_to_keep)
        result.append((raw_sequence[0], sequence))

    return result


def delete_empty_columns(raw_fed_sequences: list, verbose: bool) -> tuple[list, list]:
    """Iterates over each sequence and deletes columns
    that consist of 100% dashes. In this version, raw_feed_sequences
    is a list of tuples.

    Args:
        raw_fed_sequences (list): List of tuples containing header and sequence
        verbose (bool): Whether to print verbose output
    Returns:
        tuple[list, list]: List of tuples containing header and sequence with empty columns removed and a list of positions to keep
    """
    result = []
    positions_to_keep = []
    sequences = [x[1] for x in raw_fed_sequences]
    if sequences:
        sequence_length = len(sequences[0])
        for i in range(sequence_length):
            if any(sequence[i] != "-" for sequence in sequences):
                positions_to_keep.append(i)

        for raw_sequence in raw_fed_sequences:
            try:
                sequence = ''.join(raw_sequence[1][x] for x in positions_to_keep)
            except IndexError:
                if verbose:
                    print(f"WARNING: Sequence length is not the same as other sequences: {raw_sequence[0]}")
                continue

            result.append((raw_sequence[0], sequence))

    return result, positions_to_keep



def folder_check(output_target_path: str) -> None:
    """Checks to see if input and output directory has the necessary
    folder structure.

    if not, create missing folders.

    Args:
        output_target_path (str): Path to output directory
    Returns:
        None
    """
    output_aa_path = os.path.join(output_target_path, "aa")
    output_nt_path = os.path.join(output_target_path, "nt")
    rmtree(output_aa_path, ignore_errors=True)
    rmtree(output_nt_path, ignore_errors=True)
    os.makedirs(output_aa_path, exist_ok=True)
    os.makedirs(output_nt_path, exist_ok=True)


def make_nt(aa_file_name: str) -> str:
    """Converts AA file name to NT file name.
    
    Args:
        aa_file_name (str): AA file name
    Returns:
        str: NT file name
    """
    return aa_file_name.replace(".aa.", ".nt.")


def parse_fasta(fasta_path: str) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """Parses fasta file into header and sequences.
    
    Args:
        fasta_path (str): Path to fasta file
    Returns:
        tuple[list[tuple[str, str]], list[tuple[str, str]]]: Tuple containing list of reference record tuples and candidate record tuples
    """
    references = []
    candidates = []
    end_of_references = False

    for header, sequence in parseFasta(fasta_path):
        
        if end_of_references:
            candidates.append((header, sequence))
            continue

        if header[-1] == ".":  # Is reference
            references.append((header, sequence))
            continue
        
        candidates.append((header, sequence))
        end_of_references = True

    return references, candidates


def trim_around(
    starting_index: int,
    lower_limit: int,
    upper_limit: int,
    sequence: str,
    amt_matches: int,
    mismatches: int,
    all_dashes_by_index: list[bool],
    character_at_each_pos: list[set[str]],
    gap_present_threshold: list[bool],
) -> range:
    """Trim from starting index to the upper limit and from the starting index to the
    lower limit and return the range of the trimmed slice.
    
    Trim will stop once the amount of matches is reached or the upper limit is reached.
    
    Args:
        starting_index (int): Starting index
        lower_limit (int): Lower limit
        upper_limit (int): Upper limit
        sequence (str): Sequence
        amt_matches (int): Amount of consecutive matches required to pass
        mismatches (int): Amount of mismatches allowed in consecutive matches
        all_dashes_by_index (list[bool]): List of booleans indicating whether all references have a dash at the index
        character_at_each_pos (list[set[str]]): List of sets containing characters at each position
        gap_present_threshold (list[bool]): List of booleans indicating whether references contain gaps at the index
    Returns:
        range: Range of the trimmed slice"""
    
    offset = amt_matches - 1
    cull_end = upper_limit
    
    for i in range(starting_index, len(sequence) - 1):
        skip_first = 0
        char = sequence[i]
        mismatch = mismatches

        if i == upper_limit - offset:
            break

        if char == "-" or all_dashes_by_index[i]:
            continue

        if char not in character_at_each_pos[i]:
            skip_first = 1
            mismatch -= 1

        if mismatch < 0:
            continue

        pass_all = True
        checks = amt_matches - 1
        match_i = 1

        while checks > 0:
            if i + match_i >= upper_limit:
                pass_all = False
                break

            next_char = sequence[i + match_i]

            if next_char == "-":
                if gap_present_threshold[i + match_i]:
                    pass_all = False
                    break
                match_i += 1
            elif next_char not in character_at_each_pos[i + match_i]:
                mismatch -= 1
                if mismatch < 0:
                    pass_all = False
                    break
                match_i += 1
                checks -= 1
            else:
                match_i += 1
                checks -= 1

        if pass_all:
            cull_end = i + skip_first
            break

    cull_start = starting_index
    
    for i in range(starting_index - 1, -1, -1):
        mismatch = mismatches
        skip_last = 0

        char = sequence[i]
        
        if i < lower_limit + offset:
            break
            
        if char == "-" or all_dashes_by_index[i]:
            continue
            
        if char not in character_at_each_pos[i]:
            skip_last += 1
            mismatch -= 1

        if mismatch < 0:
            continue

        pass_all = True
        checks = amt_matches - 1
        match_i = 1

        while checks > 0:
            if i - match_i < lower_limit:
                pass_all = False
                break

            prev_char = sequence[i - match_i]

            if prev_char == "-":
                if gap_present_threshold[i - match_i]:
                    pass_all = False
                    break
                match_i += 1
            elif prev_char not in character_at_each_pos[i - match_i]:
                mismatch -= 1
                if mismatch < 0:
                    pass_all = False
                    break
                match_i += 1
                checks -= 1
            else:
                match_i += 1
                checks -= 1

        if pass_all:
            cull_start = i - skip_last + 1  # Inclusive
            break

    return range(cull_start, cull_end)


def get_data_difference(trim: int, ref: int) -> float:
    """
    Returns the difference between the trim and ref.
    """
    if trim == 0:
        return 0.0
    if ref == 0 and trim != 0:
        return 1.0
    if trim == 0 and ref == 0:
        return 1.0
    return trim / ref


def do_cull(
    sequence: str,
    sequence_length: int,
    offset: int,
    amt_matches: int,
    mismatches: int,
    all_dashes_by_index: list[bool],
    character_at_each_pos: list[set[str]],
    gap_present_threshold: list[bool],
) -> tuple[int, int, bool]:
    """
    Culls leading and trailing bp from a sequence based on if the bp is present in a reference sequence.

    Args:
        sequence (str): Sequence to cull
        sequence_length (int): Length of the sequence
        offset (int): Offset to start culling from
        amt_matches (int): Amount of consecutive matches required to pass
        mismatches (int): Amount of mismatches allowed in consecutive matches
        all_dashes_by_index (list[bool]): List of booleans indicating whether all references have a dash at the index
        character_at_each_pos (list[set[str]]): List of sets containing characters at each position
        gap_present_threshold (list[bool]): List of booleans indicating whether references contain gaps at the index
    Returns:
        tuple: Tuple containing the index of the first bp that has required consecutive matches,
        the index of the last bp that has required consecutive matches, and a boolean indicating whether the trim kicked the entire sequence
    """
    cull_start = None
    cull_end = None
    kick = False
    for i, char in enumerate(sequence):
        mismatch = mismatches
        skip_first = 0

        if i == sequence_length - offset:
            kick = True
            break

        # Don't allow cull to point of all dashes
        if char == "-":
            continue

        window_start = i

        if all_dashes_by_index[i]:
            continue
        if char not in character_at_each_pos[i]:
            skip_first = 1
            if window_start == i:
                continue
            mismatch -= 1

        if mismatch < 0:
            continue

        pass_all = True
        checks = amt_matches - 1
        match_i = 1

        while checks > 0:
            if i + match_i >= len(sequence):
                pass_all = False
                break

            if sequence[i + match_i] == "-":
                if gap_present_threshold[i + match_i]:
                    pass_all = False
                    break
                match_i += 1
            elif sequence[i + match_i] not in character_at_each_pos[i + match_i]:
                mismatch -= 1
                if mismatch < 0 or sequence[i + match_i] == "*":
                    pass_all = False
                    break
                match_i += 1
                checks -= 1
            else:
                match_i += 1
                checks -= 1

        if pass_all:
            cull_start = i + skip_first

            break

    if not kick:
        # If not kicked from Cull Start Calc. Continue
        window_end = -1
        for i in range(len(sequence) - 1, -1, -1):
            mismatch = mismatches
            skip_last = 0

            char = sequence[i]
            if i < cull_start + offset:
                kick = True
                break
            if char == "-":
                continue
            window_end = i

            if all_dashes_by_index[i]:
                # Don't allow cull to point of all dashes
                continue
            if char not in character_at_each_pos[i]:
                skip_last += 1
                if window_end == i:
                    continue
                mismatch -= 1

            if mismatch < 0:
                continue

            pass_all = True
            checks = amt_matches - 1
            match_i = 1
            while checks > 0:
                if i - match_i < 0:
                    pass_all = False
                    break

                if sequence[i - match_i] == "-":
                    if gap_present_threshold[i - match_i]:
                        pass_all = False
                        break
                    match_i += 1
                elif sequence[i - match_i] not in character_at_each_pos[i - match_i]:
                    mismatch -= 1
                    if mismatch < 0 or sequence[i - match_i] == "*":
                        pass_all = False
                        break
                    match_i += 1
                    checks -= 1
                else:
                    match_i += 1
                    checks -= 1

            if pass_all:
                cull_end = i - skip_last + 1  # Inclusive
                break

    return cull_start, cull_end, kick



def get_start_end(sequence: str) -> tuple:
    """Returns the start and end of the sequence."""
    start = next((i for i, char in enumerate(sequence) if char != "-"), 0)
    end = next((len(sequence) - i for i, char in enumerate(sequence[::-1]) if char != "-"), len(sequence))
    return start, end


def process_refs(references: list[tuple], gap_threshold:float, column_cull_percent: float, mat: dict) -> tuple:
    """
    Processes the references to generate the character_at_each_pos, gap_present_threshold , column_cull and all_dashes_by_index.

    character_at_each_pos: set of all characters aswell as the positive blosum subsitutions in each column of the references
    gap_present_threshold: boolean indicating whether the amount of data in a column is above the threshold
    column_cull: set of columns to cull
    all_dashes_by_index: boolean indicating whether all references have a dash at the index

    Args:
        references (list[tuple]): List of tuples containing the reference name and sequence
        gap_threshold (float): Threshold for the amount of gaps required in a reference for it to be a reference gap column
        column_cull_percent (float): Percent of data required in a column to allow it to pass
        mat (dict): The blosum62 subsitution matrix
    Returns:
        tuple: Tuple containing the character_at_each_pos, gap_present_threshold , column_cull and all_dashes_by_index tables

    """
    character_at_each_pos = {}
    gap_present_threshold = {}
    all_dashes_by_index = {}
    column_cull = set()

    max_ref_length = max(len(sequence) for _, sequence in references)
    all_dashes_by_index = {i: True for i in range(max_ref_length)}

    for _, sequence in references:
        for i, char in enumerate(sequence.replace("*", "-")):
            character_at_each_pos.setdefault(i, []).append(char)

    for i, chars in list(character_at_each_pos.items()):
        data_present = 1 - (chars.count("-") / len(chars))
        all_dashes_by_index[i] = data_present == 0
        gap_present_threshold[i] = data_present >= gap_threshold
        if data_present < column_cull_percent:
            column_cull.add(i * 3)

        character_at_each_pos[i] = set(chars).union(blosum_sub for char in chars for blosum_sub in mat.get(char, {}).keys())

    return character_at_each_pos, gap_present_threshold, all_dashes_by_index, column_cull



def column_cull_seqs(this_seqs: list[tuple], column_cull: set, minimum_bp: int, is_assembly: bool) -> tuple[list[tuple], list[str], set]:
    aa_out = []
    kicks = []

    this_column_cull = set()
    for nt_i in column_cull:
        i = nt_i // 3
        this_sequence = Counter()
        internal_count = 0
        for _, seq, start, end in this_seqs:
            if start <= i <= end:
                internal_count += 1
                if seq[i] != "-":
                    this_sequence[seq[i]] += 1
        # Get the count of the most occuring bp at this position
        if sum(this_sequence.values()) >= 5:
            if this_sequence.most_common()[0][1] > internal_count / 3:
                continue
        this_column_cull.add(i * 3)

    if this_column_cull:
        for header, seq, _, _ in this_seqs:
            # new_seq = "".join(
            #     [
            #         let if i * 3 not in this_column_cull else "-"
            #         for i, let in enumerate(seq)
            #     ],
            # )
            new_seq = join_with_exclusions(seq, this_column_cull)
            if len(new_seq) - new_seq.count("-") < minimum_bp:
                kicks.append(header)
                continue
            aa_out.append((header, new_seq))
    else:
        aa_out.extend([(header, sequence) for header, sequence, _, _ in this_seqs])

    return aa_out, kicks, this_column_cull


def cull_codons(out_line: str, cull_start: int, cull_end: int, amt_matches:int, mismatches: int, all_dashes_by_index: dict, character_at_each_pos: dict, gap_present_threshold: dict):
    positions_to_trim = set()
    codons = []
    for i in range(cull_start, cull_end):
        char = out_line[i]
        if char == "*":
            codons.append(i)

    non_trimmed_codons = [c for c in codons if c * 3 not in positions_to_trim]
    while non_trimmed_codons:
        i = non_trimmed_codons[int(len(non_trimmed_codons) / 2)]
        positions = trim_around(
            i,
            cull_start,
            cull_end,
            out_line,
            amt_matches,
            mismatches,
            all_dashes_by_index,
            character_at_each_pos,
            gap_present_threshold,
        )
        for x in positions:
            positions_to_trim.add(x * 3)
            out_line[x] = "-"

        left_after = out_line[cull_start:i]
        right_after = out_line[i:cull_end]

        left_side_ref_data_columns = sum(
            gap_present_threshold[x] for x in range(cull_start, i)
        )
        left_of_trim_data_columns = len(left_after) - left_after.count("-")

        right_side_ref_data_columns = sum(
            gap_present_threshold[x] for x in range(i, cull_end)
        )
        right_of_trim_data_columns = len(right_after) - right_after.count("-")

        if (
            get_data_difference(
                left_of_trim_data_columns, left_side_ref_data_columns,
            )
            < 0.55
        ):  # candidate has less than % of data columns compared to reference
            for x in range(cull_start, i):
                positions_to_trim.add(x * 3)
                out_line[x] = "-"
        if (
            get_data_difference(
                right_of_trim_data_columns, right_side_ref_data_columns,
            )
            < 0.55
        ):
            for x in range(i, cull_end):
                positions_to_trim.add(x * 3)
                out_line[x] = "-"

        non_trimmed_codons.remove(i)

    return out_line, positions_to_trim


def trim_large_gaps(aa_out: list[tuple], reference_gap_col: set, amt_matches: int, mismatches: int, post_all_dashes_by_index: dict, post_character_at_each_pos: dict, post_gap_present_threshold: dict, minimum_bp: int, debug:bool):
    """
    Trims gaps with a length of 80 or more non-reference gap column gaps.

    Args:
        aa_out (list[tuple]): List of tuples containing header and sequence
        reference_gap_col (set): Set of reference gap columns
        amt_matches (int): Amount of consecutive matches required to pass
        mismatches (int): Amount of mismatches allowed in consecutive matches
        post_all_dashes_by_index (dict): Dictionary containing boolean indicating whether all references have a dash at the index
        post_character_at_each_pos (dict): Dictionary containing set of characters at each position
        post_gap_present_threshold (dict): Dictionary containing boolean indicating whether references contain gaps at the index
        minimum_bp (int): Minimum amount of bp required to pass
        debug (bool): Whether to print debug output
    Returns:
        tuple: Tuple containing the updated aa_out, gap_pass_through to align the deletions to NT, the debug log, and kicks
    """
    log = []
    kicks = []
    gap_pass_through = {}
    for record_index, record in enumerate(aa_out):
        header, sequence = record
        if not header.endswith("."):
            gap_cull = set()
            seq_start, seq_end = get_start_end(sequence)
            change_made = False
            non_ref_gap_dash_count = 0
            raw_dash_count = 0
            out_line = list(sequence)
            for j, let in enumerate(out_line[seq_start : seq_end + 1], seq_start):
                if let == "-":
                    if j not in reference_gap_col:
                        non_ref_gap_dash_count += 1
                    raw_dash_count += 1
                else:
                    if non_ref_gap_dash_count >= 80:
                        i = j - (raw_dash_count // 2)
                        positions = trim_around(
                            i,
                            seq_start,
                            seq_end,
                            out_line,
                            amt_matches,
                            mismatches,
                            post_all_dashes_by_index,
                            post_character_at_each_pos,
                            post_gap_present_threshold,
                        )
                        for x in positions:
                            gap_cull.add(x * 3)
                            out_line[x] = "-"

                        left_after = out_line[seq_start:i]
                        right_after = out_line[i:seq_end]
                        left_side_ref_data_columns = sum(
                            post_gap_present_threshold[x] for x in range(seq_start, i)
                        )
                        left_of_trim_data_columns = len(
                            left_after,
                        ) - left_after.count("-")

                        right_side_ref_data_columns = sum(
                            post_gap_present_threshold[x] for x in range(i, seq_end)
                        )
                        right_of_trim_data_columns = len(
                            right_after,
                        ) - right_after.count("-")

                        # If both sides kicked and sequence ends up being empty keep the side with the most bp.
                        keep_left = False
                        keep_right = False

                        if (
                            get_data_difference(
                                left_of_trim_data_columns,
                                left_side_ref_data_columns,
                            )
                            < 0.55
                            and get_data_difference(
                                right_of_trim_data_columns,
                                right_side_ref_data_columns,
                            )
                            < 0.55
                        ):
                            keep_left = len(left_after) - left_after.count(
                                "-",
                            ) >= len(right_after) - right_after.count("-")
                            keep_right = not keep_left

                        if (
                            get_data_difference(
                                left_of_trim_data_columns,
                                left_side_ref_data_columns,
                            )
                            < 0.55
                            and not keep_left
                        ):  # candidate has less than % of data columns compared to reference
                            for x in range(seq_start, i):
                                gap_cull.add(x * 3)
                                out_line[x] = "-"
                        if (
                            get_data_difference(
                                right_of_trim_data_columns,
                                right_side_ref_data_columns,
                            )
                            < 0.55
                            and not keep_right
                        ):
                            for x in range(i, seq_end):
                                gap_cull.add(x * 3)
                                out_line[x] = "-"
                        change_made = True

                    non_ref_gap_dash_count = 0
                    raw_dash_count = 0
            if change_made:
                bp_after_cull = len(out_line) - out_line.count("-")
                if bp_after_cull < minimum_bp:
                    if debug:
                        removed_section = sequence[:seq_start] + sequence[seq_end:]
                        data_removed = len(removed_section) - removed_section.count(
                            "-",
                        )
                        log.append(
                            header.split("|")[0]
                            + ","
                            + header
                            + ",Not enough BP after gap cull,,"
                            + str(bp_after_cull)
                            + ","
                            + str(data_removed)
                            + "\n",
                        )
                    kicks.append(header)
                    aa_out[record_index] = None
                else:
                    aa_out[record_index] = (header, "".join(out_line))

                gap_pass_through[header] = gap_cull

    aa_out = [i for i in aa_out if i is not None]
    return aa_out, gap_pass_through, log, kicks

def do_gene(fargs: FlexcullArgs) -> None:
    """FlexCull main function. Culls input aa and nt using specified amount of matches."""
    gene_path = os.path.join(fargs.aa_input, fargs.aa_file)
    this_gene = fargs.aa_file.split(".")[0]

    printv(f"Doing: {this_gene}", fargs.verbosity, 2)

    references, candidates = parse_fasta(gene_path)

    character_at_each_pos, gap_present_threshold, all_dashes_by_index, column_cull = process_refs(
        references, fargs.gap_threshold, fargs.column_cull_percent, fargs.filtered_mat
    )

    log = []

    follow_through = {}
    offset = fargs.amt_matches - 1

    aa_out_path = os.path.join(fargs.output, "aa", fargs.aa_file.rstrip(".gz"))
    aa_out = references.copy()
    this_seqs = []

    for header, sequence in candidates:
        sequence = list(sequence)

        gene = header.split("|")[0]

        kick = False
        data_removed = 0
        data_length = 0

        sequence_length = len(sequence)

        cull_start, cull_end, kick = do_cull(
            sequence,
            sequence_length,
            offset,
            fargs.amt_matches,
            fargs.mismatches,
            all_dashes_by_index,
            character_at_each_pos,
            gap_present_threshold,
        )

        if not kick:  # If also passed Cull End Calc. Finish
            out_line = ["-"] * cull_start + sequence[cull_start:cull_end]

            characters_till_end = sequence_length - len(out_line)
            out_line += ["-"] * characters_till_end

            out_line, positions_to_trim = cull_codons(out_line, cull_start, cull_end, fargs.amt_matches, fargs.mismatches, all_dashes_by_index, character_at_each_pos, gap_present_threshold)

            out_line = "".join(out_line)


            data_length = cull_end - cull_start
            bp_after_cull = len(out_line) - out_line.count("-")

            if bp_after_cull >= fargs.bp:
                follow_through[header] = (
                    False,
                    cull_start,
                    cull_end,
                    positions_to_trim,
                )

                this_seqs.append((header, out_line, *get_start_end(out_line)))

                if fargs.debug:
                    removed_section = sequence[:cull_start] + sequence[cull_end:]
                    data_removed = len(removed_section) - removed_section.count("-")
                    log.append(
                        gene
                        + ","
                        + header
                        + ","
                        + str(cull_start)
                        + ","
                        + str(cull_end)
                        + ","
                        + str(data_length)
                        + ","
                        + str(data_removed)
                        + "\n",
                    )
            else:
                follow_through[header] = True, 0, 0, []
                if fargs.debug:
                    log.append(
                        gene
                        + ","
                        + header
                        + ",Kicked,Minimum BP Not Met,"
                        + str(bp_after_cull)
                        + ",\n",
                    )
        if kick:
            follow_through[header] = True, 0, 0, []

            if fargs.debug:
                log.append(gene + "," + header + ",Kicked,Zero Data After Cull,0,\n")

    if this_seqs:
        this_seqs, kicks, this_column_cull = column_cull_seqs(this_seqs, column_cull, fargs.gap_threshold, fargs.is_assembly)
        for header in kicks:
            follow_through[header] = True, 0, 0, []

        # remove empty columns from refs and candidates
        aa_out, aa_positions_to_keep = delete_empty_columns(aa_out + this_seqs, False)
        if len(aa_out) == len(references):
            return log  # Only refs

        # Recalcuate position based tables
        post_references = [i for i in aa_out if i[0].endswith(".")]

        post_character_at_each_pos, post_gap_present_threshold, post_all_dashes_by_index, _ = process_refs(
            post_references, fargs.gap_threshold, fargs.column_cull_percent, fargs.filtered_mat
        )

        reference_gap_col = {i for i, x in post_gap_present_threshold.items() if not x}
        
        aa_out, gap_pass_through, trim_log, kicks = trim_large_gaps(aa_out, reference_gap_col, fargs.amt_matches, fargs.mismatches, post_all_dashes_by_index, post_character_at_each_pos, post_gap_present_threshold, fargs.bp, fargs.debug)

        if fargs.debug:
            log.extend(trim_log)
        
        for kick in kicks:
            follow_through[kick] = True, 0, 0, []

        if len(aa_out) != len(references):
            writeFasta(aa_out_path, aa_out, fargs.compress)

            nt_file_name = make_nt(fargs.aa_file)
            gene_path = os.path.join(fargs.nt_input, nt_file_name)

            references, candidates = parse_fasta(gene_path)

            nt_out_path = os.path.join(fargs.output, "nt", nt_file_name.rstrip(".gz"))
            nt_out = references.copy()
            for header, sequence in candidates:
                gene = header.split("|")[0]
                kick, cull_start, cull_end, positions_to_trim = follow_through[header]

                if not kick:
                    cull_start_adjusted = cull_start * 3
                    cull_end_adjusted = cull_end * 3

                    out_line = ("-" * cull_start_adjusted) + sequence[
                        cull_start_adjusted:cull_end_adjusted
                    ]

                    characters_till_end = len(sequence) - len(out_line)
                    out_line += (
                        "-" * characters_till_end
                    )  # Add dashes till reached input distance

                    # out_line = [
                    #     out_line[i : i + 3]
                    #     if i not in positions_to_trim and i not in this_column_cull
                    #     else "---"
                    #     for i in range(0, len(out_line), 3)
                    # ]
                    # out_line = "".join(out_line)
                    out_line = join_triplets_with_exclusions(
                        out_line, positions_to_trim, this_column_cull
                    )
                    nt_out.append((header, out_line))
            nt_out = align_col_removal(nt_out, aa_positions_to_keep)
            out_nt = []
            for header, sequence in nt_out:
                gap_cull = gap_pass_through.get(header, None)
                if gap_cull:
                    out_nt.append(
                        (
                            header,
                            "".join(
                                [
                                    sequence[i : i + 3] if i not in gap_cull else "---"
                                    for i in range(0, len(sequence), 3)
                                ],
                            ),
                        ),
                    )
                else:
                    out_nt.append((header, sequence))

            writeFasta(nt_out_path, out_nt, fargs.compress)

    return log


def do_folder(folder, args: MainArgs):
    folder_time = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Processing: {folder}", args.verbose, 0)
    aa_path = os.path.join(folder, args.amino_acid)
    nt_path = os.path.join(folder, args.nucleotide)
    output_path = os.path.join(folder, args.output)
    if not os.path.exists(aa_path) or not os.path.exists(nt_path):
        printv(
            f"WARNING: Can't find aa ({aa_path}) and nt ({nt_path}) folders. Abort",
            args.verbose,
        )
        return

    nt_db_path = os.path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = wrap_rocks.RocksDB(nt_db_path)
    dbis_assembly = nt_db.get("get:isassembly")
    is_assembly = False
    if dbis_assembly and dbis_assembly == "True":
        is_assembly = True

    folder_check(output_path)
    file_inputs = [
        input_gene
        for input_gene in os.listdir(aa_path)
        if input_gene.split(".")[-1] in {"fa", "gz", "fq", "fastq", "fasta"}
    ]
    file_inputs.sort(
        key=lambda x: os.path.getsize(os.path.join(aa_path, x)), reverse=True,
    )

    blosum_mode_lower_threshold = (
        -1
        if args.blosum_strictness == "lax"
        else 0
        if args.blosum_strictness == "strict"
        else 999999
    )

    mat = bl.BLOSUM(62)
    filtered_mat = {key: {sub_key: value for sub_key, value in sub_dict.items() if value > blosum_mode_lower_threshold} for key, sub_dict in mat.items()}

    if args.processes > 1:
        arguments = []
        for input_gene in file_inputs:
            arguments.append(
                (
                    FlexcullArgs(
                        aa_path,
                        nt_path,
                        output_path,
                        args.matches,
                        input_gene,
                        args.debug,
                        args.base_pair,
                        args.verbose,
                        args.compress,
                        args.gap_threshold,
                        args.mismatches,
                        args.column_cull,
                        filtered_mat,
                        is_assembly,
                    ),
                ),
            )

        with Pool(args.processes) as pool:
            log_components = pool.starmap(do_gene, arguments, chunksize=1)
    else:
        log_components = [
            do_gene(
                FlexcullArgs(
                    aa_path,
                    nt_path,
                    output_path,
                    args.matches,
                    input_gene,
                    args.debug,
                    args.base_pair,
                    args.verbose,
                    args.compress,
                    args.gap_threshold,
                    args.mismatches,
                    args.column_cull,
                    filtered_mat,
                    is_assembly,
                ),
            )
            for input_gene in file_inputs
        ]

    if args.debug:
        log_global = []

        for component in log_components:
            log_global.extend(component)

        log_global.sort()
        log_global.insert(
            0, "Gene,Header,Cull To Start,Cull To End,Data Length,Data Removed\n",
        )
        log_out = os.path.join(output_path, "Culls.csv")
        with open(log_out, "w") as fp:
            fp.writelines(log_global)

    printv(f"Done! Took {folder_time.differential():.2f}s", args.verbose)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for folder in args.INPUT:
        do_folder(folder, args)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre FlexCull"
    raise Exception(
        msg,
    )