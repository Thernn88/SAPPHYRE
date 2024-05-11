"""FlexCull Description Goes Here.

PyLint 9.81/10
"""
from __future__ import annotations

from collections import Counter, defaultdict, namedtuple
from multiprocessing.pool import Pool
from os import listdir, makedirs, path
from shutil import rmtree
import warnings
from statistics import median, stdev
from Bio import BiopythonWarning
from blosum import BLOSUM
from sapphyre_tools import (
    find_index_pair,
    join_by_tripled_index,
    join_triplets_with_exclusions,
    join_with_exclusions,
    blosum62_distance,
    get_overlap,
)
from wrap_rocks import RocksDB

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
        "orthoset",
        "orthoset_input",
        "keep_codons",
        "blosum_max_percent",
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
        "is_assembly_or_genome",
        "is_ncg",  # non coding gene
        "keep_codons",
        "blosum_max_percent",
        "genome",
    ],
)


def align_col_removal(raw_fed_sequences: list, positions_to_keep: list) -> list:
    """Iterates over each sequence and deletes columns
    that were removed in the empty column removal.

    Args:
    ----
        raw_fed_sequences (list): List of tuples containing header and sequence
        positions_to_keep (list): List of positions to keep
    Returns:
    -------
        list: List of tuples containing header and sequence with empty columns removed
    """
    result = []
    for raw_sequence in raw_fed_sequences:
        sequence = raw_sequence[1]

        sequence = join_by_tripled_index(sequence, positions_to_keep)
        result.append((raw_sequence[0], sequence))

    return result


def delete_empty_columns(raw_fed_sequences: list, verbose: bool) -> tuple[list, list]:
    """Iterates over each sequence and deletes columns
    that consist of 100% dashes.

    Args:
    ----
        raw_fed_sequences (list): List of tuples containing header and sequence
        verbose (bool): Whether to print verbose output
    Returns:
    -------
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
                sequence = "".join(raw_sequence[1][x] for x in positions_to_keep)
            except IndexError:
                if verbose:
                    print(
                        f"WARNING: Sequence length is not the same as other sequences: {raw_sequence[0]}"
                    )
                continue

            result.append((raw_sequence[0], sequence))

    return result, positions_to_keep


def folder_check(output_target_path: str) -> None:
    """Checks to see if input and output directory has the necessary
    folder structure.

    if not, create missing folders.

    Args:
    ----
        output_target_path (str): Path to output directory
    Returns:
    -------
        None
    """
    output_aa_path = path.join(output_target_path, "aa")
    output_nt_path = path.join(output_target_path, "nt")
    rmtree(output_aa_path, ignore_errors=True)
    rmtree(output_nt_path, ignore_errors=True)
    makedirs(output_aa_path, exist_ok=True)
    makedirs(output_nt_path, exist_ok=True)


def parse_fasta(fasta_path: str) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """Parses fasta file into header and sequences.

    Args:
    ----
        fasta_path (str): Path to fasta file
    Returns:
    -------
        tuple[
            list[tuple[str, str]], list[tuple[str, str]]
        ]: Tuple containing list of reference record tuples and candidate record tuples
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


def get_internal_gaps(input_string):
    positions = []

    for i in range(*find_index_pair(input_string, "-")):
        if input_string[i] == "-":
            positions.append(i)
    
    return positions


def insert_gaps(input_string, positions):
    input_string = list(input_string)
    for coord in positions:
        input_string.insert(coord, "-")

    return ''.join(input_string)


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
    ----
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
    -------
        range: Range of the trimmed slice"""

    offset = amt_matches + 1
    cull_end = upper_limit

    # Trim to the right
    for i in range(starting_index, len(sequence) - 1):
        skip_first = 0
        char = sequence[i]
        mismatch = mismatches

        # If current index will look beyond the length of the sequence, skip
        if i == upper_limit - offset:
            break

        # If candidate is a gap or all references are dashes at this index, skip
        if char == "-" or all_dashes_by_index[i]:
            continue

        # Allow first position to be a mismatch
        if char not in character_at_each_pos[i]:
            skip_first = 1
            mismatch -= 1

        # If mismatch exhausted, skip
        if mismatch < 0:
            continue

        pass_all = True
        checks = amt_matches - 1
        match_i = 1

        # Scan until amt_matches is reached, upper limit is reached or mismatch is exhausted
        while checks > 0:
            # Scan will extend beyond the end of the sequence
            if i + match_i >= upper_limit:
                pass_all = False
                break

            next_char = sequence[i + match_i]

            if next_char == "-":
                # References don't contain gap at this index, fail
                if gap_present_threshold[i + match_i]:
                    pass_all = False
                    break
                match_i += 1
            elif next_char not in character_at_each_pos[i + match_i]:
                mismatch -= 1
                # Mismatch exhausted, fail
                if mismatch < 0:
                    pass_all = False
                    break
                match_i += 1
                checks -= 1
            else:
                # Match, move right by one and decrement checks
                match_i += 1
                checks -= 1

        if pass_all:
            cull_end = i + skip_first
            break

    cull_start = starting_index

    # Repeat cull to the left
    for i in range(starting_index - 1, -1, -1):
        mismatch = mismatches
        skip_last = 0

        char = sequence[i]

        # Scan will extend beyond the start of the sequence
        if i < lower_limit + offset:
            break

        # If candidate is a gap or all references are dashes at this index, skip
        if char == "-" or all_dashes_by_index[i]:
            continue

        # Allow first position to be a mismatch
        if char not in character_at_each_pos[i]:
            skip_last += 1
            mismatch -= 1

        # If mismatch exhausted, skip
        if mismatch < 0:
            continue

        pass_all = True
        checks = amt_matches - 1
        match_i = 1

        # Scan until amt_matches is reached, lower limit is reached or mismatch is exhausted
        while checks > 0:
            # Scan will extend beyond the start of the sequence
            if i - match_i < lower_limit:
                pass_all = False
                break

            prev_char = sequence[i - match_i]

            if prev_char == "-":
                # References don't contain gap at this index, fail
                if gap_present_threshold[i - match_i]:
                    pass_all = False
                    break
                match_i += 1
            elif prev_char not in character_at_each_pos[i - match_i]:
                mismatch -= 1
                # Mismatch exhausted, fail
                if mismatch < 0:
                    pass_all = False
                    break
                match_i += 1
                checks -= 1
            else:
                # Match, move left by one and decrement checks
                match_i += 1
                checks -= 1

        if pass_all:
            cull_start = i - skip_last + 1  # Inclusive
            break

    return range(cull_start, cull_end)


def get_data_difference(trim: int, ref: int) -> float:
    """
    Returns the difference between the trim and ref.

    Args:
    ----
        trim (int): Integer value A
        ref (int): Integer Value B
    Returns:
    -------
        float: Difference between trim and ref, q A over B
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
    blosum_max_percent,
    blosum_at_each_pos,
    gap_present_threshold: list[bool],
) -> tuple[int, int, bool]:
    """
    Culls leading and trailing bp from a sequence based on if the bp is present in a reference sequence.

    Args:
    ----
        sequence (str): Sequence to cull
        sequence_length (int): Length of the sequence
        offset (int): Offset to start culling from
        amt_matches (int): Amount of consecutive matches required to pass
        mismatches (int): Amount of mismatches allowed in consecutive matches
        all_dashes_by_index (list[bool]): List of booleans indicating whether all references have a dash at the index
        character_at_each_pos (list[set[str]]): List of sets containing characters at each position
        gap_present_threshold (list[bool]): List of booleans indicating whether references contain gaps at the index
    Returns:
    -------
        tuple: Tuple containing the index of the first bp that has required consecutive matches,
        the index of the last bp that has required consecutive matches, and a boolean indicating whether the trim kicked the entire sequence
    """
    cull_start = None
    cull_end = None
    kick = False
    for i, char in enumerate(sequence):
        mismatch = mismatches
        skip_first = 0
        # Scan will extend beyond the end of the sequence
        if i == sequence_length - offset:
            kick = True
            break

        # Don't allow cull to begin from a gap
        if char == "-":
            continue

        window_start = i

        # Don't allow cull to point of all dashes in references
        if all_dashes_by_index[i]:
            continue

        # Allow first position to be a mismatch
        if char not in character_at_each_pos[i]:
            skip_first = 1
            if window_start == i:
                continue
            mismatch -= 1

        # If mismatch exhausted, skip
        if mismatch < 0:
            continue

        blosum_matches = 0

        if char in blosum_at_each_pos[i]:
            blosum_matches += 1

        pass_all = True
        checks = amt_matches - 1
        match_i = 1

        # Scan until amt_matches is reached, end of sequence is reached or mismatch is exhausted
        while checks > 0:
            # Scan will extend beyond the end of the sequence
            if i + match_i >= len(sequence):
                pass_all = False
                break

            if sequence[i + match_i] == "-":
                # References don't contain gap at this index, fail
                if gap_present_threshold[i + match_i]:
                    pass_all = False
                    break
                # Skip if references also have gap
                match_i += 1

            # If candidate is a mismatch, decrement mismatch. If mismatch exhausted, fail.
            elif sequence[i + match_i] not in character_at_each_pos[i + match_i]:
                mismatch -= 1
                if mismatch < 0 or sequence[i + match_i] == "*":
                    pass_all = False
                    break
                match_i += 1
                checks -= 1
            else:
                # Match, move right by one and decrement checks
                if sequence[i + match_i] in blosum_at_each_pos[i + match_i]:
                    blosum_matches += 1
                match_i += 1
                checks -= 1

        if pass_all:
            # 40% of matches are blosum matches
            if blosum_max_percent != -1 and blosum_matches / match_i > blosum_max_percent:
                continue
            cull_start = i + skip_first

            break

    if not kick:
        # If not kicked from Cull Start Calc. Continue
        window_end = -1
        for i in range(len(sequence) - 1, -1, -1):
            mismatch = mismatches
            skip_last = 0

            char = sequence[i]
            # Scan will extend beyond the previous cull position
            if i < cull_start + offset:
                kick = True
                break

            # Don't allow cull to begin from a gap
            if char == "-":
                continue
            window_end = i

            # Don't allow cull to point of all dashes
            if all_dashes_by_index[i]:
                continue

            # Allow mismatch
            if char not in character_at_each_pos[i]:
                skip_last += 1
                if window_end == i:
                    continue
                mismatch -= 1

            # If mismatch exhausted, skip
            if mismatch < 0:
                continue
            
            blosum_matches = 0

            if char in blosum_at_each_pos[i]:
                blosum_matches += 1

            pass_all = True
            checks = amt_matches - 1
            match_i = 1
            # Scan until amt_matches is reached, lower limit is reached or mismatch is exhausted
            while checks > 0:
                if i - match_i < 0:
                    pass_all = False
                    break

                if sequence[i - match_i] == "-":
                    # References don't contain gap at this index, fail
                    if gap_present_threshold[i - match_i]:
                        pass_all = False
                        break
                    match_i += 1
                elif sequence[i - match_i] not in character_at_each_pos[i - match_i]:
                    mismatch -= 1
                    # Don't allow mismatch to be a stop codon
                    if mismatch < 0 or sequence[i - match_i] == "*":
                        pass_all = False
                        break
                    match_i += 1
                    checks -= 1
                else:
                    if sequence[i - match_i] in blosum_at_each_pos[i - match_i]:
                        blosum_matches += 1
                    # Match, move left by one and decrement checks
                    match_i += 1
                    checks -= 1

            if pass_all:
                if blosum_max_percent != -1 and blosum_matches / match_i > blosum_max_percent:
                    continue
                cull_end = i - skip_last + 1  # Inclusive
                break

    return cull_start, cull_end, kick


def process_refs(
    references: list[tuple], gap_threshold: float, column_cull_percent: float, mat: dict
) -> tuple:
    """
    Processes the references to generate the character_at_each_pos, gap_present_threshold , column_cull and all_dashes_by_index.

    character_at_each_pos: set of all characters aswell as the positive blosum subsitutions in each column of the references
    gap_present_threshold: boolean indicating whether the amount of data in a column is above the threshold
    column_cull: set of columns to cull
    all_dashes_by_index: boolean indicating whether all references have a dash at the index

    Args:
    ----
        references (list[tuple]): List of tuples containing the reference name and sequence
        gap_threshold (float): Threshold for the amount of gaps required in a reference for it to be a reference gap column
        column_cull_percent (float): Percent of data required in a column to allow it to pass
        mat (dict): The blosum62 subsitution matrix
    Returns:
    -------
        tuple: Tuple containing the character_at_each_pos, gap_present_threshold , column_cull and all_dashes_by_index tables

    """
    character_at_each_pos = defaultdict(list)
    blosum_at_each_pos = defaultdict(set)
    gap_present_threshold = {}
    column_cull = set()

    # Get the length of the longest reference and create
    # a dict with a default value for every index covered by the references
    max_ref_length = max(len(sequence) for _, sequence in references)
    all_dashes_by_index = {i: True for i in range(max_ref_length)}

    for _, sequence in references:
        for i, char in enumerate(sequence.replace("*", "-")):
            character_at_each_pos[i].append(char)

    for i, chars in list(character_at_each_pos.items()):
        data_present = 1 - (chars.count("-") / len(chars))
        all_dashes_by_index[i] = data_present == 0
        gap_present_threshold[i] = data_present >= gap_threshold
        if data_present < column_cull_percent:
            column_cull.add(i * 3)

        # Include all blosum subsitutions for each character in the column
        blosum_chars = {blosum_sub for char in chars for blosum_sub in mat.get(char, {}).keys()}
        chars = set(chars)
        character_at_each_pos[i] = chars.union(blosum_chars)
        blosum_at_each_pos[i] = blosum_chars - chars

    return (
        blosum_at_each_pos,
        character_at_each_pos,
        gap_present_threshold,
        all_dashes_by_index,
        column_cull,
    )


def column_cull_seqs(
    this_seqs: list[tuple], column_cull: set, minimum_bp: int
) -> tuple[list[tuple], list[str], set]:
    aa_out = []
    kicks = []

    this_column_cull = set()
    for nt_i in column_cull:
        i = nt_i // 3
        # Count the amount of each bp at this position
        this_sequence = Counter()
        internal_count = 0
        for _, seq, start, end in this_seqs:
            if start <= i <= end:
                internal_count += 1
                if seq[i] != "-":
                    this_sequence[seq[i]] += 1

        # If the most present bp has a count of 5 or greater and is present
        # in a third of sequences at this index, don't cull
        if sum(this_sequence.values()) >= 5:
            if this_sequence.most_common()[0][1] > internal_count / 3:
                continue
        this_column_cull.add(i * 3)

    # If columns require culling, use sapphyre_tools join with exclusions.
    # Check bp is still above minimum and then output
    if this_column_cull:
        for header, seq, _, _ in this_seqs:
            new_seq = join_with_exclusions(seq, this_column_cull)
            if len(new_seq) - new_seq.count("-") < minimum_bp:
                kicks.append(header)
                continue
            aa_out.append((header, new_seq))
    else:
        aa_out.extend([(header, sequence) for header, sequence, _, _ in this_seqs])

    return aa_out, kicks, this_column_cull


def cull_codons(
    out_line: str,
    cull_start: int,
    cull_end: int,
    amt_matches: int,
    mismatches: int,
    all_dashes_by_index: dict,
    character_at_each_pos: dict,
    gap_present_threshold: dict,
):
    positions_to_trim = set()
    codons = []
    # Save index of all codons
    for i in range(cull_start, cull_end):
        char = out_line[i]
        if char == "*":
            codons.append(i)
    kick = False
    non_trimmed_codons = [c for c in codons if c * 3 not in positions_to_trim]
    while non_trimmed_codons:
        # Start from the middle
        i = non_trimmed_codons[int(len(non_trimmed_codons) / 2)]

        # Use the trim around function to scan left and right of the codon
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

        # Grab the amount of columns where data is present in the refrences and
        # data is present in the candidate after trimming the codon
        left_side_ref_data_columns = sum(
            gap_present_threshold[x] for x in range(cull_start, i)
        )
        left_of_trim_data_columns = len(left_after) - left_after.count("-")

        right_side_ref_data_columns = sum(
            gap_present_threshold[x] for x in range(i, cull_end)
        )
        right_of_trim_data_columns = len(right_after) - right_after.count("-")

        left_highest_consecutive = 0
        right_highest_consecutive = 0

        this_consec = 0

        for ix, let in enumerate(left_after, cull_start):
            if let == "-":
                if not gap_present_threshold[ix]:
                    continue

                if this_consec > left_highest_consecutive:
                    left_highest_consecutive = this_consec
                this_consec = 0
            elif let in character_at_each_pos[ix]:
                this_consec += 1
        if this_consec > left_highest_consecutive:
            left_highest_consecutive = this_consec

        this_consec = 0
        for ix, let in enumerate(right_after, i):
            if let == "-":
                if not gap_present_threshold[ix]:
                    continue

                if this_consec > right_highest_consecutive:
                    right_highest_consecutive = this_consec
                this_consec = 0
            elif let in character_at_each_pos[ix]:
                this_consec += 1
        if this_consec > right_highest_consecutive:
            right_highest_consecutive = this_consec

        # If the difference between the amount of data columns in the candidate and
        # the reference is less than 55%, cull the remainder side
        cut_left, cut_right = False, False
        if (
            left_highest_consecutive < 30 and
            get_data_difference(
                left_of_trim_data_columns,
                left_side_ref_data_columns,
            )
            < 0.55
        ):
            for x in range(cull_start, i):
                positions_to_trim.add(x * 3)
                out_line[x] = "-"
            cut_left = True
        if (
            right_highest_consecutive < 30 and
            get_data_difference(
                right_of_trim_data_columns,
                right_side_ref_data_columns,
            )
            < 0.55
        ):
            for x in range(i, cull_end):
                positions_to_trim.add(x * 3)
                out_line[x] = "-"
            cut_right = True

        non_trimmed_codons.remove(i)

        if not cut_left and not cut_right:
            if gap_present_threshold[i]:
                kick = True, i
                return out_line, positions_to_trim, kick

    return out_line, positions_to_trim, kick


def trim_large_gaps(
    aa_out: list[tuple],
    reference_gap_col: set,
    amt_matches: int,
    mismatches: int,
    post_all_dashes_by_index: dict,
    post_character_at_each_pos: dict,
    post_gap_present_threshold: dict,
    minimum_bp: int,
    debug: bool,
):
    """
    Trims gaps with a length of 80 or more non-reference gap column gaps.

    Args:
    ----
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
    -------
        tuple: Tuple containing the updated aa_out, gap_pass_through to align the deletions to NT, the debug log, and kicks
    """
    log = []
    kicks = []
    gap_pass_through = {}
    for record_index, record in enumerate(aa_out):
        header, sequence = record

        if not header.endswith("."):
            gap_cull = set()
            seq_start, seq_end = find_index_pair(sequence, "-")
            change_made = False
            non_ref_gap_dash_count = 0
            raw_dash_count = 0
            out_line = list(sequence)

            for j, let in enumerate(out_line[seq_start:seq_end], seq_start):
                if let == "-":
                    if j not in reference_gap_col:
                        non_ref_gap_dash_count += 1
                    raw_dash_count += 1
                else:
                    # Consecutive gap with length of 80 or more where the majority of the
                    # gaps are not reference gaps
                    if non_ref_gap_dash_count >= 80:
                        i = j - (raw_dash_count // 2)
                        # Trim around the middle of the gap
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

                        # Grab the amount of columns where data is present in the refrences and
                        # data is present in the candidate after trimming the gap
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

                        # If both sides will end up being kicked force it to keep the side
                        # with the most bp.
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
                            ) >= len(
                                right_after
                            ) - right_after.count("-")
                            keep_right = not keep_left

                        # Kick the side that has less than 55% of the data columns compared to the reference
                        # and that isn't kept by previous logic
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
            # If change made in sequence, log the change, kick the sequence if necessary and output
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

                # Save cull data to pass onto NT sequences
                gap_pass_through[header] = gap_cull

    aa_out = [i for i in aa_out if i is not None]
    return aa_out, gap_pass_through, log, kicks


def align_to_aa_order(nt_out, aa_content):
    """
    Aligns the nt output list to the order of the aa content list.

    Args:
    ----
        nt_out (list): List of tuples containing header and sequence
        aa_content (list): List of tuples containing header and sequence
    Returns:
    -------
        list: An aligned list of tuples containing header and sequence
    """
    headers = [header for header, _ in aa_content if not header.endswith(".")]

    nt_out = dict(nt_out)
    for header in headers:
        yield (header, nt_out[header])


def do_cluster(ids, ref_coords, id_chomp_distance=100, max_distance=120):
    clusters = []
    ids.sort(key = lambda x: x[0])
    grouped_ids = defaultdict(list)
    for i, (child_index, seq_coords, start, end) in enumerate(ids):
        id = int(child_index.split("_")[0])
        grouped_ids[id].append((i, child_index, seq_coords, start, end))
        
    ids_ascending = sorted(grouped_ids.keys())
    

    req_seq_coverage = 0.5

    current_cluster = []
    
    for id in ids_ascending:
        seq_list = grouped_ids[id]
        if not current_cluster:
            current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list]
            current_index = id
            current_seqs = seq_list
            current_direction = "bi"
        else:
            passed = False
            passing_direction = None
            
            if id - current_index <= id_chomp_distance:
                for i, child_index, seq_coords, start, end in seq_list:
                    for _, _, _, current_start, current_end in current_seqs:
                        this_direction = None
                        if start == current_start and end == current_end:
                            this_direction = "bi"
                        else:
                            if start == current_start:
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
                            # distance = get_overlap(start, end, current_start, current_end, -max_distance)
                            # if distance is not None:
                            #     distance = abs(distance[1] - distance[0])
                            # if distance is not None and distance < max_distance:
                            passing_direction = this_direction
                            passed = True
                            break
                    if passed:
                        break
            
            if passed:
                current_cluster.extend([(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list])
                current_index = id
                current_seqs = seq_list
                if passing_direction != "bi":
                    current_direction = passing_direction
            else:
                if len(current_cluster) >= 2:
                    clusters.append((current_cluster[0][0], current_cluster[-1][0]))
                elif len(current_cluster) == 1:
                    if current_cluster[0][1] > req_seq_coverage:
                        clusters.append((current_cluster[0][0], current_cluster[0][0]))
                        
                current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list]
                current_index = id
                current_seqs = seq_list
                current_direction = "bi"
    
    if current_cluster:
        if len(current_cluster) >= 2:

            clusters.append((current_cluster[0][0], current_cluster[-1][0]))
        elif len(current_cluster) == 1:
            if current_cluster[0][1] > req_seq_coverage:
                cluster_coverage = current_cluster[0][1]
                clusters.append((current_cluster[0][0], current_cluster[0][0]))
                
    return clusters


def cull_reference_outliers(reference_records: list, debug: int) -> list:
    """
    Removes reference sequences which have an unusually large mean
    blosum distance. Finds the constrained blosum distance between
    each reference pull any reference with a mean 1.5x higher than
    the group mean. Returns the remaining references in a list.
    """
    distances_by_index = defaultdict(list)
    all_distances = []
    indices = {i:find_index_pair(reference_records[i][1], '-') for i in range(len(reference_records))}
    filtered = []
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
            dist = blosum62_distance(ref1[1][start:stop], ref2[1][start:stop])
            distances_by_index[i].append(dist)
            distances_by_index[j].append(dist)
            all_distances.append(dist)

    if not all_distances:
        return reference_records, filtered, 0, 0, 0

    total_median = median(all_distances)
    # ALLOWABLE_COEFFICENT = 2
    # allowable = max(total_median * ALLOWABLE_COEFFICENT, 0.3)

    std = stdev(all_distances) if len(all_distances) > 1 else 0
    allowable = total_median + (std*2)

    # if a record's mean is too high, cull it
    for index, distances in distances_by_index.items():
        this_median = median(distances)
        if this_median > allowable:# or mean > 1:
            distances_by_index[index] = None
            filtered.append( (reference_records[index], this_median, "kicked") )
        
        if debug == 2:
            filtered.append( (reference_records[index], this_median, "") )
        # else:
        #     distances_by_index[index] = mean
    # get all remaining records
    output = [reference_records[i] for i in range(len(reference_records)) if distances_by_index[i] is not None]
    return output, filtered, total_median, allowable, std


def do_gene(fargs: FlexcullArgs) -> None:
    """FlexCull main function. Culls input aa and nt using specified amount of matches."""
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    gene_path = path.join(fargs.aa_input, fargs.aa_file)
    this_gene = fargs.aa_file.split(".")[0]

    printv(f"Doing: {this_gene}", fargs.verbosity, 2)

    references, candidates = parse_fasta(gene_path)

    references, filtered_refs, total_median, allowable, std = cull_reference_outliers(references, fargs.debug)
    culled_references = []
    if filtered_refs:
        culled_references.append(f'{this_gene} total median: {total_median}\n')
        culled_references.append(f'{this_gene} threshold: {allowable}\n')
        culled_references.append(f'{this_gene} standard deviation: {std}\n')
        for ref_kick, ref_median, kick in filtered_refs:
            culled_references.append(f'{ref_kick[0]},{ref_median},{kick}\n')

    if not references:
        printv(f"No references for {this_gene} after cull", fargs.verbosity, 1)
        return [], 0, culled_references

    (
        blosum_at_each_pos,
        character_at_each_pos,
        gap_present_threshold,
        all_dashes_by_index,
        column_cull,
    ) = process_refs(
        references, fargs.gap_threshold, fargs.column_cull_percent, fargs.filtered_mat
    )

    log = []
    codon_log = []

    follow_through = {}
    offset = fargs.amt_matches + 1

    aa_out_path = path.join(fargs.output, "aa", fargs.aa_file.rstrip(".gz"))
    aa_out = references.copy()
    this_seqs = []

    get_id = lambda header: header.split("|")[3].replace("NODE_","")

    flattened_set = set()
    if fargs.genome:
        ids = []
        for header, raw_sequence in candidates:
            start, end = find_index_pair(raw_sequence, "-")
            data_cols = {i for i, let in enumerate(raw_sequence[start:end], start) if let != "-"}
            ids.append((get_id(header), data_cols, start, end))

        ref_coords = set()
        for header, sequence in references:
            start, end = find_index_pair(sequence, "-")
            for i, let in enumerate(sequence[start:end], start):
                if let != "-":
                    ref_coords.add(i)

        clusters = do_cluster(ids, ref_coords)
        if clusters:
            cluster_sets = [set(range(start, end+1)) for start, end in clusters]
            flattened_set = set.union(*cluster_sets)

    for header, raw_sequence in candidates:
        sequence = list(raw_sequence)

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
            fargs.blosum_max_percent,
            blosum_at_each_pos,
            gap_present_threshold,
        )

        # If cull didn't kick
        if not kick:
            # Pad the sequence with gaps and cull
            out_line = ["-"] * cull_start + sequence[cull_start:cull_end]

            characters_till_end = sequence_length - len(out_line)
            out_line += ["-"] * characters_till_end

            # If gene is not NCG and we don't want to keep codons, cull codons
            positions_to_trim = set()
            if not fargs.is_ncg and not fargs.keep_codons:
                out_line, positions_to_trim, kick = cull_codons(
                    out_line,
                    cull_start,
                    cull_end,
                    fargs.amt_matches,
                    fargs.mismatches,
                    all_dashes_by_index,
                    character_at_each_pos,
                    gap_present_threshold,
                )
                if kick:
                    follow_through[header] = True, 0, 0, []

                    if fargs.debug:
                        codon_log.append(f"{header},{kick[1]}\n")
                    kick = False
                    # continue

            # Join sequence and check bp after cull
            out_line = "".join(out_line)

            data_length = cull_end - cull_start
            bp_after_cull = len(out_line) - out_line.count("-")

            if get_id(header) in flattened_set or bp_after_cull >= fargs.bp:
                follow_through[header] = (
                    False,
                    cull_start,
                    cull_end,
                    positions_to_trim,
                )

                this_seqs.append((header, out_line, *find_index_pair(out_line, "-")))

                # Log removed data
                if fargs.debug:
                    rescued = " - Rescued" if get_id(header) in flattened_set and bp_after_cull < fargs.bp else ""
                    removed_section = sequence[:cull_start] + sequence[cull_end:]
                    data_removed = len(removed_section) - removed_section.count("-")
                    log.append(
                        gene
                        + ","
                        + header + rescued
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

    # If sequences still present after cull
    if this_seqs:
        # Cull columns
        this_seqs, kicks, this_column_cull = column_cull_seqs(
            this_seqs, column_cull, fargs.gap_threshold
        )
        for header in kicks:
            follow_through[header] = True, 0, 0, []

        # Remove empty columns from refs and candidates
        aa_out, aa_positions_to_keep = delete_empty_columns(aa_out + this_seqs, False)
        if len(aa_out) == len(references):
            return log, culled_references, codon_log  # Only refs

        # Recalcuate position based tables
        post_references = [i for i in aa_out if i[0].endswith(".")]

        (
            post_blosum_at_each_pos,
            post_character_at_each_pos,
            post_gap_present_threshold,
            post_all_dashes_by_index,
            _,
        ) = process_refs(
            post_references,
            fargs.gap_threshold,
            fargs.column_cull_percent,
            fargs.filtered_mat,
        )

        reference_gap_col = {i for i, x in post_gap_present_threshold.items() if not x}

        # Trim large gaps
        aa_out, gap_pass_through, trim_log, kicks = trim_large_gaps(
            aa_out,
            reference_gap_col,
            fargs.amt_matches,
            fargs.mismatches,
            post_all_dashes_by_index,
            post_character_at_each_pos,
            post_gap_present_threshold,
            fargs.bp,
            fargs.debug,
        )

        if fargs.debug:
            log.extend(trim_log)

        for kick in kicks:
            follow_through[kick] = True, 0, 0, []

        # If sequences still present after column cull and gap cull, output
        if len(aa_out) != len(references):
            writeFasta(aa_out_path, aa_out, fargs.compress)

            nt_file_name = fargs.aa_file.replace(".aa.", ".nt.")
            gene_path = path.join(fargs.nt_input, nt_file_name)

            references, candidates = parse_fasta(gene_path)

            # Align culls to NT
            nt_out_path = path.join(fargs.output, "nt", nt_file_name.rstrip(".gz"))
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

                    out_line = "".join(out_line)
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

            # Align order
            out_nt = align_to_aa_order(out_nt, aa_out)

            writeFasta(nt_out_path, out_nt, fargs.compress)

    return log, culled_references, codon_log


def do_folder(folder, args: MainArgs, non_coding_gene: set):
    folder_time = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Processing: {folder}", args.verbose, 0)
    aa_path = path.join(folder, args.amino_acid)
    nt_path = path.join(folder, args.nucleotide)
    output_path = path.join(folder, args.output)
    if not path.exists(aa_path) or not path.exists(nt_path):
        printv(
            f"WARNING: Can't find aa ({aa_path}) and nt ({nt_path}) folders. Abort",
            args.verbose,
        )
        return

    nt_db_path = path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(nt_db_path)
    dbis_assembly = nt_db.get("get:isassembly") == "True"
    dbis_genome = nt_db.get("get:isgenome") == "True"
    is_assembly_or_genome = dbis_assembly or dbis_genome

    folder_check(output_path)
    file_inputs = [
        input_gene
        for input_gene in listdir(aa_path)
        if input_gene.split(".")[-1] in {"fa", "gz", "fq", "fastq", "fasta"}
    ]
    file_inputs.sort(
        key=lambda x: path.getsize(path.join(aa_path, x)),
        reverse=True,
    )

    blosum_mode_lower_threshold = (
        -1
        if args.blosum_strictness == "lax"
        else 0
        if args.blosum_strictness == "strict"
        else 999999
    )

    if blosum_mode_lower_threshold == -1:
        blosum_max_percent = -1

    mat = BLOSUM(62)
    filtered_mat = {
        key: {
            sub_key: value
            for sub_key, value in sub_dict.items()
            if value > blosum_mode_lower_threshold
        }
        for key, sub_dict in mat.items()
    }

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
                    is_assembly_or_genome,
                    input_gene in non_coding_gene,
                    args.keep_codons,
                    args.blosum_max_percent,
                    dbis_genome,
                    
                ),
            ),
        )
    if args.processes > 1:
        with Pool(args.processes) as pool:
            log_components = pool.starmap(do_gene, arguments, chunksize=1)
    else:
        log_components = [
            do_gene(
                arg[0]
            )
            for arg in arguments
        ]

    if args.debug:
        log_global = []
        ref_log_global = []
        codon_log_global = []

        for component, ref_kicks, codon_kicks in log_components:
            log_global.extend(component)
            ref_log_global.extend(ref_kicks)
            codon_log_global.extend(codon_kicks)

        log_global.sort()
        log_global.insert(
            0,
            "Gene,Header,Cull To Start,Cull To End,Data Length,Data Removed\n",
        )
        log_out = path.join(output_path, "Culls.csv")
        with open(log_out, "w") as fp:
            fp.writelines(log_global)

        ref_log_global.insert(0, "Header,Median\n")

        ref_log_out = path.join(output_path, "Reference_culls.csv")
        with open(ref_log_out, "w") as fp:
            fp.writelines(ref_log_global)

        codon_log_global.insert(0, "Header,Index\n")

        codon_log_out = path.join(output_path, "Codon_culls.csv")
        with open(codon_log_out, "w") as fp:
            fp.writelines(codon_log_global)
            
    printv(
        f"Done! Took {folder_time.differential():.2f}s.",
        args.verbose,
    )


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input
    orthoset_db_path = path.join(orthosets_dir, orthoset, "rocksdb")
    orthoset_db = RocksDB(orthoset_db_path)
    orthoset_non_coding_genes = orthoset_db.get("get:nc_genes")
    orthoset_non_coding_genes = (
        set(orthoset_non_coding_genes.split(","))
        if orthoset_non_coding_genes
        else set()
    )

    for folder in args.INPUT:
        do_folder(folder, args, orthoset_non_coding_genes)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre FlexCull"
    raise Exception(
        msg,
    )
