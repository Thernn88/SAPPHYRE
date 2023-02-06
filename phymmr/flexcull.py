"""
FlexCull Description Goes Here

PyLint 9.81/10
"""
from __future__ import annotations

import os
from collections import namedtuple
from itertools import chain
from multiprocessing.pool import Pool
from shutil import rmtree

from .utils import printv, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode

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
        "match_percent",
        "compress",
        "gap_threshold",
        "mismatches",
        "column_cull"
    ],
)


def align_col_removal(raw_fed_sequences: list, positions_to_keep: list) -> list:
    """
    Iterates over each sequence and deletes columns
    that were removed in the empty column removal.
    """
    # raw_sequences = [
    #     i.replace("\n", "") for i in raw_fed_sequences if i.replace("\n", "") != ""
    # ]

    result = []
    raw_sequences = [*chain.from_iterable(raw_fed_sequences)]
    for i in range(0, len(raw_sequences), 2):
        # result.append(raw_sequences[i])

        sequence = raw_sequences[i + 1]

        sequence = [sequence[i * 3 : (i * 3) + 3] for i in positions_to_keep]
        x = ''.join(sequence)
        result.append((raw_sequences[i], "".join(sequence)))

    return result

def delete_empty_columns(raw_fed_sequences: list, verbose: bool) -> tuple[list, list]:
    """
    Iterates over each sequence and deletes columns
    that consist of 100% dashes. In this version, raw_feed_sequences
    is a list of tuples.
    """
    result = []
    # sequences = []
    # raw_sequences = [
    #     i.replace("\n", "") for i in raw_fed_sequences if i.replace("\n", "") != ""
    # ]
    # for i in range(0, len(raw_sequences), 2):
    #     sequences.append(raw_sequences[i + 1])
    raw_sequences = [*chain.from_iterable(raw_fed_sequences)]
    sequences = [x[1] for x in raw_fed_sequences]
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
                # result.append(raw_sequences[i])
            except IndexError:
                printv(
                    f"WARNING: Sequence length is not the same as other sequences: {raw_sequences[i]}",
                    verbose,
                    0,
                )
                continue
            sequence = "".join(sequence)

            result.append((raw_sequences[i],sequence))

    return result, positions_to_keep
def folder_check(output_target_path: str, input_target_path: str) -> str:
    """
    Checks to see if input and output directory has the necessary
    folder structure

    if not, create missing folders.

    Finally, checks if linux memory temp folder exists
    else creates tmp folder in the input

    Returns path to tmp folder
    """
    output_aa_path = os.path.join(output_target_path, "aa")
    output_nt_path = os.path.join(output_target_path, "nt")
    rmtree(output_aa_path, ignore_errors=True)
    rmtree(output_nt_path, ignore_errors=True)
    os.makedirs(output_aa_path, exist_ok=True)
    os.makedirs(output_nt_path, exist_ok=True)

    target_tmp_path = "/run/shm"

    if not os.path.exists(target_tmp_path):
        target_tmp_path = "/dev/shm"

    if not os.path.exists(target_tmp_path):
        target_tmp_path = os.path.join(input_target_path, "tmp")
        os.makedirs(target_tmp_path, exist_ok=True)

    return target_tmp_path


def make_nt(aa_file_name: str) -> str:
    """
    Converts AA file name to NT file name
    """
    return aa_file_name.replace(".aa.", ".nt.")


def parse_fasta(fasta_path: str) -> tuple:
    """
    Parses fasta file into header and sequences
    """
    references = []
    candidates = []

    for header, sequence in parseFasta(fasta_path):
        if header[-1] == ".":  # Is reference
            references.append((header, sequence))

        else:
            candidates.append((header, sequence))

    return references, candidates


def trim_around(starting_index, lower_limit, upper_limit, sequence, amt_matches, mismatches, match_percent, all_dashes_by_index, character_at_each_pos, gap_present_threshold) -> None:
    """
    Trim around a given position in a sequence
    """
    offset = amt_matches - 1
    cull_end = upper_limit
    for i in range(starting_index, len(sequence)-1):
        skip_first = 0
        char = sequence[i]
        mismatch = mismatches

        if i == upper_limit - offset:
            break
        
        if char == "-":
            continue
        if all_dashes_by_index[i]:
            continue
        if (
            not character_at_each_pos[i].count(char) / len(character_at_each_pos[i])
            >= match_percent
        ):
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
            
            if sequence[i + match_i] == "-":
                if gap_present_threshold[i + match_i]:
                    pass_all = False
                    break
                match_i += 1
            elif (
                not character_at_each_pos[i + match_i].count(sequence[i + match_i])
                / len(character_at_each_pos[i + match_i])
                >= match_percent
            ):
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
    for i in range(starting_index-1, -1, -1):
        mismatch = mismatches
        skip_last = 0

        char = sequence[i]
        if i < lower_limit + offset:
            break
        if char == "-":
            continue

        all_dashes_at_position = all_dashes_by_index[i]
        if all_dashes_at_position:
            # Don't allow cull to point of all dashes
            continue
        if (
            not character_at_each_pos[i].count(char)
            / len(character_at_each_pos[i])
            >= match_percent
        ):
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

            if sequence[i - match_i] == "-":
                if gap_present_threshold[i - match_i]:
                    pass_all = False
                    break
                match_i += 1
            elif (
                not character_at_each_pos[i - match_i].count(
                    sequence[i - match_i]
                )
                / len(character_at_each_pos[i - match_i])
                >= match_percent
            ):
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


def get_data_difference(trim, ref):
    if trim == 0:
        return 0
    if ref == 0 and trim != 0:
        return 1
    if trim == 0 and ref == 0:
        return 1
    return trim / ref


def get_cull_start_end(sequence, mismatches, sequence_length, match_percent, amt_matches, all_dashes_by_index, character_at_each_pos, gap_present_threshold):
    offset = amt_matches - 1
    cull_start = None
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

        if all_dashes_by_index[i]:
            continue
        if (
            not character_at_each_pos[i].count(char) / len(character_at_each_pos[i])
            >= match_percent
        ):
            skip_first = 1
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
            elif (
                not character_at_each_pos[i + match_i].count(sequence[i + match_i])
                / len(character_at_each_pos[i + match_i])
                >= match_percent
            ):
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
            cull_start = i + skip_first
            break

    if not kick:
        # If not kicked from Cull Start Calc. Continue
        cull_end = None
        for i in range(len(sequence)-1, -1, -1):
            mismatch = mismatches
            skip_last = 0

            char = sequence[i]
            
            if i < cull_start + offset:
                kick = True
                break
            if char == "-":
                continue

            if all_dashes_by_index[i]:
                # Don't allow cull to point of all dashes
                continue
            if (
                not character_at_each_pos[i].count(char)
                / len(character_at_each_pos[i])
                >= match_percent
            ):
                skip_last += 1
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
                elif (
                    not character_at_each_pos[i - match_i].count(
                        sequence[i - match_i]
                    )
                    / len(character_at_each_pos[i - match_i])
                    >= match_percent
                ):
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
                cull_end = i - skip_last + 1 # Inclusive
                break
    
    return cull_start, cull_end, kick

def cull_internal_gaps(org_sequence, gap_present_threshold, sequence_length):
    # Cull internal gaps
    bp_encountered = False
    gap_start = None
    gap_end = None
    cull_this_gap = False
    sequence = org_sequence.copy()
    for i, let in enumerate(sequence):
        if let != "-":
            bp_encountered = True
        
        if bp_encountered:
            if let == "-":
                if not gap_start:
                    gap_start = i

                if gap_present_threshold[i]:
                    cull_this_gap = True
            else:
                if cull_this_gap:
                    gap_end = i
                    for j in range(gap_start, gap_end):
                        sequence[j] = None

                gap_start = None
                gap_end = None
                cull_this_gap = False

    result = [i for i in sequence if i is not None]
    return result + ["-"] * (sequence_length-len(result))

def do_gene(
    aa_input: str,
    nt_input: str,
    output: str,
    amt_matches: int,
    aa_file: str,
    match_percent: bool,
    debug: bool,
    bp: int,
    verbosity: int,
    compress: bool,
    gap_threshold: float,
    mismatches: int,
    column_cull_percent: float,
) -> None:
    """
    FlexCull main function. Culls input aa and nt using specified amount of matches
    """
    gene_path = os.path.join(aa_input, aa_file)
    this_gene = aa_file.split(".")[0]

    printv(f"Doing: {this_gene}", verbosity, 2)

    references, candidates = parse_fasta(gene_path)

    character_at_each_pos = {}
    gap_present_threshold = {}

    #  make a list, each index will match an index in the reference strings
    #  check this later instead of counting hyphens
    #  convert this to an int if you need to break on more than one data character
    max_ref_length = 0
    for _, sequence in references:
        max_ref_length = max(max_ref_length, len(sequence))
    all_dashes_by_index = [True] * max_ref_length

    column_cull = set()

    for header, sequence in references:
        for i, char in enumerate(sequence):
            if char == "*":
                char = "-"
            if i not in character_at_each_pos:
                character_at_each_pos[i] = [char]
            else:
                character_at_each_pos[i].append(char)
            #  if char isnt a hyphen, this positions can't be all dashes
            if char != "-":
                all_dashes_by_index[i] = False
    
    for i, chars in character_at_each_pos.items():
        data_present = 1 - (chars.count("-") / len(chars))
        gap_present_threshold[i] = data_present >= gap_threshold
        if data_present < column_cull_percent:
            column_cull.add(i*3)

    log = []

    follow_through = {}

    aa_out_path = os.path.join(output, "aa", aa_file.rstrip(".gz"))
    aa_out = references.copy()

    for header, sequence in candidates:
        sequence = list(sequence)
        
        gene = header.split("|")[0]

        if gene not in follow_through:
            follow_through[gene] = {}

        data_removed = 0
        data_length = 0

        sequence_length = len(sequence)

        #Internal Gap Cull
        internal_culled_sequence = cull_internal_gaps(sequence, gap_present_threshold, sequence_length)
        cull_start, cull_end, kick = get_cull_start_end(internal_culled_sequence, mismatches, sequence_length, match_percent, amt_matches, all_dashes_by_index, character_at_each_pos, gap_present_threshold)
        if not kick:
            internal_culled = internal_culled_sequence[cull_start:cull_end]
            internal_cull_bp = len(internal_culled) - internal_culled.count("-")
        else:
            internal_cull_bp = -1

        # cull_start, cull_end, kick = get_cull_start_end(sequence, mismatches, sequence_length, match_percent, amt_matches, all_dashes_by_index, character_at_each_pos, gap_present_threshold)

        if not kick:  # If also passed Cull End Calc. Finish
            candidate_trimmed_sequence = sequence[cull_start:cull_end]
            candidate_trimmed_bp = len(candidate_trimmed_sequence) - candidate_trimmed_sequence.count("-")

            trimmed_sequence = candidate_trimmed_sequence if candidate_trimmed_bp >= internal_cull_bp else internal_culled_sequence

            data_before = sequence_length - sequence.count("-")
            out_line = ["-"] * cull_start + trimmed_sequence

            characters_till_end = sequence_length - len(out_line)
            out_line += ["-"] * characters_till_end

            positions_to_trim = set()
            codons = []
            if (len(out_line) - out_line.count("-")) != data_before:
                for i in range(cull_start, cull_end):
                    char = out_line[i]
                    if char == "*":
                        codons.append(i)
            
            non_trimmed_codons = [c for c in codons if c*3 not in positions_to_trim]
            while non_trimmed_codons:
                i = non_trimmed_codons[int(len(non_trimmed_codons)/2)]
                positions = trim_around(i, cull_start, cull_end, out_line, amt_matches, mismatches, match_percent, all_dashes_by_index, character_at_each_pos, gap_present_threshold)
                for x in positions:
                    positions_to_trim.add(x*3)
                    out_line[x] = "-"
                
                left_after = out_line[cull_start:i]
                right_after = out_line[i:cull_end]

                left_side_ref_data_columns = sum([gap_present_threshold[x] for x in range(cull_start, i)])
                left_of_trim_data_columns = len(left_after) - left_after.count("-")

                right_side_ref_data_columns = sum([gap_present_threshold[x] for x in range(i, cull_end)])
                right_of_trim_data_columns = len(right_after) - right_after.count("-")

                if get_data_difference(left_of_trim_data_columns, left_side_ref_data_columns) < 0.55: # candidate has less than % of data columns compared to reference
                    for x in range(cull_start, i):
                        positions_to_trim.add(x*3)
                        out_line[x] = "-"
                if get_data_difference(right_of_trim_data_columns, right_side_ref_data_columns) < 0.55:
                    for x in range(i, cull_end):
                        positions_to_trim.add(x*3)
                        out_line[x] = "-"

                non_trimmed_codons = [c for c in codons if c*3 not in positions_to_trim]

            out_line = "".join([let if i*3 not in column_cull else "-" for i, let in enumerate(out_line)])  
            if kick:
                follow_through[gene][header] = True, 0, 0, []
                if debug:
                    log.append(
                        gene
                        + ","
                        + header
                        + ",Kicked,Codon cull found no match,0,\n"
                    )
                continue

            # The cull replaces data positions with dashes to maintain the same alignment
            # while removing the bad data

            data_length = cull_end - cull_start
            bp_after_cull = len(out_line) - out_line.count("-")

            if bp_after_cull >= bp:
                follow_through[gene][header] = False, cull_start, cull_end, positions_to_trim

                aa_out.append((header, out_line))

                if debug:
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
                        + "\n"
                    )
            else:
                follow_through[gene][header] = True, 0, 0, []
                if debug:
                    log.append(
                        gene
                        + ","
                        + header
                        + ",Kicked,Minimum BP Not Met,"
                        + str(bp_after_cull)
                        + ",\n"
                    )
        if kick:
            follow_through[gene][header] = True, 0, 0, []

            if debug:
                log.append(gene + "," + header + ",Kicked,Zero Data After Cull,0,\n")

    # remove empty columns from refs and
    aa_out, aa_positions_to_keep = delete_empty_columns(aa_out, False)
    if len(aa_out) == len(references):
        return log  # Only refs

    writeFasta(aa_out_path, aa_out, compress)

    nt_file_name = make_nt(aa_file)
    gene_path = os.path.join(nt_input, nt_file_name)

    references, candidates = parse_fasta(gene_path)

    nt_out_path = os.path.join(output, "nt", nt_file_name.rstrip(".gz"))
    nt_out = references.copy()
    for header, sequence in candidates:
        gene = header.split("|")[0]
        kick, cull_start, cull_end, positions_to_trim = follow_through[gene][header]

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

            out_line = [out_line[i:i+3] if i not in positions_to_trim and i not in column_cull else "---" for i in range(0, len(out_line), 3)]
            out_line = "".join(out_line)

            nt_out.append((header, out_line))
    nt_out = align_col_removal(nt_out, aa_positions_to_keep)
    writeFasta(nt_out_path, nt_out, compress)

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
    folder_check(output_path, folder)
    file_inputs = [
        input_gene
        for input_gene in os.listdir(aa_path)
        if input_gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]
    file_inputs.sort(
        key=lambda x: os.path.getsize(os.path.join(aa_path, x)), reverse=True
    )

    if args.processes > 1:
        arguments = []
        for input_gene in file_inputs:
            arguments.append(
                (
                    aa_path,
                    nt_path,
                    output_path,
                    args.matches,
                    input_gene,
                    args.match_percent,
                    args.debug,
                    args.base_pair,
                    args.verbose,
                    args.compress,
                    args.gap_threshold,
                    args.mismatches,
                    args.column_cull
                )
            )

        with Pool(args.processes) as pool:
            log_components = pool.starmap(do_gene, arguments, chunksize=1)
    else:
        log_components = [
            do_gene(
                aa_path,
                nt_path,
                output_path,
                args.matches,
                input_gene,
                args.match_percent,
                args.debug,
                args.base_pair,
                args.verbose,
                args.compress,
                args.gap_threshold,
                args.mismatches,
                args.column_cull
            )
            for input_gene in file_inputs
        ]

    if args.debug:
        log_global = []

        for component in log_components:
            log_global.extend(component)

        log_global.sort()
        log_global.insert(
            0, "Gene,Header,Cull To Start,Cull To End,Data Length,Data Removed\n"
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
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr FlexCull"
    )
