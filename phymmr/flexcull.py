"""
FlexCull Description Goes Here

PyLint 9.81/10
"""
from __future__ import annotations

import os
from collections import namedtuple
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
        "mismatches"
    ],
)


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


def trim_around(starting_index, sequence, amt_matches, mismatches, match_percent, all_dashes_by_index, character_at_each_pos, gap_present_threshold) -> None:
    """
    Trim around a given position in a sequence
    """
    offset = amt_matches - 1
    cull_end = None
    for i in range(starting_index, len(sequence)-1):
        skip_first = 0
        char = sequence[i]
        mismatch = mismatches

        if i == len(sequence) - offset:
            kick = True
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
            cull_end = i + skip_first
            break

    if not cull_end:
        return True, None

    cull_start = None
    kick = False
    for i in range(starting_index-1, -1, -1):
        mismatch = mismatches
        skip_last = 0

        char = sequence[i]
        if i < 0 + offset:
            kick = True
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
            cull_start = i - skip_last + 1  # Inclusive
            break
    if kick or not cull_start:
        return True, None
    return kick, range(cull_start, cull_end)


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
    mismatches: int
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

    for header, sequence in references:
        for i, char in enumerate(sequence):
            if i not in character_at_each_pos:
                character_at_each_pos[i] = [char]
            else:
                character_at_each_pos[i].append(char)
            #  if char isnt a hyphen, this positions can't be all dashes
            if char != "-":
                all_dashes_by_index[i] = False
    
    for i, chars in character_at_each_pos.items():
        gap_present_threshold[i] = chars.count("-") / len(chars) < gap_threshold

    log = []

    follow_through = {}
    offset = amt_matches - 1

    aa_out_path = os.path.join(output, "aa", aa_file.rstrip(".gz"))
    aa_out = references.copy()

    for header, sequence in candidates:
        sequence = list(sequence)
        
        gene = header.split("|")[0]

        if gene not in follow_through:
            follow_through[gene] = {}

        cull_start = None
        kick = False
        data_removed = 0
        data_length = 0

        sequence_length = len(sequence)

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
                    cull_end = i - skip_last + 1  # Inclusive
                    break

        if not kick:  # If also passed Cull End Calc. Finish
            out_line = ["-"] * cull_start + sequence[cull_start:cull_end]

            characters_till_end = sequence_length - len(out_line)
            out_line += ["-"] * characters_till_end

            positions_to_trim = []
            for i, char in enumerate(out_line):
                if char == "*":
                    kick, positions = trim_around(i, out_line, amt_matches, mismatches, match_percent, all_dashes_by_index, character_at_each_pos, gap_present_threshold)
                    if kick:
                        break
                    for i in positions:
                        positions_to_trim.append(i)
                        out_line[i] = "-"

            out_line = "".join(out_line)     
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

            out_line = [out_line[i:i+3] for i in range(0, len(out_line), 3)]
            for pos in positions_to_trim:
                out_line[pos] = "---"
            out_line = "".join(out_line)

            nt_out.append((header, out_line))

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
                    args.mismatches
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
                args.mismatches
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