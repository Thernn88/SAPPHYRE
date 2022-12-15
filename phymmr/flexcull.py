"""
FlexCull Description Goes Here

PyLint 9.81/10
"""
from __future__ import annotations

import os
from collections import namedtuple
from multiprocessing.pool import Pool
from time import time

from Bio import AlignIO

from .utils import printv
from .timekeeper import TimeKeeper, KeeperMode

MainArgs = namedtuple(
    "MainArgs",
    [
        'verbose',
        'processes',
        'debug',
        'INPUT',
        'output',
        'amino_acid',
        'nucleotide',
        'matches',
        'base_pair',
        'match_percent'
    ]
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
    os.makedirs(output_aa_path, exist_ok=True)
    os.makedirs(output_nt_path, exist_ok=True)

    target_tmp_path = "/run/shm"

    if not os.path.exists(target_tmp_path):
        target_tmp_path = '/dev/shm'

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
    raw_references = []

    with open(fasta_path, encoding="UTF-8") as fasta_io:
        fasta_file = AlignIO.parse(fasta_io, "fasta")
        for seq_record in fasta_file:
            for seq in seq_record:
                header = seq.name
                sequence = str(seq.seq)

                if header[-1] == '.':  #Is reference
                    references.append((header, sequence))

                    raw_references.append(">"+header+"\n"+sequence+"\n")

                else:
                    candidates.append((header, sequence))

    raw_references = "".join(raw_references)

    return references, candidates, raw_references


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
) -> None:
    """
    FlexCull main function. Culls input aa and nt using specified amount of matches
    """
    gene_path = os.path.join(aa_input, aa_file)
    this_gene = aa_file.split(".")[0]

    printv("Doing: {this_gene}", verbosity, 2)

    references, candidates, raw_references = parse_fasta(gene_path)

    character_at_each_pos = {}

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
                
    log = []

    follow_through = {}
    offset = amt_matches - 1

    aa_out_path = os.path.join(output, "aa", aa_file)
    aa_out = [raw_references]

    for header, sequence in candidates:
        gene = header.split("|")[0].replace(">", "")

        if gene not in follow_through:
            follow_through[gene] = {}

        cull_start = None
        kick = False
        data_removed = 0

        data_length = 0

        sequence_length = len(sequence)

        for i, char in enumerate(sequence):
            if i == sequence_length - offset:
                kick = True
                break

            # Don't allow cull to point of all dashes
            if char == "-":
                continue
            
            all_dashes_at_position = all_dashes_by_index[i]
            if all_dashes_at_position:
                continue
            if not character_at_each_pos[i].count(char) / len(character_at_each_pos[i]) >= match_percent:
                continue

            pass_all = True
            for match_i in range(1, amt_matches):
                    if sequence[i + match_i] == "-" or (not character_at_each_pos[i + match_i].count(sequence[i + match_i]) / len(character_at_each_pos[i + match_i]) >= match_percent):
                        pass_all = False
                        break

            if pass_all:
                cull_start = i
                break

        if not kick:
            # If not kicked from Cull Start Calc. Continue
            cull_end = None
            for i_raw in range(len(sequence)):
                i = sequence_length - 1 - i_raw  # Start from end

                char = sequence[i]

                if i < cull_start + offset:
                    kick = True
                    break
                if char == "-":
                    continue

                all_dashes_at_position = all_dashes_by_index[i]
                if all_dashes_at_position:
                    # Don't allow cull to point of all dashes
                    continue
                if not character_at_each_pos[i].count(char) / len(character_at_each_pos[i]) >= match_percent:
                    continue

                pass_all = True
                for match_i in range(1, amt_matches):
                    if sequence[i - match_i] == "-" or (not character_at_each_pos[i - match_i].count(sequence[i - match_i]) / len(character_at_each_pos[i - match_i]) >= match_percent):
                        pass_all = False
                        break

                if pass_all:
                    cull_end = i + 1 #Inclusive
                    break

        if not kick:  # If also passed Cull End Calc. Finish
            out_line = ("-" * cull_start) + sequence[cull_start:cull_end] 
            
            characters_till_end = sequence_length - len(out_line)
            out_line += ("-" * characters_till_end)

            # The cull replaces data positions with dashes to maintain the same alignment
            # while removing the bad data

            data_length = cull_end - cull_start
            bp_after_cull = len(out_line) - out_line.count("-")

            if bp_after_cull >= bp:
                follow_through[gene][header] = False, cull_start, cull_end

                aa_out.append(">" + header + "\n")
                aa_out.append(out_line + "\n")

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
                follow_through[gene][header] = True, 0, 0
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
            follow_through[gene][header] = True, 0, 0

            if debug:
                log.append(gene + "," + header + ",Kicked,Zero Data After Cull,0,\n")

    if len(aa_out) == 1:
        return log #Only refs

    with open(aa_out_path, "w", encoding="UTF-8") as fp:
        fp.write("".join(aa_out))

    nt_file_name = make_nt(aa_file)
    gene_path = os.path.join(nt_input, nt_file_name)

    references, candidates, raw_references = parse_fasta(gene_path)
    nt_out_path = os.path.join(output, "nt", nt_file_name)
    with open(nt_out_path, "w", encoding="UTF-8") as nt_out:
        nt_out.write(raw_references)
        for header, sequence in candidates:
            gene = header.split("|")[0].replace(">", "")

            kick, cull_start, cull_end = follow_through[gene][header]

            if not kick:
                cull_start_adjusted = cull_start * 3
                cull_end_adjusted = cull_end * 3

                out_line = ("-" * cull_start_adjusted) + sequence[cull_start_adjusted:cull_end_adjusted]

                characters_till_end = len(sequence) - len(out_line)
                out_line += (
                    "-" * characters_till_end
                )  # Add dashes till reached input distance

                nt_out.write(">" + header + "\n")
                nt_out.write(out_line + "\n")

    return log

def do_folder(folder, args: MainArgs):
    folder_time = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Processing: {folder}", args.verbose, 0)
    aa_path = os.path.join(folder, args.amino_acid)
    nt_path = os.path.join(folder, args.nucleotide)
    output_path = os.path.join(folder, args.output)
    if not os.path.exists(aa_path) or not os.path.exists(nt_path):
        printv(f"WARNING: Can't find aa ({aa_path}) and nt ({nt_path}) folders. Abort", args.verbose)
        return
    available_tmp_path = folder_check(output_path, folder)
    file_inputs = [
        input_gene
        for input_gene in os.listdir(aa_path)
        if ".aa" in input_gene
    ]
    file_inputs.sort(
        key=lambda x : os.path.getsize(os.path.join(aa_path, x)), reverse=True
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
                )
            )

        with Pool(args.processes) as pool:
            log_components = pool.starmap(do_gene, arguments, chunksize=1)
    else:
        log_components = [do_gene(
                            aa_path,
                            nt_path,
                            output_path,
                            args.matches,
                            input_gene,
                            args.match_percent,
                            args.debug,
                            args.base_pair,
                            args.verbose,
                        ) for input_gene in file_inputs]

    if args.debug:
        log_global = []

        for component in log_components:
            log_global.extend(component)

        log_global.sort()
        log_global.insert(0, "Gene,Header,Cull To Start,Cull To End,Data Length,Data Removed\n")
        log_out = os.path.join(output_path, 'Culls.csv')
        with open(log_out,'w') as fp:
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

