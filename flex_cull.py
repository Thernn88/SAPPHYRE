"""
FlexCull Description Goes Here

PyLint 9.81/10
"""
import argparse
import os
from multiprocessing.pool import Pool
from time import time


def folder_check(output_target_path: str, input_target_path: str) -> str:
    """
    Checks to see if input and output directory has the necessary
    folder structure

    if not, create missing folders.

    Finally, checks if linux memory temp folder exists
    else creates tmp folder in the input

    Returns path to tmp folder
    """
    if not os.path.exists(output_target_path):
        os.mkdir(output_target_path)

    output_aa_path = os.path.join(output_target_path, "aa")
    if not os.path.exists(output_aa_path):
        os.mkdir(output_aa_path)

    output_nt_path = os.path.join(output_target_path, "nt")
    if not os.path.exists(output_nt_path):
        os.mkdir(output_nt_path)

    target_tmp_path = "/run/shm"

    if not os.path.exists(target_tmp_path):
        target_tmp_path = '/dev/shm'

    if not os.path.exists(target_tmp_path):
        target_tmp_path = os.path.join(input_target_path, "tmp")
        if not os.path.exists(target_tmp_path):
            os.mkdir(target_tmp_path)

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
    is_reference_header = lambda header: header.count("|") == 2

    lines = []
    with open(fasta_path, encoding="UTF-8") as fasta_io:
        lines = list(fasta_io)

    for i in range(0, len(lines), 2):
        header = lines[i]
        seq = lines[i + 1]

        if is_reference_header(header):
            references.append((header.rstrip(), seq.rstrip()))

            raw_references.extend([header + "\n", seq + "\n"])

        else:
            candidates.append((header.strip(), seq.rstrip()))

    return references, candidates, raw_references


def main(
    aa_input: str,
    nt_input: str,
    output: str,
    amt_matches: int,
    aa_file: str,
    tmp_path: str,
    debug: bool,
) -> None:
    """
    FlexCull main function. Culls input aa and nt using specified amount of matches
    """
    gene_path = os.path.join(aa_input, aa_file)

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
                character_at_each_pos[i] = []
            character_at_each_pos[i].append(char)

            #  if char isnt a hyphen, this positions can't be all dashes
            if char != "-":
                all_dashes_by_index[i] = False

    if debug:
        log = ["Gene,Header,Cull To Start,Cull To End,Data Length,Data Removed\n"]

    out_lines = raw_references.copy()
    follow_through = {}
    offset = amt_matches - 1
    print(len(candidates))
    for header, sequence in candidates:
        gene = header.split("|")[0].replace(">", "")

        if gene not in follow_through:
            follow_through[gene] = {}

        cull_start = None
        kick = False
        data_removed = 0

        data_length = 0

        for i, char in enumerate(sequence):
            all_dashes_at_position = not all_dashes_by_index[i]
            if i == len(sequence) - offset:
                kick = True
                break

            if all_dashes_at_position and char != "-":
                # Don't allow cull to point of all dashes
                this_matches = char in character_at_each_pos[i]
            else:
                this_matches = False

            matches = [this_matches]
            for match_i in range(1, amt_matches):
                character_matches_ref = (
                    sequence[i + match_i] in character_at_each_pos[i + match_i]
                )
                matches.append(character_matches_ref)

            if sum(matches) == len(matches):
                # If all points checked match
                cull_start = i
                break

        if not kick:
            # If not kicked from Cull Start Calc. Continue
            cull_end = None
            for i_raw, char in enumerate(sequence):
                all_dashes_at_position = not all_dashes_by_index[i]
                i = len(sequence) - 1 - i_raw  # Start from end
                if i < cull_start + offset:
                    kick = True
                    break

                if all_dashes_at_position and char != "-":
                    # Don't allow cull to point of all dashes
                    this_matches = char in character_at_each_pos[i]
                else:
                    this_matches = False

                matches = [this_matches]
                for match_i in range(1, amt_matches):
                    character_matches_ref = (
                        sequence[i - match_i] in character_at_each_pos[i - match_i]
                    )

                    matches.append(character_matches_ref)

                if sum(matches) == len(matches):
                    # If all points checked match
                    cull_end = i + 1
                    break

        if not kick:  # If also passed Cull End Calc. Finish
            out_line = "-" * cull_start  # Cull start
            out_line += sequence[cull_start:cull_end]  # Add Culled Sequence

            characters_till_end = len(sequence) - len(out_line)
            out_line += (
                "-" * characters_till_end
            )  # Add dashes till reached input distance

            # The cull replaces data positions with dashes to maintain the same alignment
            # while removing the bad data

            removed_section = sequence[:cull_start] + sequence[cull_end:]
            data_removed = len(removed_section) - removed_section.count("-")

            data_length = cull_end - cull_start
            bp_after_cull = len(out_line) - out_line.count("-")

            if bp_after_cull >= args.bp:
                follow_through[gene][header] = False, cull_start, cull_end

                out_lines.append(header + "\n")
                out_lines.append(out_line + "\n")

                if debug:
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

    aa_out_path = os.path.join(output, "aa", aa_file)
    with open(aa_out_path, "w", encoding="UTF-8") as aa_out:
        aa_out.writelines(out_lines)

    this_gene = aa_file.split(".")[0]

    if debug:
        log_out_path = os.path.join(tmp_path, this_gene + ".csv")
        with open(log_out_path, "w", encoding="UTF-8") as log_out:
            log_out.writelines(log)

    nt_file_name = make_nt(aa_file)
    gene_path = os.path.join(nt_input, nt_file_name)
    # gene_content = open(gene_path).read()

    references, candidates, raw_references = parse_fasta(gene_path)
    out_lines = raw_references.copy()
    for header, sequence in candidates:
        gene = header.split("|")[0].replace(">", "")

        kick, cull_start, cull_end = follow_through[gene][header]

        if not kick:
            cull_start_adjusted = cull_start * 3
            cull_end_adjusted = cull_end * 3

            out_line = "-" * cull_start_adjusted  # Cull start
            out_line += sequence[
                cull_start_adjusted:cull_end_adjusted
            ]  # Add Culled Sequence

            characters_till_end = len(sequence) - len(out_line)
            out_line += (
                "-" * characters_till_end
            )  # Add dashes till reached input distance

            out_lines.append(header + "\n")
            out_lines.append(out_line + "\n")

    nt_out_path = os.path.join(output, "nt", nt_file_name)
    with open(nt_out_path, "w", encoding="UTF-8") as nt_out:
        nt_out.writelines(out_lines)


def consolidate(log_paths: list) -> str:
    """Consolidates each individual gene log to
    global log"""

    consolidated_log_out = []

    for i, log_path in enumerate(log_paths):
        with open(log_path, encoding="UTF-8") as log_in:
            log_lines = log_in.readlines()

        if i != 0:
            log_lines = log_lines[1:]

        consolidated_log_out.extend(log_lines)

    return consolidated_log_out


def run_command(arg_tuple: tuple) -> None:
    """
    Calls the main() function parallel in each thread
    """
    aa_input, nt_input, output, matches, aa_file, tmp_path, debug = arg_tuple
    main(aa_input, nt_input, output, matches, aa_file, tmp_path, debug)


if __name__ == "__main__":
    start = time()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, default="parent", help="Parent input path."
    )
    parser.add_argument(
        "-o", "--output", type=str, default="trimmed", help="Output Directory."
    )
    parser.add_argument(
        "-aa", "--aa", type=str, default="mafft", help="AA Folder Name."
    )
    parser.add_argument(
        "-nt", "--nt", type=str, default="nt_aligned", help="NT Folder Name."
    )
    parser.add_argument(
        "-m",
        "--matches",
        type=int,
        default=3,
        help="Amount of nucleotides that have to match reference.",
    )
    parser.add_argument(
        "-bp", "--bp", type=int, default=15, help="Minimum bp after cull."
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=2,
        help="Number of threads used to call processes.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled Output log of culls.",
    )
    args = parser.parse_args()
    allowed_extensions = ["fa", "fas", "fasta"]

    for taxa in os.listdir(args.input):
        print(f"Doing taxa {taxa}")
        aa_path = os.path.join(args.input, taxa, args.aa)
        nt_path = os.path.join(args.input, taxa, args.nt)
        output_path = os.path.join(args.input, taxa, args.output)

        if os.path.exists(aa_path) and os.path.exists(nt_path):
            available_tmp_path = folder_check(output_path, os.path.join(args.input, taxa))

            file_inputs = [input_gene for input_gene in os.listdir(aa_path) if ".aa" in input_gene]

            if args.processes:
                arguments = []
                for input_gene in file_inputs:
                    arguments.append(
                        (
                            aa_path,
                            nt_path,
                            output_path,
                            args.matches,
                            input_gene,
                            available_tmp_path,
                            args.debug,
                        )
                    )

                with Pool(args.processes) as pool:
                    pool.map(run_command, arguments, chunksize=1)
            else:
                for input_gene in file_inputs:
                    main(
                        aa_path,
                        nt_path,
                        output_path,
                        args.matches,
                        input_gene,
                        available_tmp_path,
                        args.debug,
                    )

            if args.debug:
                logs = [os.path.join(available_tmp_path, log_present) for log_present in os.listdir(available_tmp_path)]
                log_global = consolidate(logs)
                log_global.sort()
                log_out = os.path.join(output_path,'Culls.csv')
                open(log_out,'w').writelines(log_global)

            time_taken = time()
            time_taken = time_taken - start

            print(f"Finished in {time_taken} seconds")
        else:
            print(f"Can't find aa and nt folder for taxa {taxa}")
