import argparse
import os
from multiprocessing.pool import Pool
from time import time

def folder_check(output_path, input_path): 
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    output_aa_path = os.path.join(output_path, "aa")
    if not os.path.exists(output_aa_path):
        os.mkdir(output_aa_path)

    output_nt_path = os.path.join(output_path, "nt")
    if not os.path.exists(output_nt_path):
        os.mkdir(output_nt_path)

    tmp_path = "/run/shm"

    if not os.path.exists(tmp_path):
        tmp_path = os.path.join(input_path, "tmp")
        if not os.path.exists(tmp_path):
            os.mkdir(tmp_path)

    return tmp_path


def deinterleave(fasta_lines: list) -> str:
    result = []
    this_out = list()
    for line in fasta_lines:
        if line[0] == ">":
            if this_out:
                result.append("".join(this_out))
            result.append(line)
            this_out = list()
        else:
            this_out.append(line.strip())
    if this_out != list():
        result.append("".join(this_out))
    return result


def make_nt(aa_file_name: str) -> str:
    """Converts AA file name to NT file name"""

    file_parts = aa_file_name.split(".")
    file_parts[file_parts.index("aa")] = "nt"
    file_name = ".".join(file_parts)

    return file_name


def is_reference_header(header: str) -> bool:
    """'Returns true if header has only three fields"""
    fields = header.split("|")

    if len(fields) == 3:
        return True
    else:
        return False


def parse_fasta(fasta_path: str) -> tuple:
    """Parses fasta file into header and sequences"""
    references = []
    candidates = []
    raw_references = []

    lines = list()
    with open(fasta_path) as f:
        lines = [line for line in f]

    # lines = fasta_content.split('\n')
    # while '' in lines:
    #     lines.remove('')

    lines = deinterleave(lines)

    for i in range(0, len(lines), 2):
        header = lines[i]
        seq = lines[i + 1]

        if is_reference_header(header):
            references.append((header.rstrip(), seq.rstrip()))

            # raw_references.append(header+'\n')
            raw_references.append(header)
            # raw_references.append(seq+'\n')
            raw_references.append(seq + "\n")

        else:
            candidates.append((header.strip(), seq.rstrip()))

    return references, candidates, raw_references


def main(aa_input, nt_input, output, amt_matches, aa_file, tmp_path):
    gene_path = os.path.join(aa_input, aa_file)
    # gene_content = open(gene_path).read()

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

    log = ["Gene,Header,Cull To Start,Cull To End,Data Length,Data Removed\n"]
    out_lines = raw_references.copy()
    follow_through = {}
    offset = amt_matches - 1
    print(len(candidates))
    for header, sequence in candidates:
        gene, taxa, taxa_id, seq_id = header.split("|")
        new_gene = gene[1:]
        gene = gene.replace(">", "")
        if gene not in follow_through:
            follow_through[gene] = {}

        cull_start = None
        kick = False
        data_removed = 0
        data_length = 0

        for i, char in enumerate(sequence):
            # all_dashes_at_position = character_at_each_pos[i].count('-') != len(character_at_each_pos[i])
            all_dashes_at_position = not all_dashes_by_index[i]
            if i == len(sequence) - offset:
                kick = True
                break

            if (
                all_dashes_at_position and sequence[i] != "-"
            ):  # Don't allow cull to point of all dashes
                this_matches = sequence[i] in character_at_each_pos[i]
            else:
                this_matches = False

            matches = [this_matches]
            for x in range(1, amt_matches):
                character_matches_ref = sequence[i + x] in character_at_each_pos[i + x]
                matches.append(character_matches_ref)
                ### this should probably be a boolean check###
            if sum(matches) == len(matches):  # If all points checked match
                cull_start = i
                break

        if not kick:  # If not kicked from Cull Start Calc. Continue
            cull_end = None
            for i_raw, char in enumerate(sequence):
                #  this boolean is the opposite of what the name implies, so negate the array value
                # all_dashes_at_position = character_at_each_pos[i].count('-') != len(character_at_each_pos[i])
                all_dashes_at_position = not all_dashes_by_index[i]
                i = len(sequence) - 1 - i_raw  # Start from end
                if i < cull_start + offset:
                    kick = True
                    break

                # if sequence[i] != '-' and all_dashes_at_position: #Don't allow cull to point of all dashes
                if (
                    all_dashes_at_position and sequence[i] != "-"
                ):  # Don't allow cull to point of all dashes
                    this_matches = sequence[i] in character_at_each_pos[i]
                else:
                    this_matches = False

                matches = [this_matches]
                for x in range(1, amt_matches):
                    character_matches_ref = (
                        sequence[i - x] in character_at_each_pos[i - x]
                    )
                    matches.append(character_matches_ref)

                if sum(matches) == len(matches):  # If all points checked match
                    cull_end = i + 1
                    break

        if not kick:  # If also passed Cull End Calc. Finish
            out_line = "-" * cull_start  # Cull start
            out_line += sequence[cull_start:cull_end]  # Add Culled Sequence

            characters_till_end = len(sequence) - len(out_line)
            out_line += (
                "-" * characters_till_end
            )  # Add dashes till reached input distance

            """
            The cull replaces data positions with dashes to maintain the same alignment
            but remove the bad data
            """

            removed_section = sequence[:cull_start] + sequence[cull_end:]
            data_removed = len(removed_section) - removed_section.count("-")

            data_length = cull_end - cull_start
            bp_after_cull = len(out_line) - out_line.count("-")

            if bp_after_cull >= args.bp:
                follow_through[gene][header] = False, cull_start, cull_end

                out_lines.append(header + "\n")
                out_lines.append(out_line + "\n")

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
            log.append(gene + "," + header + ",Kicked,Zero Data After Cull,0,\n")

    aa_out = os.path.join(output, "aa", aa_file)
    open(aa_out, "w").writelines(out_lines)

    this_gene = aa_file.split(".")[0]
    log_out = os.path.join(tmp_path, this_gene + ".csv")
    open(log_out, "w").writelines(log)

    nt_file_name = make_nt(aa_file)
    gene_path = os.path.join(nt_input, nt_file_name)
    # gene_content = open(gene_path).read()

    references, candidates, raw_references = parse_fasta(gene_path)
    out_lines = raw_references.copy()
    for header, sequence in candidates:
        gene, taxa, taxa_id, seq_id = header.split("|")
        gene = gene.replace(">", "")
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

    nt_out = os.path.join(output, "nt", nt_file_name)
    open(nt_out, "w").writelines(out_lines)


def consolidate(log_paths: list) -> str:
    """Consolidates each individual gene log to
    global log"""

    consolidated_log_out = []

    for i, log_path in enumerate(log_paths):
        log_lines = open(log_path).readlines()
        if i != 0:
            log_lines = log_lines[1:]

        consolidated_log_out.extend(log_lines)

    return consolidated_log_out


def run_command(arg_tuple: tuple) -> None:
    aa_input, nt_input, output, matches, aa_file, tmp_path = arg_tuple
    main(aa_input, nt_input, output, matches, aa_file, tmp_path)


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
    args = parser.parse_args()
    allowed_extensions = ["fa", "fas", "fasta"]

    for taxa in os.listdir(args.input):
        print("Doing taxa {}".format(taxa))
        aa_path = os.path.join(args.input, taxa, args.aa)
        nt_path = os.path.join(args.input, taxa, args.nt)
        output_path = os.path.join(args.input, taxa, args.output)

        if os.path.exists(aa_path) and os.path.exists(nt_path):
            tmp_path = folder_check(output_path, os.path.join(args.input, taxa))

            file_inputs = [gene for gene in os.listdir(aa_path) if ".aa" in gene]

            if args.processes:
                arguments = list()
                for gene in file_inputs:
                    arguments.append(
                        (aa_path, nt_path, output_path, args.matches, gene, tmp_path)
                    )

                with Pool(args.processes) as pool:
                    pool.map(run_command, arguments, chunksize=1)
            else:
                for gene in file_inputs:
                    main(aa_path, nt_path, output_path, args.matches, gene, tmp_path)

            # logs = [os.path.join(tmp_path,log) for log in os.listdir(tmp_path)]
            # log_global = consolidate(logs)
            # log_global.sort()
            # log_out = os.path.join(output_path,'Culls.csv')
            # time1 = time()
            # with open(log_out,'w') as f:
            #     f.writelines(log_global)
            # time2 = time()
            # print(f'writelines: {time2-time1}')
            # time1 = time()
            # with open(log_out, 'w', newline="") as f:
            #     csv_write = csv.writer(f)
            #     csv_write.writerows(log_global)
            # time2 = time()
            # print(f'csv: {time2-time1}')
            time_taken = time()
            time_taken = time_taken - start

            print("Finished in {} seconds".format(time_taken))
            # open('runtime.txt','w').writelines('Finished in {} seconds'.format(round(time_taken)))
        else:
            print("Can't find aa and nt folder for taxa {}".format(taxa))
