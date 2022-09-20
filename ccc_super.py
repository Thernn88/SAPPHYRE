"""
    Recursively cleans and clones AA and NT
    files found in taxa directories.

    9.71/10 PyLint
"""
import argparse
import os


def folder_check(folder: str) -> None:
    """
    Creates target directory if target directory does not exist
    """
    if not os.path.exists(folder):
        os.mkdir(folder)


def get_filenames(folder: str) -> list:
    """
    Returns all file names within a folder excluding files with prefix "." or "$"
    """
    result = [item for item in os.listdir(folder) if item[0] not in [".", "$"]]
    return result


def sequence_is_reference(header):
    """
    Returns True if the reference header identifier is present in the header
    """
    return header[-1] == "."


def clone_and_clean_files(
    aa_fastas: list, out_dir_clean_aa: str, nt_fastas: str, out_dir_clean_nt: str
) -> None:
    """
    Clones aa and nt files into clean directory with formatted headers
    """
    rev_f_tags = {}
    anti_dupe_tags = {}
    anti_dupe_counts = {}

    aa_content = {}

    for fasta in aa_fastas:
        this_gene = os.path.basename(fasta)

        with open(fasta, encoding="UTF-8") as fasta_io:
            lines = [line for line in fasta_io if line != "\n"]

            if this_gene in aa_content:
                aa_content[this_gene].extend(lines)
            else:
                aa_content[this_gene] = lines

    for fasta, lines in aa_content.items():

        aa_new_file_path = os.path.join(out_dir_clean_aa, fasta)
        gene = fasta.split(".")[0]
        if gene not in rev_f_tags:
            rev_f_tags[gene] = {}
            anti_dupe_tags[gene] = {}

        headers = [line.split("|")[2] for line in lines if line[0] == ">"]
        end_of_references = False
        with open(aa_new_file_path, "w", encoding="UTF-8") as output:
            already_written = set()
            for i in range(0, len(lines), 2):
                header = lines[i].strip()
                sequence = lines[i + 1].strip()

                if end_of_references is False:
                    if sequence_is_reference(header):
                        header_gene, taxa_name, taxa_id, identifier = header.split("|")
                        new_header = "|".join(
                            [header_gene, taxa_name, taxa_id, identifier]
                        )  # its faster to change a list
                        if new_header in already_written:
                            continue
                        already_written.add(new_header)
                    else:
                        end_of_references = True
                
                if end_of_references is True:
                    header_gene, taxa, taxa_id, node, _, frame = header.split("|")
                    if not taxa_id[-1].isnumeric() and taxa_id[-2] == "_":
                        taxa_id = taxa_id[:-2]

                    new_header = [header_gene, "|", taxa, "|", taxa_id, "|", node]

                    if headers.count(node) > 1:
                        if node not in rev_f_tags[gene]:
                            rev_f_tags[gene][node] = {}

                        if "revcomp" in frame:
                            rev_f_tags[gene][node][i] = "_R"
                            new_header.append("_R")
                        else:
                            rev_f_tags[gene][node][i] = "_F"

                            new_header.append("_F")

                        if not node in anti_dupe_tags[gene]:
                            anti_dupe_tags[gene][node] = {}

                        if node in anti_dupe_counts:
                            anti_dupe_counts[node] += 1
                            count = anti_dupe_counts[node]
                        else:
                            anti_dupe_counts[node] = 1
                            count = 1

                        header_addition = f"_{count}"
                        anti_dupe_tags[gene][node][i] = header_addition

                        new_header.append(header_addition)

                    new_header = "".join(new_header)

                    output.write(new_header)
                    output.write("\n")
                    output.write(sequence)
                    output.write("\n")

    nt_content = {}

    for fasta in nt_fastas:
        this_gene = os.path.basename(fasta)

        with open(fasta, encoding="UTF-8") as fasta_io:
            lines = [line for line in fasta_io if line != "\n"]

            if this_gene in nt_content:
                nt_content[this_gene].extend(lines)
            else:
                nt_content[this_gene] = lines

    for fasta, lines in nt_content.items():
        gene = fasta.split(".")[0]

        nt_new_file_path = os.path.join(out_dir_clean_nt, fasta)

        headers = [line.split("|")[2] for line in lines if ">" in line]

        end_of_references = False
        with open(nt_new_file_path, "w", encoding="UTF-8") as output:
            already_written = set()
            for i in range(0, len(lines), 2):
                header = lines[i].strip()
                sequence = lines[i + 1].strip()

                if end_of_references is False:
                    if sequence_is_reference(header):
                        header_gene, taxa_name, taxa_id, identifier = header.split("|")
                        new_header = "|".join(
                            [header_gene, taxa_name, taxa_id, identifier]
                        )  # Format reference header
                        if new_header in already_written:
                            continue
                        already_written.add(new_header)
                    else:
                        end_of_references = True
                
                if end_of_references is True:
                    header_gene, taxa, taxa_id, node, _, frame = header.split("|")
                    if not taxa_id[-1].isnumeric() and taxa_id[-2] == "_":
                        taxa_id = taxa_id[:-2]

                    new_header = [
                        header_gene,
                        "|",
                        taxa,
                        "|",
                        taxa_id,
                        "|",
                        node,
                    ]  # Format candidate header

                    if headers.count(node) > 1:
                        new_header.append(rev_f_tags[gene][node][i])
                        new_header.append(anti_dupe_tags[gene][node][i])
                    new_header = "".join(new_header)
                    
                output.write(new_header)
                output.write("\n")
                output.write(sequence)
                output.write("\n")


def main():
    """
    Crawls super dir input for taxa folders to clean and clone
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Super folder")
    parser.add_argument(
        "-o", "--output", required=True, help="Location of cloned directory"
    )
    args = parser.parse_args()

    folder_check(args.output)

    taxas = [
        taxa
        for taxa in os.listdir(args.input)
        if os.path.isdir(os.path.join(args.input, taxa))
    ]
    taxa_groups = {}

    for taxa in taxas:
        taxa_name = taxa.split(".")[0]
        case_a = taxa_name[-1].isnumeric() and taxa_name[-2] == "_"
        case_b = (
            taxa_name[-1].isnumeric() and taxa_name[-2] == "R" and taxa_name[-3] == "_"
        )
        case_c = (
            not taxa_name[-1].isnumeric()
            and taxa_name[-2].isnumeric()
            and taxa_name[-3] == "R"
            and taxa_name[-4] == "_"
        )
        if case_a or case_b or case_c:  # If _# or _R# or _R#X
            taxa_name = "_".join(taxa_name.split("_")[:-1]) + ".fa"

        if taxa_name not in taxa_groups:
            taxa_groups[taxa_name] = []

        taxa_groups[taxa_name].append(taxa)

    for taxa, to_combine_taxa in taxa_groups.items():
        aa_fastas = []
        nt_fasta = []

        for combine_taxa in to_combine_taxa:
            aa_in = os.path.join(os.path.join(args.input, combine_taxa), "aa")
            aa_fastas.extend([os.path.join(aa_in, aa) for aa in get_filenames(aa_in)])

            nt_in = os.path.join(os.path.join(args.input, combine_taxa), "nt")
            nt_fasta.extend([os.path.join(nt_in, nt) for nt in get_filenames(nt_in)])

        clone = os.path.join(args.output, taxa)
        folder_check(clone)

        aa_out = os.path.join(clone, "aa")
        nt_out = os.path.join(clone, "nt")

        if not os.path.exists(aa_out):
            os.mkdir(aa_out)
        if not os.path.exists(nt_out):
            os.mkdir(nt_out)

        clone_and_clean_files(aa_fastas, aa_out, nt_fasta, nt_out)


if __name__ == "__main__":
    main()
