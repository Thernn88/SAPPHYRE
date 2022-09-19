"""
Aligns references and runs mafft align on all genes of each taxa within a folder of taxa folders.

PyLint 8.38/10
"""

import argparse
import os
from multiprocessing.pool import ThreadPool
from threading import Lock


def run_command(arg_tuple: tuple) -> None:
    """
    Aligns the references and calls mafft on the target gene
    """
    string, gene, taxa, parent, aa_path, nt_path, lock = arg_tuple

    align_references(aa_path, nt_path)
    with lock:
        print(gene, aa_path, nt_path)
    COMMAND = string.format(gene, taxa, parent)
    os.system(COMMAND)


def align_references(aa_path: str, nt_path: str) -> None:
    """
    Aligns the nt references with the same order as the aa references
    """
    order = []
    with open(aa_path, encoding = "UTF-8") as aa_file:
        for line in aa_file:
            if line != "":
                if ">" in line:
                    fields = line.split("|")
                    if len(fields) <= 3:
                        order.append(fields[1]) # Save order of reference taxa name

                    else:
                        break

    nt_references = {}
    nt_out_lines = []

    with open(nt_path, encoding = "UTF-8") as nt_file:
        content = nt_file.read()
        lines = content.split("\n")

        for i in range(0, len(lines), 2):
            if lines[i] != "":
                header = lines[i]
                seq = lines[i + 1]

                fields = header.split("|")
                if len(fields) > 3:
                    nt_out_lines.append(header)
                    nt_out_lines.append(seq)

                else:
                    nt_references[fields[1]] = (header, seq)

    to_add = []

    for taxa_name in order:
        header, seq = nt_references[taxa_name]
        to_add.append(header)
        to_add.append(seq)

    nt_out_lines = to_add + nt_out_lines

    with open(nt_path, "w", encoding = "UTF-8") as nt_out:
        nt_out.write("\n".join(nt_out_lines))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, default="Parent", help="Parent input path."
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=0,
        help="Number of threads used to call processes.",
    )
    args = parser.parse_args()

    COMMAND = "perl pal2nal.mod.pl -output fasta {2}/{1}/mafft/{0}.aa.fa {2}/{1}/nt/{0}.nt.fa > {2}/{1}/nt_aligned/{0}.nt.fa"

    for taxa in os.listdir(args.input):
        print(f"Doing taxa {taxa}")
        taxa_path = os.path.join(args.input, taxa)
        nt_aligned_path = os.path.join(taxa_path, "nt_aligned")
        mafft_path = os.path.join(taxa_path, "mafft")

        if os.path.exists(mafft_path):
            if not os.path.exists(nt_aligned_path):
                os.mkdir(nt_aligned_path)

            genes = [gene.split(".")[0] for gene in os.listdir(mafft_path) if ".aa" in gene]

            if args.processes:
                arguments = []
                lock = Lock()
                for gene in genes:

                    aa_path = os.path.join(args.input, taxa, "mafft", gene + ".aa.fa")
                    nt_path = os.path.join(args.input, taxa, "nt", gene + ".nt.fa")

                    arguments.append(
                        (COMMAND, gene, taxa, args.input, aa_path, nt_path, lock)
                    )
                with ThreadPool(args.processes) as pool:
                    pool.map(run_command, arguments, chunksize=1)
            else:
                for gene in genes:
                    print(gene)

                    aa_path = os.path.join(args.input, taxa, "mafft", gene + ".aa.fa")
                    nt_path = os.path.join(args.input, taxa, "nt", gene + ".nt.fa")

                    align_references(aa_path, nt_path)
                    os.system(COMMAND.format(gene, taxa, args.input))
        else:
            print("Can't find mafft folder for taxa {}".format(taxa))
