"""
Placeholder module docstring

Aligns references and runs mafft on all genes of a taxa -i
"""
import argparse
import os
from multiprocessing.pool import ThreadPool
from threading import Lock


def run_command(arg_tuple: tuple) -> None:
    """
    Aligns the references and calls mafft on the target gene
    """
    string, parent, taxa, gene, aa_path, nt_path, lock = arg_tuple

    align_references(aa_path, nt_path)
    with lock:
        print(gene)
    COMMAND = string.format( 
        parent = parent,
        taxa = taxa,
        gene = gene
        )
    os.system(COMMAND)


def align_references(aa_path: str, nt_path: str) -> None:
    """
    Aligns the nt references with the same order as the aa references
    """
    order = []
    with open(aa_path, encoding = "UTF-8") as aa_file:
        for line in aa_file:
            line = line.strip()
            if line != "":
                if ">" in line:
                    if line[-1] == '.':
                        fields = line.split("|")
                        order.append(fields[1]) # Save order of reference taxa name
                    else:
                        break

    nt_references = {}
    nt_out_lines = []

    with open(nt_path, 'r+', encoding = "UTF-8") as fp:
        content = fp.read()
        lines = content.split("\n")

        end_of_references = False
        for i in range(0, len(lines), 2):
            if lines[i] != "":
                header = lines[i]
                seq = lines[i + 1]

                if end_of_references is False:
                    if header[-1] == '.':
                        fields = header.split("|")
                        nt_references[fields[1]] = (header, seq)
                    else:
                        end_of_references = True

                if end_of_references is True:
                    nt_out_lines.append(header+'\n')
                    nt_out_lines.append(seq+'\n')

        

        to_add = []

        for taxa_name in order:
            header, seq = nt_references[taxa_name]
            to_add.append(header+'\n')
            to_add.append(seq+'\n')

        nt_out_lines = to_add + nt_out_lines

        fp.seek(0)
        fp.writelines(nt_out_lines)
        fp.truncate()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', 
        '--input', 
        type=str, 
        default='Parent',
        help='Parent input path.'
    )

    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=0,
        help="Number of threads used to call processes.",
    )
    args = parser.parse_args()

    for taxa in os.listdir(args.input):
        taxa_path = os.path.join(args.input, taxa)
        if os.path.isdir(taxa_path):
            mafft_path = os.path.join(taxa_path, "mafft")
            if os.path.isdir(mafft_path):
                nt_aligned_path = os.path.join(taxa_path, "nt_aligned")
                if not os.path.exists(nt_aligned_path):
                    os.mkdir(nt_aligned_path)

                genes = [gene.split(".")[0] for gene in os.listdir(mafft_path) if ".aa" in gene]

                COMMAND = "perl pal2nal.mod.pl -output fasta {parent}/{taxa}/mafft/{gene}.aa.fa {parent}/{taxa}/nt/{gene}.nt.fa > {parent}/{taxa}/nt_aligned/{gene}.nt.fa"

                if args.processes:
                    arguments = []
                    lock = Lock()
                    for target_gene in genes:
                        target_aa_path = os.path.join(mafft_path, target_gene + ".aa.fa")
                        target_nt_path = os.path.join(nt_aligned_path, target_gene + ".nt.fa")

                        arguments.append((COMMAND, args.input, taxa, target_gene, target_aa_path, target_nt_path, lock))
                    with ThreadPool(args.processes) as pool:
                        pool.map(run_command, arguments, chunksize=1)
                else:
                    for target_gene in genes:
                        print(target_gene)
                        target_aa_path = os.path.join(mafft_path, target_gene + ".aa.fa")
                        target_nt_path = os.path.join(nt_aligned_path, target_gene + ".nt.fa")

                        align_references(target_aa_path, target_nt_path)
                        os.system(
                            COMMAND.format( 
                                parent = args.input,
                                taxa = taxa,
                                gene = target_gene
                                )
                            )

if __name__ == "__main__":
    main()