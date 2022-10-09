import argparse
import os
from multiprocessing.pool import ThreadPool
from threading import Lock
from time import time

def run_command(arg_tuple: tuple) -> None:
    string, gene_file, result_file, gene, temp_folder, lock = arg_tuple
    if lock is not None:
        with lock:
            print(gene)
    else:
        print(gene)

    ref_og_hashmap = {}

    tmp_gene_file = os.path.join(temp_folder,gene+'.tmp')
    with open(tmp_gene_file, "w") as fp_out, open(gene_file) as fp_in:
        for line in fp_in:
            if line[0] == ">" and line[-2] == ".":
                ref_og_hashmap[line.split("|")[2]] = line
                next(fp_in)  # Skip this line and next line
            else:
                fp_out.write(line)

    command = string.format(tmp_gene_file, result_file, gene)
    os.system(command)

    os.remove(tmp_gene_file)

    # Overwrite reference headers with original headers
    non_ref_found = False
    with open(result_file, "r+") as fp_out:
        out = []
        for line in fp_out:
            if line[0] == ">" and not non_ref_found:
                if "|" in line:
                    non_ref_found = True
                    out.append(line)
                else:
                    out.append(ref_og_hashmap[line[1:].strip()])
            else:
                out.append(line)

        fp_out.seek(0)
        fp_out.writelines(out)
        fp_out.truncate()


parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input", type=str, default="Parent", help="Path to parent folder for input"
)
parser.add_argument(
    "-p",
    "--processes",
    type=int,
    default=0,
    help="Number of threads used to call processes.",
)
parser.add_argument(
    "-oi",
    "--orthoset_input",
    type=str,
    default="orthosets",
    help="Path to directory of Orthosets folder",
)
parser.add_argument(
    "-o", "--orthoset", type=str, default="Ortholog_set_Mecopterida_v4", help="Orthoset"
)
args = parser.parse_args()

mafft_folder = "mafft"
aa_folder = "aa"
aln_folder = "aln"

aln_path = os.path.join(args.orthoset_input, args.orthoset, aln_folder)

if os.path.exists("/run/shm"):
    tmp_dir = "/run/shm"
elif os.path.exists("/dev/shm"):
    tmp_dir = "/dev/shm"
else:
    tmp_dir = args.input

temp_folder = os.path.join(tmp_dir, "tmp")
delete_on_exit = False
if not os.path.exists(temp_folder):
    delete_on_exit = True
    os.makedirs(temp_folder, exist_ok=True)


for taxa in os.listdir(args.input):
    start = time()
    print("Doing taxa {}".format(taxa))
    mafft_path = os.path.join(args.input, taxa, mafft_folder)
    aa_path = os.path.join(args.input, taxa, aa_folder)
    if os.path.exists(aa_path):
        if not os.path.exists(mafft_path):
            os.mkdir(mafft_path)

        genes = [gene.split(".")[0] for gene in os.listdir(aa_path) if ".aa" in gene]

        # command = 'mafft --anysymbol --auto --quiet --thread -1  --addfragments {0} --thread -1 '+aln_path+'/{2}.aln.fa > {1}'
        command = (
            "mafft-linsi --anysymbol --quiet --linelength -1 --addfragments {0} --thread -1 "
            + aln_path
            + "/{2}.aln.fa > {1}"
        )

        if args.processes:
            arguments = list()
            lock = Lock()
            for gene in genes:
                gene_file = os.path.join(aa_path, gene + ".aa.fa")
                result_file = os.path.join(mafft_path, gene + ".aa.fa")
                arguments.append((command, gene_file, result_file, gene, temp_folder, lock))
            with ThreadPool(args.processes) as pool:
                pool.map(run_command, arguments, chunksize=1)
        else:
            for gene in genes:
                gene_file = os.path.join(aa_path, gene + ".aa.fa")
                result_file = os.path.join(mafft_path, gene + ".aa.fa")
                run_command((command, gene_file, result_file, gene, temp_folder, None))

        print("Took {:.2f}s".format(time() - start))
    else:
        print("Can't find aa folder for taxa {}".format(taxa))

if delete_on_exit:
    os.remove(temp_folder)