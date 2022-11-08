from __future__ import annotations

import os
from collections import namedtuple
from multiprocessing.pool import ThreadPool
from tempfile import TemporaryDirectory, NamedTemporaryFile
from threading import Lock
from time import time

from .utils import printv, gettempdir

MAFFT_FOLDER = "mafft"
AA_FOLDER = "aa"
ALN_FOLDER = "aln"


def process_genefile(filewrite, fileread):
    ref_og_hashmap = {}
    with open(fileread) as fp:
        for line in fp:
            if line[0] == ">" and line[-2] == ".":
                ref_og_hashmap[line.split("|")[2]] = line
                next(fp)  # Skip this line and next line
            else:
                filewrite.write(line)
    return ref_og_hashmap


CmdArgs = namedtuple(
    "CmdArgs", ["string", "gene_file", "result_file", "gene", "lock", "verbose"]
)


def run_command(args: CmdArgs) -> None:
    if args.lock is not None:
        with args.lock:
            printv(args.gene, args.verbose)
    else:
        printv(args.gene, args.verbose)

    with TemporaryDirectory(dir=gettempdir()) as tmpdir, NamedTemporaryFile(mode="w+", dir=tmpdir) as tmpfile:
        ref_og_hashmap = process_genefile(tmpfile, args.gene_file)
        tmpfile.file.flush()
        command = args.string.format(
            tmpfile=tmpfile.name,
            resultfile=args.result_file,
            gene=args.gene
        )
        printv(f"Executing command: {command}", args.verbose, 2)
        os.system(command)

    # Overwrite reference headers with original headers
    non_ref_found = False
    with open(args.result_file, "r+") as fp_out:
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


def do_folder(folder, args):
    print(f"### Processing folder {folder}")
    start = time()
    mafft_path = os.path.join(folder, MAFFT_FOLDER)
    aa_path = os.path.join(folder, AA_FOLDER)
    if not os.path.exists(aa_path):
        print(f"Can't find aa ({aa_path}) folder. Abort")
        return
    os.makedirs(mafft_path, exist_ok=True)

    genes = [gene.split(".")[0] for gene in os.listdir(aa_path) if ".aa" in gene]
    aln_path = os.path.join(args.orthoset_input, args.orthoset, ALN_FOLDER)
    # command = 'mafft --anysymbol --auto --quiet --thread -1  --addfragments {0} --thread -1 '+aln_path+'/{2}.aln.fa > {1}'
    command = "mafft --anysymbol --quiet --linelength -1 --addfragments {tmpfile} --thread -1 %s/{gene}.aln.fa > {resultfile}" % aln_path

    if args.processes > 1:
        arguments = []
        func = arguments.append
        lock = Lock()
    else:
        func = run_command
        lock = None

    for gene in genes:
        gene_file = os.path.join(aa_path, gene + ".aa.fa")
        result_file = os.path.join(mafft_path, gene + ".aa.fa")
        func(
            CmdArgs(command, gene_file, result_file, gene, lock, args.verbose)
        )

    if args.processes > 1:
        with ThreadPool(args.processes) as pool:
            pool.map(run_command, arguments, chunksize=1)

    print("Took {:.2f}s".format(time() - start))


def main(args):
    if not all(os.path.exists(i) for i in args.INPUT):
        print("ERROR: All folders passed as argument must exists.")
        return False
    for folder in args.INPUT:
        do_folder(folder, args)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr mafft"
    )
