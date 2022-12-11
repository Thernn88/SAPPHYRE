from __future__ import annotations

import os
from collections import namedtuple
from multiprocessing.pool import ThreadPool
from tempfile import TemporaryDirectory, NamedTemporaryFile
from threading import Lock
from time import time

from .utils import printv, gettempdir
from .timekeeper import TimeKeeper, KeeperMode

MAFFT_FOLDER = "mafft"
AA_FOLDER = "aa"
ALN_FOLDER = "aln"


def process_genefile(filewrite, fileread):
    ref_og_hashmap = {}
    cand_og_hashmap = {}
    with open(fileread) as fp:
        for line in fp:
            if line[0] == ">":
                if line[-2] == ".":
                    ref_og_hashmap[line.split("|")[2]] = line
                    next(fp)  # Skip this line and next line
                else:
                    cand_og_hashmap[line[:244]] = line
                    filewrite.write(line)
            else:
                filewrite.write(line)
    return ref_og_hashmap, cand_og_hashmap


CmdArgs = namedtuple(
    "CmdArgs", ["string", "gene_file", "result_file", "gene", "lock", "verbose"]
)


def run_command(args: CmdArgs) -> None:
    if args.lock is not None:
        with args.lock:
            printv(f"Doing: {args.gene}", args.verbose, 2)
    else:
        printv(f"Doing: {args.gene}", args.verbose, 2)

    with TemporaryDirectory(dir=gettempdir()) as tmpdir, NamedTemporaryFile(mode="w+", dir=tmpdir) as tmpfile:
        ref_og_hashmap, cand_og_hashmap = process_genefile(tmpfile, args.gene_file)
        tmpfile.file.flush()
        command = args.string.format(
            tmpfile=tmpfile.name,
            resultfile=args.result_file,
            gene=args.gene
        )
        printv(f"Executing command: {command}", args.verbose, 2)
        os.system(command)

    # Overwrite reference headers with original headers
    with open(args.result_file, "r+") as fp_out:
        out = []
        for line in fp_out:
            if line[0] == ">":
                if "|" in line:
                    out.append(cand_og_hashmap[line[:244]])
                else:
                    out.append(ref_og_hashmap[line[1:].strip()])
            else:
                out.append(line)

        fp_out.seek(0)
        fp_out.writelines(out)
        fp_out.truncate()


def do_folder(folder, args):
    printv(f"Processing: {os.path.basename(folder)}", args.verbose)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    mafft_path = os.path.join(folder, MAFFT_FOLDER)
    aa_path = os.path.join(folder, AA_FOLDER)
    if not os.path.exists(aa_path):
        printv(f"ERROR: Can't find aa ({aa_path}) folder. Abort", args.verbose, 0)
        return
    os.makedirs(mafft_path, exist_ok=True)

    genes = [gene.split(".")[0] for gene in os.listdir(aa_path) if ".aa" in gene]
    #genes.sort(key = lambda x : os.path.getsize(os.path.join(aa_path, x + ".aa.fa")), reverse=True)
    orthoset_path = os.path.join(args.orthoset_input, args.orthoset)
    aln_path = os.path.join(orthoset_path, ALN_FOLDER)
    if not os.path.exists(orthoset_path):
        printv('ERROR: Orthoset path not found.', args.verbose, 0)
        return False
    if not os.path.exists(aln_path):
        printv("ERROR: Aln folder not found.", args.verbose, 0)
        return False
    cmd = "mafft"
    if args.linsi:
        cmd = "mafft-linsi"
    command = (
        "%s --anysymbol --quiet --addfragments {tmpfile} --thread -1 %s/{gene}.aln.fa > {resultfile}"
        % (cmd, aln_path)
    )

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

    printv(f"Done! Took {time_keeper.differential():.2f}s", args.verbose)
    return True


def main(args):
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    for folder in args.INPUT:
        success = do_folder(folder, args)
        if not success:
            return False
    
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {time_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr mafft"
    )
