from __future__ import annotations

import os
from collections import namedtuple
from multiprocessing.pool import ThreadPool
from shutil import rmtree
from tempfile import TemporaryDirectory, NamedTemporaryFile
from threading import Lock

from .utils import printv, gettempdir, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode

MAFFT_FOLDER = "mafft"
AA_FOLDER = "aa"
ALN_FOLDER = "aln"


def process_genefile(filewrite, fileread):
    ref_og_hashmap = {}
    cand_og_hashmap = {}
    for header, sequence in parseFasta(fileread):

        if header.endswith("."):
            ref_og_hashmap[header.split("|")[2]] = header
        else:
             cand_og_hashmap[header[:242]] = header
             filewrite.write(f">{header}\n")
             filewrite.write(sequence+"\n")
    return ref_og_hashmap, cand_og_hashmap


CmdArgs = namedtuple(
    "CmdArgs",
    ["string", "gene_file", "result_file", "gene", "lock", "verbose", "compress"],
)


def run_command(args: CmdArgs) -> None:
    if args.lock is not None:
        with args.lock:
            printv(f"Doing: {args.gene}", args.verbose, 2)
    else:
        printv(f"Doing: {args.gene}", args.verbose, 2)

    with TemporaryDirectory(dir=gettempdir()) as tmpdir, NamedTemporaryFile(
        mode="w+", dir=tmpdir
    ) as tmpfile:
        ref_og_hashmap, cand_og_hashmap = process_genefile(tmpfile, args.gene_file)
        tmpfile.file.flush()
        command = args.string.format(
            tmpfile=tmpfile.name, resultfile=args.result_file, gene=args.gene
        )
        printv(f"Executing command: {command}", args.verbose, 2)
        os.system(command)

    # Overwrite reference headers with original headers
    out = []
    for header, sequence in parseFasta(args.result_file):
        if "|" in header:
            out.append((cand_og_hashmap[header[:242]],sequence))

        else:
            out.append((ref_og_hashmap[header.strip()], sequence))

    if args.compress:
        os.unlink(args.result_file)
    writeFasta(args.result_file, out, args.compress)


def do_folder(folder, args):
    printv(f"Processing: {os.path.basename(folder)}", args.verbose)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    mafft_path = os.path.join(folder, MAFFT_FOLDER)
    aa_path = os.path.join(folder, AA_FOLDER)
    if not os.path.exists(aa_path):
        printv(f"ERROR: Can't find aa ({aa_path}) folder. Abort", args.verbose, 0)
        return
    rmtree(mafft_path, ignore_errors=True)
    os.mkdir(mafft_path)

    genes = [
        gene
        for gene in os.listdir(aa_path)
        if gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]
    # genes.sort(key = lambda x : os.path.getsize(os.path.join(aa_path, x + ".aa.fa")), reverse=True)
    orthoset_path = os.path.join(args.orthoset_input, args.orthoset)
    aln_path = os.path.join(orthoset_path, ALN_FOLDER)
    if not os.path.exists(orthoset_path):
        printv("ERROR: Orthoset path not found.", args.verbose, 0)
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

    for file in genes:
        gene = file.split(".")[0]
        gene_file = os.path.join(aa_path, file)
        result_file = os.path.join(mafft_path, file.rstrip(".gz"))
        func(
            CmdArgs(
                command, gene_file, result_file, gene, lock, args.verbose, args.compress
            )
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
    raise Exception("Cannot be called directly, please use the module:\nphymmr mafft")
