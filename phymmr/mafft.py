from __future__ import annotations
import argparse
import gzip
from itertools import combinations

import os
from collections import namedtuple
from multiprocessing.pool import ThreadPool
from shutil import rmtree
from tempfile import TemporaryDirectory, NamedTemporaryFile
from threading import Lock
from .utils import printv, gettempdir, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode
from Bio.Seq import Seq
KMER_LEN = 21   

def find_kmers(fasta):
    kmers = {}
    for header, sequence in parseFasta(fasta):
        kmers[header] = set()
        for i in range(0, len(sequence)-KMER_LEN):
            kmer = sequence[i:i+KMER_LEN]
            if "*" not in kmer:
                kmers[header].add(kmer)
    return kmers

MAFFT_FOLDER = "mafft"
AA_FOLDER = "aa"
ALN_FOLDER = "aln"


def process_genefile(fileread):
    data = []
    targets = {}
    for header, sequence in parseFasta(fileread):
        if not header.endswith("."):
            data.append((header, sequence))
        else:
            targets[header.split("|")[2]] = header
            
    return len(data), data, targets


CmdArgs = namedtuple(
    "CmdArgs",
    ["string", "gene_file", "result_file", "gene", "lock", "verbose", "compress", "aln_path", "debug"],
)


def run_command(args: CmdArgs) -> None:
    keeper = TimeKeeper(KeeperMode.DIRECT)
    debug = args.debug
    if args.lock is not None:
        with args.lock:
            printv(f"Doing: {args.gene}", args.verbose, 2)
    else:
        printv(f"Doing: {args.gene}", args.verbose, 2)
    temp_dir = gettempdir()

    aligned_ingredients = []

    if debug:
        intermediates = "intermediates"
        this_intermediates = os.path.join(intermediates, args.gene)
        if not os.path.exists(this_intermediates):
            os.mkdir(this_intermediates)

    aln_file = os.path.join(args.aln_path, args.gene + ".aln.fa")
    to_merge = []

    with TemporaryDirectory(dir=temp_dir) as parent_tmpdir, TemporaryDirectory(dir=parent_tmpdir) as raw_files_tmp, TemporaryDirectory(dir=parent_tmpdir) as aligned_files_tmp:
        printv(f"Cleaning NT file. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
        seq_count, data, targets = process_genefile(args.gene_file)

        cluster_time = 0
        align_time = 0
        merge_time = 0

        if seq_count == 1:
            printv(f"Outputting singleton alignment. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
            aligned_file = os.path.join(aligned_files_tmp,f"{args.gene}_cluster0")
            data = [(header, str(Seq(sequence).translate())) for header, sequence in data]
            writeFasta(aligned_file, data)
            aligned_ingredients.append(aligned_file)
        else:
            printv(f"Generating Cluster. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
            clusters = []
            cluster_children = {header: [] for header, _ in data}
            with NamedTemporaryFile(mode="w+", dir=parent_tmpdir) as tmpfile:
                writeFasta(tmpfile.name, data)
                tmpfile.flush()

                kmers = find_kmers(tmpfile.name)
                master_headers = list(kmers.keys())
                for i in range(len(master_headers) -1, -1, -1):  # reverse iteration
                    master_header = master_headers[i]
                    master = kmers[master_header]
                    if master:
                        for key in master_headers[:i]:
                            candidate = kmers[key]
                            if candidate:
                                if not master.isdisjoint(candidate):
                                    candidate.update(master)
                                    children = [master_header.strip(">")] + cluster_children[master_header.strip(">")]
                                    cluster_children[key.strip(">")].extend(children)
                                    cluster_children[master_header.strip(">")] = None
                                    kmers[master_header] = None
                                    break
                            

            for master, children in cluster_children.items():
                if children is not None:
                    clusters.append([master, *children])
            cluster_time = keeper.differential()
            printv(f"Found {seq_count} sequences over {len(clusters)} clusters. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
            printv(f"Aligning Clusters. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
            
            items = os.listdir(raw_files_tmp)
            items.sort(key = lambda x: int(x.split("_cluster")[-1]))

            for i, cluster in enumerate(clusters):
                printv(f"Aligning cluster {i}. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
                aligned_cluster = os.path.join(aligned_files_tmp, f"{args.gene}_cluster{len(aligned_ingredients)}")
                
                cluster_seqs =  [(header, sequence) for header, sequence in data if header in cluster]

                if len(cluster_seqs) == 1:
                    to_merge.extend(cluster_seqs)
                    
                else:
                    aligned_ingredients.append(aligned_cluster)
                    raw_cluster = os.path.join(raw_files_tmp, f"{args.gene}_cluster{len(aligned_ingredients)}")
                    writeFasta(raw_cluster, cluster_seqs)

                    command = args.string.format(
                            in_file=raw_cluster, out_file=aligned_cluster
                        )
                    
                    os.system(command)
                    if debug:
                        printv(command, args.verbose, 3)
                        writeFasta(os.path.join(this_intermediates, f"{args.gene}_cluster{len(aligned_ingredients)}_aligned"), parseFasta(aligned_cluster))
        align_time = keeper.differential() - cluster_time
        if to_merge:
            merged_singleton_final = os.path.join(aligned_files_tmp, f"{args.gene}_cluster{len(aligned_ingredients)}")
            writeFasta(merged_singleton_final, to_merge)
        if aligned_ingredients:
            if not to_merge:
                out_file = args.result_file
            else:
                out_file = os.path.join(parent_tmpdir, f"{args.gene}_mafft_merged")
            aligned_ingredients.sort(key = lambda x: int(x.split("_cluster")[-1]))
            printv(f"Merging Alignments. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
            with NamedTemporaryFile(
                    dir = parent_tmpdir, mode="w+"
                ) as tmp, NamedTemporaryFile(
                    dir = parent_tmpdir, mode="w+"
                ) as tmp_special, NamedTemporaryFile(
                    dir = parent_tmpdir, mode="w+"
                ) as tmp_merged:
                sequences = []
                for header, sequence in parseFasta(aln_file):
                    if header in targets:
                        sequences.append((targets[header], sequence))

                empty_columns = None
                for header, sequence in sequences:
                    if empty_columns is None:
                        empty_columns = [True] * len(sequence)

                    for col, let in enumerate(sequence):
                        if let != "-":
                            empty_columns[col] = False
                
                to_write=  []
                for header, sequence in sequences:
                    to_write.append((header, "".join([let for col, let in enumerate(sequence) if not empty_columns[col]])))

                writeFasta(tmp.name, to_write)
                if debug:
                    writeFasta(os.path.join(this_intermediates, "references.fa"), to_write)
                tmp.flush()

                sequence_done = set()

                lines = []
                total = 0 
                aligned_to_write = []
                for item in aligned_ingredients:
                    file = os.path.basename(item)
                    lines.append(file)
                    to_write = []
                    for header, sequence in parseFasta(item):
                        if sequence not in sequence_done:
                            to_write.append((header, sequence))
                            sequence_done.add(sequence)
                            
                    seq_count = len(to_write)
                    aligned_to_write.extend(to_write)
                    line = " " + " ".join(map(str, range(1+total, seq_count+1+total))) + f" # {file}"
                    total += seq_count
                    lines.append(line)
                
                tmp_special.write("\n".join(lines))
                tmp_special.flush()

                writeFasta(tmp_merged.name, aligned_to_write)
                tmp_merged.flush()

                os.system(f"mafft --quiet --merge {tmp_special.name} {tmp_merged.name} > {out_file}")
                if args.debug:
                    print(f"mafft --quiet --merge {tmp_special.name} {tmp_merged.name} > {out_file}")

        if to_merge:
            if aligned_ingredients:
                os.system(f"mafft --anysymbol --quiet --jtt 1 --addfragments {merged_singleton_final} --thread 1 {out_file} > {args.result_file}")
                if debug:
                    printv(f"mafft --anysymbol --quiet --jtt 1 --addfragments {merged_singleton_final} --thread 1 {out_file} > {args.result_file}", args.verbose, 3)
                #Deinterleave
                try:
                    to_write = [(header, sequence) for header, sequence in parseFasta(args.result_file)]
                except:
                    print("SINGLETON MERGE FAILED:",args.gene,merged_singleton_final, out_file, args.result_file)
                writeFasta(args.result_file, to_write)
                if debug:
                    writeFasta(os.path.join(this_intermediates, os.path.basename(args.result_file)), to_write)
            else:
                with NamedTemporaryFile(dir = parent_tmpdir, mode="w+") as tmp:
                    sequences = []
                    for header, sequence in parseFasta(aln_file):
                        if header in targets:
                            sequences.append((targets[header], sequence))

                    empty_columns = None
                    for header, sequence in sequences:
                        if empty_columns is None:
                            empty_columns = [True] * len(sequence)

                        for col, let in enumerate(sequence):
                            if let != "-":
                                empty_columns[col] = False
                    
                    to_write=  []
                    for header, sequence in sequences:
                        to_write.append((header, "".join([let for col, let in enumerate(sequence) if not empty_columns[col]])))

                    writeFasta(tmp.name, to_write)
                    tmp.flush()
                
                    out_file = args.result_file
                    os.system(f"mafft --anysymbol --quiet --jtt 1 --addfragments {merged_singleton_final} --thread 1 {tmp.name} > {out_file}")
                    if debug:
                        printv(f"mafft --anysymbol --quiet --jtt 1 --addfragments {merged_singleton_final} --thread 1 {tmp.name} > {out_file}", args.verbose, 3)
                    empty_columns = None
                    for header, sequence in parseFasta(out_file):
                        if not empty_columns:
                            empty_columns = [True] * len(sequence)

                        for col, let in enumerate(sequence):
                            if let != "-":
                                empty_columns[col] = False
                    
                    to_write=  []
                    for header, sequence in parseFasta(out_file):
                        to_write.append((header, "".join([let for col, let in enumerate(sequence) if not empty_columns[col]])))

                    writeFasta(out_file, to_write)

                    if debug:
                        writeFasta(os.path.join(this_intermediates, os.path.basename(out_file)), to_write)
        merge_time = keeper.differential() - cluster_time - align_time
        printv(f"Done. {args.gene} took {keeper.differential():.2f}s", args.verbose, 2) # Debug

    return args.gene, cluster_time, align_time, merge_time, keeper.differential()

def do_folder(folder, args):
    printv(f"Processing: {os.path.basename(folder)}", args.verbose)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    mafft_path = os.path.join(folder, MAFFT_FOLDER)
    aa_path = os.path.join(folder, AA_FOLDER)
    if not os.path.exists(aa_path):
        printv(f"ERROR: Can't find aa ({aa_path}) folder. Abort", args.verbose, 0)
        return False
    rmtree(mafft_path, ignore_errors=True)
    os.mkdir(mafft_path)

    genes = [
        gene
        for gene in os.listdir(aa_path)
        if gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]
    orthoset_path = os.path.join(args.orthoset_input, args.orthoset)
    aln_path = os.path.join(orthoset_path, ALN_FOLDER)
    if not os.path.exists(orthoset_path):
        printv("ERROR: Orthoset path not found.", args.verbose, 0)
        return False
    if not os.path.exists(aln_path):
        printv("ERROR: Aln folder not found.", args.verbose, 0)
        return False

    command = f"clustalo -i {{in_file}} -o {{out_file}} --threads=1 --iter=2 --full --full-iter --force"

    intermediates = "intermediates"
    if not os.path.exists(intermediates):
        os.mkdir(intermediates)

    if args.processes > 1:
        arguments = []
        func = arguments.append
        lock = Lock()
    else:
        func = run_command
        lock = None

    times = []
    for file in genes:
        gene = file.split(".")[0]
        gene_file = os.path.join(aa_path, file)
        result_file = os.path.join(mafft_path, file.rstrip(".gz"))
        if func == run_command:
            times.append(
                run_command(
                        CmdArgs(
                        command, gene_file, result_file, gene, lock, args.verbose, args.compress, aln_path, args.debug
                    )
                )
            )
        else:
            func(
                (CmdArgs(
                    command, gene_file, result_file, gene, lock, args.verbose, args.compress, aln_path, args.debug
                ),)
            )

    if args.processes > 1:
        with ThreadPool(args.processes) as pool:
            times = pool.starmap(run_command, arguments, chunksize=1)

    if args.debug:
        with open("mafft_times.csv", "w") as fp:
            fp.write("\n".join(["Gene,Clustering,Aligning,Merging,Total"]+[",".join(map(str, i)) for i in times]))

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
