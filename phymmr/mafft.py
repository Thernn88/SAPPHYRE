from __future__ import annotations
import argparse
import gzip

import os
from collections import namedtuple
from multiprocessing.pool import ThreadPool
from shutil import rmtree
from tempfile import TemporaryDirectory, NamedTemporaryFile
from threading import Lock
from .utils import printv, gettempdir, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode
from Bio.Seq import Seq
KMER_LEN = 45   

### Retrofitted from suchapalaver/fastas2kmers

def indexfasta(filename):
    infile = open(filename, 'rb') # opens and reads the fasta file in binary
    chunksize = 1024*1024 # it reads in these chunks
    filepos = 0
    headstart = list()
    headend = list()
    while True: # Exit loop when, chunk by chunk, we've gone through the whole file
        content = infile.read(chunksize)
        if len(content) == 0:
            break
        chunkpos = 0 # chunks away!
        while chunkpos != -1: # exit this loop when we're at the file's end
            chunkpos = content.find(b'>', chunkpos) # chunk from 1st identifier after current chunkpos
            if chunkpos != -1: # i.e. when there are no more left, find will return -1
                headstart.append(chunkpos + filepos) # headstart from '>' in record
                chunkpos += 1 # chunking beyond previous '>' to get to next record
        for i in range(len(headend), len(headstart)): # how many records we're looking for
            chunkpos = max(0, headstart[i] - filepos)
            chunkpos = content.find(b'\n', chunkpos)
            if chunkpos != -1:
                headend.append(chunkpos + filepos)
        filepos += len(content)
    infile.close()
    # Eliminating wrong headers due to extra > in header line
    for i in range(len(headstart)-1, 0, -1):
        if headend[i] == headend[i-1]:
            del headstart[i]
            del headend[i]            
    headstart.append(filepos)
    fastaindex = list()
    with open(filename, 'rb') as fh:
        for i in range(len(headend)):
            seq_start = headend[i]+1
            seq_end = headstart[i+1] - 1
            fh.seek(headstart[i])
            identifier = str(fh.read((headend[i]) - headstart[i])).strip("b'\n ")
            fastaindex.append((identifier, (seq_start, seq_end, seq_end-seq_start)))
    return fastaindex

def indexsequence(seq: str) -> int:
    """
    Given a sequence, find the first and last characters, then
    return the values as a tuple in a list
    """
    pointer = 0
    seqindex = list()
    while len(seq) > pointer:
        potenstart = [seq.find(b'a', pointer), seq.find(b't', pointer), seq.find(b'c', pointer), seq.find(b'g', pointer)]
        realstart = min(potenstart)
        if realstart == -1:
            # happens rarely, so slow code is ok, apparently
            potenstart = [i for i in potenstart if i > -1]
            if len(potenstart) == 0:
                break
            realstart = min(potenstart)
        realend = seq.find(b'N', realstart)
        if realend == -1:
            realend = len(seq)
        seqindex.append((realstart, realend))
        pointer = realend
    return seqindex

def find_kmers(fasta, i):
    transtable = bytes.maketrans(b'ATCGMRYKVHDBWmrykvhdbxnsw', b'atcgNNNNNNNNNNNNNNNNNNNNN')
    kmer_len = KMER_LEN
    infile = open(fasta, 'rb')
    infile.seek(i[1][0])     
    identifier = i[0]
    seq = infile.read(i[1][1] - i[1][0]+1).translate(transtable, b'\r\n\t ')
    infile.close()
    # Index sequence
    seqindex = indexsequence(seq)
    subset = set()
    seqdict = dict()
    for start, stop in seqindex:
        for i in range(start, stop-kmer_len+1):
            kmer = str(seq[i:i+kmer_len]).strip("b'").upper()
            subset.add(kmer)
        seqdict[identifier] = subset
        yield seqdict

###

MAFFT_FOLDER = "mafft"
AA_FOLDER = "aa"
ALN_FOLDER = "aln"


def process_genefile(fileread):
    data = []
    for header, sequence in parseFasta(fileread):
        if not header.endswith("."):
            data.append((header, sequence))
            
    return len(data), data


CmdArgs = namedtuple(
    "CmdArgs",
    ["string", "gene_file", "result_file", "gene", "lock", "verbose", "compress", "aln_path", "debug"],
)


def get_targets(gene_file: str) -> list[str]:
    targets = {}
    with open(gene_file, "r") as fp:
        for line in fp:
            if line.startswith(">"):
                if line.strip().endswith("."):
                    targets[line.split("|")[2]] = line.strip()[1:]
    return targets

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
        if not os.path.exists(intermediates):
            os.mkdir(intermediates)
        this_intermediates = os.path.join(intermediates, args.gene)
        if not os.path.exists(this_intermediates):
            os.mkdir(this_intermediates)

    aln_file = os.path.join(args.aln_path, args.gene + ".aln.fa")
    to_merge = []

    nt_file = args.gene_file.replace("/aa/", "/nt/").replace(".aa.", ".nt.")

    with TemporaryDirectory(dir=temp_dir) as parent_tmpdir, TemporaryDirectory(dir=parent_tmpdir) as raw_files_tmp, TemporaryDirectory(dir=parent_tmpdir) as aligned_files_tmp:
        printv(f"Cleaning NT file. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
        seq_count, data = process_genefile(nt_file)
        targets = get_targets(args.gene_file)

        if seq_count == 1:
            printv(f"Outputting singleton alignment. Elapsed time: {keeper.differential():.2f}, args.verbose, 3") # Debug
            aligned_file = os.path.join(aligned_files_tmp,f"{args.gene}_cluster0")
            data = [(header, str(Seq(sequence).translate())) for header, sequence in data]
            writeFasta(aligned_file, data)
            aligned_ingredients.append(aligned_file)
        else:
            printv(f"Generating Cluster. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
            clusters = []
            cluster_children = {header: [] for header, _ in data}
            kmers = {}
            with NamedTemporaryFile(mode="w+", dir=parent_tmpdir) as tmpfile:
                writeFasta(tmpfile.name, data)
                tmpfile.flush()

                my_indexes = indexfasta(tmpfile.name)

                
                for i in my_indexes:
                    for seq_dict in find_kmers(tmpfile.name, i):
                        for key, this_kmers in seq_dict.items():  # for header in seq_dicts
                            kmers[key] = this_kmers
                master_headers = list(kmers.keys())
                for i in range(len(master_headers) -1, -1, -1):  # reverse iteration
                    master_header = master_headers[i]
                    master = kmers[master_header]
                    if master:
                        # for key in list(kmers.keys())[:i]:
                        for key in master_headers[:i]:
                            candidate = kmers[key]
                            if candidate:
                                # intersecting_kmers = len(master.intersection(candidate))

                                # if intersecting_kmers > 0:
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

            data = [(header, str(Seq(sequence).translate())) for header, sequence in data]
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

        if to_merge:
            merged_singleton_final = os.path.join(aligned_files_tmp, f"{args.gene}_cluster{len(aligned_ingredients)}")
            writeFasta(merged_singleton_final, to_merge)
        if aligned_ingredients:
            aligned_ingredients.sort(key = lambda x: int(x.split("_cluster")[-1]))
            printv(f"Merging Alignments. Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
            printv("Merging 1 of", len(aligned_ingredients),f"Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
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
                if debug:
                    writeFasta(os.path.join(this_intermediates, "references.fa"), to_write)
                tmp.flush()
                
                out_file = os.path.join(parent_tmpdir,"part_0.fa")
                if not to_merge and len(aligned_ingredients) == 1:
                    out_file = args.result_file
                
                os.system(f"clustalo --p1 {tmp.name} --p2 {aligned_ingredients[0]} -o {out_file} --threads=1 --is-profile --force")
                if debug:
                    printv(f"clustalo --p1 {tmp.name} --p2 {aligned_ingredients[0]} -o {out_file} --threads=1 --is-profile --force", args.verbose, 3)
                    writeFasta(os.path.join(this_intermediates, os.path.basename(out_file)), parseFasta(out_file))

            for i, file in enumerate(aligned_ingredients[1:]):
                printv("Merging", i+2, "of", len(aligned_ingredients),f"Elapsed time: {keeper.differential():.2f}", args.verbose, 3) # Debug
                out_file = os.path.join(parent_tmpdir,f"part_{i+1}.fa")
                in_file = os.path.join(parent_tmpdir,f"part_{i}.fa")
                if not to_merge and i == len(aligned_ingredients) - 2:
                    out_file = args.result_file

               
                os.system(f"clustalo --p1 {in_file} --p2 {file} -o {out_file} --threads=1 --is-profile --force")
                if debug:
                    printv(f"clustalo --p1 {in_file} --p2 {file} -o {out_file} --threads=1 --is-profile --force", args.verbose, 3)
                    writeFasta(os.path.join(this_intermediates, os.path.basename(out_file)), parseFasta(out_file))

        if to_merge:
            if aligned_ingredients:
                os.system(f"mafft --anysymbol --quiet --jtt 1 --addfragments {merged_singleton_final} --thread 1 {out_file} > {args.result_file}")
                if debug:
                    printv(f"mafft --anysymbol --quiet --jtt 1 --addfragments {merged_singleton_final} --thread 1 {out_file} > {args.result_file}", args.verbose, 3)

                #Deinterleave
                to_write = [(header, sequence) for header, sequence in parseFasta(args.result_file)]
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
                        print(f"mafft --anysymbol --quiet --jtt 1 --addfragments {merged_singleton_final} --thread 1 {tmp.name} > {out_file}")
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

        printv(f"Done. {args.gene} took {keeper.differential():.2f}", args.verbose, 2) # Debug
            

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
                command, gene_file, result_file, gene, lock, args.verbose, args.compress, aln_path, args.debug
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
