from __future__ import annotations

import json
import math
import os
from collections import Counter
from dataclasses import dataclass
from itertools import count
from multiprocessing.pool import Pool
from shutil import rmtree
import sys
from tempfile import TemporaryDirectory, NamedTemporaryFile
from time import time

import wrap_rocks
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord

from .utils import printv, gettempdir
from .timekeeper import TimeKeeper, KeeperMode


class Result:
    __slots__ = (
        "header",
        "target", 
        "evalue", 
        "score", 
        "blast_start", 
        "blast_end",
        "reftaxon",
        "hmmsearch_id",
        "gene"
    )

    def __init__(self, gene, hmmsearch_id, header, reftaxon, target, evalue, bit_score, q_start, q_end):
        self.gene = gene
        self.target = target
        self.header = header
        self.hmmsearch_id = hmmsearch_id
        self.reftaxon = reftaxon
        self.evalue = float(evalue)
        self.score = float(bit_score)
        self.blast_start = int(q_start)
        self.blast_end = int(q_end)

    def to_json(self):
        return {
                    "hmmId": self.hmmsearch_id,
                    "gene": self.gene,
                    "target": self.target,
                    "refTaxon": self.reftaxon,
                }

def run_process(args, input_path) -> None:
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    taxa = os.path.basename(input_path)
    printv(f'Processing: {taxa}', args.verbose, 0)
    printv("Grabbing reference data from Orthoset DB.", args.verbose)
    # make dirs
    blast_path = os.path.join(input_path, "blast")
    if args.overwrite:
        if os.path.exists(blast_path):
            rmtree(blast_path)
    os.makedirs(blast_path, exist_ok=True)

    num_threads = args.processes

    db_path = os.path.join(orthosets_dir, orthoset, "rocksdb")
    #Spud check
    if not os.path.exists(db_path):
        if input(f"Could not find orthoset DB at {db_path}. Would you like to generate it? Y/N: ").lower() == "y":
            print("Attempting to generate DB")
            os.system(f"python3 -m phymmr -p {num_threads} Makeref {orthoset}.sqlite -s {orthoset}")
        else:
            print("Aborting")
            sys.exit(1)
    orthoset_db = wrap_rocks.RocksDB(db_path)

    target_to_taxon = json.loads(orthoset_db.get("getall:targetreference"))

    del orthoset_db

    diamond_db_path = os.path.join(orthosets_dir, orthoset, "diamond", orthoset+'.dmnd')

    printv(f"Done! Took {time_keeper.lap():.2f}s. Writing reference data to DB", args.verbose)

    db_path = os.path.join(input_path, "rocksdb", "hits")
    db = wrap_rocks.RocksDB(db_path)

    printv(f"Done! Took {time_keeper.lap():.2f}s. Grabbing HMM data from DB", args.verbose)

    global_hmm_object_raw = db.get("hmmbatch:all")
    global_hmm_batches = global_hmm_object_raw.split(",")

    with TemporaryDirectory(dir=gettempdir()) as tmpdir, NamedTemporaryFile(dir=tmpdir, mode="w+") as tmpfile: #seq_hits
        out_lines = []
        for batch_i in global_hmm_batches:
            key = f"diamondbatch:{batch_i}"
            seq_json = db.get(key)
            seq_hits = json.loads(seq_json)
            for hit in seq_hits:
                out_lines.append(hit)
        tmpfile.writelines(out_lines)
        tmpfile.flush()

        os.system(f"diamond blastp --db {diamond_db_path} --query {tmpfile.name} --outfmt 6 qseqid sseqid evalue bitscore qstart qend --very-sensitive --masking 0 --threads {num_threads} --out {os.path.join(blast_path, 'blast.tsv')}")

    to_write = []
    with open(os.path.join(blast_path, 'blast.tsv'), "r") as fp:
        for line in fp:
            header, reference, evalue, score, start, end = line.strip().split("\t")
            
            gene, reftaxon = target_to_taxon[reference]

            header,hmmsearch_id = header.split("_hmmid")

            this_result = Result(gene, hmmsearch_id, header, reftaxon, reference, evalue, score, start, end)

            to_write.append(this_result.to_json())

    total = count()
    counter = count()
    this_batch = []
    recipe = []
    batch_i = 1
    for result in to_write:
        this_batch.append(result)
        if next(counter) == args.max_blast_batch_size:
            recipe.append(batch_i)
            db.put(f'blastbatch:{batch_i}', json.dumps(this_batch))
            batch_i += 1
            counter = count()
            this_batch = []
        next(total)
    
    if this_batch:
        recipe.append(batch_i)
        db.put(f'blastbatch:{batch_i}', json.dumps(this_batch))
        batch_i += 1

    db.put('blastbatch:all', ','.join(map(str, recipe)))

    printv(f"Writing {next(total)} results over {batch_i} batches took {time_keeper.lap():.2f}s", args.verbose)
    printv(f"Done! Took {time_keeper.differential():.2f}s overall", args.verbose)

def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for input_path in args.INPUT:
        run_process(args, input_path)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr BlastPal"
    )
