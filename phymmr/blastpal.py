from __future__ import annotations

import json
import math
import os
from collections import Counter
from dataclasses import dataclass
from itertools import count
from multiprocessing.pool import Pool
from shutil import rmtree
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
        "query_id",
        "target", 
        "evalue", 
        "log_evalue", 
        "score", 
        "blast_start", 
        "blast_end",
        "reftaxon",
        "hmmsearch_id",
        "gene"
    )

    def __init__(self, gene, query_id, subject_id, evalue, bit_score, q_start, q_end):
        self.gene = gene
        self.query_id, self.hmmsearch_id = str(query_id).split('_hmmid')
        self.target = str(subject_id)
        self.evalue = float(evalue)
        self.log_evalue = math.log(self.evalue) if self.evalue != 0 else -999.0
        self.score = float(bit_score)
        self.blast_start = int(q_start)
        self.blast_end = int(q_end)
        self.reftaxon = None

    def to_json(self):
        return {
                    "hmmId": self.hmmsearch_id,
                    "gene": self.gene,
                    "target": self.target,
                    "refTaxon": self.reftaxon,
                }


@dataclass
class GeneConfig:  # FIXME: I am not certain about types.
    gene: str
    gene_sequences: list
    ref_names: dict
    blast_path: str
    blast_db_path: str
    blast_minimum_score: float
    blast_minimum_evalue: float

    def __post_init__(self):
        self.blast_file_path = os.path.join(self.blast_path, f"{self.gene}.blast")
        self.blast_file_done = f"{self.blast_file_path}.done"

def blast(
    gene_conf: GeneConfig,
    verbose: bool,
) -> None:
    with TemporaryDirectory(dir=gettempdir()) as tmpdir, NamedTemporaryFile(dir=tmpdir, mode="w+") as tmpfile:
            tmpfile.write("".join([f">{header}\n{sequence}\n" for header, sequence in gene_conf.gene_sequences])) # Write sequences
            tmpfile.flush() # Flush the internal buffer so it can be read by blastp

            cmd = (
                "blastp -outfmt '7 qseqid sseqid evalue bitscore qstart qend' "
                "-evalue '{evalue_threshold}' -threshold '{score_threshold}' "
                "-num_threads '{num_threads}' -db '{db}' -query '{queryfile}' "
                "-out '{outfile}'".format(
                    evalue_threshold=gene_conf.blast_minimum_evalue,
                    score_threshold=gene_conf.blast_minimum_score,
                    num_threads=2,
                    db=gene_conf.blast_db_path,
                    queryfile=tmpfile.name,
                    outfile=gene_conf.blast_file_path,
                )
            )
            os.system(cmd)
            os.rename(gene_conf.blast_file_path, gene_conf.blast_file_done)
            printv(f"Blasted: {gene_conf.gene}", verbose, 2)

def do(
    gene_conf: GeneConfig,
    verbose: bool,
):
    if (
        not os.path.exists(gene_conf.blast_file_done)
        or os.path.getsize(gene_conf.blast_file_done) == 0
    ):
        blast(gene_conf, verbose)

    gene_out = []

    with open(gene_conf.blast_file_done, "r") as fp:
        for line in fp:
            if line == "\n" or line[0] == "#":
                continue
            fields = line.split("\t")

            this_result = Result(gene_conf.gene, *fields)

            if this_result.target in gene_conf.ref_names: # Hit target not valid
                this_result.reftaxon = gene_conf.ref_names[this_result.target]

                # Although we have a threshold in the Blast call. Some still get through.
                if (
                    this_result.score >= gene_conf.blast_minimum_score
                ):
                    gene_out.append(this_result.to_json())

    return gene_out


def get_set_id(orthoset_db_con, orthoset):
    """
    Retrieves orthoset id from orthoset db
    """
    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute(f'SELECT id FROM orthograph_set_details WHERE name = "{orthoset}";')
    for row in rows:
        return row[0]

    raise Exception("Orthoset {} id cant be retrieved".format(orthoset))


def run_process(args, input_path) -> None:
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    taxa = os.path.basename(input_path)
    printv(f'Processing: {taxa}', args.verbose, 0)
    printv("Grabbing Reference data from SQL.", args.verbose)
    # make dirs
    blast_path = os.path.join(input_path, "blast")
    if args.overwrite:
        if os.path.exists(blast_path):
            rmtree(blast_path)
    os.makedirs(blast_path, exist_ok=True)

    db_path = os.path.join(orthosets_dir, orthoset, "rocksdb")
    orthoset_db = wrap_rocks.RocksDB(db_path)

    target_to_taxon = json.loads(orthoset_db.get("getall:targetreference"))

    del orthoset_db

    blast_db_path = os.path.join(orthosets_dir, orthoset, "blast", orthoset)

    printv(f"Done! Took {time_keeper.lap():.2f}s. Writing reference data to DB", args.verbose)

    db_path = os.path.join(input_path, "rocksdb", "hits")
    db = wrap_rocks.RocksDB(db_path)

    printv(f"Done! Took {time_keeper.lap():.2f}s. Grabbing HMM data from DB", args.verbose)

    gene_to_hits = {}

    global_hmm_object_raw = db.get("hmmbatch:all")
    global_hmm_batches = global_hmm_object_raw.split(",")
    genes = Counter()

    for batch_i in global_hmm_batches:
        key = f"hmmbatch:{batch_i}"
        hmm_json = db.get(key)
        hmm_hits = json.loads(hmm_json)

        for hmm_object in hmm_hits:
            hmm_id = hmm_object["hmm_id"]
            header = hmm_object["header"].strip()
            gene = hmm_object["gene"]
            genes[gene] = genes.get(gene, 0) + 1

            gene_to_hits.setdefault(gene, [])
            gene_to_hits[gene].append((header + f"_hmmid{hmm_id}", hmm_object["hmm_sequence"]))

    printv(f"Grabbed HMM Data. Took: {time_keeper.lap():.2f}s. Found {sum(genes.values())} hits", args.verbose)

    del global_hmm_object_raw

    num_threads = args.processes
    # Run
    if num_threads <= 1:
        to_write = [
            do(
                GeneConfig(
                    gene=gene,
                    gene_sequences=gene_to_hits[gene],
                    ref_names=target_to_taxon[gene],
                    blast_path=blast_path,
                    blast_db_path=blast_db_path,
                    blast_minimum_score=args.blast_minimum_score,
                    blast_minimum_evalue=args.blast_minimum_evalue,
                ),
                args.verbose
            )
            for gene, _ in genes.most_common()
        ]
    else:
        arguments = [
            (
                GeneConfig(
                    gene=gene,
                    gene_sequences=gene_to_hits[gene],
                    ref_names=target_to_taxon[gene],
                    blast_path=blast_path,
                    blast_db_path=blast_db_path,
                    blast_minimum_score=args.blast_minimum_score,
                    blast_minimum_evalue=args.blast_minimum_evalue,
                ),
                args.verbose,
            )
            for gene, _ in genes.most_common()
        ]

        with Pool(num_threads) as pool:
            to_write = pool.starmap(do, arguments, chunksize=1)

    printv(f"Got Blast Results. Took {time_keeper.lap():.2f}s. Writing to DB", args.verbose)

    total = count()
    counter = count()
    this_batch = []
    recipe = []
    batch_i = 1
    for batch in to_write:
        for result in batch:
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
