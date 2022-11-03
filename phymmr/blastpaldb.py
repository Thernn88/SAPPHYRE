import argparse
import json
import math
import os
import sqlite3
import wrap_rocks
from dataclasses import dataclass
from multiprocessing.pool import Pool
from shutil import rmtree
from time import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter

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
        "ref_sequence"
    )

    def __init__(self, query_id, subject_id, evalue, bit_score, q_start, q_end):
        self.query_id = str(query_id)
        self.target = int(subject_id)
        self.evalue = float(evalue)
        self.log_evalue = math.log(self.evalue) if self.evalue != 0 else -999.0
        self.score = float(bit_score)
        self.blast_start = int(q_start)
        self.blast_end = int(q_end)
        self.reftaxon = None
        self.ref_sequence = None

    def to_json(self):
        return {
                    "target": self.target,
                    "reftaxon": self.reftaxon,
                    "ref_sequence": self.ref_sequence
                }

@dataclass
class GeneConfig:  # FIXME: I am not certain about types.
    gene: str
    tmp_path: str
    gene_sequences: list
    ref_names: dict
    blast_path: str
    blast_db_path: str
    blast_minimum_score: float
    blast_minimum_evalue: float

    def __post_init__(self):
        self.blast_file_path = os.path.join(self.blast_path, f"{self.gene}.blast")
        self.blast_file_done = f"{self.blast_file_path}.done"
        self.fa_file = os.path.join(self.tmp_path, f"{self.gene}.fa")


def do(
    gene_conf: GeneConfig,
    verbose: bool,
    prog="blastp",
):
    if (
        not os.path.exists(gene_conf.blast_file_done)
        or os.path.getsize(gene_conf.blast_file_done) == 0
    ):
        printv("Blasted: {}".format(gene_conf.gene), verbose)
        with open(gene_conf.fa_file, "w") as fp:
            fw = FastaWriter(fp)
            fw.write_file(gene_conf.gene_sequences)

        cmd = (
            "{prog} -outfmt '7 qseqid sseqid evalue bitscore qstart qend' "
            "-evalue '{evalue_threshold}' -threshold '{score_threshold}' "
            "-num_threads '{num_threads}' -db '{db}' -query '{queryfile}' "
            "-out '{outfile}'".format(
                prog=prog,
                evalue_threshold=gene_conf.blast_minimum_evalue,
                score_threshold=gene_conf.blast_minimum_score,
                num_threads=2,
                db=gene_conf.blast_db_path,
                queryfile=gene_conf.fa_file,
                outfile=gene_conf.blast_file_path,
            )
        )
        os.system(cmd)
        os.remove(gene_conf.fa_file)

        os.rename(gene_conf.blast_file_path, gene_conf.blast_file_done)

    gene_out = {}
    this_return = []

    with open(gene_conf.blast_file_done, "r") as fp:
        for line in fp:
            if line == "\n" or line[0] == "#":
                continue
            fields = line.split("\t")

            this_result = Result(*fields)

            if this_result.target in gene_conf.ref_names: # Hit target not valid
                this_result.reftaxon, this_result.ref_sequence = gene_conf.ref_names[this_result.target]

                # Although we have a threshold in the Blast call. Some still get through.
                if (
                    this_result.score >= gene_conf.blast_minimum_score
                    and this_result.evalue <= gene_conf.blast_minimum_evalue
                ):
                    _, hmmsearch_id = this_result.query_id.split("_hmmid")

                    gene_out.setdefault(hmmsearch_id, [])
                    gene_out[hmmsearch_id].append(this_result.to_json())

    for hmmsearch_id in gene_out:
        this_out_results = gene_out[hmmsearch_id]
        key = "blastfor:{}".format(hmmsearch_id)

        data = json.dumps(this_out_results)

        this_return.append((key, data, len(this_out_results)))
    
    return this_return

def get_set_id(orthoset_db_con, orthoset):
    """
    Retrieves orthoset id from orthoset db
    """

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute(f'SELECT id FROM orthograph_set_details WHERE name = "{orthoset}";')

    for row in rows:
        return row[0]

    raise Exception("Orthoset {} id cant be retrieved".format(orthoset))

def get_ref_taxon_for_genes(set_id, orthoset_db_con):
    data = {}
    orthoset_db_cur = orthoset_db_con.cursor()
    query = 'SELECT id, sequence FROM orthograph_aaseqs'
    rows = orthoset_db_cur.execute(query)

    for row in rows:
        id, sequence = row
        data[id] = sequence

    result = {}
    query = f'''SELECT DISTINCT
        orthograph_taxa.name,
        orthograph_orthologs.ortholog_gene_id,
        orthograph_aaseqs.id
        FROM orthograph_orthologs
        INNER JOIN orthograph_sequence_pairs
            ON orthograph_orthologs.sequence_pair = orthograph_sequence_pairs.id
        INNER JOIN orthograph_aaseqs
            ON orthograph_sequence_pairs.aa_seq = orthograph_aaseqs.id
        INNER JOIN orthograph_set_details
            ON orthograph_orthologs.setid = orthograph_set_details.id
        INNER JOIN orthograph_taxa
            ON orthograph_aaseqs.taxid = orthograph_taxa.id
        WHERE orthograph_set_details.id = "{set_id}"'''

    rows = orthoset_db_cur.execute(query)

    for row in rows:
        name, gene, id = row
        
        if gene not in result:
            result[gene] = {id:(name, data[id])}
        else:
            result[gene][id] = (name, data[id])

    return result


def printv(msg, verbosity) -> None:
    if verbosity:
        print(msg)


def run_process(args, input_path) -> None:
    start = time()
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    taxa = os.path.basename(input_path)
    print(f'Begin BlastPal for {taxa}')
    printv("Grabbing Reference data from SQL.", args.verbose)
    # make dirs
    blast_path = os.path.join(input_path, "blast")
    if args.overwrite:
        if os.path.exists(blast_path):
            rmtree(blast_path)
    os.makedirs(blast_path, exist_ok=True)

    if os.path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif os.path.exists("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = os.path.join(input_path, "tmp")
        os.makedirs(tmp_path, exist_ok=True)

    # grab gene reftaxon
    orthoset_db_path = os.path.join(orthosets_dir, orthoset + ".sqlite")
    orthoset_db_con = sqlite3.connect(orthoset_db_path)

    orthoset_id = get_set_id(orthoset_db_con, orthoset)

    sql_start = time()
    ref_taxon = get_ref_taxon_for_genes(orthoset_id, orthoset_db_con)

    blast_db_path = os.path.join(orthosets_dir, orthoset, "blast", orthoset)

    printv("Done! Took {:.2f}s. Grabbing HMM data from DB".format(time() - sql_start), args.verbose)

    db_path = os.path.join(input_path, "rocksdb", "hits")
    db = wrap_rocks.RocksDB(db_path)

    gene_to_hits = {}

    grab_hmm_start = time()

    global_hmm_object_raw = db.get("hmmbatch:all")
    global_hmm_batches = global_hmm_object_raw.split(",")
    hit_count = 0

    for batch_i in global_hmm_batches:
        key = f"hmmbatch:{batch_i}"
        hmm_json = db.get(key)
        hmm_hits = json.loads(hmm_json)

        for hmm_object in hmm_hits:
            hit_count += 1
            hmm_id = hmm_object["hmm_id"]
            header = hmm_object["header"].strip()
            gene = hmm_object["gene"]

            gene_to_hits.setdefault(gene, [])
            gene_to_hits[gene].append(
                SeqRecord(Seq(hmm_object["hmm_sequence"]), id=header + f"_hmmid{hmm_id}", description=""))

    genes = list(gene_to_hits.keys())

    printv(
        "Grabbed HMM Data. Took: {:.2f}s. Found {} hits".format(
            time() - grab_hmm_start, hit_count
        ), args.verbose
    )

    del global_hmm_object_raw

    blast_start = time()
    num_threads = args.processes
    # Run
    if num_threads <= 1:
        to_write = [
            do(
                GeneConfig(
                    gene=gene,
                    tmp_path=tmp_path,
                    gene_sequences=gene_to_hits[gene],
                    ref_names=ref_taxon[gene],
                    blast_path=blast_path,
                    blast_db_path=blast_db_path,
                    blast_minimum_score=args.blast_minimum_score,
                    blast_minimum_evalue=args.blast_minimum_evalue,
                ),
                args.verbose
            )
            for gene in genes
        ]
    else:
        arguments = [
            (GeneConfig(
                gene=gene,
                tmp_path=tmp_path,
                gene_sequences=gene_to_hits[gene],
                ref_names=ref_taxon[gene],
                blast_path=blast_path,
                blast_db_path=blast_db_path,
                blast_minimum_score=args.blast_minimum_score,
                blast_minimum_evalue=args.blast_minimum_evalue,
            ),
             args.verbose,)
            for gene in genes
        ]

        with Pool(num_threads) as pool:
            to_write = pool.starmap(do, arguments, chunksize=1)

    printv(
        "Got Blast Results. Took {:.2f}s. Writing to DB".format(
            time() - blast_start
        ), args.verbose
    )

    write_start = time()

    i = 0
    for batch in to_write:
        for key, data, count in batch:
            db.put(key, data)
            i += count

    printv("Writing {} results took {:.2f}s".format(i, time() - write_start), args.verbose)

    print("Done. Took {:.2f}s overall".format(time() - start))


def main(args):
    if not all(os.path.exists(i) for i in args.INPUT):
        print("ERROR: All folders passed as argument must exists.")
        return False
    for input_path in args.INPUT:
        run_process(args, input_path)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr BlastPal"
    )
