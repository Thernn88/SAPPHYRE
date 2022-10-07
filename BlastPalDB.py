import argparse
import json
import math
import os
from dataclasses import dataclass
from multiprocessing.pool import Pool
from shutil import rmtree
from time import time

import wrap_rocks


@dataclass
class GeneConfig:  # FIXME: I am not certain about types.
    gene: str
    tmp_path: str
    gene_sequence: str
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
    prog="blastp",
    evalue_threshold=0.00001,
    score_threshold=40,
    num_threads=1,
):
    if (
        not os.path.exists(gene_conf.blast_file_done)
        or os.path.getsize(gene_conf.blast_file_done) == 0
    ):
        print("Blasted:", gene_conf.gene)
        with open(gene_conf.fa_file, "w") as target_handle:
            for target, sequence in gene_conf.gene_sequences:
                target_handle.write(">" + target + "\n" + sequence + "\n")

        cmd = (
            "{prog} -outfmt '7 qseqid sseqid evalue bitscore qstart qend' "
            "-evalue '{evalue_threshold}' -threshold '{score_threshold}' "
            "-num_threads '{num_threads}' -db '{db}' -query '{queryfile}' "
            "-out '{outfile}'".format(
                prog=prog,
                evalue_threshold=evalue_threshold,
                score_threshold=score_threshold,
                num_threads=2,
                db=gene_conf.blast_db_path,
                queryfile=gene_conf.fa_file,
                outfile=gene_conf.blast_path,
            )
        )
        os.system(cmd)
        os.remove(gene_conf.fa_file)

        os.rename(gene_conf.blast_path, gene_conf.blast_file_done)

    gene_out = {}
    this_return = []

    with open(gene_conf.blast_file_done, "r") as fp:
        for line in fp:
            if line == "\n" or line[0] == "#":
                continue
            fields = line.split("\t")

            query_id, subject_id, evalue, bit_score, q_start, q_end = fields
            if (
                float(bit_score) >= gene_conf.blast_minimum_score
                and float(evalue) <= gene_conf.blast_minimum_evalue
            ):
                try:
                    log_evalue = str(math.log(float(evalue)))
                except ValueError:  # evaluate is 0, math domain error
                    log_evalue = "-999"

                query_id, hmmsearch_id = query_id.split("_hmmid")

                this_out = {
                    "target": int(subject_id),
                    "score": float(bit_score),
                    "evalue": float(evalue),
                    "log_evalue": float(log_evalue),
                    "blast_start": int(q_start),
                    "blast_end": int(q_end),
                }

                gene_out.setdefault(hmmsearch_id, [])
                gene_out[hmmsearch_id].append(this_out)

    for hmmsearch_id in gene_out:
        this_out_results = gene_out[hmmsearch_id]
        key = "blastfor:{}".format(hmmsearch_id)

        data = json.dumps(this_out_results)

        this_return.append((key, data, len(this_out_results)))

    

    return this_return


def main():
    start = time()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="PhyMMR/Acroceridae/SRR6453524.fa",
        help="Path to directory of Input folder",
    )
    parser.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="PhyMMR/orthosets",
        help="Path to directory of Orthosets folder",
    )
    parser.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    parser.add_argument(
        "-bs",
        "--blast_minimum_score",
        type=float,
        default=40.0,
        help="Minimum score filter in blast.",
    )
    parser.add_argument(
        "-be",
        "--blast_minimum_evalue",
        type=float,
        default=0.00001,
        help="Minimum evalue filter in blast.",
    )
    parser.add_argument(
        "-ovw",
        "--overwrite",
        action="store_true",
        help="Overwrite existing blast results.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Number of threads used to call processes.",
    )
    args = parser.parse_args()

    print("Grabbing HMM data from db.")

    input_path = args.input
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    # make dirs
    blast_path = os.path.join(input_path, "blast")
    if args.overwrite:
        rmtree(blast_path)
    os.makedirs(blast_path, exist_ok=True)

    if os.path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif os.path.exists("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = os.path.join(input_path, "tmp")
        os.makedirs(tmp_path, exist_ok=True)

    blast_db_path = os.path.join(orthosets_dir, orthoset, "blast", orthoset)

    db_path = os.path.join(input_path, "rocksdb")
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
                (header + f"_hmmid{hmm_id}", hmm_object["hmm_sequence"])
            )

    genes = list(gene_to_hits.keys())

    print(
        "Grabbed HMM Data. Took: {:.2f}s. Grabbed {} hits.".format(
            time() - grab_hmm_start, hit_count
        )
    )

    del global_hmm_object_raw

    num_threads = args.processes
    # Run
    if num_threads == 1:
        to_write = [
            do(
                GeneConfig(
                    gene=gene,
                    tmp_path=tmp_path,
                    gene_sequence=gene_to_hits[gene],
                    blast_path=blast_path,
                    blast_db_path=blast_db_path,
                    blast_minimum_score=args.blast_minimum_score,
                    blast_minimum_evalue=args.blast_minimum_evalue,
                )
            )
            for gene in genes
        ]
    else:
        arguments = [
            (GeneConfig(
                gene=gene,
                tmp_path=tmp_path,
                gene_sequence=gene_to_hits[gene],
                blast_path=blast_path,
                blast_db_path=blast_db_path,
                blast_minimum_score=args.blast_minimum_score,
                blast_minimum_evalue=args.blast_minimum_evalue,
            ),)
            for gene in genes
        ]

        with Pool(num_threads) as pool:
            to_write = pool.starmap(do, arguments, chunksize=1)

    print("Writing to DB.")

    write_start = time()

    i = 0
    for batch in to_write:
        for key, data, count in batch:
            db.put(key, data)
            i += count  # REVIEW: why do we count?

    print(
        "Done. Took {:.2f}s overall. Writing {} results took {:.2f}s.".format(
            time() - start, i, time() - write_start
        )
    )


if __name__ == "__main__":
    main()
