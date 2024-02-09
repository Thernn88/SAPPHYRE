from argparse import ArgumentParser
from collections import defaultdict
from itertools import count
from pathlib import Path
from tempfile import TemporaryDirectory
from .timekeeper import TimeKeeper, KeeperMode
from wrap_rocks import RocksDB
import os
from msgspec import json
from .diamond import ReporterHit as Hit
from .prepare import group_taxa_in_glob
from .utils import gettempdir, printv
from pyfastx import Fastq
from sapphyre_tools import bio_revcomp

def get_diamondhits(
    rocks_hits_db: RocksDB,
) -> dict[str, list[Hit]]:
    """Returns a dictionary of gene to corresponding hits.

    Args:
    ----
        rocks_hits_db (RocksDB): RocksDB instance
    Returns:
        dict[str, list[Hit]]: Dictionary of gene to corresponding hits
    """
    present_genes = rocks_hits_db.get("getall:presentgenes")

    genes_to_process = present_genes.split(",")

    gene_based_results = []
    for gene in genes_to_process:
        gene_result = rocks_hits_db.get_bytes(f"gethits:{gene}")
        this_hits = json.decode(gene_result, type=list[Hit])

        gene_based_results.append((gene, this_hits))

    return gene_based_results


def correct_folder(folder, args):
    if not args.fastq:
        msg = "FastQ file not provided, please provide one using -fq (path)."
        raise FileNotFoundError(msg)
    nt_db_path = os.path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(nt_db_path)

    if nt_db.get("get:isassembly") == "True" or nt_db.get("get:isgenome") == "True":
        printv("Skipping Error Connection", args.verbose, 0)
        return True


    original_posiitons = json.decode(nt_db.get("getall:original_positions"), type=dict[str, int])

    hits_db = RocksDB(os.path.join(folder, "rocksdb", "hits"))


    diamond_hits = get_diamondhits(hits_db)

    positions_to_keep = {}
    gene_to_hit = defaultdict(list)

    for gene, hits in diamond_hits:
        for i, hit in enumerate(hits):
            positions_to_keep[original_posiitons[hit.node]] = (gene, i)
            gene_to_hit[gene].append(hit)

    keep_set = set(positions_to_keep.keys())
    header_to_old_pos = {}
    printv(f"Found {len(positions_to_keep)} hits from diamond.", args.verbose, 1)

    with TemporaryDirectory(dir=gettempdir()) as tempdir:
        output = os.path.join(tempdir, "output.fastq")
        result = os.path.join(tempdir, "result")
        os.makedirs(result)
        with open(output, "w") as out:
            for i, (header, seq, qual) in enumerate(Fastq(args.fastq, build_index = False)):
                if i in keep_set:
                    header_to_old_pos[header] = i
                    out.write(f"@{header}\n{seq}\n+{header}\n{qual}\n")

        printv("Correcting reads using spades.", args.verbose, 0)
        os.system(f"spades --only-error-correction --disable-gzip-output -s '{output}' -o '{result}' > /dev/null")

        printv("Reading result.", args.verbose, 1)
        expected_outupt_location = os.path.join(result, "corrected")
        result_file = None
        for item in os.listdir(expected_outupt_location):
            if item.endswith(".fastq"):
                result_file = os.path.join(expected_outupt_location, item)
                break
        
        if not result_file:
            msg = "Could not find the result file"
            raise FileNotFoundError(msg)

        gene_output = defaultdict(list)
        for i, (header, seq, qual) in enumerate(Fastq(result_file, build_index = False)):   
            gene, hit_index = positions_to_keep[header_to_old_pos[header]]
            hit = gene_to_hit[gene][hit_index]
            hit.seq = seq[hit.qstart - 1 : hit.qend]
            if hit.frame < 0:
                hit.seq = bio_revcomp(hit.seq)
            
            gene_output[gene].append(hit)

        printv("Writing corrected reads.", args.verbose, 1)
        encoder = json.Encoder()
        results = 0
        for gene, hits in gene_output.items():
            if hits:
                results += len(hits)
                hits_db.put_bytes(f"gethits:{gene}", encoder.encode(hits))

        printv(f"Stored {results} corrected reads", args.verbose, 1)

def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for folder in args.INPUT:
        correct_folder(folder, args)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True