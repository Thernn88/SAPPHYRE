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

    hits_result = rocks_hits_db.get_bytes(f"getall:hits")
    this_hits = json.decode(hits_result, type=list[Hit])


    return this_hits


def N_trim(parent_sequence: str):
    """
    Breaks sequence into chunks at Ns, and yields chunks that are longer than the minimum sequence length.
    """
    if "N" in parent_sequence:
        # Get N indices and start and end of sequence
        indices = (
            [0]
            + [i for i, ltr in enumerate(parent_sequence) if ltr == "N"]
            + [len(parent_sequence)]
        )

        for i in range(0, len(indices) - 1):
            start = indices[i]
            end = indices[i + 1]

            raw_seq = parent_sequence[start + 1 : end]
            yield raw_seq
    else:
        yield parent_sequence

def correct_folder(folder, args):
    nt_db_path = os.path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(nt_db_path)

    if nt_db.get("get:isassembly") == "True" or nt_db.get("get:isgenome") == "True":
        printv("Skipping Error Connection", args.verbose, 0)
        return True


    original_posiitons = json.decode(nt_db.get("getall:original_positions"), type=dict[str, tuple[int, int, int|None]])
    original_inputs = json.decode(nt_db.get("getall:original_inputs"), type=list[str])

    hits_db = RocksDB(os.path.join(folder, "rocksdb", "hits"))


    diamond_hits = get_diamondhits(hits_db)

    positions_to_keep = defaultdict(dict)
    keep_set = defaultdict(set)

    hit_count = 0
    for i, hit in enumerate(diamond_hits):
        hit_count += 1
        file_index, line_index, n_index = original_posiitons[hit.node]

        keep_set[file_index].add(line_index)
        positions_to_keep[file_index].setdefault(line_index, []).append((hit.gene, i, n_index))

    
    header_to_old_pos = {}
    old_seq = {}
    printv(f"Found {hit_count} hits from diamond.", args.verbose, 1)

    with TemporaryDirectory(dir=gettempdir()) as tempdir:
        output = os.path.join(tempdir, "output.fastq")
        result = os.path.join(tempdir, "result")
        os.makedirs(result)
        with open(output, "w") as out:
            for file_index, file in enumerate(original_inputs):
                for i, (header, seq, qual) in enumerate(Fastq(file, build_index = False)):
                    if i in keep_set[file_index]:
                        header_to_old_pos[header] = file_index, i
                        old_seq[i] = seq
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
        output_indices = set()
        for i, (header, seq, qual) in enumerate(Fastq(result_file, build_index = False)):   
            from_file, from_index = header_to_old_pos[header]

            for gene, hit_index, n_index in positions_to_keep[from_file][from_index]:
                hit = diamond_hits[hit_index]
                output_indices.add(hit_index)

                if n_index is not None:
                    if "N" in seq:
                        seq = list(N_trim(seq))[n_index]
                    else:
                        #N was deleted
                        # print(header)
                        continue

                hit.seq = seq[hit.qstart - 1 : hit.qend]
                if hit.frame < 0:
                    hit.seq = bio_revcomp(hit.seq)

                gene_output[gene].append(hit)
        
        to_check = set(range(len(diamond_hits))) - output_indices
        for i in to_check:
            hit = diamond_hits[i]
            seq = old_seq[original_posiitons[hit.node][1]]
            
            gene_output[hit.gene].append(hit)

        printv("Writing corrected reads.", args.verbose, 1)
        encoder = json.Encoder()
        results = 0
        for gene, hits in gene_output.items():
            if hits:
                results += len(hits)
                hits_db.put_bytes(f"gethits:{gene}", encoder.encode(hits))

        hits_db.put("getall:presentgenes", ",".join(list(gene_output.keys())))

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