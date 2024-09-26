from collections import defaultdict
from itertools import combinations
from math import ceil
from pathlib import Path
from shutil import rmtree
import subprocess
from tempfile import NamedTemporaryFile
from time import time
from wrap_rocks import RocksDB
from os import listdir, makedirs, mkdir, path
from .utils import gettempdir, parseFasta, printv, writeFasta
from sapphyre_tools import (
    find_index_pair,
    get_overlap,
    entropy_filter,
)
from Bio.Seq import Seq
from multiprocessing import Pool
from msgspec import Struct, json
from .pal2nal import worker
import xxhash
from .hmmsearch import HmmHit as Hit

class ExonerateHit(Struct):
    node: str
    frame: int
    coding_coords: str
    score: float
    query_coords: str
    query_name: str
    seq: str

def get_diamondhits(
    rocks_hits_db: RocksDB,
) -> dict[str, list[Hit]]:
    """Returns a dictionary of gene to corresponding hits.

    Args:
    ----
        rocks_hits_db (RocksDB): RocksDB instance
        list_of_wanted_genes (set): Set of genes to filter by
    Returns:
        dict[str, list[Hit]]: Dictionary of gene to corresponding hits
    """
    present_genes = rocks_hits_db.get("getall:presentgenes")
    if not present_genes:
        printv("ERROR: No genes found in hits database", 0)
        printv("Please make sure Diamond completed successfully", 0)
        return None
    genes_to_process = present_genes.split(",")
    
    gene_based_results = []
    for gene in genes_to_process:
        gene_result = rocks_hits_db.get_bytes(f"gethmmhits:{gene}")
        if not gene_result:
            printv(
                f"WARNING: No hits found for {gene}. If you are using a gene list file this may be a non-issue",
                0,
            )
            continue
        
        gene_based_results.append((gene, gene_result))

    return gene_based_results



def do_gene(out_path, orthoset_raw_path, gene_name, this_hits, verbose):
    printv(f"Processing {gene_name}", verbose, 2)
    raw_path = path.join(orthoset_raw_path, gene_name+".fa")
    hits = json.decode(this_hits, type=list[Hit])
    
    gene_seqs = [
        (f"{hit.node}_{hit.frame}",hit.seq) for hit in hits
    ]
    
    if not gene_seqs:
        return gene_name, []
    
    this_out = path.join(out_path, f"{gene_name}.txt")
    if not path.exists(this_out) or path.getsize(this_out) == 0:
        with NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()) as f, open(this_out, "w") as f2:
            writeFasta(f.name, gene_seqs)
            command = [
                "exonerate",
                "--geneticcode", "1",
                "--ryo", '>%ti|%tcb/%tce|%s|%qcb/%qce|%qi\n%tcs\n',
                "--score", "50",
                #"--model", "protein2genome",
                #"--showcigar", "no",
                "--showvulgar", "no",
                "--showalignment", "no",
                "--verbose", "0",
                raw_path,
                f.name
            ]
            
            subprocess.run(command, stdout=f2)
            
            if path.getsize(this_out) == 0:
                return gene_name, []
        
    out = defaultdict(list)
    for header, seq in parseFasta(this_out, True):
        node, coding_coords, score, query_coords, query_name = header.split("|")
        node, frame = map(int, node.split("_"))
        score = float(score)
        out[(node, frame)].append(ExonerateHit(node,frame,coding_coords, score, query_coords, query_name, seq))
    
    result = []
    for hit in hits:
        this_exonerate_hits = out.get((hit.node, hit.frame), [])
        if not this_exonerate_hits:
            continue
        this_best_hit = max(this_exonerate_hits, key=lambda x: x.score)
        
        coding_start, coding_end = map(int, this_best_hit.coding_coords.split("/"))
        
        # if frame > 0:
        
        frame = 1
        if coding_start > coding_end:
            coding_start, coding_end = coding_end, coding_start
            frame = -1
        
        qstart = hit.qstart + coding_start
        qend = hit.qstart + coding_end
        frame = (qstart % 3 + 1) * frame
        
        new_hit = Hit(
            node = this_best_hit.node,
            score = this_best_hit.score,
            frame = frame,
            evalue = hit.evalue,
            qstart = qstart,
            qend = qend,
            gene = gene_name,
            query = this_best_hit.query_name,
            uid = hit.uid,
            refs = [],
            seq = this_best_hit.seq
        )
        result.append(new_hit)
            
    return gene_name, result
        
def do_folder(folder, args):
    hits_db_path = path.join(folder, "rocksdb", "hits")
    hits_db = RocksDB(hits_db_path)
    
    rocks_db_path = path.join(folder, "rocksdb", "sequences", "nt")
    rocksdb_db = RocksDB(str(rocks_db_path))
    is_genome = rocksdb_db.get("get:isgenome")
    if is_genome == "True":
        printv("Running exonerate.", args.verbose, 1)
    else:
        return True
    
    orthoset_path = path.join(args.orthoset_input, args.orthoset)
    orthoset_raw_path = path.join(orthoset_path, "raw")
    
    transcripts_mapped_to = get_diamondhits(
        hits_db,
    )
    
    exonerate_folder = path.join(folder, "exonerate")
    makedirs(exonerate_folder, exist_ok=True)
    
    # result = [do_gene(orthoset_raw_path, gene.split(".")[0], hits) for gene, hits in transcripts_mapped_to]
    
    if args.processes <= 1:
        result = [do_gene(exonerate_folder, orthoset_raw_path, gene.split(".")[0], hits, args.verbose) for gene, hits in transcripts_mapped_to]
    else:
        with Pool(args.processes) as p:
            result = p.starmap(do_gene, [(exonerate_folder, orthoset_raw_path, gene.split(".")[0], hits, args.verbose) for gene, hits in transcripts_mapped_to])
    
    for gene, hits in result:
        hits_db.put_bytes(f"exonerate:{gene}", json.encode(hits))
        
    return True
        
def main(args):
    # global_time = TimeKeeper(KeeperMode.DIRECT)
    # print(f"Processing: {path.basename(args.INPUT)}")
    for folder in args.INPUT:
        success = do_folder(folder, args)
    # printv(f"Done! Took {global_time.differential():.2f} seconds", args.verbose, 0)
    return success