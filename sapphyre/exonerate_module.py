from collections import defaultdict
from itertools import combinations
from math import ceil
from pathlib import Path
from shutil import rmtree
import subprocess
from tempfile import NamedTemporaryFile
from time import time
from wrap_rocks import RocksDB
from os import mkdir, path
from .utils import gettempdir, parseFasta, printv, writeFasta
from sapphyre_tools import (
    entropy_filter,
)
from multiprocessing import Pool
from msgspec import json
from .diamond import ReporterHit as Hit

def get_head_to_seq(nt_db):
    """Get a dictionary of headers to sequences.

    Args:
    ----
        nt_db (RocksDB): The NT rocksdb database.
        recipe (list[str]): A list of each batch index.
    Returns:
    -------
        dict[str, str]: A dictionary of headers to sequences.
    """
    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
        
    head_to_seq = {}
    for i in recipe:
        lines = nt_db.get_bytes(f"ntbatch:{i}").decode().splitlines()
        head_to_seq.update(
            {
                int(lines[i][1:]): lines[i + 1]
                for i in range(0, len(lines), 2)
                if lines[i] != ""
            },
        )

    return head_to_seq


class Node:
    def __init__(self, head, seq, start, end, score, frame, ref_start, ref_end, ref_name) -> None:
        self.head = head
        self.seq = seq
        self.start = start
        self.end = end
        self.score = score
        self.frame = frame
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ref_name = ref_name



def calculate_new_frame(start, end, original_frame):
    """
    Calculate the new frame based on the original frame and the positions on the sequence.
    
    Parameters:
    start (int): The start position of the new sequence (inclusive).
    end (int): The end position of the new sequence (inclusive).
    original_frame (int): The original reading frame (can be positive or negative).
    
    Returns:
    int: The new frame adjusted for the given positions.
    """
    # Ensure the start is always less than or equal to end
    if start > end:
        start, end = end, start
    
    # Determine the shift in frame caused by moving from the start to end
    frame_shift = (end - start) % 3
    
    # Adjust the original frame
    if original_frame < 0:
        # Negative strand, we reverse the frame shift
        new_frame = original_frame - frame_shift
        if new_frame < -3:
            new_frame = -((abs(new_frame) % 3) + 1)
    else:
        # Positive strand, normal frame shift
        new_frame = original_frame + frame_shift
        if new_frame > 3:
            new_frame = (new_frame % 3) or 3  # Reset to 1, 2, or 3 if it goes beyond 3
    
    return new_frame


class exonerate:
    def __init__(self, folder, chomp_max_distance, orthoset_raw_path, exonerate_path, max_extend, target_to_taxon, debug, entropy_percent) -> None:
        self.folder = folder
        self.chomp_max_distance = chomp_max_distance
        self.orthoset_raw_path = orthoset_raw_path
        self.exonerate_path = exonerate_path
        self.max_extend = max_extend
        self.target_to_taxon = target_to_taxon
        self.debug = debug
        self.entropy_percent = entropy_percent
        
    def run(self, batches):
        batch_result = []
        additions = 0
        for gene, hits in batches:
            diamond_hits = json.decode(hits, type=list[Hit])
            gene_name = gene.split(".")[0]
            print(gene)
            
            if self.debug:
                f = open(path.join(self.exonerate_path, f"{gene_name}.fa"), "w")
            else:
                f = NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir())
                
            this_seqs = [(f"{hit.node}_{hit.frame}",hit.seq) for hit in diamond_hits]
            raw_path = path.join(self.orthoset_raw_path, gene_name+".fa")
            final_output = []
            with open(f"{gene_name}.txt", "w") as result:
                writeFasta(f.name, this_seqs)
                f.flush()
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
                                    
                subprocess.run(command, stdout=result)
                if path.getsize(result.name) == 0:
                    continue

                this_results = defaultdict(list)
                for header, sequence in parseFasta(result.name, True):
                    header, coords, score, ref_coords, ref_id = header.split("|")
                    
                    header, original_frame = header.split("_")
                    original_frame = int(original_frame)
                    start, end = map(int, coords.split("/"))
                    ref_start, ref_end = map(int, ref_coords.split("/"))
                    
                    frame = calculate_new_frame(start, end, original_frame)
                      
                    this_results[header].append(Node(int(header), sequence, start, end, int(score), frame, ref_start, ref_end, ref_id))
                    
                this_nodes = []
                for header, nodes in this_results.items():
                    for frame in range(-3, 4):
                        if frame == 0:
                            continue
                        
                        frame_nodes = [node for node in nodes if abs(node.frame) == frame]
                        if frame_nodes:
                            this_nodes.append(max(frame_nodes, key=lambda x: x.score))
        
                # kicked_ids = set()
                # overlap_min = 0.1
                # for node_a, node_b in combinations(this_nodes, 2):
                #     if node_a.head in kicked_ids or node_b.head in kicked_ids:
                #         continue
                #     overlap = get_overlap(node_a.start, node_a.end, node_b.start, node_b.end, 1)
                #     if overlap:
                #         amount = overlap[1] - overlap[0]
                #         percent = amount / min((node_a.end - node_a.start), (node_b.end - node_b.start))

                        
                #         if percent > overlap_min:
                #             if node_a.score > node_b.score:
                #                 kicked_ids.add(node_b.head)
                #             else:
                #                 kicked_ids.add(node_a.head)
                
                # this_nodes = [node for node in this_nodes if node.head not in kicked_ids]
                
                sequence_template = ">{}|{}\n{}\n"
                as_fasta = [sequence_template.format(node.head, node.frame, node.seq) for node in this_nodes]
                passing_fasta = entropy_filter(as_fasta, self.entropy_percent)
                passing_headers = set()
                for header, _ in passing_fasta:
                    passing_headers.add(header.strip()[1:])
                        

                for node in this_nodes:
                    if not f"{node.head}|{node.frame}" in passing_headers:
                        continue
                    
                    target = node.ref_name
                    key = f"{gene_name}|{target}"
                    
                    _, ref, _ = self.target_to_taxon[key]
                    additions += 1
                    final_output.append(Hit(
                        node=node.head,
                        frame=node.frame,
                        evalue=0,
                        qstart=node.start,
                        qend=node.end,
                        gene=gene_name,
                        query=ref,
                        uid=round(time()),
                        seq = node.seq
                    ))
            del f
            batch_result.append((gene, final_output))
        return batch_result, additions

       
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
    if not present_genes:
        printv("ERROR: No genes found in hits database", 0)
        printv("Please make sure Diamond completed successfully", 0)
        return None
    genes_to_process = present_genes.split(",")

    gene_based_results = []
    for gene in genes_to_process:
        gene_result = rocks_hits_db.get_bytes(f"gethits:{gene}")
        if not gene_result:
            printv(
                f"WARNING: No hits found for {gene}. If you are using a gene list file this may be a non-issue",
                0,
            )
            continue
        gene_based_results.append((gene, gene_result))

    return genes_to_process, gene_based_results
     

def do_folder(folder, args):
    exonerate_path = path.join(folder, "exonerate")
    if path.exists(exonerate_path):
        rmtree(exonerate_path)
    mkdir(exonerate_path)
        
    # exonerate_aa_path = path.join(exonerate_path, "aa")
    # mkdir(exonerate_aa_path)
        
    # exonerate_nt_path = path.join(exonerate_path, "nt")
    # mkdir(exonerate_nt_path)
    
    # aa_input = path.join(folder, "align")# "outlier", "excise", "aa")
    # nt_input = path.join(folder, "nt")#"outlier", "excise", "nt")

    # genes = []
    # for gene in listdir(aa_input):
    #     genes.append(gene)
    

    orthoset_path = path.join(args.orthoset_input, args.orthoset)
    orthoset_raw_path = path.join(orthoset_path, "raw")
    
    hits_db = RocksDB(path.join(folder, "rocksdb", "hits"))
    diamond_genes, transcripts_mapped_to = get_diamondhits(
        hits_db
    )
    
    per_batch = ceil(len(diamond_genes) / args.processes)
    batches = [transcripts_mapped_to[i : i + per_batch] for i in range(0, len(diamond_genes), per_batch)]
    
    orthoset_db_path = path.join(args.orthoset_input, args.orthoset, "rocksdb")
    orthoset_db = RocksDB(orthoset_db_path)
    
    target_to_taxon_raw = orthoset_db.get_bytes("getall:targetreference")
    target_to_taxon = json.decode(
        target_to_taxon_raw,
        type=dict[str, list[str | int]],
    )
    batch_result = []
    if args.processes <= 1:
        exonerate_obj = exonerate(folder, args.chomp_max_distance, orthoset_raw_path, exonerate_path, args.max_extend, target_to_taxon, args.debug, args.entropy_percent)
        for batch in batches:
            batch_result.append(exonerate_obj.run(batch))
    else:
        with Pool(args.processes) as pool:
            batch_result.extend(pool.starmap(
                exonerate(folder, args.chomp_max_distance, orthoset_raw_path, exonerate_path, args.max_extend, target_to_taxon, args.debug, args.entropy_percent).run,
                batches,
            ))
            
    encoder = json.Encoder()
    total_new_seqs = 0
    for batchr, additions in batch_result:
        total_new_seqs += additions
        for gene, hits in batchr:
            hits_db.put_bytes(f"get_exoneratehits:{gene}", encoder.encode(hits))

    printv(f"Found {total_new_seqs} sequences", args.verbose, 1)


    return True    

def main(args):
    success = False
    if isinstance(args.INPUT, list):
        success = all(do_folder(folder, args) for folder in args.INPUT)
    elif isinstance(args.INPUT, str):
        success = do_folder(args.INPUT, args)
    return success

if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
