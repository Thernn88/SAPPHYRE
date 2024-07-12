from collections import Counter, defaultdict
from functools import cached_property
from itertools import combinations, product
from math import ceil
from multiprocessing import Pool
from os import listdir, makedirs, path
from pathlib import Path
from shutil import move, rmtree
import copy
from tempfile import NamedTemporaryFile
import warnings
from msgspec import Struct, json
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    constrained_distance,
    bio_revcomp,
)
from .directional_cluster import cluster_ids, within_distance, node_to_ids, quick_rec
from wrap_rocks import RocksDB
from Bio.Seq import Seq
from Bio import BiopythonWarning

from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta

class NODE(Struct):
    header: str
    frame: int
    sequence: str
    nt_sequence: str
    start: int
    end: int

def insert_gaps(input_string, positions, offset):
    input_string = list(input_string)
    
    for coord in positions:
        input_string.insert(offset + coord, "-")

    return "".join(input_string)
  
def int_first_id(x):
    return int(x.split("_")[0])
  
def generate_sequence(ids, frame, head_to_seq):
    id_to_coords = {}
    
    prev_og = head_to_seq[ids[0]]
    start, end = 0, len(prev_og)
    id_to_coords[ids[0]] = (start, end)
    for i, child in enumerate(ids):
        if i == 0:
            continue
        prev_og += head_to_seq[child][250:]
        start = end
        end = len(prev_og)
        id_to_coords[child] = (start, end)

    if frame < 0:
        prev_og = bio_revcomp(prev_og)
        
    return id_to_coords, prev_og
  
            
def reverse_pwm_splice(aa_nodes, cluster_sets, ref_consensus, head_to_seq, log_output, minimum_gap = 30, max_gap = 180, ref_gap_thresh = 0.5, majority_gaps = 0.33):
    new_aa = []
    new_nt = []
    id_count = Counter()
    
    for node in aa_nodes:
        for id in node_to_ids(node.header.split("|")[3]):
            id_count[id] += 1

    for cluster_set in cluster_sets:
        aa_subset = [node for node in aa_nodes if (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
        aa_subset.sort(key=lambda x: x.start)
            
        for i in range(1, len(aa_subset)):
            node_a = aa_subset[i - 1]
            node_b = aa_subset[i]
            
            overlap = get_overlap(node_a.start * 3, node_a.end * 3, node_b.start * 3, node_b.end * 3, -max_gap)
            if not overlap:
                continue
            
            amount = overlap[1] - overlap[0]
            if amount >= 0:
                continue
            
            if abs(amount) < minimum_gap:
                continue
            
            log_output.append("Gap between {} and {} of size {}".format(node_a.header, node_b.header, abs(amount)))
            
            ids = set(map(int_first_id, node_a.header.split("|")[3].replace("NODE_", "").split("&&"))).union(set(map(int_first_id, node_b.header.split("|")[3].replace("NODE_", "").split("&&"))))
            genomic_range = list(range(min(ids), max(ids) + 1))
            
            gap_start = overlap[1]
            gap_end = overlap[0]
            
            ref_gaps = []
            for x, y in enumerate(range(gap_start//3, gap_end//3)):
                if ref_consensus[y].count("-") / len(ref_consensus[y]) >= ref_gap_thresh:
                    ref_gaps.append(x)
                    
            # if len(ref_gaps) / (gap_end//3 - gap_start//3) >= majority_gaps:
            #     log_output.append("Too many gaps in reference")
            #     continue
            # if "NODE_521670&&521671" in node_a.header:
            #     print(ref_gaps)
                
            # x = ["N"] * (gap_end//3 - gap_start//3)
            # for gap in ref_gaps:
            #     x[gap] = "-"
                
            # if "NODE_521670&&521671" in node_a.header:
            #     print("".join(x))

            id_to_coords, genomic_sequence = generate_sequence(genomic_range, node_a.frame, head_to_seq)
                
            node_a_kmer = node_a.nt_sequence[node_a.start * 3: node_a.end * 3]
            node_b_kmer = node_b.nt_sequence[node_b.start * 3: node_b.end * 3]
            
            node_a_internal_gap = [i for i, let in enumerate(node_a_kmer) if let == "-"]
            node_b_internal_gap = [i for i, let in enumerate(node_b_kmer) if let == "-"]
            
            node_a_len = len(node_a_kmer)
            
            node_a_kmer = node_a_kmer.replace("-", "")
            node_b_kmer = node_b_kmer.replace("-", "")
            
            node_a_og_start = genomic_sequence.find(node_a_kmer)
            node_b_og_start = genomic_sequence.find(node_b_kmer)
            
            if node_a_og_start == -1 or node_b_og_start == -1:
                log_output.append("Unable to find genomic coords\n")
                continue
            
            genomic_sequence = insert_gaps(genomic_sequence, node_a_internal_gap, node_a_og_start)
            genomic_sequence = insert_gaps(genomic_sequence, node_b_internal_gap, node_b_og_start)

            end_of_a = node_a_og_start + node_a_len
            start_of_b = node_b_og_start
            
            splice_region = genomic_sequence[end_of_a: start_of_b + 3]
            if splice_region == "":
                log_output.append("Splice region empty\n")
                continue
             
            kmer_size = (abs(amount) // 3) - len(ref_gaps)
            # input(kmer_size)
            
            log_output.append("Reference seqs:")
            for x in range(len(ref_consensus[0])):
                this_line = ""
                for y in range(gap_start//3, gap_end//3):
                    this_line += ref_consensus[y][x]
                log_output.append(this_line)
            log_output.append("")
            log_output.append("Raw splice region:")
            log_output.append(splice_region)
            log_output.append("")
            log_output.append("Translated Splice Sequences:")
            best_kmer = None
            best_score = None
            results = []
            for frame in range(3):
                protein_seq = str(Seq(splice_region[frame:]).translate())
                log_output.append(protein_seq)
                for i in range(0, len(protein_seq) - kmer_size):
                    kmer = protein_seq[i: i + kmer_size]
                    kmer = insert_gaps(kmer, ref_gaps, 0)
                    kmer_score = sum(ref_consensus[i].count(let) for i, let in enumerate(kmer, gap_start//3))
                    
                    if kmer_score == 0:
                        continue
                    
                    best_qstart = (i * 3) + frame
                    best_qend = best_qstart + (kmer_size * 3)
                    results.append((kmer, kmer_score, best_qstart, best_qend, frame))
            
            log_output.append("")
                       
            if results:
                best_kmer, best_score, best_qstart, best_qend, best_frame = max(results, key=lambda x: x[1])
                best_frame += 1
                if node_a.frame < 0:
                    best_frame = -best_frame
                target_score = best_score * 0.9
                
                other = [f"{res[0]} - {res[1]}" for res in results if res[1] >= target_score and res[0] != best_kmer]
                
                ids_in_qstart = [id for id, range in id_to_coords.items() if best_qstart >= range[0] and best_qstart <= range[1] or best_qend >= range[0] and best_qend <= range[1]]
                new_header_fields = node_a.header.split("|")
                final_ids = []
                for id in ids_in_qstart:
                    
                    count = id_count[id]
                    if count == 0:
                        final_ids.append(str(id))
                    else:
                        final_ids.append(f"{id}_{count}")
                        
                    id_count[id] += 1
                    
                log_output.append("Best match: {} - Score: {} - Other possible matches within 10% of score: {}".format(best_kmer, best_score, len(other)))
                if other:
                    log_output.append("Other matches:")
                    for o in other:
                        log_output.append(o)
                    
                this_id = "&&".join(final_ids)
                new_header_fields[3] = f"NODE_{this_id}"
                new_header_fields[4] = str(best_frame)
                new_header_fields[5] = "1"
                new_header = "|".join(new_header_fields)
                
                new_aa_sequence = ("-" * (gap_start//3)) + best_kmer
                new_aa_sequence += ("-" * (len(node_a.sequence) - len(new_aa_sequence)))
                
                nt_seq = splice_region[best_qstart: best_qend]
                new_nt_seq = "-" * gap_start + nt_seq 
                new_nt_seq += "-" * (len(node_a.nt_sequence) - len(nt_seq))
                
                new_aa.append((new_header, new_aa_sequence))
                new_nt.append((new_header, new_nt_seq))
            else:
                log_output.append("No suitable kmer found")
            log_output.append("")
                
                
    return new_nt, new_aa


def do_genes(genes, input_aa, input_nt, seq_source, out_aa_path, out_nt_path):
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    batch_result = []
    head_to_seq = {int(head): seq for head, seq in parseFasta(seq_source)}
    for gene in genes:
        print(gene)
        batch_result.extend(do_gene(gene, input_aa, input_nt, head_to_seq, out_aa_path, out_nt_path))
        
    return batch_result

def do_gene(gene, input_aa, input_nt, head_to_seq, out_aa_path, out_nt_path):
    aa_nodes = []
    reference_cluster_data = set()
    ref_consensus = defaultdict(list)
    
    raw_aa = []
    raw_nt = []
    
    for header, seq in parseFasta(path.join(input_aa, gene)):
        raw_aa.append((header, seq))
        if header.endswith('.'):
            start, end = find_index_pair(seq, "-")
            for i, bp in enumerate(seq):
                if i >= start and i <= end:
                    if bp != "-":
                        reference_cluster_data.add(i)
                    ref_consensus[i].append(bp)
                else:
                    ref_consensus[i].append(" ")
            continue

        frame = int(header.split("|")[4])
        aa_nodes.append(NODE(header, frame, seq, None, *find_index_pair(seq, "-")))
        
        
    nt_sequences = {}
    for header, seq in parseFasta(path.join(input_nt, gene.replace(".aa.", ".nt."))):
        nt_sequences[header] = seq
        raw_nt.append((header, seq))
        
    for node in aa_nodes:
        node.nt_sequence = nt_sequences[node.header]
        
    log_output = []
        
    cluster_sets = [None]
    ids = []
    for node in aa_nodes:
        start, end = find_index_pair(node.sequence, "-")
        ids.append(quick_rec(node.header.split("|")[3], node.frame, node.sequence, start, end))

    max_gap_size = round(len(aa_nodes[0].sequence) * 0.3) # Half MSA length

    clusters, _ = cluster_ids(ids, 100, max_gap_size, reference_cluster_data, req_seq_coverage=0) #TODO: Make distance an arg
    
    if clusters:
        cluster_sets = [set(range(a, b+1)) for a, b, _ in clusters]
       
    new_nt, new_aa = reverse_pwm_splice(aa_nodes, cluster_sets, ref_consensus, head_to_seq, log_output)
    
    aa_seqs = raw_aa+new_aa
    
    aa_references = [i for i in aa_seqs if i[0].endswith('.')]
    aa_candidates = [i for i in aa_seqs if not i[0].endswith('.')]
    
    aa_candidates.sort(key=lambda x: int(x[0].split("|")[3].split("&&")[0].split("_")[1]))
    
    aa_header_order = [i[0] for i in aa_candidates]
    nt_sequences = {header: seq for header, seq in raw_nt+new_nt}
    
    nt_seqs = [(header, nt_sequences[header]) for header in aa_header_order]

    writeFasta(path.join(out_aa_path, gene), aa_references+aa_candidates)
    writeFasta(path.join(out_nt_path, gene.replace(".aa.", ".nt.")), nt_seqs)
                    
    return log_output    
                        
def do_folder(folder, args):
    print(folder)
    tk = TimeKeeper(KeeperMode.DIRECT)
    rocksdb_path = path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(rocksdb_path)
    
    head_to_seq = get_head_to_seq(nt_db)
    
    input_aa_path = path.join(folder, "trimmed", "aa")
    input_nt_path = path.join(folder, "trimmed", "nt")
    
    output_folder = path.join(folder, "motif")
    out_aa_path = path.join(output_folder, "aa")
    out_nt_path = path.join(output_folder, "nt")
    
    if path.exists(output_folder):
        rmtree(output_folder)
    makedirs(out_aa_path, exist_ok=True)
    makedirs(out_nt_path, exist_ok=True)
    
    head_to_seq_source = NamedTemporaryFile(dir=gettempdir(), prefix="seqs_")
    writeFasta(head_to_seq_source.name, head_to_seq.items())
    
    genes = [i for i in listdir(input_aa_path) if ".fa" in i]
    per_batch = ceil(len(genes) / args.processes)
    arguments = [(genes[i: i+per_batch], input_aa_path, input_nt_path, head_to_seq_source.name, out_aa_path, out_nt_path) for i in range(0, len(genes), per_batch)]
    new_seqs = []
    if args.processes == 1:
        for arg in arguments:
            new_seqs.extend(do_genes(*arg))
    else:
        with Pool(args.processes) as pool:
            batches = pool.starmap(do_genes, arguments)
        for batch in batches:
            new_seqs.extend(batch)
            
    motif_log = path.join(output_folder, "motif.log")
    with open(motif_log, "w") as f:
        f.write("\n".join(new_seqs))
        
    del head_to_seq_source
        
    print(f"Done! Took {tk.differential():.2f}s")
    
    return True

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


def main(args):
    success = False
    if isinstance(args.INPUT, list):
        success = all(do_folder(folder, args) for folder in args.INPUT)
    elif isinstance(args.INPUT, str):
        success = do_folder(args.INPUT, args)
    return success


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
