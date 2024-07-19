from collections import Counter, defaultdict
from functools import cached_property
from itertools import combinations, product
from math import ceil, floor
from multiprocessing import Pool
from os import listdir, makedirs, path
from pathlib import Path
from shutil import move, rmtree
import copy
from statistics import median
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

def insert_gaps(input_string, positions, offset, is_nt = False):
    if is_nt:
        input_triplets = [input_string[i:i+3] for i in range(0, len(input_string), 3)]
        
        for coord in positions:
            input_triplets.insert(offset + coord, "---")
            
        return "".join(input_triplets)
    else:
        input_string = list(input_string)
        
        for coord in positions:
            input_string.insert(offset + coord, "-")

        return "".join(input_string)
  
def int_first_id(x):
    return int(x.split("_")[0])

def scan_sequence(starting_id, left, frame, head_to_seq):
    coords = []
    ids = [starting_id]
    this_seq = head_to_seq[starting_id]
    start, end = 0, len(this_seq)
    coords.append((start, end))

    i = 0
    while True:
        if left:
            i -= 1
            if head_to_seq[starting_id + i][-250:] != this_seq[:250]:
                break
            
            this_seq = head_to_seq[starting_id + i][:-250] + this_seq
        else:
            i += 1
            
            if head_to_seq[starting_id + i][:250] != this_seq[-250:]:
                break
            
            this_seq += head_to_seq[starting_id + i][250:]
            
        ids.append(starting_id + i)
        start = end
        end = len(this_seq)
        coords.append((start, end))
        if end >= 15000:
            break
        
    if frame < 0:
        this_seq = bio_revcomp(this_seq)
        
    if left:
        coords = coords[::-1]
        
    id_to_coords = {id: coord for id, coord in zip(ids, coords)}
    
    return id_to_coords, this_seq

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
  
  
def filter_ref_consensus(log_output, ref_count, ref_consensus, ref_gaps, gap_start, gap_end, coverage_thresh):
    log_output.append("Reference seqs:")
    ref_seqs = []
    for y in range(ref_count):
        this_seq = "".join(ref_consensus[x][y] for x in range(gap_start, gap_end) if x not in ref_gaps)
        
        this_seq_coverage = 1 - ((this_seq.count("-") + this_seq.count(" ")) / len(this_seq))
        ref_seqs.append((this_seq, this_seq_coverage))
        
        
    target_coverage = coverage_thresh * max(ref_seq[1] for ref_seq in ref_seqs)
    ref_indices = [i for i, ref_seq in enumerate(ref_seqs) if ref_seq[1] >= target_coverage]
    ref_seqs = [seq for seq, cov in ref_seqs if cov >= target_coverage]
    
    this_consensus = {}
    for x in range(gap_start, gap_end):
        if x in ref_gaps:
            continue
        this_consensus[x] = [ref_consensus[x][y] for y in ref_indices if ref_consensus[x][y] != " " and ref_consensus[x][y] != "-"]
  
    return this_consensus, ref_seqs


def count_ref_gaps(gap_start, gap_end, ref_consensus, ref_gap_thresh):
    ref_gaps = set()
    insert_at = []
    longest_consecutive = None
    consecutive_non_gap = 0
    for i, x in enumerate(range(gap_start, gap_end)):
        if ref_consensus[x].count("-") / len(ref_consensus[x]) >= ref_gap_thresh:
            ref_gaps.add(x)
            insert_at.append(i)
            if longest_consecutive is None or consecutive_non_gap > longest_consecutive:
                longest_consecutive = consecutive_non_gap
            consecutive_non_gap = 0
        else:
            consecutive_non_gap += 1
            
    if longest_consecutive is None or consecutive_non_gap > longest_consecutive:
        longest_consecutive = consecutive_non_gap
        
    if longest_consecutive is None:
        longest_consecutive = 0
        
    return ref_gaps, insert_at, longest_consecutive

            
def scan_kmer(amount, log_output, splice_region, ref_gaps, flex, this_consensus, gap_start, gap_end, max_score, stop_penalty):
    kmer_size = (abs(amount) // 3) - len(ref_gaps) - flex
    # input(kmer_size)
    
    log_output.append("Raw splice region:")
    log_output.append(splice_region)
    log_output.append("")
    log_output.append("Translated Splice Sequences:")
    results = []
    
    all_posibilities = set()
    for i, cons in this_consensus.items():
        all_posibilities.update(cons)
        
    rows = [""] * (len(all_posibilities) + 1)  
    rows[0] += "N\t"
    for i, let in enumerate(all_posibilities):
        rows[i+1] += f"{let}\t"
    for x in range(gap_start, gap_end):
        if x in ref_gaps:
            continue
        rows[0] += str(x) + "\t"   
        for i, let in enumerate(all_posibilities):
            rows[i+1] += f"{this_consensus[x].count(let)}\t"
            
    ref_cols = [x for x in range(gap_start, gap_end) if x not in ref_gaps]
    
    highest_possible_score = 0

    for x in range(gap_start, gap_end):
        if x in ref_gaps:
            continue
        max_count = max(this_consensus[x].count(let) for let in all_posibilities)
        highest_possible_score += min(max_count, max_score)
    
    for frame in range(3):
        protein_seq = str(Seq(splice_region[frame:]).translate())
        log_output.append(protein_seq)
        for i in range(0, len(protein_seq) - kmer_size):
            kmer = protein_seq[i: i + kmer_size]
            for offset in range(0, flex + 1):
                kmer_score = 0
                for kmer_i in range(kmer_size):
                    kmer_score += min(this_consensus[ref_cols[kmer_i + offset]].count(kmer[kmer_i]), max_score)
                
                kmer_score -= (stop_penalty * kmer.count("*"))
                
                best_qstart = (i * 3) + frame
                best_qend = best_qstart + (kmer_size * 3)
                results.append((kmer, kmer_score, best_qstart, best_qend, frame))
    
    log_output.append("")
            
    return rows, highest_possible_score, results


def finalise_seq(node, rows, highest_possible_score, insert_at, results, gap_start, last_node, seq, ids_to_coords, id_count, log_output, minimum_aa):
    new_header = None
    new_aa_sequence = None
    new_nt_seq = None
    best_kmer, best_score, best_qstart, best_qend, best_frame = max(results, key=lambda x: x[1])
    if len(best_kmer) >= minimum_aa:
        best_frame += 1
        if node.frame < 0:
            best_frame = -best_frame
        target_score = best_score * 0.8
        
        other = [f"{res[0]} - {res[1]}" for res in results if res[1] >= target_score and res[0] != best_kmer]
        
        ids_in_qstart = [id for id, range in ids_to_coords.items() if best_qstart >= range[0] and best_qstart <= range[1] or best_qend >= range[0] and best_qend <= range[1]]
        new_header_fields = node.header.split("|")
        final_ids = []
        for id in ids_in_qstart:
                
                count = id_count[id]
                if count == 0:
                    final_ids.append(str(id))
                else:
                    final_ids.append(f"{id}_{count}")
                    
                id_count[id] += 1
        if len(best_kmer) >= 10:
            log_output.append("Threshold: {}".format(round(highest_possible_score * 0.5)))
        else:
            log_output.append("Threshold: {}".format(round(highest_possible_score * 0.70)))
        log_output.append("Best match: {} - Score: {} - Highest possible score: {} - Other possible matches within 20% of score: {}".format(best_kmer, best_score, highest_possible_score, len(other)))
        log_output.append("\n".join(rows))
        if other:
            log_output.append("Other matches:")
            for o in other:
                log_output.append(o)
        

        if len(best_kmer) >= 10:
            if best_score >= round(highest_possible_score * 0.5):
                this_id = "&&".join(final_ids)
                new_header_fields[3] = f"NODE_{this_id}"
                new_header_fields[4] = str(best_frame)
                new_header_fields[5] = "1"
                new_header = "|".join(new_header_fields)
                
                best_kmer = insert_gaps(best_kmer, insert_at, 0)
                
                new_aa_sequence = ("-" * (gap_start)) + best_kmer
                new_aa_sequence += ("-" * (len(last_node.sequence) - len(new_aa_sequence)))
                
                nt_seq = seq[best_qstart: best_qend]
                nt_seq = insert_gaps(nt_seq, insert_at, 0, True)
                
                new_nt_seq = ("-" * (gap_start * 3)) + nt_seq
                new_nt_seq += "-" * (len(last_node.nt_sequence) - len(new_nt_seq))
            else:
                log_output.append("Failed score threshold")
        else:
            if best_score >= round(highest_possible_score * 0.7):
                this_id = "&&".join(final_ids)
                new_header_fields[3] = f"NODE_{this_id}"
                new_header_fields[4] = str(best_frame)
                new_header_fields[5] = "1"
                new_header = "|".join(new_header_fields)
                
                best_kmer = insert_gaps(best_kmer, insert_at, 0)
                
                new_aa_sequence = ("-" * (gap_start)) + best_kmer
                new_aa_sequence += ("-" * (len(last_node.sequence) - len(new_aa_sequence)))
                
                nt_seq = seq[best_qstart: best_qend]
                nt_seq = insert_gaps(nt_seq, insert_at, 0, True)
                
                new_nt_seq = ("-" * (gap_start * 3)) + nt_seq
                new_nt_seq += "-" * (len(last_node.nt_sequence) - len(new_nt_seq))
            else:
                log_output.append("Failed score threshold")
    else:
        log_output.append("Kmer does not have enough bp: {}/{}".format(len(best_kmer), minimum_aa)) 
            
    return new_header, new_aa_sequence, new_nt_seq

            
def scan_last_node(gap_start, gap_end, minimum_gap_bp, max_gap_bp, last_node, log_output, ref_consensus, head_to_seq, ref_count, ref_gap_thresh, leftright_ref_coverage, min_consec_char, id_count, new_aa, new_nt, max_score, stop_penalty, flex, minimum_aa):
    amount = (gap_end - gap_start) * 3
    log_output.append("Right trailing gap of {} with size {}".format(last_node.header, abs(amount)))
    
    if amount < minimum_gap_bp or amount >= max_gap_bp:
        log_output.append("Gap too small or too large\n")
        return
        
    last_id = list(map(int_first_id, last_node.header.split("|")[3].replace("NODE_", "").split("&&")))[-1]
    ids_to_coords, seq = scan_sequence(last_id, False, last_node.frame, head_to_seq)
    
    last_node_kmer = last_node.nt_sequence[last_node.start * 3: last_node.end * 3].replace("-", "")
    
    last_node_og_start = seq.find(last_node_kmer)
    if last_node_og_start == -1:
        log_output.append("Could not find exon on genome\n")
        return
    
    last_node_og_end = last_node_og_start + len(last_node_kmer)
    
    seq = seq[last_node_og_end:]
    
    ref_gaps, insert_at, longest_consecutive = count_ref_gaps(gap_start, gap_end, ref_consensus, ref_gap_thresh)
    if longest_consecutive < min_consec_char:
        # log_output.append("No non-gap region with 5 char\n")
        return
    
    this_consensus, ref_seqs = filter_ref_consensus(log_output, ref_count, ref_consensus, ref_gaps, gap_start, gap_end, leftright_ref_coverage)
    
    log_output.append("\n".join(ref_seqs))
    log_output.append("")
    
    rows, highest_possible_score, results = scan_kmer(amount, log_output, seq, ref_gaps, flex, this_consensus, gap_start, gap_end, max_score, stop_penalty)
    
    if results:
        new_header, new_aa_sequence, new_nt_seq = finalise_seq(last_node, rows, highest_possible_score, insert_at, results, gap_start, last_node, seq, ids_to_coords, id_count, log_output, minimum_aa)
        if new_header:
            new_aa.append((new_header, new_aa_sequence))
            new_nt.append((new_header, new_nt_seq))  
    else:
        log_output.append("No suitable kmer found")
        
    log_output.append("")   
            
def scan_first_node(gap_start, gap_end, minimum_gap_bp, max_gap_bp, first_node, log_output, ref_consensus, head_to_seq, ref_count, ref_gap_thresh, leftright_ref_coverage, min_consec_char, id_count, new_aa, new_nt, max_score, stop_penalty, flex, minimum_aa):
    amount = (gap_end - gap_start) * 3
    log_output.append("Left leading gap of {} with size {}".format(first_node.header, abs(amount)))
    if amount < minimum_gap_bp or amount >= max_gap_bp:
        log_output.append("Gap too small or too large\n")
        return
    
    first_id = list(map(int_first_id, first_node.header.split("|")[3].replace("NODE_", "").split("&&")))[0]
    ids_to_coords, seq = scan_sequence(first_id, True, first_node.frame, head_to_seq)
    first_node_kmer = first_node.nt_sequence[first_node.start * 3: first_node.end * 3].replace("-", "")
    
    first_node_og_start = seq.find(first_node_kmer)
    
    if first_node_og_start == -1:
        log_output.append("Could not find exon on genome\n")
        return
    
    seq = seq[:first_node_og_start]
    
    ref_gaps, insert_at, longest_consecutive = count_ref_gaps(gap_start, gap_end, ref_consensus, ref_gap_thresh)
        
    if longest_consecutive is None or longest_consecutive < min_consec_char:
        # log_output.append("No non-gap region with 5 char\n")
        return
    this_consensus, ref_seqs = filter_ref_consensus(log_output, ref_count, ref_consensus, ref_gaps, gap_start, gap_end, leftright_ref_coverage)
    
    log_output.append("\n".join(ref_seqs))
    log_output.append("")
    
    rows, highest_possible_score, results = scan_kmer(amount, log_output, seq, ref_gaps, flex, this_consensus, gap_start, gap_end, max_score, stop_penalty)
    
    if results:
        new_header, new_aa_sequence, new_nt_seq = finalise_seq(first_node, rows, highest_possible_score, insert_at, results, gap_start, first_node, seq, ids_to_coords, id_count, log_output, minimum_aa)
        if new_header:
            new_aa.append((new_header, new_aa_sequence))
            new_nt.append((new_header, new_nt_seq))
    else:
        log_output.append("No suitable kmer found")
        
    log_output.append("")
            
            
def align_and_trim_seq(node_a, node_b, genomic_sequence):
    node_a_kmer = node_a.nt_sequence[node_a.start * 3: node_a.end * 3]
    node_b_kmer = node_b.nt_sequence[node_b.start * 3: node_b.end * 3]

    node_a_internal_gap = [i for i, let in enumerate(node_a_kmer) if let == "-"]
    node_b_internal_gap = [i for i, let in enumerate(node_b_kmer) if let == "-"]
    node_a_len = len(node_a_kmer)
    node_b_len = len(node_b_kmer)
    
    node_a_kmer = node_a_kmer.replace("-", "")
    node_b_kmer = node_b_kmer.replace("-", "")

    node_a_og_start = genomic_sequence.find(node_a_kmer)
    node_b_og_start = genomic_sequence.find(node_b_kmer)
    
    if node_a_og_start == -1 or node_b_og_start == -1:
        return None
    
    genomic_sequence = insert_gaps(genomic_sequence, node_a_internal_gap, node_a_og_start)
    genomic_sequence = insert_gaps(genomic_sequence, node_b_internal_gap, node_b_og_start+len(node_a_internal_gap))

    start_of_a = node_a_og_start    
    end_of_a = start_of_a + node_a_len
    start_of_b = node_b_og_start + len(node_a_internal_gap)
    end_of_b = start_of_b + node_b_len
    
    if end_of_a < start_of_b:
        splice_region = genomic_sequence[end_of_a: start_of_b]
    else:
        splice_region = genomic_sequence[end_of_b: start_of_a]
          
    return splice_region  
            
def reverse_pwm_splice(aa_nodes, cluster_sets, ref_consensus, head_to_seq, log_output, ref_count, ref_median_start, ref_median_end):
    new_aa = []
    new_nt = []
    id_count = Counter()
    
    flex = 1
    max_score = 100
    stop_penalty = 10
    ref_coverage_thresh = 0.7
    leftright_ref_coverage = 0.8
    minimum_gap_bp = 15
    max_gap_bp = 180
    ref_gap_thresh = 0.7
    min_consec_char = 5
    minimum_aa = 5
    
    for node in aa_nodes:
        for id in node_to_ids(node.header.split("|")[3]):
            id_count[id] += 1

    for cluster_set in cluster_sets:
        aa_subset = [node for node in aa_nodes if (cluster_set is None or within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0))]
        if aa_subset[0].frame < 0:
            aa_subset.sort(key=lambda x: x.start, reverse=True)
        else:
            aa_subset.sort(key=lambda x: x.start)
            
        if aa_subset[0].frame < 0:
            first_node = aa_subset[-1]
            last_node = aa_subset[0]
        else:
            first_node = aa_subset[0]
            last_node = aa_subset[-1]
            
        gap_start = ref_median_start
        gap_end = first_node.start
        
        scan_first_node(gap_start, gap_end, minimum_gap_bp, max_gap_bp, first_node, log_output, ref_consensus, head_to_seq, ref_count, ref_gap_thresh, leftright_ref_coverage, min_consec_char, id_count, new_aa, new_nt, max_score, stop_penalty, flex, minimum_aa)
        
        gap_start = last_node.end
        gap_end = ref_median_end
        scan_last_node(gap_start, gap_end, minimum_gap_bp, max_gap_bp, last_node, log_output, ref_consensus, head_to_seq, ref_count, ref_gap_thresh, leftright_ref_coverage, min_consec_char, id_count, new_aa, new_nt, max_score, stop_penalty, flex, minimum_aa)        
            
        for i in range(1, len(aa_subset)):
            node_a = aa_subset[i - 1]
            node_b = aa_subset[i]
            
            overlap = get_overlap(node_a.start * 3, node_a.end * 3, node_b.start * 3, node_b.end * 3, -max_gap_bp)
            if not overlap:
                continue
            
            amount = overlap[1] - overlap[0]
            if amount >= 0:
                continue
            
            if abs(amount) < minimum_gap_bp:
                continue
            
            ids = set(map(int_first_id, node_a.header.split("|")[3].replace("NODE_", "").split("&&"))).union(set(map(int_first_id, node_b.header.split("|")[3].replace("NODE_", "").split("&&"))))
            genomic_range = list(range(min(ids), max(ids) + 1))
            
            gap_start = overlap[1] // 3
            gap_end = overlap[0] // 3
            
            ref_gaps, insert_at, longest_consecutive = count_ref_gaps(gap_start, gap_end, ref_consensus, ref_gap_thresh)
                
            if longest_consecutive < min_consec_char:
                # log_output.append("No non-gap region with 5 char\n")
                continue
            
            log_output.append("Gap between {} and {} of size {}".format(node_a.header, node_b.header, abs(amount)))
                    
            this_consensus, ref_seqs = filter_ref_consensus(log_output, ref_count, ref_consensus, ref_gaps, gap_start, gap_end, ref_coverage_thresh)
            
            log_output.append("\n".join(ref_seqs))
            log_output.append("")

            id_to_coords, genomic_sequence = generate_sequence(genomic_range, node_a.frame, head_to_seq)
                
            splice_region = align_and_trim_seq(node_a, node_b, genomic_sequence)
            
            if splice_region is None:
                log_output.append("Could not find coords of mapped sequences")
                continue
            
            if splice_region == "":
                log_output.append("Splice region empty\n")
                continue
            
            rows, highest_possible_score, results = scan_kmer(amount, log_output, splice_region, ref_gaps, flex, this_consensus, gap_start, gap_end, max_score, stop_penalty)
                       
            if results:
                new_header, new_aa_sequence, new_nt_seq = finalise_seq(node_a, rows, highest_possible_score, insert_at, results, gap_start, node_b, splice_region, id_to_coords, id_count, log_output, minimum_aa)
                if new_header:
                    new_aa.append((new_header, new_aa_sequence))
                    new_nt.append((new_header, new_nt_seq))
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
    ref_starts, ref_ends = [], []
    
    raw_aa = []
    raw_nt = []
    ref_count = 0
    
    for header, seq in parseFasta(path.join(input_aa, gene)):
        raw_aa.append((header, seq))
        if header.endswith('.'):
            ref_count += 1
            start, end = find_index_pair(seq, "-")
            ref_starts.append(start)
            ref_ends.append(end)
                
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
        
    ref_median_start = floor(median(ref_starts))
    ref_median_end = ceil(median(ref_ends))
        
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

    max_gap_bp_size = round(len(aa_nodes[0].sequence) * 0.3) # Half MSA length

    clusters, _ = cluster_ids(ids, 100, max_gap_bp_size, reference_cluster_data, req_seq_coverage=0) #TODO: Make distance an arg
    
    if clusters:
        cluster_sets = [set(range(a, b+1)) for a, b, _ in clusters]
       
    new_nt, new_aa = reverse_pwm_splice(aa_nodes, cluster_sets, ref_consensus, head_to_seq, log_output, ref_count, ref_median_start, ref_median_end)
    
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
    
    genes = [i for i in listdir(input_aa_path) if ".fa"]
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
