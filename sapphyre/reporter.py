from __future__ import annotations

from argparse import Namespace
from collections import Counter, defaultdict, namedtuple
from itertools import combinations, groupby, product
from multiprocessing.pool import Pool
from operator import itemgetter
from os import makedirs, path, system
from shutil import rmtree
import subprocess
from tempfile import NamedTemporaryFile

from msgspec import json
from wrap_rocks import RocksDB
from xxhash import xxh3_64
from sapphyre_tools import (
    find_index_pair,
    get_overlap
)
from . import rocky
from .hmmsearch import HmmHit
from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta
from Bio.Seq import Seq
from parasail import blosum62, sw_trace_scan_profile_16, profile_create_16
import pyfamsa

MainArgs = namedtuple(
    "MainArgs",
    [
        "verbose",
        "processes",
        "debug",
        "INPUT",
        "orthoset_input",
        "orthoset",
        "compress",
        "minimum_bp",
        "gene_list_file",
        "keep_output",
    ],
)


# Extend hit with new functions
class Hit(HmmHit):#, frozen=True):
    raw_node: int = None
    header: str = None
    coords: list[tuple[int, int]] = None
    base_header: str = None
    aa_sequence: str = None
    parent: str = None
    strand: str = None
    parent_start: int = None
    chomp_start: int = None
    chomp_end: int = None
    children: list = []
    raw_children: list = []
    
    def get_merge_header(self):
        return "&&".join(map(str, [self.node] + self.children))

def get_hmmresults(
    rocks_hits_db: RocksDB,
    list_of_wanted_genes: list,
    is_genome,
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
        printv("Please make sure Hmmsearch completed successfully", 0)
        return None
    genes_to_process = list_of_wanted_genes or present_genes.split(",")
    
    if is_genome:
        decoder = json.Decoder(type=list[Hit])

    gene_based_results = []
    for gene in genes_to_process:
        gene_result = rocks_hits_db.get_bytes(f"gethmmhits:{gene}")
        if not gene_result:
            printv(
                f"WARNING: No hits found for {gene}. If you are using a gene list file this may be a non-issue",
                0,
            )
            continue
        if is_genome:
            gene_based_results.append((gene, decoder.decode(gene_result)))
        else:
            gene_based_results.append((gene, gene_result))

    return gene_based_results


def get_gene_variants(rocks_hits_db: RocksDB) -> dict[str, list[str]]:
    """Grabs the target variants from the hits database.

    Args:
    ----
        rocks_hits_db (RocksDB): RocksDB instance
    Returns:
        dict: Dictionary of gene to corresponding target variants
    """
    return json.decode(
        rocks_hits_db.get("getall:target_variants"),
        type=dict[str, list[str]],
    )


def get_toprefs(rocks_nt_db: RocksDB) -> list[str]:
    """Grabs the top references from the nt database.

    Args:
    ----
        rocks_nt_db (RocksDB): RocksDB instance
    Returns:
        list: List of top references
    """
    return (
        json.decode(rocks_nt_db.get("getall:valid_refs"), type=dict[str, list]),
        rocks_nt_db.get("get:isassembly") == "True",
        rocks_nt_db.get("get:isgenome") == "True",
    )


def translate_cdna(cdna_seq):
    """Translates Nucleotide sequence to Amino Acid sequence.

    Args:
    ----
        cdna_seq (str): Nucleotide sequence
    Returns:
        str: Amino Acid sequence
    """
    if len(cdna_seq) % 3 != 0:
        printv("WARNING: NT Sequence length is not divisable by 3", 0)

    # try:
    #     return translate(cdna_seq)
    # except:
    return str(Seq(cdna_seq).translate())

def get_core_sequences(
    gene: str,
    orthoset_db: RocksDB,
) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    """Returns the core reference sequences for a given gene.

    Args:
    ----
        gene (str): Gene name
        orthoset_db (RocksDB): RocksDB instance
    Returns:
        tuple: Tuple of core AA and NT sequences
    """
    core_seqs = json.decode(
        orthoset_db.get(f"getcore:{gene}"),
        type=dict[str, list[tuple[str, str, str]]],
    )
    return core_seqs["aa"]#, core_seqs["nt"]


def print_core_sequences(
    gene,
    core_sequences,
    target_taxon,
    top_refs,
):
    """Returns a filtered list of headers and sequences for the core sequences.

    Args:
    ----
        gene (str): Gene name
        core_sequences (list): List of core sequences
        target_taxon (str): Target taxa

    """
    result = []
    for taxon, taxa_id, seq in sorted(core_sequences):
        # Filter out non-target hits and variants
        if taxon not in top_refs:
            continue
        if target_taxon:
            if f"{gene}|{taxa_id}" not in target_taxon and taxa_id not in target_taxon: #TODO: Slow to handle both 
                continue  
        

        header = gene + "|" + taxon + "|" + taxa_id + "|."
        result.append((header, seq))

    return result


def print_unmerged_sequences(
    hits: list,
    is_assembly_or_genome: bool,
) -> tuple[dict[str, list], list[tuple[str, str]], list[tuple[str, str]]]:
    
    aa_result = []
    nt_result = []
    header_to_score = {}
    for hit in hits:
        # Write unique sequence
        aa_result.append((hit.header, hit.aa_sequence))
        nt_result.append((hit.header, hit.seq))

        if is_assembly_or_genome:
            header_to_score[hit.node] = max(header_to_score.get(hit.node, 0), hit.score)  

    return aa_result, nt_result, header_to_score


OutputArgs = namedtuple(
    "OutputArgs",
    [
        "gene",
        "list_of_hits",
        "aa_out_path",
        "debug_path",
        "taxa_id",
        "nt_out_path",
        "verbose",
        "compress",
        "target_taxon",
        "prepare_dupes",
        "top_refs",
        "minimum_bp",
        "debug",
        "is_assembly_or_genome",
        "is_genome",
    ],
)


def tag(sequences: list[tuple[str, str]], prepare_dupes: dict[str, dict[str, int]], reporter_dupes: dict[str, list[str]]) -> list[tuple[str, str]]:
    """Tags the output with the prepare dupes and the dupes.

    Args:
    ----
        sequences (list): List of sequences
        prepare_dupes (dict): Prepare dupes
        dupes (dict): Dupes
    Returns:
        list: Tagged output
    """
    output = []
    for header, sequence in sequences:
        this_node = header.split("|")[3].replace("NODE_","")
        nodes = list(map(lambda x: int(x.split("_")[0]), this_node.split("&&")))
        dupes = 1 + sum(prepare_dupes.get(node, 0) for node in nodes) + sum(prepare_dupes.get(child, 1) for child in reporter_dupes.get(this_node, []))
        header = f"{header}|{dupes}"
        output.append((header, sequence))

    return output


def merge_hits(hits: list[Hit], minimum_bp_overlap = 30) -> tuple[list[Hit], list[str]]:
    hits.sort(key = lambda x: x.raw_node)
    log = ["Merges for: "+hits[0].gene]
    
    merge_occured = True
    max_recursion = 15
    CLUSTER_DISTANCE = 300
    while merge_occured and max_recursion:
        max_recursion -= 1
        merge_occured = False
        
        clusters = []
        current_cluster = []
        current_index = None
        for i, hit in enumerate(hits):
            if hit is None:
                continue
            if current_index is not None:
                if hit.raw_node - current_index >= CLUSTER_DISTANCE:
                    clusters.append(current_cluster)
                    current_cluster = []
                
            current_cluster.append(i)
            current_index = hit.raw_node

        if current_cluster:
            clusters.append(current_cluster)
            
        for indices in clusters:
            for i, j in combinations(indices, 2):
                hit_a = hits[i]
                hit_b = hits[j]
                if hit_a is None or hit_b is None:
                    continue
                
                if hit_a.strand != hit_b.strand:
                    continue
                
                if get_overlap(hit_a.chomp_start, hit_a.chomp_end, hit_b.chomp_start, hit_b.chomp_end, minimum_bp_overlap) is None:
                    continue
                
                if not any(b - a <= 1 for a, b in product([hit_a.raw_node] + hit_a.raw_children, [hit_b.raw_node] + hit_b.raw_children)):
                    continue
                
                if abs(hit_a.chomp_start - hit_b.chomp_start) % 3 != 0:
                    continue # different frame same seq
                
                if hit_a.strand == "-":
                    gaps_a_start = max(hit_b.chomp_end - hit_a.chomp_end, 0)
                    gaps_b_start = max(hit_a.chomp_end - hit_b.chomp_end, 0)
                else:
                    gaps_a_start = max(hit_a.chomp_start - hit_b.chomp_start, 0)
                    gaps_b_start = max(hit_b.chomp_start - hit_a.chomp_start, 0)

                a_align_seq = ("-" * gaps_a_start) + hit_a.seq
                b_align_seq = "-" * gaps_b_start + hit_b.seq

                alignment_length = max(len(b_align_seq),len(a_align_seq))
                a_align_seq += "-" * (alignment_length - len(a_align_seq))
                b_align_seq += "-" * (alignment_length - len(b_align_seq))
                
                a_align_start, a_align_end = find_index_pair(a_align_seq, "-")
                b_align_start, b_align_end = find_index_pair(b_align_seq, "-")
                
                overlap = get_overlap(a_align_start, a_align_end, b_align_start, b_align_end, 1)
                if overlap is None:
                    continue
                
                a_kmer = a_align_seq[overlap[0]:overlap[1]]
                b_kmer = b_align_seq[overlap[0]:overlap[1]]
                
                if translate_cdna(a_kmer) != translate_cdna(b_kmer):
                    continue
                
                merged_seq = "".join([a_align_seq[i] if a_align_seq[i] != "-" else b_align_seq[i] for i in range(len(a_align_seq))])

                if hit_a.raw_node == hit_b.raw_node:
                    log.append("WARNING: {} and {} same base node merge".format(hit_a.node, hit_b.node))
                
                if hit_a.evalue != 0 and hit_b.evalue != 0:
                    hit_a.evalue = min(hit_a.evalue, hit_b.evalue)
                    hit_a.query = min((hit_a.query, hit_a.evalue), (hit_b.query, hit_b.evalue), key = lambda x: x[1])[0]
                elif hit_a.evalue == 0:
                    hit_a.evalue = hit_b.evalue
                    hit_a.query = hit_b.query
                
                hit_a.children.append(hit_b.node)
                hit_a.raw_children.append(hit_b.raw_node)
                
                hit_a.children.extend(hit_b.children)
                hit_a.raw_children.extend(hit_b.raw_children)
                
                hit_a.seq = merged_seq
                
                hit_a.parent_start = min(hit_a.parent_start, hit_b.parent_start)
                
                # Update coords
                hit_a.chomp_start = min(hit_a.chomp_start, hit_b.chomp_start)
                hit_a.chomp_end = max(hit_a.chomp_end, hit_b.chomp_end)
                
                hits[j] = None
                merge_occured = True

    for hit in hits:
        if hit is not None and hit.children:
            log.append(hit.get_merge_header())
    return [i for i in hits if i is not None], log

def translate_sequences(hits):
    for hit in hits:
        hit.aa_sequence = translate_cdna(hit.seq)


def check_minimum_bp(hits, minimum_bp):
    out_hits = []
    for hit in hits:
        if len(hit.seq) >= minimum_bp:
            out_hits.append(hit)

    return out_hits


def do_dupe_check(hits, header_template, is_assembly_or_genome, taxa_id):
    base_header_template = "{}_{}"
    header_mapped_x_times = Counter()
    base_header_mapped_already = {}
    seq_mapped_already = {}
    exact_hit_mapped_already = set()
    dupes = defaultdict(list)
    
    for i, hit in enumerate(hits):
        base_header = str(hit.node)
        hit.header = header_template.format(hit.gene, hit.query, taxa_id, hit.node, hit.frame)
        unique_hit = None
        if not is_assembly_or_genome:
            # Hash the NT sequence and the AA sequence + base header
            unique_hit = xxh3_64(base_header + hit.aa_sequence).hexdigest()
            nt_seq_hash = xxh3_64(hit.seq).hexdigest()
            # Filter and save NT dupes
            if nt_seq_hash in seq_mapped_already:
                mapped_to = seq_mapped_already[nt_seq_hash]
                dupes.setdefault(mapped_to, []).append(base_header)
                hits[i] = None
                continue
            seq_mapped_already[nt_seq_hash] = base_header

        # If the sequence is unique
        if is_assembly_or_genome or unique_hit not in exact_hit_mapped_already:
            # Remove subsequence dupes from same read
            if base_header in base_header_mapped_already:
                (
                    already_mapped_index,
                    already_mapped_hit,
                ) = base_header_mapped_already[base_header]
                # Dont kick if assembly
                if not is_assembly_or_genome:
                    already_mapped_sequence = already_mapped_hit.aa_sequence
                    if len(hit.aa_sequence) > len(already_mapped_sequence):
                        if already_mapped_sequence in hit.aa_sequence:
                            hits[already_mapped_index] = None
                            continue
                    else:
                        if hit.aa_sequence in already_mapped_sequence:
                            hits[i] = None
                            continue

                if base_header in header_mapped_x_times:
                    # Make header unique
                    current_count = header_mapped_x_times[base_header]
                    header_mapped_x_times[base_header] += 1
                    modified_header = base_header_template.format(base_header, current_count)
                    hit.node = modified_header
                    hit.header = header_template.format(hit.gene, hit.query, taxa_id, hit.node, hit.frame)

                    header_mapped_x_times[base_header] += 1
            else:
                base_header_mapped_already[base_header] = i, hit


            if base_header not in header_mapped_x_times:
                header_mapped_x_times[base_header] = 1
            exact_hit_mapped_already.add(unique_hit)
            
    return [i for i in hits if i is not None], dupes


def pairwise_sequences(hits, debug_fp, ref_seqs, min_gaps=10):
    ref_dict = {taxa: seq for taxa, _, seq in ref_seqs}
    internal_introns_removed = []
    intron_coordinates = defaultdict(list)
    
    # aligner = pyfamsa.Aligner(threads=1)
    for hit in hits:
        ref_seq = ref_dict[hit.query]
        # this_aa = hit.aa_sequence
        # profile = profile_create_16(this_aa, blosum62)
        # result = sw_trace_scan_profile_16(
        #     profile,
        #     ref_seq,
        #     2,
        #     1,
        # )
        # if not hasattr(result, "traceback"):
        #     continue #alignment failed, to investigate
        # this_aa, ref_seq = result.traceback.query, result.traceback.ref
        to_remove_final = set() 
        for alignment_method in ["clustalo", "mafft"]:
            with NamedTemporaryFile("w", dir=gettempdir()) as in_file,  NamedTemporaryFile("w", dir=gettempdir()) as out_file:
                in_file.write(f">hit\n{hit.aa_sequence}\n>query\n{ref_seq}\n")
                in_file.flush()
                
                if alignment_method == "clustalo":
                    cmd = ["clustalo", "-i", in_file.name, "-o", out_file.name, "--thread=1", "--force"]
                    subprocess.run(cmd, stdout=subprocess.DEVNULL)
                else:
                    system(
                        f"mafft --localpair --quiet --thread 1 --anysymbol '{in_file.name}' > '{out_file.name}'"
                    )
                
                aligned_seqs = dict(parseFasta(out_file.name, True))
        
            ref_seq = aligned_seqs["query"]
            this_aa = aligned_seqs["hit"]

            if debug_fp:
                debug_fp.write(f">{alignment_method}: {hit.header}\n{this_aa}\n>{hit.query}\n{ref_seq}\n\n")
            
            ref_start, ref_end = find_index_pair(ref_seq, "-")
            internal_ref_gaps = [i for i in range(ref_start, ref_end) if ref_seq[i] == "-"]
            
            # group consecutive gaps
            to_remove = {}
            
            groups = [(id, list(map(int,map(itemgetter(1),g)))) for id, g in groupby(enumerate(internal_ref_gaps),lambda x:x[0]-x[1])]
            merge_occured = True
            while merge_occured:
                
                merge_occured = False
                for i in range(len(groups)-1):
                    if groups[i] is None or groups[i+1] is None:
                        continue
                    
                    key_a, group_a = groups[i]
                    group_b = groups[i+1][1]
                    if min(group_b) - max(group_a) > 10:
                        continue
                    
                    groups[i] = (key_a, list(range(min(group_a), max(group_b))))
                    groups[i+1] = None
                    merge_occured = True
                    break
                    
                groups = [i for i in groups if i is not None]
            for k, group in groups:
                if len(group) >= min_gaps:
                    for i in group:
                        to_remove[i] = k

            current_non_aligned = 0
            to_remove_unaligned = defaultdict(list)
            for i, let in enumerate(this_aa):
                if let != "-":
                    current_non_aligned += 1
        
                if i in to_remove:
                    group = to_remove[i]
                    to_remove_unaligned[group].append(current_non_aligned)
                    
            for group in to_remove_unaligned.values():
                
                left_stops = []
                right_stops = []
                for i in range(0, max(group) + 1):
                    if hit.aa_sequence[i] == "*":
                        left_stops.append(i * 3)
                        
                for i in range(min(group) - 1, len(hit.aa_sequence)):
                    if hit.aa_sequence[i] == "*":
                        right_stops.append(i * 3)
                    
                middle_of_gap = group[len(group) // 2] * 3
                left_most_codon = [i for i in left_stops if max(i, middle_of_gap) - min(i, middle_of_gap) <= 30 + len(group)] # Within 30 bp of middle
                
                start_left_scan = left_most_codon[0] if left_most_codon else max(group) * 3
                
                right_most_codon = [i for i in right_stops if max(i, middle_of_gap) - min(i, middle_of_gap) <= 30 + len(group)] # Within 30 bp of middle
                
                start_right_scan = right_most_codon[-1] if right_most_codon else min(group) * 3

                gt_coords = []
                ag_coords = []
                for i in range(0, start_left_scan):
                    if hit.seq[i:i+2] == "GT":
                        gt_coords.append(i)
                
                for i in range(start_right_scan, len(hit.seq) - 2):
                    if hit.seq[i:i+2] == "AG":
                        ag_coords.append(i + 2)
                        
                gt_coords.reverse()
                
                intron = ""
                for gt_coord, ag_coord in product(gt_coords, ag_coords):
                    intron = hit.seq[gt_coord:ag_coord]
                    if gt_coord >= ag_coord:
                        continue
                    
                    if (ag_coord - gt_coord) % 3 != 0:
                        continue
                    
                    if ag_coord - gt_coord < 30:
                        continue

                    this_to_remove = set(range(gt_coord, ag_coord))
                    if this_to_remove.intersection(to_remove_final):
                        continue
                    
                    to_remove_final.update(this_to_remove)
                    break
                
                internal_introns_removed.append(f"{hit.header}\n{hit.query}\nRemoved group of size {len(group)} at {group[0]}-{group[-1]} on pairwise alignment\n{intron}\n")
            
            nt_seq = "".join([let for i, let in enumerate(hit.seq) if i not in to_remove_final])
            
            triplets = {nt_seq[i:i+3] for i in range(0, len(nt_seq), 3)}
            if not "TAG" in triplets and not "TAA" in triplets  and not "TGA" in triplets:
                break
            
         # TODO instead of traversing the whole list just check the first and last element of the group
        exons = []
        exon_start = None

        for i in range(hit.chomp_start, hit.chomp_end + 1):
            relative_index = i - hit.chomp_start  # Convert to zero-based index relative to the start coordinate

            if relative_index not in to_remove_final:  # Check if this index is not in the gaps
                if exon_start is None:
                    exon_start = i  # Start a new exon
            else:
                if exon_start is not None:
                    # If we're in a gap and an exon was ongoing, close it
                    exons.append((exon_start, i - 1))  # End the current exon
                    exon_start = None

        # Step 3: If the last exon continues to the end, close it
        if exon_start is not None:
            exons.append((exon_start, hit.chomp_end))
        
        hit.coords = exons
            
        to_remove_parent = set()
        for i in to_remove_final:
            to_remove_parent.add(i + hit.chomp_start - hit.parent_start)
        
        if to_remove_final:
            intron_coordinates[hit.header] = (hit.chomp_start - hit.parent_start - 1, hit.chomp_end - hit.parent_start, to_remove_final)
        
        hit.seq = nt_seq
    
    return internal_introns_removed, intron_coordinates


def detect_repeat_nucleotides(hits: list[Hit], min_bp = 9) -> list[Hit]:
    repeats_removed = []
    for hit in hits:
        before = hit.seq
        before = (hit.chomp_start, hit.chomp_end)
        change_made = False
        # single bp
        #start
        possible = hit.seq[0]
        repeat_bp = 1
        for i in range(1, len(hit.seq)):
            if hit.seq[i] == possible:
                repeat_bp += 1
            else:
                break
            
        if repeat_bp >= min_bp:
            repeat_bp = repeat_bp - (repeat_bp % 3)
            hit.seq = hit.seq[repeat_bp:]
            hit.aa_sequence = hit.aa_sequence[repeat_bp//3:]
            if hit.frame < 0: # Todo fix coords in splice or do this prior to coord flip to stop this jank
                hit.chomp_end -= repeat_bp
            else:
                hit.chomp_start += repeat_bp
            change_made = True
            
            
        # end
        possible = hit.seq[-1]
        repeat_bp = 1
        for i in range(len(hit.seq) - 2, -1, -1):
            if hit.seq[i] == possible:
                repeat_bp += 1
            else:
                break
        
        if repeat_bp >= min_bp:
            repeat_bp = repeat_bp - (repeat_bp % 3)
            hit.seq = hit.seq[:-repeat_bp]
            hit.aa_sequence = hit.aa_sequence[:-repeat_bp//3]
            if hit.frame < 0: # Todo fix coords in splice or do this prior to coord flip to stop this jank
                hit.chomp_start += repeat_bp
            else:
                hit.chomp_end -= repeat_bp
            change_made = True
        
        # two bp
        
        # start
        possible = hit.seq[:2]
        repeat_bp = 2
        for i in range(2, len(hit.seq) - 1, 2):
            if hit.seq[i:i+2] == possible:
                repeat_bp += 2
            else:
                break
            
        if repeat_bp >= min_bp:
            repeat_bp = repeat_bp - (repeat_bp % 3)
            hit.seq = hit.seq[repeat_bp:]
            hit.aa_sequence = hit.aa_sequence[repeat_bp//3:]
            if hit.frame < 0: # Todo fix coords in splice or do this prior to coord flip to stop this jank
                hit.chomp_end -= repeat_bp
            else:
                hit.chomp_start += repeat_bp
            change_made = True
            
        # end
        possible = hit.seq[-2:]
        repeat_bp = 2
        for i in range(len(hit.seq) - 4, -1, -2):
            if hit.seq[i:i+2] == possible:
                repeat_bp += 2
            else:
                break
            
        if repeat_bp >= min_bp:
            repeat_bp = repeat_bp - (repeat_bp % 3)
            hit.seq = hit.seq[:-repeat_bp]
            hit.aa_sequence = hit.aa_sequence[:-repeat_bp//3]
            if hit.frame < 0: # Todo fix coords in splice or do this prior to coord flip to stop this jank
                hit.chomp_start += repeat_bp
            else:
                hit.chomp_end -= repeat_bp
            change_made = True
            
        if change_made:
            repeats_removed.append(f"{hit.header}\n{before}\n{hit.seq}\n")

    return hits, repeats_removed

def merge_and_write(oargs: OutputArgs) -> tuple[str, dict, int]:
    """Merges, dedupes and writes the output for a given gene.

    Args:
    ----
        oargs (OutputArgs): Output arguments
    Returns:
        tuple:
            Tuple containing the gene name,
            a dict of removed duplicates and the number of sequences written
    """
    t_gene_start = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Doing output for: {oargs.gene}", oargs.verbose, 2)

    # Unpack the hits
    this_hits = oargs.list_of_hits if oargs.is_genome else json.decode(oargs.list_of_hits, type=list[Hit])
    gene_nodes = [hit.node for hit in this_hits]

    # Get reference sequences
    core_sequences = get_core_sequences(
        oargs.gene,
        rocky.get_rock("rocks_orthoset_db"),
    )
   
    this_aa_path = path.join(oargs.aa_out_path, oargs.gene + ".aa.fa")
    header_template = "{}|{}|{}|NODE_{}|{}"
        
    translate_sequences(this_hits)
    
    this_hits = check_minimum_bp(this_hits, oargs.minimum_bp)
        
    this_hits, this_gene_dupes = do_dupe_check(this_hits, header_template, oargs.is_assembly_or_genome, oargs.taxa_id)
    before_merge_count = len(this_hits)
        
    merge_log = []
    removed_introns = []
    repeats_removed = []
    intron_coordinates = {}
    if oargs.is_genome or False: # Set False to disable
        this_hits, repeats_removed = detect_repeat_nucleotides(this_hits)
        
        this_hits, merge_log = merge_hits(this_hits)
        for hit in this_hits:
            hit.header = header_template.format(hit.gene, hit.query, oargs.taxa_id, hit.get_merge_header(), hit.frame)

        # Refresh translation
        translate_sequences(this_hits)
        
        for hit in this_hits:
            hit.coords = [(hit.chomp_start, hit.chomp_end)]
        
        if oargs.debug:
            debug_fp = open(path.join(oargs.debug_path, oargs.gene + ".debug"), "w")
        else:
            debug_fp = None
            
        removed_introns, intron_coordinates = pairwise_sequences(this_hits, debug_fp, core_sequences)
        
    translate_sequences(this_hits)
        
    # Trim and save the sequences
    aa_output, nt_output, header_to_score = print_unmerged_sequences(
        this_hits,
        oargs.is_assembly_or_genome,
    )

    aa_output = tag(aa_output, oargs.prepare_dupes, this_gene_dupes)
    nt_output = tag(nt_output, oargs.prepare_dupes, this_gene_dupes)

    if aa_output:
        # If valid sequences were found, insert the present references
        aa_core_sequences = print_core_sequences(
            oargs.gene,
            core_sequences,
            oargs.target_taxon,
            oargs.top_refs,
        )
        # Write the output
        writeFasta(this_aa_path, aa_core_sequences + aa_output, oargs.compress)

        this_nt_path = path.join(oargs.nt_out_path, oargs.gene + ".nt.fa")
        writeFasta(this_nt_path, nt_output, oargs.compress)

    printv(
        f"{oargs.gene} took {t_gene_start.differential():.2f}s. Had {len(aa_output)} sequences",
        oargs.verbose,
        2,
    )
    return (
        oargs.gene,
        removed_introns,
        before_merge_count,
        header_to_score,
        gene_nodes,
        [(hit.parent, hit.get_merge_header(), hit.coords, hit.strand, hit.frame) for hit in this_hits],
        merge_log,
        intron_coordinates,
        repeats_removed
    )


def get_prepare_dupes(rocks_nt_db: RocksDB) -> dict[str, dict[str, int]]:
    prepare_dupe_counts = json.decode(
        rocks_nt_db.get("getall:gene_dupes"), type=dict[str, dict[str, int]]
    )
    return prepare_dupe_counts


def do_taxa(taxa_path: str, taxa_id: str, args: Namespace):
    """Main function for processing a given taxa.

    Args:
    ----
        path (str): Path to the taxa directory
        taxa_id (str): Taxa ID
        args (Namespace): Reporter arguments
    Returns:
        bool: True if the taxa was processed successfully, False otherwise
    """
    printv(f"Processing: {taxa_id}", args.verbose, 0)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    num_threads = args.processes

    # Grab gene list file if present
    if args.gene_list_file:
        with open(args.gene_list_file) as fp:
            list_of_wanted_genes = fp.read().splitlines()
    else:
        list_of_wanted_genes = []

    aa_out = "aa"
    nt_out = "nt"

    aa_out_path = path.join(taxa_path, aa_out)
    nt_out_path = path.join(taxa_path, nt_out)

    if not args.keep_output:
        if path.exists(aa_out_path):
            rmtree(aa_out_path)
        if path.exists(nt_out_path):
            rmtree(nt_out_path)

    makedirs(aa_out_path, exist_ok=True)
    makedirs(nt_out_path, exist_ok=True)

    printv(
        f"Initialized databases. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reciprocal Hmmer hits.",
        args.verbose,
    )
    
    top_refs, is_assembly, is_genome = get_toprefs(rocky.get_rock("rocks_nt_db"))

    transcripts_mapped_to = get_hmmresults(
        rocky.get_rock("rocks_hits_db"),
        list_of_wanted_genes,
        is_genome,
    )

    printv(
        f"Got hits. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reference data.",
        args.verbose,
    )

    target_taxon = get_gene_variants(rocky.get_rock("rocks_hits_db"))
    gene_dupes = get_prepare_dupes(rocky.get_rock("rocks_nt_db"))

    coords_path = path.join(taxa_path, "coords")
    if path.exists(coords_path):
        rmtree(coords_path)

    makedirs(coords_path, exist_ok=True)
    
    debug_path = path.join(taxa_path, "pairwise_debug")
    if path.exists(debug_path):
        rmtree(debug_path)
    if args.debug:
        makedirs(debug_path, exist_ok=True)

    printv(
        f"Got reference data. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Writing sequences.",
        args.verbose,
    )

    arguments: list[OutputArgs | None] = []
    original_coords = {}
    if is_genome:
        raw_data = rocky.get_rock("rocks_nt_db").get("getall:original_coords")
        
        if raw_data:
            original_coords = json.decode(raw_data, type = dict[str, tuple[str, int, int, int, int]])
        
    for gene, transcript_hits in transcripts_mapped_to:
        if transcript_hits:
            if is_genome:
                for hit in transcript_hits:
                    parent, chomp_start, chomp_end, _, chomp_len = original_coords.get(str(hit.node), (None, None, None, None, None))
                    hit.parent = parent
                    hit.raw_node = hit.node
                    if hit.frame < 0:
                        hit.strand = "-"
                        hit.parent_start = chomp_start
                        hit.chomp_start = (chomp_len - hit.qend) + chomp_start
                        hit.chomp_end = (chomp_len - hit.qstart) + chomp_start - 1
                    else:
                        hit.strand = "+"
                        hit.parent_start = chomp_start
                        hit.chomp_start = hit.qstart + chomp_start
                        hit.chomp_end = hit.qend + chomp_start - 1
                    
            arguments.append(
                (
                    OutputArgs(
                        gene,
                        transcript_hits,
                        aa_out_path,
                        debug_path,
                        taxa_id,
                        nt_out_path,
                        args.verbose,
                        args.compress,
                        set(target_taxon.get(gene, [])),
                        gene_dupes.get(gene, {}),
                        top_refs.get(gene, []),
                        args.minimum_bp,
                        args.debug,
                        is_assembly or is_genome,
                        is_genome
                    ),
                ),
            )
    if num_threads > 1:
        with Pool(num_threads) as pool:
            recovered = pool.starmap(merge_and_write, arguments, chunksize=1)
    else:
        recovered = [merge_and_write(i[0]) for i in arguments]
        
    
    printv(
        f"Done! Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Writing coord data.",
        args.verbose,
    )    
        
    final_count = 0
    removed_total = 0
    intron_removal_log = []
    this_gene_based_scores = {}
    this_gene_based_intron_removals = {}
    global_out = []
    parent_gff_output = defaultdict(list)
    end_bp = {}
    gff_output = ["##gff-version\t3"]
    global_merge_log = []
    global_repeats_log = []
    
    for gene, remove_introns, amount, scores, nodes, gff, merge_log, gene_intron_coordinates, gene_repeats in recovered:
        this_gene_based_intron_removals[gene] = gene_intron_coordinates
        intron_removal_log.extend(remove_introns)
        removed_total += len(remove_introns)
        global_merge_log.extend(merge_log)
        global_repeats_log.extend(gene_repeats)
        out_data = defaultdict(list)
        if original_coords:
            for node in nodes:
                tup = original_coords.get(str(node), None)
                if tup:
                    parent, chomp_start, chomp_end, input_len, chomp_len = tup
                    if parent not in end_bp:
                        end_bp[parent] = input_len
                    out_data[parent].append((node, chomp_start, chomp_end))
                
            for parent, node, coords, strand, frame in gff:
                if len(coords) > 1:
                    for i, (act_start, act_end) in enumerate(coords):
                        let = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")[i]
                        parent_gff_output[parent].append(((act_start), f"{parent}\tSapphyre\texon\t{act_start}\t{act_end}\t.\t{strand}\t.\tName={node} Exon {let};Parent={node};Note={frame};"))
                else:
                    act_start, act_end = coords[0]
                    parent_gff_output[parent].append(((act_start), f"{parent}\tSapphyre\texon\t{act_start}\t{act_end}\t.\t{strand}\t.\tParent={node};Note={frame};"))
                    
        
            if out_data:
                global_out.append(f"### {gene} ###")
                for og, data in out_data.items():
                    data.sort(key = lambda x: (x[1]))
                    global_out.append(f">{og}\n")
                    for node, start, end in data:
                        global_out.append(f"{node}\t{start}-{end}\n")
        
        final_count += amount
        this_gene_based_scores[gene] = scores
        
    with open(path.join(coords_path,"Diamond_merges.txt"), "w") as fp:
        fp.write("\n".join(global_merge_log))
    
    if args.debug:
        with open(path.join(taxa_path,"Intron_removal.txt"), "w") as fp:
            fp.write("\n".join(intron_removal_log))
        with open(path.join(taxa_path,"Repeats_removed.txt"), "w") as fp:
            fp.write("\n".join(global_repeats_log))
        
    if is_genome:
        
        for parent, rows in parent_gff_output.items():
            end = end_bp[parent]
            gff_output.append(f"##sequence-region\t{parent}\t{1}\t{end}")
            rows.sort(key = lambda x: (x[0]))
            gff_output.extend(i[1] for i in rows)
        
        if gff_output:
            with open(path.join(coords_path, "coords.gff"), "w") as fp:
                fp.write("\n".join(gff_output))
        if args.debug:
            with open(path.join(taxa_path, "coords.txt"), "w") as fp:
                fp.write("\n".join(global_out))

    key = "getall:hmm_gene_scores"
    data = json.encode(this_gene_based_scores)
    rocky.get_rock("rocks_nt_db").put_bytes(key, data)
    
    key = "getall:intron_removals"
    data = json.encode(this_gene_based_intron_removals)
    rocky.get_rock("rocks_nt_db").put_bytes(key, data)

    printv(
        f"Done! Took {time_keeper.differential():.2f}s overall. Coords took {time_keeper.lap():.2f}s. Found {final_count} sequences. Removed {removed_total} introns.",
        args.verbose,
    )

    return True


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exist.", args.verbose, 0)
        return False
    rocky.create_pointer(
        "rocks_orthoset_db",
        path.join(args.orthoset_input, args.orthoset, "rocksdb"),
    )
    result = []
    for input_path in args.INPUT:
        rocks_db_path = path.join(input_path, "rocksdb")
        rocky.create_pointer(
            "rocks_nt_db",
            path.join(rocks_db_path, "sequences", "nt"),
        )

        rocky.create_pointer("rocks_hits_db", path.join(rocks_db_path, "hits"))
        result.append(
            do_taxa(
                input_path,
                path.basename(input_path).split(".")[0],
                args,
            ),
        )
        rocky.close_pointer("rocks_nt_db")
        rocky.close_pointer("rocks_hits_db")
    rocky.close_pointer("rocks_orthoset_db")
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return all(result)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Reporter"
    raise Exception(
        msg,
    )
