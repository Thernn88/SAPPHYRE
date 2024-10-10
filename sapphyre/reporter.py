from __future__ import annotations

from argparse import Namespace
from collections import Counter, defaultdict, namedtuple
from itertools import combinations, product
from multiprocessing.pool import Pool
from os import makedirs, path
from shutil import rmtree

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
from .utils import printv, writeFasta
from Bio.Seq import Seq

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
    base_header: str = None
    aa_sequence: str = None
    parent: str = None
    strand: str = None
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
        if target_taxon:
            if f"{gene}|{taxa_id}" not in target_taxon and taxa_id not in target_taxon: #TODO: Slow to handle both 
                continue
        else:
            if taxon not in top_refs:
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
                if hits[i] is None or hits[j] is None:
                    continue
                
                if hits[i].strand != hits[j].strand:
                    continue
                
                if get_overlap(hits[i].chomp_start, hits[i].chomp_end, hits[j].chomp_start, hits[j].chomp_end, minimum_bp_overlap) is None:
                    continue
                
                if not any(b - a <= 1 for a, b in product([hits[i].raw_node] + hits[i].raw_children, [hits[j].raw_node] + hits[j].raw_children)):
                    continue
                
                if abs(hits[i].chomp_start - hits[j].chomp_start) % 3 != 0:
                    continue # different frame same seq
                
                if hits[i].strand == "-":
                    gaps_a_start = max(hits[j].chomp_end - hits[i].chomp_end, 0)
                    gaps_b_start = max(hits[i].chomp_end - hits[j].chomp_end, 0)
                else:
                    gaps_a_start = max(hits[i].chomp_start - hits[j].chomp_start, 0)
                    gaps_b_start = max(hits[j].chomp_start - hits[i].chomp_start, 0)

                a_align_seq = ("-" * gaps_a_start) + hits[i].seq
                b_align_seq = "-" * gaps_b_start + hits[j].seq

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

                if hits[i].raw_node == hits[j].raw_node:
                    log.append("WARNING: {} and {} same base node merge".format(hits[i].node, hits[j].node))
                
                hits[i].children.append(hits[j].node)
                hits[i].raw_children.append(hits[j].raw_node)
                
                hits[i].children.extend(hits[j].children)
                hits[i].raw_children.extend(hits[j].raw_children)
                
                hits[i].seq = merged_seq
                
                # Update coords
                hits[i].chomp_start = min(hits[i].chomp_start, hits[j].chomp_start)
                hits[i].chomp_end = max(hits[i].chomp_end, hits[j].chomp_end)
                
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


def do_dupe_check(hits, header_template, is_assembly_or_genome, dupe_debug_fp, taxa_id):
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
                if dupe_debug_fp:
                    dupe_debug_fp.write(
                        f"{hit.header}\n{hit.seq}\nis an nt dupe of\n{mapped_to}\n\n",
                    )
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
                    already_mapped_header = already_mapped_hit.header
                    already_mapped_sequence = already_mapped_hit.aa_sequence
                    if len(hit.aa_sequence) > len(already_mapped_sequence):
                        if already_mapped_sequence in hit.aa_sequence:
                            hits[already_mapped_index] = None
                            continue
                    else:
                        if hit.aa_sequence in already_mapped_sequence:
                            if dupe_debug_fp:
                                dupe_debug_fp.write(
                                    f"{hit.header}\n{hit.aa_sequence}\nis an aa dupe of\n{already_mapped_header}\n\n",
                                )
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
    debug_dupes = None
    if oargs.debug:
        makedirs(f"align_debug/{oargs.gene}", exist_ok=True)
        debug_dupes = open(f"align_debug/{oargs.gene}/{oargs.taxa_id}.dupes", "w")

    header_template = "{}|{}|{}|NODE_{}|{}"
        
    translate_sequences(this_hits)
    
    this_hits = check_minimum_bp(this_hits, oargs.minimum_bp)
        
    this_hits, this_gene_dupes = do_dupe_check(this_hits, header_template, oargs.is_assembly_or_genome, debug_dupes, oargs.taxa_id)
        
    merge_log = []
    if oargs.is_genome or False: # Set False to disable
        this_hits, merge_log = merge_hits(this_hits)
        for hit in this_hits:
            hit.header = header_template.format(hit.gene, hit.query, oargs.taxa_id, hit.get_merge_header(), hit.frame)

        # Refresh translation
        translate_sequences(this_hits)
        
    # Trim and save the sequences
    aa_output, nt_output, header_to_score = print_unmerged_sequences(
        this_hits,
        oargs.is_assembly_or_genome,
    )
    if debug_dupes:
        debug_dupes.close()

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
    return oargs.gene, len(aa_output), header_to_score, gene_nodes, [(hit.parent, hit.get_merge_header(), hit.chomp_start, hit.chomp_end, hit.strand, hit.frame) for hit in this_hits], merge_log


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
                        hit.chomp_start = (chomp_len - hit.qend) + chomp_start
                        hit.chomp_end = (chomp_len - hit.qstart) + chomp_start - 1
                    else:
                        hit.strand = "+"
                        hit.chomp_start = hit.qstart + chomp_start
                        hit.chomp_end = hit.qend + chomp_start - 1
                    
            arguments.append(
                (
                    OutputArgs(
                        gene,
                        transcript_hits,
                        aa_out_path,
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
    if args.debug:
        makedirs("align_debug", exist_ok=True)

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
    this_gene_based_scores = {}
    global_out = []
    parent_gff_output = defaultdict(list)
    end_bp = {}
    gff_output = ["##gff-version\t3"]
    global_merge_log = []
    
    for gene, amount, scores, nodes, gff, merge_log in recovered:
        global_merge_log.extend(merge_log)
        out_data = defaultdict(list)
        if original_coords:
            for node in nodes:
                tup = original_coords.get(str(node), None)
                if tup:
                    parent, chomp_start, chomp_end, input_len, chomp_len = tup
                    if parent not in end_bp:
                        end_bp[parent] = input_len
                    out_data[parent].append((node, chomp_start, chomp_end))
                
            for parent, node, act_start, act_end, strand, frame in gff:
                parent_gff_output[parent].append(((act_start), f"{parent}\tSapphyre\texon\t{act_start}\t{act_end}\t.\t{strand}\t.\tID={gene};Name={gene};Description={node};Note={frame};"))
                    
        
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

    printv(
        f"Done! Took {time_keeper.differential():.2f}s overall. Coords took {time_keeper.lap():.2f}s. Found {final_count} sequences.",
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
