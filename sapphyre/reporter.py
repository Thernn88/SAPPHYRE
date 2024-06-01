from __future__ import annotations

from argparse import Namespace
from collections import Counter, defaultdict, namedtuple
from itertools import combinations, product
from multiprocessing.pool import Pool
from os import makedirs, path
from shutil import rmtree
from typing import TextIO

from blosum import BLOSUM
from msgspec import json
from parasail import blosum62, sw_trace_scan_profile_16, profile_create_16
#from sapphyre_tools import translate
from wrap_rocks import RocksDB
from xxhash import xxh3_64
from Bio.Seq import Seq
from sapphyre_tools import (
    find_index_pair,
    get_overlap
)

from . import rocky
from .hmmsearch import HmmHit
from .timekeeper import KeeperMode, TimeKeeper
from .utils import printv, writeFasta

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
        "matches",
        "blosum_strictness",
        "minimum_bp",
        "gene_list_file",
        "keep_output",
    ],
)


# Extend hit with new functions
class Hit(HmmHit):#, frozen=True):
    parent: str = None
    strand: str = None
    chomp_start: int = None
    chomp_end: int = None
    children: list = []
    
    def get_bp_trim(
        self,
        this_aa: str,
        references: dict[str, str],
        matches: int,
        is_positive_match: callable,
        debug_fp: TextIO,
        header: str,
        mat: dict,
        exact_match_amount: int,
        top_refs: set,
    ) -> tuple[int, int]:
        """Get the bp to trim from each end so that the alignment matches for 'matches' bp.

        BP Trim has 3 different modes. Exact, Strict, and Lax. Exact is the most strict, and
        will only match if the reference and query are identical. Strict will match if the
        reference and query are identical, or if the blosum substitution matrix score is
        greater than 0. Lax will match if the reference and query are identical, or if the
        blosum substitution matrix score is greater than or equal to 0.

        Args:
        ----
            this_aa (str): The amino acid sequence to trim.
            references (dict[str, str]): A dictionary of reference sequences.
            matches (int): The number of matches to look for.
            mode (str): The mode to use. Can be "exact", "strict", or "lax".
            debug_fp (TextIO): A file pointer to write debug information to.
            header (str): The header of the sequence.

        Returns:
        -------
            tuple[int, int]: The number of bp to trim from the start and end of the sequence.
        """
        MISMATCH_AMOUNT = 1
        GAP_PENALTY = 2
        EXTEND_PENALTY = 1

        # Create blosum pairwise aligner profile
        profile = profile_create_16(this_aa, blosum62)

        if debug_fp:
            debug_fp.write(f">{header}\n{this_aa}\n")

        # For each reference sequence
        best_alignment = None
        best_alignment_score = 0

        if not self.refs:
            return 0, 0

        for number, ref in enumerate(self.refs):
            if not ref.ref in top_refs:
                continue
            ref_seq = references[ref.query]
            ref_seq = ref_seq[ref.start - 1 : ref.end]

            # Pairwise align the reference and query
            result = nw_trace_scan_profile_16(
                profile,
                ref_seq,
                GAP_PENALTY,
                EXTEND_PENALTY,
            )

            # Get the aligned sequence with the highest score
            if not hasattr(result, "traceback"):
                continue #alignment failed, to investigate
            this_aa, ref_seq = result.traceback.query, result.traceback.ref
            if result.score > best_alignment_score:
                best_alignment = (this_aa, ref_seq)
                best_alignment_score = result.score

            this_aa_len = len(this_aa)
            if debug_fp:
                debug_fp.write(
                    f">ref_{number} = {result.score}\n{ref_seq}\n>aln\n{this_aa}\n"
                )

        if debug_fp:
            debug_fp.write("\n")

        if best_alignment:
            reg_start = None
            reg_end = None
            # Find which start position matches for 'matches' bases
            skip_l = 0
            for i, let in enumerate(this_aa):
                this_pass = True
                if let == "-":
                    skip_l += 1
                    continue

                l_mismatch = MISMATCH_AMOUNT
                l_exact_matches = 0
                for j in range(0, matches):
                    if (i + j) > (this_aa_len - 1):
                        this_pass = False
                        break
                    if this_aa[i + j] == ref_seq[i + j]:
                        l_exact_matches += 1

                    if not is_positive_match(ref_seq[i + j], this_aa[i + j], mat):
                        if j == 0 or this_aa[i + j] == "*":
                            this_pass = False
                            break
                        l_mismatch -= 1
                        if l_mismatch < 0:
                            this_pass = False
                            break

                if this_pass and l_exact_matches >= exact_match_amount:
                    reg_start = i - skip_l
                    break

            # Find which end position matches for 'matches' bases
            skip_r = 0
            for i in range(this_aa_len - 1, -1, -1):
                this_pass = True
                if this_aa[i] == "-":
                    skip_r += 1
                    continue

                r_mismatch = MISMATCH_AMOUNT
                r_exact_matches = 0
                for j in range(0, matches):
                    if i - j < 0:
                        this_pass = False
                        break

                    if this_aa[i - j] == ref_seq[i - j]:
                        r_exact_matches += 1

                    if not is_positive_match(ref_seq[i - j], this_aa[i - j], mat):
                        if j == 0 or this_aa[i - j] == "*":
                            this_pass = False
                            break

                        r_mismatch -= 1
                        if r_mismatch < 0:
                            this_pass = False
                            break

                if this_pass and r_exact_matches >= exact_match_amount:
                    reg_end = this_aa_len - i - (1 + skip_r)
                    break
            return reg_start, reg_end

        # If no trims were found, return None
        return None, None

    def get_merge_header(self):
        return "&&".join(map(str, [self.node] + self.children))

def get_diamondhits(
    rocks_hits_db: RocksDB,
    list_of_wanted_genes: list,
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
    genes_to_process = list_of_wanted_genes or present_genes.split(",")

    gene_based_results = []
    for gene in genes_to_process:
        gene_result = rocks_hits_db.get_bytes(f"gethmmhits:{gene}")
        if not gene_result:
            printv(
                f"WARNING: No hits found for {gene}. If you are using a gene list file this may be a non-issue",
                0,
            )
            continue
        gene_based_results.append((gene, json.decode(gene_result, type=list[Hit])))

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
    return core_seqs["aa"], core_seqs["nt"]


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
    gene: str,
    taxa_id: str,
    core_aa_seqs: list,
    trim_matches: int,
    is_positive_match: callable,
    minimum_bp: int,
    debug_fp: TextIO,
    dupe_debug_fp: TextIO,
    verbose: int,
    mat: dict,
    is_assembly_or_genome: bool,
    exact_match_amount: int,
    top_refs: set,
) -> tuple[dict[str, list], list[tuple[str, str]], list[tuple[str, str]]]:
    """Returns a list of unique trimmed sequences for a given gene with formatted headers.

    For every hit in the given hits list the header is formatted, sequence is trimmed and
    translated to AA and the results are deduplicated.

    Args:
    ----
        hits (list): List of hits
        gene (str): Gene name
        taxa_id (str): Taxa ID
        core_aa_seqs (list): List of core AA sequences
        trim_matches (int): Number of matches to trim
        blosum_strictness (str): Blosum mode
        minimum_bp (int): Minimum number of bp after trim
        debug_fp (TextIO): Debug file pointer
    Returns:
        tuple:
            Tuple containg a dict of removed duplicates,
            a list containing AA and a list containing NT sequences
    """
    header_template = "{}|{}|{}|NODE_{}|{}"
    base_header_template = "{}_{}"
    aa_result = []
    nt_result = []
    header_maps_to_where = {}
    header_mapped_x_times = Counter()
    base_header_mapped_already = {}
    seq_mapped_already = {}
    exact_hit_mapped_already = set()
    dupes = defaultdict(list)
    header_to_score = {}
    for hit in hits:
        base_header = hit.get_merge_header()
        reference_frame = str(hit.frame)

        # Format header to gene|taxa_name|taxa_id|sequence_id|frame
        header = header_template.format(gene, hit.query, taxa_id, base_header, reference_frame)
        # Translate to AA
        nt_seq = hit.seq
        aa_seq = translate_cdna(nt_seq)

        # Trim to match reference
        r_start = 0
        if not is_assembly_or_genome:
            r_start, r_end = hit.get_bp_trim(
                aa_seq,
                core_aa_seqs,
                trim_matches,
                is_positive_match,
                debug_fp,
                header,
                mat,
                exact_match_amount,
                top_refs,
            )
            if r_start is None or r_end is None:
                printv(f"WARNING: Trim kicked: {hit.node}|{hit.frame}", verbose, 2)
                continue

            if r_end == 0:
                nt_seq = nt_seq[(r_start * 3) :]
                aa_seq = aa_seq[r_start:]
            else:
                nt_seq = nt_seq[(r_start * 3) : -(r_end * 3)]
                aa_seq = aa_seq[r_start:-r_end]

            if debug_fp:
                debug_fp.write(f">{header}\n{aa_seq}\n\n")

        # Check if new seq is over bp minimum
        data_after = len(aa_seq)

        if data_after >= minimum_bp:
            unique_hit = None

            if not is_assembly_or_genome:
                # Hash the NT sequence and the AA sequence + base header
                unique_hit = xxh3_64(base_header + aa_seq).hexdigest()
                nt_seq_hash = xxh3_64(nt_seq).hexdigest()
                # Filter and save NT dupes
                if nt_seq_hash in seq_mapped_already:
                    mapped_to = seq_mapped_already[nt_seq_hash]
                    dupes.setdefault(mapped_to, []).append(base_header)
                    if dupe_debug_fp:
                        dupe_debug_fp.write(
                            f"{header}\n{nt_seq}\nis an nt dupe of\n{mapped_to}\n\n",
                        )
                    continue
                seq_mapped_already[nt_seq_hash] = base_header

            # If the sequence is unique
            if is_assembly_or_genome or unique_hit not in exact_hit_mapped_already:
                # Remove subsequence dupes from same read
                if base_header in base_header_mapped_already:
                    (
                        already_mapped_header,
                        already_mapped_sequence,
                    ) = base_header_mapped_already[base_header]
                    # Dont kick if assembly
                    if not is_assembly_or_genome:
                        if len(aa_seq) > len(already_mapped_sequence):
                            if already_mapped_sequence in aa_seq:
                                aa_result[header_maps_to_where[already_mapped_header]] = (
                                    header,
                                    aa_seq,
                                )
                                nt_result[header_maps_to_where[already_mapped_header]] = (
                                    header,
                                    nt_seq,
                                )
                                continue
                        else:
                            if aa_seq in already_mapped_sequence:
                                if dupe_debug_fp:
                                    dupe_debug_fp.write(
                                        f"{header}\n{aa_seq}\nis an aa dupe of\n{already_mapped_header}\n\n",
                                    )
                                continue

                    if base_header in header_mapped_x_times:
                        # Make header unique
                        current_count = header_mapped_x_times[base_header]
                        header_mapped_x_times[base_header] += 1
                        modified_header = base_header_template.format(base_header, current_count)
                        hit.node = modified_header
                        header = header_template.format(
                            gene,
                            hit.query,
                            taxa_id,
                            modified_header,
                            reference_frame,
                        )

                        header_mapped_x_times[base_header] += 1
                else:
                    base_header_mapped_already[base_header] = header, aa_seq

                # Save the index of the sequence output
                header_maps_to_where[header] = len(
                    aa_result,
                )

                # Write unique sequence
                aa_result.append((header, aa_seq))
                nt_result.append((header, nt_seq))


                header_to_score[hit.node] = max(header_to_score.get(hit.node, 0), hit.score)

                if base_header not in header_mapped_x_times:
                    header_mapped_x_times[base_header] = 1
                exact_hit_mapped_already.add(unique_hit)

    return dupes, aa_result, nt_result, header_to_score, hits


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
        "matches",
        "blosum_strictness",
        "EXACT_MATCH_AMOUNT",
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
        node = header.split("|")[3]
        dupes = prepare_dupes.get(node, 1) + sum(
            prepare_dupes.get(node, 1)
            for node in reporter_dupes.get(node, [])
        )
        header = f"{header}|{dupes}"
        output.append((header, sequence))

    return output


def merge_hits(hits: list[Hit]) -> tuple[list[Hit], list[str]]:
    GAP_PENALTY = 1
    EXTEND_PENALTY = 1
    hits.sort(key = lambda x: (x.node))
    log = ["Merges for: "+hits[0].gene]
    indices = range(len(hits))
    merge_occured = True
    while merge_occured:
        merge_occured = False
        for i, j in combinations(indices, 2):
            if hits[i] is None or hits[j] is None:
                continue
            
            if not any(b - a <= 1 for a, b in product(([hits[i].node] + hits[i].children), ([hits[j].node] + hits[j].children))):
                continue
            
            if get_overlap(hits[i].chomp_start, hits[i].chomp_end, hits[j].chomp_start, hits[j].chomp_end, 1) is None:
                continue
            
            a_seq = translate_cdna(hits[i].seq)
            b_seq = translate_cdna(hits[j].seq)
            
            profile = profile_create_16(a_seq, blosum62)
            result = sw_trace_scan_profile_16(
                    profile,
                    b_seq,
                    GAP_PENALTY,
                    EXTEND_PENALTY,
                )
            
            a_align, b_align = result.traceback.query, result.traceback.ref

            if a_align != b_align:
                continue

            a_coord = result.end_query * 3
            b_coord = result.end_ref * 3
            
            # Same strand
            if (hits[i].frame / abs(hits[i].frame)) != (hits[j].frame / abs(hits[j].frame)):
                continue
            
            # Same frame
            # if hits[i].frame != hits[j].frame:
            #     continue
            
            if hits[i].node == hits[j].node:
                log.append("WARNING: {} and {} same node merge".format(hits[i].node, hits[j].node))
            
            hits[i].children.append(hits[j].node)
            hits[i].seq = hits[i].seq[:a_coord] + hits[j].seq[b_coord:]
            
            hits[i].chomp_end = hits[j].chomp_end
            hits[j] = None
            merge_occured = True
        
    for hit in hits:
        if hit is not None and hit.children:
            log.append(hit.get_merge_header())
    return [i for i in hits if i is not None], log


def trim_and_write(oargs: OutputArgs) -> tuple[str, dict, int]:
    """Trims, dedupes and writes the output for a given gene.

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
    this_hits = oargs.list_of_hits
    gene_nodes = [hit.node for hit in this_hits]
    
    if oargs.is_genome:
        this_hits, merge_log = merge_hits(this_hits)

    # Get reference sequences
    core_sequences, core_sequences_nt = get_core_sequences(
        oargs.gene,
        rocky.get_rock("rocks_orthoset_db"),
    )

    core_seq_aa_dict = {target: seq for _, target, seq in core_sequences}
    this_aa_path = path.join(oargs.aa_out_path, oargs.gene + ".aa.fa")
    debug_alignments = None
    debug_dupes = None
    if oargs.debug:
        makedirs(f"align_debug/{oargs.gene}", exist_ok=True)
        debug_alignments = open(
            f"align_debug/{oargs.gene}/{oargs.taxa_id}.alignments",
            "w",
        )
        debug_alignments.write(
            "GAP_PENALTY: 2\nEXTEND_PENALTY: 1\n",
        )
        debug_dupes = open(f"align_debug/{oargs.gene}/{oargs.taxa_id}.dupes", "w")

    # Initialize blosum matrix and distance function
    mat = BLOSUM(62)
    if oargs.blosum_strictness == "exact":

        def dist(bp_a, bp_b, _):
            return bp_a == bp_b and bp_a != "-" and bp_b != "-"

    elif oargs.blosum_strictness == "strict":

        def dist(bp_a, bp_b, mat):
            return mat[bp_a][bp_b] > 0.0 and bp_a != "-" and bp_b != "-"

    else:

        def dist(bp_a, bp_b, mat):
            return mat[bp_a][bp_b] >= 0.0 and bp_a != "-" and bp_b != "-"

    # Trim and save the sequences
    this_gene_dupes, aa_output, nt_output, header_to_score, this_hits = print_unmerged_sequences(
        this_hits,
        oargs.gene,
        oargs.taxa_id,
        core_seq_aa_dict,
        oargs.matches,
        dist,
        oargs.minimum_bp,
        debug_alignments,
        debug_dupes,
        oargs.verbose,
        mat,
        oargs.is_assembly_or_genome,
        oargs.EXACT_MATCH_AMOUNT,
        oargs.top_refs,
    )
    if debug_alignments:
        debug_alignments.close()
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
    return oargs.gene, this_gene_dupes, len(aa_output), header_to_score, gene_nodes, [(hit.parent, hit.get_merge_header(), hit.chomp_start, hit.chomp_end, hit.strand, hit.frame) for hit in this_hits], merge_log


def get_prepare_dupes(rocks_nt_db: RocksDB) -> dict[str, dict[str, int]]:
    prepare_dupe_counts = json.decode(
        rocks_nt_db.get("getall:gene_dupes"), type=dict[str, dict[str, int]]
    )
    return prepare_dupe_counts


def do_taxa(taxa_path: str, taxa_id: str, args: Namespace, EXACT_MATCH_AMOUNT: int):
    """Main function for processing a given taxa.

    Args:
    ----
        path (str): Path to the taxa directory
        taxa_id (str): Taxa ID
        args (Namespace): Reporter arguments
        EXACT_MATCH_AMOUNT (int): Number of exact matches to look for
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
        f"Initialized databases. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reciprocal diamond hits.",
        args.verbose,
    )

    transcripts_mapped_to = get_diamondhits(
        rocky.get_rock("rocks_hits_db"),
        list_of_wanted_genes,
    )

    printv(
        f"Got hits. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reference data.",
        args.verbose,
    )

    target_taxon = get_gene_variants(rocky.get_rock("rocks_hits_db"))
    top_refs, is_assembly, is_genome = get_toprefs(rocky.get_rock("rocks_nt_db"))
    gene_dupes = get_prepare_dupes(rocky.get_rock("rocks_nt_db"))

    coords_path = path.join(taxa_path, "coords")
    if path.exists(coords_path):
        rmtree(coords_path)

    makedirs(coords_path, exist_ok=True)

    printv(
        f"Got reference data. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Trimming hits to alignment coords.",
        args.verbose,
    )

    arguments: list[OutputArgs | None] = []
    original_coords = {}
    if is_genome:
        raw_data = rocky.get_rock("rocks_nt_db").get("getall:original_coords")
        
        if raw_data:
            original_coords = json.decode(raw_data, type = dict[str, tuple[str, int, int, int, int]])
        
    for gene, transcript_hits in transcripts_mapped_to:
        
        for hit in transcript_hits:
            parent, chomp_start, chomp_end, _, chomp_len = original_coords.get(str(hit.node), (None, None, None, None, None))
            hit.parent = parent
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
                    args.matches,
                    args.blosum_strictness,
                    EXACT_MATCH_AMOUNT,
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
            recovered = pool.starmap(trim_and_write, arguments, chunksize=1)
    else:
        recovered = [trim_and_write(i[0]) for i in arguments]
        
    
    printv(
        f"Done! Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Writing coord data.",
        args.verbose,
    )    
        
    final_count = 0
    this_gene_based_dupes = {}
    this_gene_based_scores = {}
    global_out = []
    parent_gff_output = defaultdict(list)
    end_bp = {}
    gff_output = ["##gff-version\t3"]
    global_merge_log = []
    
    for gene, dupes, amount, scores, nodes, gff, merge_log in recovered:
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
                parent_gff_output[parent].append(((act_start), f"{parent}\tSapphyre\texon\t{act_start}\t{act_end}\t.\t{strand}\t.\tID={node};Parent={gene};Note={frame};"))
                    
        
            if out_data:
                with open(path.join(coords_path,gene+".txt"), "w") as fp:
                    global_out.append(f"### {gene} ###")
                    for og, data in out_data.items():
                        data.sort(key = lambda x: (x[1]))
                        fp.write(f">{og}\n")
                        global_out.append(f">{og}\n")
                        for node, start, end in data:
                            fp.write(f"{node}\t{start}-{end}\n")
                            global_out.append(f"{node}\t{start}-{end}\n")
        
        final_count += amount
        this_gene_based_dupes[gene] = dupes
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


    key = "getall:reporter_dupes"
    data = json.encode(this_gene_based_dupes)
    rocky.get_rock("rocks_nt_db").put_bytes(key, data)

    key = "getall:hmm_gene_scores"
    data = json.encode(this_gene_based_scores)
    rocky.get_rock("rocks_nt_db").put_bytes(key, data)

    printv(
        f"Done! Took {time_keeper.differential():.2f}s overall. Trim took {time_keeper.lap():.2f}s and found {final_count} sequences.",
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
    EXACT_MATCH_AMOUNT = 2
    if args.matches < EXACT_MATCH_AMOUNT:
        printv(
            f"WARNING: Impossible match paramaters. {EXACT_MATCH_AMOUNT} exact matches required whereas only {args.matches} matches are checked.",
            args.verbose,
            0,
        )
        printv(
            f"Setting exact matches to {args.matches}.\n",
            args.verbose,
            0,
        )
        EXACT_MATCH_AMOUNT = args.matches
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
                EXACT_MATCH_AMOUNT,
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
