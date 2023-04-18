from __future__ import annotations
from argparse import Namespace

import os
import shutil
from collections import Counter, defaultdict, namedtuple
from multiprocessing.pool import Pool
from typing import Optional, TextIO, Union
from msgspec import json, Struct
import parasail as ps
import blosum as bl
import phymmr_tools
from Bio.Seq import Seq
import xxhash
from wrap_rocks import RocksDB
from . import rocky
from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv, writeFasta
from .diamond import ReporterHit

MISMATCH_AMOUNT = 1
EXACT_MATCH_AMOUNT = 4
GAP_PENALTY = 2
EXTEND_PENALTY = 1

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
        "blosum_mode",
        "minimum_bp",
        "gene_list_file",
        "clear_output",
    ],
)

# Extend hit with new functions
class Hit(ReporterHit):
    def get_bp_trim(
        self,
        this_aa: str,
        references: dict[str, str],
        matches: int,
        mode: str,
        debug_fp: TextIO,
        header: str,
    ) -> tuple[int, int]:
        """
        Get the bp to trim from each end so that the alignment matches for 'matches' bp.

        BP Trim has 3 different modes. Exact, Strict, and Lax. Exact is the most strict, and
        will only match if the reference and query are identical. Strict will match if the
        reference and query are identical, or if the blosum substitution matrix score is
        greater than 0. Lax will match if the reference and query are identical, or if the
        blosum substitution matrix score is greater than or equal to 0.

        Args:
            this_aa (str): The amino acid sequence to trim.
            references (dict[str, str]): A dictionary of reference sequences.
            matches (int): The number of matches to look for.
            mode (str): The mode to use. Can be "exact", "strict", or "lax".
            debug_fp (TextIO): A file pointer to write debug information to.
            header (str): The header of the sequence.
        Returns:
            tuple[int, int]: The number of bp to trim from the start and end of the sequence.
        """
        # Create the distance function based on the current mode
        if mode == "exact":

            def dist(bp_a, bp_b, _):
                return bp_a == bp_b and bp_a != "-" and bp_b != "-"

        elif mode == "strict":

            def dist(bp_a, bp_b, mat):
                return mat[bp_a][bp_b] > 0.0 and bp_a != "-" and bp_b != "-"

        else:

            def dist(bp_a, bp_b, mat):
                return mat[bp_a][bp_b] >= 0.0 and bp_a != "-" and bp_b != "-"

        mat = bl.BLOSUM(62)
        reg_starts = []
        reg_ends = []

        # Create blosum pairwise aligner profile
        profile = ps.profile_create_16(this_aa, ps.blosum62)

        if debug_fp:
            debug_fp.write(f">{header}\n{this_aa}\n")

        # For each reference sequence
        for number, ref in enumerate(self.ref_hits):
            # Trim to candidate alignment coords
            ref_seq = references[ref.target]
            ref_seq = ref_seq[ref.sstart - 1 : ref.send]

            # Pairwise align the reference and query
            result = ps.nw_trace_scan_profile_16(
                profile, ref_seq, GAP_PENALTY, EXTEND_PENALTY
            )

            # Get the aligned sequences
            this_aa, ref_seq = result.traceback.query, result.traceback.ref
            this_aa_len = len(this_aa)
            if debug_fp:
                debug_fp.write(f">ref_{number}\n{ref_seq}\n")  # DEBUG

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

                    if not dist(ref_seq[i + j], this_aa[i + j], mat):
                        if j == 0 or this_aa[i + j] == "*":
                            this_pass = False
                            break
                        l_mismatch -= 1
                        if l_mismatch < 0:
                            this_pass = False
                            break

                if this_pass and l_exact_matches >= EXACT_MATCH_AMOUNT:
                    reg_starts.append((i - skip_l))
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

                    if not dist(ref_seq[i - j], this_aa[i - j], mat):
                        if j == 0 or this_aa[i - j] == "*":
                            this_pass = False
                            break

                        r_mismatch -= 1
                        if r_mismatch < 0:
                            this_pass = False
                            break

                if this_pass and r_exact_matches >= EXACT_MATCH_AMOUNT:
                    reg_ends.append(this_aa_len - i - (1 + skip_r))
                    break
        if debug_fp:
            debug_fp.write("\n")

        # Return smallest trims
        if reg_starts and reg_ends:
            return min(reg_starts), min(reg_ends)

        # If no trims were found, return None
        return None, None

    def trim_to_coords(self):
        """
        Trims the hit's est_seq to the alignment coords
        """
        self.est_seq = self.est_seq[self.qstart - 1 : self.qend]
        if "revcomp" in self.header:
            self.est_seq = phymmr_tools.bio_revcomp(self.est_seq)


def get_diamondhits(
    rocks_hits_db: RocksDB, list_of_wanted_genes: list
) -> dict[str, list[Hit]]:
    """
    Returns a dictionary of gene to corresponding hits.

    Args:
        rocks_hits_db (RocksDB): RocksDB instance
        list_of_wanted_genes (set): Set of genes to filter by
    Returns:
        dict[str, list[Hit]]: Dictionary of gene to corresponding hits
    """
    present_genes = rocks_hits_db.get("getall:presentgenes").split(",")
    genes_to_process = list_of_wanted_genes or present_genes

    decoder = json.Decoder(list[Hit]) 

    gene_based_results = defaultdict(list)
    for gene in genes_to_process:
        gene_based_results[gene] = decoder.decode(rocks_hits_db.get_bytes(f"gethits:{gene}"))

    return gene_based_results


def get_gene_variants(rocks_hits_db: RocksDB) -> dict[str, list[str]]:
    """
    Grabs the target variants from the hits database.

    Args:
        rocks_hits_db (RocksDB): RocksDB instance
    Returns:
        dict: Dictionary of gene to corresponding target variants
    """
    return json.decode(rocks_hits_db.get("getall:target_variants"), type=dict[str, list[str]])


def get_toprefs(rocks_nt_db: RocksDB) -> list[str]:
    """
    Grabs the top references from the nt database.

    Args:
        rocks_nt_db (RocksDB): RocksDB instance
    Returns:
        list: List of top references
    """
    return rocks_nt_db.get("getall:valid_refs").split(",")


def translate_cdna(cdna_seq):
    """
    Translates Nucleotide sequence to Amino Acid sequence.

    Args:
        cdna_seq (str): Nucleotide sequence
    Returns:
        str: Amino Acid sequence
    """
    if len(cdna_seq) % 3 != 0:
        printv("WARNING: NT Sequence length is not divisable by 3", 0)

    return str(Seq(cdna_seq).translate())


def get_core_sequences(
    gene: str, orthoset_db: RocksDB
) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    """
    Returns the core reference sequences for a given gene.

    Args:
        gene (str): Gene name
        orthoset_db (RocksDB): RocksDB instance
    Returns:
        tuple: Tuple of core AA and NT sequences
    """
    core_seqs = json.decode(orthoset_db.get(f"getcore:{gene}"), type=dict[str, list[tuple[str, str, str]]])
    return core_seqs["aa"], core_seqs["nt"]


def print_core_sequences(
    gene, core_sequences, target_taxon, top_refs, header_seperator="|", identifier="."
):
    """
    Returns a filtered list of headers and sequences for the core sequences.

    Args:
        gene (str): Gene name
        core_sequences (list): List of core sequences
        target_taxon (str): Target taxa

    """
    result = []
    for taxon, taxa_id, seq in sorted(core_sequences):
        # Filter out non-target hits and variants
        if target_taxon:
            if not taxa_id in target_taxon:
                continue
        else:
            if not taxon in top_refs:
                continue

        header = (
            gene
            + header_seperator
            + taxon
            + header_seperator
            + taxa_id
            + header_seperator
            + identifier
        )
        result.append((header, seq))

    return result


def print_unmerged_sequences(
    hits: list,
    gene: str,
    taxa_id: str,
    core_aa_seqs: list,
    trim_matches: int,
    blosum_mode: str,
    minimum_bp: int,
    debug_fp: TextIO,
) -> tuple[dict[str, list], list[tuple[str, str]], list[tuple[str, str]]]:
    """
    Returns a list of unique trimmed sequences for a given gene with formatted headers.

    For every hit in the given hits list the header is formatted, sequence is trimmed and
    translated to AA and the results are deduplicated.

    Args:
        hits (list): List of hits
        gene (str): Gene name
        taxa_id (str): Taxa ID
        core_aa_seqs (list): List of core AA sequences
        trim_matches (int): Number of matches to trim
        blosum_mode (str): Blosum mode
        minimum_bp (int): Minimum number of bp after trim
        debug_fp (TextIO): Debug file pointer
    Returns:
        tuple:
            Tuple containg a dict of removed duplicates,
            a list containing AA and a list containing NT sequences
    """
    aa_result = []
    nt_result = []
    header_maps_to_where = {}
    header_mapped_x_times = Counter()
    base_header_mapped_already = {}
    seq_mapped_already = {}
    exact_hit_mapped_already = set()
    dupes = defaultdict(list)
    header_seperator = "|"

    for hit in hits:
        base_header = hit.header
        reference_frame = str(hit.frame)

        if hit.frame < 0:
            hit.header = (
                hit.header + "|[revcomp]:[translate("+str(abs(hit.frame))+")]"
            )
        else:
            hit.header = hit.header + "|[translate("+str(hit.frame)+")]"

        # Format header to gene|taxa_name|taxa_id|sequence_id|frame
        header = (
            gene
            + header_seperator
            + hit.reftaxon
            + header_seperator
            + taxa_id
            + header_seperator
            + base_header
            + header_seperator
            + reference_frame
        )

        # Trim to alignment coords
        hit.trim_to_coords()

        # Translate to AA
        nt_seq = hit.est_seq
        aa_seq = translate_cdna(nt_seq)

        # Trim to match reference
        r_start, r_end = hit.get_bp_trim(
            aa_seq, core_aa_seqs, trim_matches, blosum_mode, debug_fp, header
        )
        if r_start is None or r_end is None:
            print(f"WARNING: Trim kicked: {hit.header}")
            continue

        if r_end == 0:
            nt_seq = nt_seq[(r_start * 3) :]
            aa_seq = aa_seq[r_start:]
        else:
            nt_seq = nt_seq[(r_start * 3) : -(r_end * 3)]
            aa_seq = aa_seq[r_start:-r_end]

        if debug_fp:
            debug_fp.write(f">{header}\n{aa_seq}\n")

        # Check if new seq is over bp minimum
        data_after = len(aa_seq)

        if data_after >= minimum_bp:
            # Hash the NT sequence and the AA sequence + base header
            unique_hit = xxhash.xxh3_64(base_header + aa_seq).hexdigest()
            nt_seq_hash = xxhash.xxh3_64(nt_seq).hexdigest()
            # Filter and save NT dupes
            if nt_seq_hash in seq_mapped_already:
                mapped_to = seq_mapped_already[nt_seq_hash]
                dupes.setdefault(mapped_to, []).append(base_header)
                continue
            seq_mapped_already[nt_seq_hash] = base_header

            # If the sequence is unique
            if unique_hit not in exact_hit_mapped_already:
                # Remove subsequence dupes from same read
                if base_header in base_header_mapped_already:
                    (
                        already_mapped_header,
                        already_mapped_sequence,
                    ) = base_header_mapped_already[base_header]

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
                            continue

                    if base_header in header_mapped_x_times:
                        # Make header unique
                        old_header = base_header
                        header = (
                            gene
                            + header_seperator
                            + hit.reftaxon
                            + header_seperator
                            + taxa_id
                            + header_seperator
                            + base_header
                            + f"_{header_mapped_x_times[old_header]}"
                            + header_seperator
                            + reference_frame
                        )

                        header_mapped_x_times[base_header] += 1
                else:
                    base_header_mapped_already[base_header] = header, aa_seq

                header_maps_to_where[header] = len(
                    aa_result
                )  # Save the index of the sequence output

                # Write unique sequence
                aa_result.append((header, aa_seq))
                nt_result.append((header, nt_seq))

                header_mapped_x_times.setdefault(base_header, 1)
                exact_hit_mapped_already.add(unique_hit)

    return dupes, aa_result, nt_result


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
        "top_refs",
        "matches",
        "blosum_mode",
        "minimum_bp",
        "debug",
    ],
)


def trim_and_write(oargs: OutputArgs) -> tuple[str, dict, int]:
    """
    Trims, dedupes and writes the output for a given gene.

    Args:
        oargs (OutputArgs): Output arguments
    Returns:
        tuple: 
            Tuple containing the gene name,
            a dict of removed duplicates and the number of sequences written
    """

    t_gene_start = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Doing output for: {oargs.gene}", oargs.verbose, 2)

    core_sequences, core_sequences_nt = get_core_sequences(
        oargs.gene, rocky.get_rock("rocks_orthoset_db")
    )

    core_seq_aa_dict = {target: seq for _, target, seq in core_sequences}
    this_aa_path = os.path.join(oargs.aa_out_path, oargs.gene + ".aa.fa")
    debug_alignments = None
    if oargs.debug:
        os.makedirs(f"align_debug/{oargs.gene}", exist_ok=True)  # DEBUG
        debug_alignments = open(
            f"align_debug/{oargs.gene}/{oargs.taxa_id}.fa", "w"
        )  # DEBUG
        debug_alignments.write(
            f"GAP_PENALTY: {GAP_PENALTY}\nEXTEND_PENALTY: {EXTEND_PENALTY}\n"
        )  # DEBUG

    this_gene_dupes, aa_output, nt_output = print_unmerged_sequences(
        oargs.list_of_hits,
        oargs.gene,
        oargs.taxa_id,
        core_seq_aa_dict,
        oargs.matches,
        oargs.blosum_mode,
        oargs.minimum_bp,
        debug_alignments,
    )
    if debug_alignments:
        debug_alignments.close()
    if aa_output:
        aa_core_sequences = print_core_sequences(
            oargs.gene, core_sequences, oargs.target_taxon, oargs.top_refs
        )
        writeFasta(this_aa_path, aa_core_sequences + aa_output, oargs.compress)

        this_nt_path = os.path.join(oargs.nt_out_path, oargs.gene + ".nt.fa")

        nt_core_sequences = print_core_sequences(
            oargs.gene, core_sequences_nt, oargs.target_taxon, oargs.top_refs
        )
        writeFasta(this_nt_path, nt_core_sequences + nt_output, oargs.compress)

    printv(
        f"{oargs.gene} took {t_gene_start.differential():.2f}s. Had {len(aa_output)} sequences",
        oargs.verbose,
        2,
    )
    return oargs.gene, this_gene_dupes, len(aa_output)


def do_taxa(path: str, taxa_id: str, args: Namespace):
    """
    Main function for processing a given taxa.

    Args:
        path (str): Path to the taxa directory
        taxa_id (str): Taxa ID
        args (Namespace): Reporter arguments
    Returns:
        bool: True if the taxa was processed successfully, False otherwise
    """

    printv(f"Processing: {taxa_id}", args.verbose, 0)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    num_threads = args.processes
    if not isinstance(num_threads, int) or num_threads < 1:
        num_threads = 1

    if os.path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif os.path.exists("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = os.path.join(path, "tmp")
        os.makedirs(tmp_path, exist_ok=True)

    if args.gene_list_file:
        with open(args.gene_list_file) as fp:
            list_of_wanted_genes = fp.read().split("\n")
    else:
        list_of_wanted_genes = []

    aa_out = "aa"
    nt_out = "nt"

    aa_out_path = os.path.join(path, aa_out)
    nt_out_path = os.path.join(path, nt_out)

    if args.clear_output:
        if os.path.exists(aa_out_path):
            shutil.rmtree(aa_out_path)
        if os.path.exists(nt_out_path):
            shutil.rmtree(nt_out_path)

    os.makedirs(aa_out_path, exist_ok=True)
    os.makedirs(nt_out_path, exist_ok=True)

    printv(
        f"Initialized databases. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reciprocal diamond hits.",
        args.verbose,
    )

    transcripts_mapped_to = get_diamondhits(
        rocky.get_rock("rocks_hits_db"), list_of_wanted_genes
    )

    printv(
        f"Got hits. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reference data.",
        args.verbose,
    )

    target_taxon = get_gene_variants(rocky.get_rock("rocks_hits_db"))
    top_refs = get_toprefs(rocky.get_rock("rocks_nt_db"))

    printv(
        f"Got reference data. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Trimming hits to alignment coords.",
        args.verbose,
    )

    arguments: list[Optional[OutputArgs]] = []
    for gene in sorted(
        transcripts_mapped_to, key=lambda k: len(transcripts_mapped_to[k]), reverse=True
    ):
        arguments.append(
            (
                OutputArgs(
                    gene,
                    transcripts_mapped_to[gene],
                    aa_out_path,
                    taxa_id,
                    nt_out_path,
                    args.verbose,
                    args.compress,
                    set(target_taxon.get(gene, [])),
                    top_refs,
                    args.matches,
                    args.blosum_mode,
                    args.minimum_bp,
                    args.debug,
                ),
            )
        )
    if args.debug:
        os.makedirs("align_debug", exist_ok=True)
    # this sorting the list so that the ones with the most hits are first
    if num_threads > 1:
        if num_threads > 1:
            with Pool(num_threads) as pool:
                recovered = pool.starmap(trim_and_write, arguments, chunksize=1)

    else:
        recovered = [trim_and_write(i[0]) for i in arguments]

    final_count = 0
    this_gene_based_dupes = {}

    for gene, dupes, amount in recovered:
        final_count += amount
        this_gene_based_dupes[gene] = dupes

    key = "getall:reporter_dupes"
    data = json.encode(this_gene_based_dupes)
    rocky.get_rock("rocks_nt_db").put_bytes(key, data)

    printv(
        f"Done! Took {time_keeper.differential():.2f}s overall. Trim took {time_keeper.lap():.2f}s and found {final_count} sequences.",
        args.verbose,
    )

    return True


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exist.", args.verbose, 0)
        return False
    rocky.create_pointer(
        "rocks_orthoset_db", os.path.join(args.orthoset_input, args.orthoset, "rocksdb")
    )
    result = []
    if args.matches < EXACT_MATCH_AMOUNT:
        printv(f"ERROR: Impossible match paramaters. {EXACT_MATCH_AMOUNT} exact matches required whereas only {args.matches} matches are checked.", args.verbose, 0)
        printv(f"Please increase the number of matches to at least {EXACT_MATCH_AMOUNT} or change the minimum amount of exact matches.\n", args.verbose, 0)
        return False
    for input_path in args.INPUT:
        rocks_db_path = os.path.join(input_path, "rocksdb")
        rocky.create_pointer(
            "rocks_nt_db", os.path.join(rocks_db_path, "sequences", "nt")
        )
        rocky.create_pointer("rocks_hits_db", os.path.join(rocks_db_path, "hits"))
        result.append(
            do_taxa(
                path=input_path,
                taxa_id=os.path.basename(input_path).split(".")[0],
                args=args,
            )
        )
        rocky.close_pointer("rocks_nt_db")
        rocky.close_pointer("rocks_hits_db")
    rocky.close_pointer("rocks_orthoset_db")
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return all(result)


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nsapphyre Reporter"
    )
