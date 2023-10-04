from __future__ import annotations

from os import makedirs, path, system
from shutil import rmtree
from argparse import Namespace
from collections import Counter, defaultdict, namedtuple
from multiprocessing.pool import Pool
from tempfile import NamedTemporaryFile, TemporaryDirectory
from xxhash import xxh3_64
from Bio.Seq import Seq
from msgspec import json
from wrap_rocks import RocksDB

from . import rocky
from .diamond import ReporterHit
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
        "clear_output",
    ],
)

class Hit(ReporterHit):
    aa_seq: str = None
    aa_start: int = None
    aa_end: int = None

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
    present_genes = rocks_hits_db.get("getall:presentgenes").split(",")
    genes_to_process = list_of_wanted_genes or present_genes

    gene_based_results = []
    for gene in genes_to_process:
        gene_based_results.append((gene,rocks_hits_db.get_bytes(f"gethits:{gene}")))

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
    return rocks_nt_db.get("getall:valid_refs").split(","), rocks_nt_db.get("get:isassembly") == "True"


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
            if taxa_id not in target_taxon:
                continue
        else:
            if taxon not in top_refs:
                continue

        header = (
            gene
            + "|"
            + taxon
            + "|"
            + taxa_id
            + "|."
        )
        result.append((header, seq))

    return result


def print_unmerged_sequences(
    hits: list,
    gene: str,
    taxa_id: str,
    minimum_bp: int,
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
    aa_result = []
    nt_result = []
    header_maps_to_where = {}
    header_mapped_x_times = Counter()
    base_header_mapped_already = {}
    seq_mapped_already = {}
    exact_hit_mapped_already = set()
    dupes = defaultdict(list)
    for hit in hits:
        base_header = hit.node
        reference_frame = str(hit.frame)

        # Format header to gene|taxa_name|taxa_id|sequence_id|frame
        header = (
            gene
            + "|"
            + hit.query
            + "|"
            + taxa_id
            + "|"
            + base_header
            + "|"
            + reference_frame
        )

        nt_seq = hit.seq
        aa_seq = hit.aa_seq

        # Check if new seq is over bp minimum
        data_after = len(aa_seq)

        if data_after >= minimum_bp:
            # Hash the NT sequence and the AA sequence + base header
            unique_hit = xxh3_64(base_header + aa_seq).hexdigest()
            nt_seq_hash = xxh3_64(nt_seq).hexdigest()
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
                            + "|"
                            + hit.query
                            + "|"
                            + taxa_id
                            + "|"
                            + base_header
                            + f"_{header_mapped_x_times[old_header]}"
                            + "|"
                            + reference_frame
                        )

                        header_mapped_x_times[base_header] += 1
                else:
                    base_header_mapped_already[base_header] = header, aa_seq

                header_maps_to_where[header] = len(
                    aa_result,
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
        "minimum_bp",
        "debug",
        "is_assembly",
        "tmp_path",
    ],
)


def exonerate_queries(
        hits: list[Hit],
        hit_indices: list[int],
        core_sequence: tuple[str, str],
        tmp_path: str,
    ):
    head_to_coords = defaultdict(list)
    with TemporaryDirectory(dir=tmp_path) as tmp_dir:
        with NamedTemporaryFile(dir=tmp_dir, suffix="_query.fa") as tmp_query, NamedTemporaryFile(dir=tmp_dir, suffix="_target.fa") as tmp_target, NamedTemporaryFile(dir=tmp_dir, suffix="_output.fa") as tmp_output:
            tmp_query.write(f">{core_sequence[0]}\n{core_sequence[1]}\n".encode())
            tmp_query.flush()

            for i in hit_indices:
                hit = hits[i]
                tmp_target.write(f">{hit.node}_{hit.frame}\n{hit.aa_seq}\n".encode())

            tmp_target.flush()

            system(f'exonerate --score 20 --ryo "%ti > %tab %tae\n" --subopt 0 --geneticcode 1 --model affine:local --querytype protein --targettype protein --verbose 0 --showvulgar no --showalignment no --query {tmp_query.name} --target {tmp_target.name} > {tmp_output.name}')
    
            for line in tmp_output.read().decode().split("\n"):
                if line.strip():
                    head_to_coords[line.split(" > ")[0]].append(tuple(map(int, line.split(" > ")[1].split(" "))))

            for i in hit_indices:
                hit = hits[i]
                results = head_to_coords[f"{hit.node}_{hit.frame}"]
                if results:
                    for result in results:
                        if hit.aa_start is None:
                            hit.aa_start, hit.aa_end = result
                            continue

                        result_length = result[1] - result[0]
                        if hit.aa_end - hit.aa_start < result_length:
                            hit.aa_start, hit.aa_end = result

                        
def exonerate_and_trim(oargs: OutputArgs) -> tuple[str, dict, int]:
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

    core_sequences, core_sequences_nt = get_core_sequences(
        oargs.gene,
        rocky.get_rock("rocks_orthoset_db"),
    )

    core_seq_dict = defaultdict(list)

    for _, taxa_id, seq in core_sequences:
        core_seq_dict[taxa_id] = (taxa_id, seq)

    this_aa_path = path.join(oargs.aa_out_path, oargs.gene + ".aa.fa")
    this_hits = json.decode(oargs.list_of_hits, type = list[Hit])

    # translate sequences
    for hit in this_hits:
        hit.aa_seq = translate_cdna(hit.seq)

    candidate_seq_dict = defaultdict(list)
    for i, hit in enumerate(this_hits):
        for ref in hit.refs:
            candidate_seq_dict[ref.query].append(i)

    for query, hit_indices in candidate_seq_dict.items():
        exonerate_queries(this_hits, hit_indices, core_seq_dict.get(query), oargs.tmp_path)

    this_hits = [hit for hit in this_hits if hit.aa_start is not None and hit.aa_end is not None]

    for hit in this_hits:
        hit.aa_seq = hit.aa_seq[hit.aa_start:hit.aa_end]
        hit.seq = hit.seq[hit.aa_start*3:hit.aa_end*3]

    this_gene_dupes, aa_output, nt_output = print_unmerged_sequences(
        this_hits,
        oargs.gene,
        oargs.taxa_id,
        oargs.minimum_bp,
    )

    if aa_output:
        aa_core_sequences = print_core_sequences(
            oargs.gene,
            core_sequences,
            oargs.target_taxon,
            oargs.top_refs,
        )
        writeFasta(this_aa_path, aa_core_sequences + aa_output, oargs.compress)

        this_nt_path = path.join(oargs.nt_out_path, oargs.gene + ".nt.fa")

        nt_core_sequences = print_core_sequences(
            oargs.gene,
            core_sequences_nt,
            oargs.target_taxon,
            oargs.top_refs,
        )
        writeFasta(this_nt_path, nt_core_sequences + nt_output, oargs.compress)

    printv(
        f"{oargs.gene} took {t_gene_start.differential():.2f}s. Had {len(aa_output)} sequences",
        oargs.verbose,
        2,
    )
    return oargs.gene, this_gene_dupes, len(aa_output)


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
    if not isinstance(num_threads, int) or num_threads < 1:
        num_threads = 1

    if path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif path.exists("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = taxa_path.join(taxa_path, "tmp")
        makedirs(tmp_path, exist_ok=True)

    if args.gene_list_file:
        with open(args.gene_list_file) as fp:
            list_of_wanted_genes = fp.read().split("\n")
    else:
        list_of_wanted_genes = []

    aa_out = "aa"
    nt_out = "nt"

    aa_out_path = path.join(taxa_path, aa_out)
    nt_out_path = path.join(taxa_path, nt_out)

    if args.clear_output:
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
    top_refs, is_assembly = get_toprefs(rocky.get_rock("rocks_nt_db"))

    printv(
        f"Got reference data. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Trimming hits to alignment coords.",
        args.verbose,
    )

    arguments: list[OutputArgs | None] = []
    for gene, transcript_hits in transcripts_mapped_to:
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
                    top_refs,
                    args.minimum_bp,
                    args.debug,
                    is_assembly,
                    tmp_path,
                ),
            ),
        )
    if args.debug:
        makedirs("align_debug", exist_ok=True)
    # this sorting the list so that the ones with the most hits are first
    if num_threads > 1:
        with Pool(num_threads) as pool:
            recovered = pool.starmap(exonerate_and_trim, arguments, chunksize=1)

    else:
        recovered = [exonerate_and_trim(i[0]) for i in arguments]

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
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exist.", args.verbose, 0)
        return False
    rocky.create_pointer(
        "rocks_orthoset_db",
        path.join(args.orthoset_input, args.orthoset, "rocksdb"),
    )
    result = []
    EXACT_MATCH_AMOUNT = 4
    if args.matches < EXACT_MATCH_AMOUNT:
        printv(
            f"ERROR: Impossible match paramaters. {EXACT_MATCH_AMOUNT} exact matches required whereas only {args.matches} matches are checked.",
            args.verbose,
            0,
        )
        printv(
            f"Please increase the number of matches to at least {EXACT_MATCH_AMOUNT} or change the minimum amount of exact matches.\n",
            args.verbose,
            0,
        )
        return False
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
