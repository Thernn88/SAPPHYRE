from argparse import Namespace
from collections import Counter, defaultdict
from decimal import Decimal
from itertools import combinations, count
from math import ceil
from multiprocessing.pool import Pool
from os import makedirs, path, remove, stat, system
from shutil import rmtree
from tempfile import NamedTemporaryFile, TemporaryDirectory
from time import time
from typing import Union

from msgspec import Struct, json
from numpy import float32, float64, int8, uint16, where
from pandas import DataFrame, read_csv
from sapphyre_tools import bio_revcomp, get_overlap
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta


class ReferenceHit(Struct, frozen=True):
    query: str
    ref: Union[str, None]
    start: int
    end: int


class Hit(Struct, frozen=True):
    node: str
    target: str
    frame: int
    evalue: float
    score: float
    qstart: int
    qend: int
    sstart: int
    send: int
    gene: str
    uid: int
    ref: str
    refs: list[ReferenceHit]


class ReporterHit(Struct):
    node: str
    frame: int
    qstart: int
    qend: int
    gene: str
    query: str
    uid: int
    refs: list[ReferenceHit]
    seq: str = None


class ProcessingArgs(Struct, frozen=True):
    i: int
    grouped_data: DataFrame
    target_to_taxon: dict[str, tuple[str, str, int]]
    debug: bool
    result_fp: str
    pairwise_refs: set[str]
    is_assembly_or_genome: bool
    evalue_threshold: float


class InternalArgs(Struct, frozen=True):
    gene: str
    this_hits: list[Hit]
    debug: bool
    internal_percent: float


class InternalReturn(Struct, frozen=True):
    gene: str
    kick_count: int
    log: list[str]
    this_kicks: set


class MultiArgs(Struct, frozen=True):
    hits: list[Hit]
    debug: bool


class MultiReturn(Struct, frozen=True):
    kick_count: int
    log: list[str]
    this_kicks: set


class ConvertArgs(Struct, frozen=True):
    hits: list[Hit]
    gene: str


class ConvertReturn(Struct, frozen=True):
    gene: str
    hits: list[ReporterHit]


def get_overlap_amount(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Get the overlap between two ranges.

    Args:
    ----
        a_start (int): The starting position of range A.
        a_end (int): The ending position of range A.
        b_start (int): The starting position of range B.
        b_end (int): The ending position of range B.

    Returns:
    -------
        int: The amount of overlap between the two ranges.
    """
    overlap_coords = get_overlap(
        a_start,
        a_end,
        b_start,
        b_end,
        1,
    )
    if overlap_coords == None:
        return 0

    return overlap_coords[1] - overlap_coords[0]


def get_score_difference(score_a: float, score_b: float) -> float:
    """Get the decimal difference between two scores.

    Args:
    ----
        score_a (float): The first score.
        score_b (float): The second score.

    Returns:
    -------
        float: The decimal difference between the two scores.
    """
    # If either score is zero return zero.
    if score_a == 0.0 or score_b == 0.0:
        return 0.0

    # Return the decimal difference between the largest score and the smallest score.
    return max(score_a, score_b) / min(score_a, score_b)


def multi_filter(this_args: MultiArgs) -> tuple[list, int, list]:
    """Filter out hits that imprecisely map to multiple genes.

    Args:
    ----
        this_args (MultiArgs): The arguments for the multi_filter function.

    Returns:
    -------
        tuple[list, int, list]: A tuple containing a list of filtered hits, the number of hits filtered, and a list of debug rows.
    """
    # Global variables for the multi_filter function.
    MULTI_PERCENTAGE_OF_OVERLAP = 0.3
    MULTI_SCORE_DIFFERENCE = 1.1

    log = []

    # Assign the highest scoring hit as the 'master'
    master = this_args.hits[0]

    # Assign all the lower scoring hits as 'candidates'
    candidates = this_args.hits[1:]
    kicks = set()
    passes = len(this_args.hits)

    for i, candidate in enumerate(candidates, 1):
        # Skip if the candidate has been kicked already or if master and candidate are internal hits:
        if not candidate or master.gene == candidate.gene:
            continue

        # Get the distance between the master and candidate
        distance = master.qend - master.qstart

        # Get the amount and percentage overlap between the master and candidate
        amount_of_overlap = get_overlap_amount(
            master.qstart,
            master.qend,
            candidate.qstart,
            candidate.qend,
        )
        percentage_of_overlap = amount_of_overlap / distance

        # If the overlap is greater than 30% and the score difference is greater than 5%
        if percentage_of_overlap >= MULTI_PERCENTAGE_OF_OVERLAP:
            score_difference = get_score_difference(master.score, candidate.score)
            if score_difference >= MULTI_SCORE_DIFFERENCE:
                # Kick the candidate
                passes -= 1
                kicks.add(this_args.hits[i].uid)
                kicks.add(candidates[i - 1].uid)
                if this_args.debug:
                    log.append(
                        (
                            candidate.gene,
                            candidate.node,
                            candidate.ref,
                            candidate.score,
                            candidate.qstart,
                            candidate.qend,
                            "Kicked out by",
                            master.gene,
                            master.node,
                            master.ref,
                            master.score,
                            master.qstart,
                            master.qend,
                        ),
                    )
            else:
                # If the score difference is not greater than 5% trigger miniscule score.
                # Miniscule score means hits map to loosely to multiple genes thus
                # we can't determine a viable hit on a single gene and must kick all.
                if this_args.debug:
                    log.extend(
                        [
                            (
                                hit.gene,
                                hit.node,
                                hit.ref,
                                hit.score,
                                hit.qstart,
                                hit.qend,
                                "Kicked due to miniscule score",
                                master.gene,
                                master.node,
                                master.ref,
                                master.score,
                                master.qstart,
                                master.qend,
                            )
                            for hit in this_args.hits
                            if hit
                        ],
                    )

                kicks = set(hit.uid for hit in this_args.hits)

                return MultiReturn(len(this_args.hits), log, kicks)

    return MultiReturn(len(this_args.hits) - passes, log, kicks)


def internal_filter(this_args: InternalArgs) -> tuple[set, list, int]:
    """Filters out overlapping hits who map to the same gene.

    Args:
    ----
        this_args (InternalArgs): The arguments for the internal filter.

    Returns:
    -------
        tuple[set, list, int]:
            A tuple containing the set of kicked hit headers,
            the log of kicked hits, and the number of hits kicked.
    """
    log = []
    kicks = count()
    this_kicks = set()

    for header_hits in this_args.this_hits:
        for hit_a, hit_b in combinations(header_hits, 2):
            if hit_a.uid in this_kicks or hit_b.uid in this_kicks:
                continue

            # Calculate overlap percent over the internal overlap's length
            overlap_amount = get_overlap_amount(
                hit_a.qstart,
                hit_a.qend,
                hit_b.qstart,
                hit_b.qend,
            )
            internal_start = min(hit_a.qstart, hit_b.qstart)
            internal_end = max(hit_a.qend, hit_b.qend)
            internal_length = internal_end - internal_start
            percent = overlap_amount / internal_length

            # If bp overlap and the percent of overlap is greater than or equal to the internal percent arg
            if overlap_amount > 0 and percent >= this_args.internal_percent:
                # Iterate the kicks counter and set the kick attribute to True
                next(kicks)

                # Add the hit to the kicked set
                this_kicks.add(hit_b.uid)

                if this_args.debug:
                    log.append(
                        "".join(
                            [
                                f"{this_args.gene}, {hit_b.node}, {round(hit_b.score, 2)}, {hit_b.qstart}, ",
                                f"{hit_b.qend} Internal kicked out by {this_args.gene}, {hit_a.node}, ",
                                f"{round(hit_a.score, 2)}, {hit_a.qstart}, {hit_a.qend}, {round(percent, 3)}",
                            ],
                        )
                    )

    return InternalReturn(this_args.gene, next(kicks), log, this_kicks)


def convert_and_cull(this_args: ConvertArgs) -> ConvertReturn:
    """Converts a list of the larger Hit class to a smaller ReporterHit class and culls
    left over data

    Args:
    ----
        this_args (ConvertArgs): The arguments for the convert_and_cull function.
    Returns:
    -------
        ConvertReturn: A struct containing the gene and the list of ReporterHits.
    """
    output = []
    for hit in this_args.hits:
        output.append(
            ReporterHit(
                hit.node,
                hit.frame,
                hit.qstart,
                hit.qend,
                hit.gene,
                hit.ref,
                hit.uid,
                hit.refs,
            ),
        )
    return ConvertReturn(this_args.gene, output)


def process_lines(pargs: ProcessingArgs) -> tuple[dict[str, Hit], int, list[str]]:
    """Process a subset of lines from the Diamond tsv result.

    Args:
    ----
        pargs (ProcessingArgs): The arguments for processing the lines.

    Returns:
    -------
        tuple[dict[str, Hit], int, list[str]]:
            A tuple containing the dictionary of hits grouped by header,
            the number of hits kicked, and the list of kicked hits.
    """
    output = defaultdict(list)
    multi_kicks = 0
    this_log = []

    # Grab groups of data based on base header
    for _, header_df in pargs.grouped_data.groupby("header"):
        frame_to_hits = defaultdict(list)
        target_to_hits = defaultdict(list)
        pool_scores = {}

        hits = []
        has_gene_key = False
        for i, row in enumerate(header_df.values):
            # These values are referenced multiple times
            key = row[1]
            sstart = row[7]
            send = row[8]
            frame = row[2]
            qstart = row[5]
            qend = row[6]
            evalue = row[3]

            if i == 0:
                has_gene_key = "|" in key

            # Skip if evalue is greater than threshold
            if evalue > pargs.evalue_threshold:
                continue

            # Negative frames coords are in reverse order.
            if frame < 0:
                qend, qstart = qstart, qend

            # Get the gene and taxon from the target orthoset data
            if has_gene_key:
                target = key.split("|",1)[1]
            else:
                target = key

            gene, ref, _ = pargs.target_to_taxon[key]

            # Create a list of ReferenceHit objects starting with the
            # reference hit for the current row
            refs = [ReferenceHit(target, ref, sstart, send)]

            this_hit = Hit(
                row[0],
                target,
                row[2],
                row[3],
                row[4],
                qstart,
                qend,
                sstart,
                send,
                gene,
                abs(hash(time())) // (pargs.i + 1),
                ref,
                refs,
            )

            # Used to calculate the sum of each targets' individual scores
            if target not in pool_scores:
                pool_scores[target] = this_hit.score
            else:
                pool_scores[target] += this_hit.score

            # Add the hit to a dict based on target
            target_to_hits[target].append(this_hit)

        if pool_scores:
            top_score = max(pool_scores.values())
            top_score *= 0.9

            # Filter hits based on target where the sum of scores is not greater than 90% of
            # the highest targets' sum.
            for target, hits in target_to_hits.items():
                if pool_scores[target] >= top_score:
                    for this_hit in hits:
                        frame_to_hits[this_hit.frame].append(this_hit)

        for hits in frame_to_hits.values():
            hits.sort(key=lambda x: x.score, reverse=True)
            genes_present = {hit.gene for hit in hits}

            # Run multi filter if hits map to multiple genes
            kicks = set()
            if len(genes_present) > 1:
                result = multi_filter(MultiArgs(hits, pargs.debug))
                kicks = result.this_kicks
                multi_kicks += result.kick_count
                if pargs.debug:
                    this_log.extend(result.log)

            if any(hits):
                # Delegate hits into a gene based dict. If multi was ran then apply kicks
                if len(genes_present) > 1:
                    gene_hits = defaultdict(list)
                    for hit in hits:
                        if hit.uid not in kicks:
                            gene_hits[hit.gene].append(hit)
                else:
                    gene_hits = {
                        hits[0].gene: [hit for hit in hits if hit.uid not in kicks]
                    }

                # Output the top hit for each gene with the remaining hits as references
                for hits in gene_hits.values():
                    top_hit = hits[0]

                    if pargs.is_assembly_or_genome:
                        ref_seqs = [
                            ReferenceHit(hit.target, hit.ref, hit.sstart, hit.send)
                            for hit in hits[1:]
                        ]
                    else:
                        ref_seqs = [
                            ReferenceHit(hit.target, None, hit.sstart, hit.send)
                            for hit in hits[1:]
                            if hit.ref in pargs.pairwise_refs
                        ]
                    top_hit.refs.extend(ref_seqs)

                    output[top_hit.gene].append(top_hit)

    result = (
        "\n".join([",".join([str(i) for i in line]) for line in this_log]),
        multi_kicks,
        output,
    )

    # Write the results as a compressed msgspec Struct JSON object to the temporary file
    # allocated to this thread
    with open(pargs.result_fp, "wb") as result_file:
        result_file.write(json.encode(result))


def get_head_to_seq(nt_db, recipe):
    """Get a dictionary of headers to sequences.

    Args:
    ----
        nt_db (RocksDB): The NT rocksdb database.
        recipe (list[str]): A list of each batch index.
    Returns:
    -------
        dict[str, str]: A dictionary of headers to sequences.
    """

    head_to_seq = {}
    for i in recipe:
        lines = nt_db.get_bytes(f"ntbatch:{i}").decode().split("\n")
        head_to_seq.update(
            {
                lines[i][1:]: lines[i + 1]
                for i in range(0, len(lines), 2)
                if lines[i] != ""
            },
        )

    return head_to_seq


def top_reference_realign(orthoset_raw_path, orthoset_aln_path, top_refs, target_to_taxon, top_path, gene):
    out = []

    gene_path = path.join(orthoset_aln_path, gene+".aln.fa")
    if not path.exists(gene_path):
        gene_path = path.join(orthoset_raw_path, gene+".fa")
        source = parseFasta(gene_path, True)
    else:
        source = parseFasta(gene_path, True)
        
    for header, seq in source:
        key = f"{gene}|{header}"
        if target_to_taxon.get(header, set()) in top_refs or target_to_taxon.get(key, set()) in top_refs:
            out.append((header, seq.replace("-", "")))        
        
    out_path = path.join(top_path, gene+".aln.fa")

    if len(out) == 1:
        writeFasta(out_path, out)
        return
    
    if path.exists(out_path):
        if len(out) == len(list(parseFasta(out_path, True))):
            return
    
    with NamedTemporaryFile(dir=gettempdir(), prefix=f"{gene}_") as tmp_prealign, NamedTemporaryFile(dir=gettempdir(), prefix=f"{gene}_") as tmp_result:
        tmp_prealign.write("\n".join([f">{i}\n{j}" for i, j in out]).encode())
        tmp_prealign.flush()

        system(
            f"clustalo -i '{tmp_prealign.name}' -o '{tmp_result.name}' --thread=1 --full --force"
        )

        recs = list(parseFasta(tmp_result.name, True))
        
        del_columns = set()
        for i in range(len(recs[0][1])):
            if all(j[1][i] == "-" or j[1][i] == "X" for j in recs):
                del_columns.add(i)

        out = [(header, "".join([let for i, let in enumerate(seq) if i not in del_columns])) for header, seq in recs]

        writeFasta(out_path, out)


def run_process(args: Namespace, input_path: str) -> bool:
    """Run the main process on the input path.

    Args:
    ----
        args (Namespace): The arguments for the program.
        input_path (str): The path to the input directory.

    Returns:
    -------
        bool: Whether the process was successful or not.
    """
    # Due to the thread bottleneck of chunking a ceiling is set on the threads post reporter
    THREAD_CAP = 32
    # Amount of overshoot in estimating end
    OVERSHOOT_AMOUNT = 1.5
    # Minimum headers to try to estimate thread distrubution
    MINIMUM_HEADERS = 32000
    # Minimum amount of hits to delegate to a process
    MINIMUM_CHUNKSIZE = 50
    # Hard coded values for assembly datasets
    ASSEMBLY_EVALUE = 5
    ASSEMBLY_MIN_ORF = 15
    # Default min orf value
    NORMAL_MIN_ORF = 15

    json_encoder = json.Encoder()

    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    sensitivity = args.sensitivity
    top_amount = args.top

    taxa = path.basename(input_path)
    printv(f"Processing: {taxa}", args.verbose, 0)
    printv("Grabbing necessary directories.", args.verbose)
    # make dirs
    diamond_path = path.join(input_path, "diamond")
    if args.overwrite and path.exists(diamond_path):
        rmtree(diamond_path)
    makedirs(diamond_path, exist_ok=True)

    num_threads = args.processes
    post_threads = args.processes if args.processes < THREAD_CAP else THREAD_CAP
    if post_threads > 1:
        pool = Pool(post_threads)
    orthoset_db_path = path.join(orthosets_dir, orthoset, "rocksdb")
    diamond_db_path = path.join(
        orthosets_dir,
        orthoset,
        "diamond",
        orthoset + ".dmnd",
    )
    # Legacy db check
    if not path.exists(orthoset_db_path) or not path.exists(diamond_db_path):
        gen = False
        print("Could not find orthoset directory or diamond database.")
        if path.exists(f"{orthoset}.sqlite"):
            if (
                input(
                    f"Could not find orthoset DB at {orthoset_db_path}.\nDetected legacy orthoset: {orthoset}.sqlite\nWould you like to generate new Orthoset DB? Y/N: ",
                ).lower()
                == "y"
            ):
                print("Attempting to generate DB")
                system(
                    f"python3 -m sapphyre -p {num_threads} Makeref {orthoset}.sqlite -s {orthoset} -all",
                )
                gen = True
        if not gen:
            print("Aborting")
            return False
    printv(
        f"Done! Took: {time_keeper.lap():.2f}s. Grabbing reference data from Orthoset DB",
        args.verbose,
    )
    orthoset_db = RocksDB(orthoset_db_path)

    target_to_taxon_raw = orthoset_db.get_bytes("getall:targetreference")
    if not target_to_taxon_raw:
        print(
            "Diamond DB is missing data. Try regenerating your diamond DB using 'Makeref --diamond'."
        )
        return False

    target_to_taxon = json.decode(
        target_to_taxon_raw,
        type=dict[str, list[str | int]],
    )

    del orthoset_db

    printv(
        f"Got reference data. Took: {time_keeper.lap():.2f}s. Elapsed: {time_keeper.differential():.2f}s. Grabbing NT sequences",
        args.verbose,
    )

    nt_db_path = path.join(input_path, "rocksdb", "sequences", "nt")

    if not path.exists(nt_db_path):
        msg = "Could not find NT database. Either your dataset is corrupt or Prepare has failed to run."
        raise ValueError(msg)

    nt_db = RocksDB(nt_db_path)

    # Check if sequence is assembly
    dbis_assembly = nt_db.get("get:isassembly") == "True"
    dbis_genome = nt_db.get("get:isgenome") == "True"
    evalue = args.evalue
    min_orf = NORMAL_MIN_ORF
    is_assembly_or_genome = dbis_assembly or dbis_genome
    if is_assembly_or_genome:
        evalue = ASSEMBLY_EVALUE
        min_orf = ASSEMBLY_MIN_ORF

    precision = Decimal("5") * Decimal("0.1") ** Decimal(
        evalue,
    )

    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
    else:
        msg = "Nucleotide sequences not found in database. Did Prepare succesfully finish?"
        raise ValueError(
            msg,
        )
    dupe_counts = json.decode(nt_db.get_bytes("getall:dupes"), type=dict[str, int])

    out_path = path.join(diamond_path, f"{sensitivity}")
    extension_found = False
    for extension in [".tsv.gz", ".tsv.tar.gz", ".tsv"]:
        possible = out_path + extension
        if path.exists(possible):
            out_path = possible
            extension_found = True
            break

    printv(
        f"Got NT sequences. Took: {time_keeper.lap():.2f}s. Elapsed: {time_keeper.differential():.2f}s. Running Diamond",
        args.verbose,
    )

    if not path.exists(out_path) or stat(out_path).st_size == 0:
        with TemporaryDirectory(dir=gettempdir()) as dir, NamedTemporaryFile(
            dir=dir,
        ) as input_file:
            for i in recipe:
                input_file.write(nt_db.get_bytes(f"ntbatch:{i}"))
            input_file.flush()

            quiet = "--quiet" if args.verbose == 0 else ""

            if not extension_found:
                out_path += ".tsv"

            time_keeper.lap()  # Reset timer
            system(
                f"diamond blastx -d {diamond_db_path} -q {input_file.name} -o {out_path} --{sensitivity}-sensitive --masking 0 -e {precision} --compress 1 --outfmt 6 qseqid sseqid qframe evalue bitscore qstart qend sstart send {quiet} --top {top_amount} --min-orf {min_orf} --max-hsps 0 -p {num_threads}",
            )
            if not path.exists(path.join(out_path)) and path.exists(path.join(out_path+".gz")):
                out_path += ".gz"
                
            input_file.seek(0)

        printv(
            f"Diamond completed successfully. Took: {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Reading file to memory",
            args.verbose,
        )
    else:
        printv(
            f"Found existing Diamond output. Elapsed time {time_keeper.differential():.2f}s. Reading file to memory",
            args.verbose,
        )

    db = RocksDB(path.join(input_path, "rocksdb", "hits"))
    output = defaultdict(list)
    multi_kicks = 0

    global_log = []
    dupe_divy_headers = defaultdict(set)

    if stat(out_path).st_size == 0:
        printv("Diamond returned zero hits.", args.verbose, 0)
        return True

    df = read_csv(
        out_path,
        engine="pyarrow",
        delimiter="\t",
        header=None,
        names=[
            "header",
            "target",
            "frame",
            "evalue",
            "score",
            "qstart",
            "qend",
            "sstart",
            "send",
        ],
        dtype={
            "header": str,
            "target": str,
            "frame": int8,
            "evalue": float64,
            "score": float32,
            "qstart": uint16,
            "qend": uint16,
            "sstart": uint16,
            "send": uint16,
        },
    )
    printv(
        f"Done! Took: {time_keeper.lap():.2f}s. Elapsed: {time_keeper.differential():.2f}s. Calculating targets.",
        args.verbose,
    )
    target_counts = df["target"].value_counts()
    combined_count = Counter()
    taxon_to_targets = defaultdict(list)
    target_to_gene = {}
    for target, count in target_counts.to_dict().items():
        gene, ref_taxa, _ = target_to_taxon[target]
        taxon_to_targets[ref_taxa].append(target)
        target_to_gene[target] = gene
        combined_count[ref_taxa] += count

    top_refs = set()
    pairwise_refs = set()
    top_targets = set()
    most_common = combined_count.most_common()
    target_count = min(most_common[0:args.top_ref], key=lambda x: x[1])[1]
    #target_count = min_count * (1 - args.top_ref)
    for taxa, count in most_common:
        if count >= target_count:
            top_refs.add(taxa)
            top_targets.update(taxon_to_targets[taxa])

    gene_target_to_taxa = defaultdict(dict)
    for target, (gene, taxa, _) in target_to_taxon.items():
        gene_target_to_taxa[gene][target] = taxa

    target_has_hit = set(df["target"].unique())
    df = df[(df["target"].isin(top_targets))]
    headers = df["header"].unique()
    if len(headers) > 0:
        per_thread = ceil(len(headers) / post_threads)
        if per_thread <= MINIMUM_CHUNKSIZE:
            per_thread = MINIMUM_CHUNKSIZE
            chunks = ceil(len(headers) / per_thread)
        else:
            chunks = post_threads

        arguments = []
        if len(headers) < MINIMUM_HEADERS:
            indices = []
            for x, i in enumerate(range(0, len(headers), per_thread), 1):
                start_index = 0 if x == 1 else end_index + 1

                if x != chunks:
                    last_header = headers[i + per_thread - 1]
                    end_index = (
                        where(
                            df[start_index:]["header"].values == last_header,
                        )[
                            0
                        ][-1]
                        + start_index
                    )

                else:
                    end_index = len(df) - 1

                indices.append((start_index, end_index))
        else:
            estimated_end = ceil((len(df) / chunks) * OVERSHOOT_AMOUNT)

            indices = []
            for x, i in enumerate(range(0, len(headers), per_thread), 1):
                start_index = 0 if x == 1 else end_index + 1

                if x != chunks:
                    last_header = headers[i + per_thread - 1]
                    try:
                        end_index = (
                            where(
                                df[start_index : (start_index + estimated_end)][
                                    "header"
                                ].values
                                == last_header,
                            )[0][-1]
                            + start_index
                        )
                    except IndexError:
                        # Outliers cause our estimation to go whack.
                        end_index = (
                            where(
                                df[start_index:]["header"].values == last_header,
                            )[
                                0
                            ][-1]
                            + start_index
                        )

                else:
                    end_index = len(df) - 1

                indices.append((start_index, end_index))
        temp_files = [NamedTemporaryFile(dir=gettempdir()) for _ in range(chunks)]

        arguments = (
            ProcessingArgs(
                i,
                df.iloc[index[0] : index[1] + 1],
                target_to_taxon,
                args.debug,
                temp_files[i].name,
                pairwise_refs,
                is_assembly_or_genome,
                precision,
            )
            for i, index in enumerate(indices)
        )

        printv(
            f"Got targets. Took: {time_keeper.lap():.2f}s. Elapsed: {time_keeper.differential():.2f}s. Processing data.",
            args.verbose,
        )
        if post_threads > 1:
            pool.map(process_lines, arguments)
        else:
            for arg in arguments:
                process_lines(arg)

        printv(
            f"Data processed. Took: {time_keeper.lap():.2f}s. Elapsed: {time_keeper.differential():.2f}s. Consolidating threads",
            args.verbose,
        )

        del arguments
        decoder = json.Decoder(tuple[str, int, dict[str, list[Hit]]])
        for temp_file in temp_files:
            with open(temp_file.name, "rb") as fp:
                this_log, mkicks, this_output = decoder.decode(fp.read())

            for gene, hits in this_output.items():
                output[gene].extend(hits)

            multi_kicks += mkicks

            if args.debug:
                global_log.append(this_log)
        del (
            this_output,
            mkicks,
            this_log,
        )
        printv(
            f"Done! Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Doing internal filters",
            args.verbose,
        )

        requires_internal = defaultdict(dict)
        internal_order = []
        for gene, hits in output.items():
            this_counter = Counter([i.node for i in hits]).most_common()
            if this_counter[0][1] > 1:
                this_hits = sum(i[1] for i in this_counter if i[1] > 1)
                this_common = {i[0] for i in this_counter if i[1] > 1}
                for hit in [i for i in hits if i.node in this_common]:
                    requires_internal[gene].setdefault(hit.node, []).append(hit)

                internal_order.append((gene, this_hits))

        internal_order.sort(key=lambda x: x[1], reverse=True)

        if post_threads > 1:
            # with Pool(post_threads) as pool:
            internal_results = pool.map(
                internal_filter,
                [
                    InternalArgs(
                        gene,
                        list(requires_internal[gene].values()),
                        args.debug,
                        args.internal_percent,
                    )
                    for gene, _ in internal_order
                ],
            )
        else:
            internal_results = [
                internal_filter(
                    InternalArgs(
                        gene,
                        list(requires_internal[gene].values()),
                        args.debug,
                        args.internal_percent,
                    )
                )
                for gene, _ in internal_order
            ]

        internal_kicks = 0
        for result in internal_results:
            internal_kicks += result.kick_count
            if args.debug:
                global_log.extend(result.log)

            output[gene] = [i for i in output[gene] if i.uid not in result.this_kicks]

        printv(
            f"Filtering done. Took {time_keeper.lap():.2f}s."
            + f" Elapsed time {time_keeper.differential():.2f}s. Doing variant filter",
            args.verbose,
        )
        arguments = []
        present_genes = []
        present_genes = list(output.keys())
        output = output.items()

        # DOING VARIANT FILTER
        variant_filter = defaultdict(list)

        for target, ref_tuple in target_to_taxon.items():
            gene, ref_taxon, data_length = ref_tuple
            if ref_taxon in top_refs:
                variant_filter[gene].append((ref_taxon, target, data_length))

        dict_items = list(variant_filter.items())
        for gene, targets in dict_items:
            target_taxons = [i[0] for i in targets]
            if len(target_taxons) != len(list(set(target_taxons))):
                this_counts = Counter(target_taxons)
                out_targets = [i[1] for i in targets if this_counts[i[0]] == 1]
                for target, count in this_counts.most_common():
                    if count == 1:
                        continue

                    this_targets = [i for i in targets if i[0] == target]
                    variants_with_hits = sum(
                        i[1] in target_has_hit for i in this_targets
                    )
                    all_variants_kicked = variants_with_hits == 0
                    if all_variants_kicked:
                        reintroduce = max(this_targets, key=lambda x: x[2])
                        out_targets.append(reintroduce[1])
                        continue

                    out_targets.extend(
                        [i[1] for i in this_targets if i[1] in target_has_hit],
                    )

                variant_filter[gene] = out_targets
            else:
                variant_filter.pop(gene, -1)

        variant_filter = {k: list(v) for k, v in variant_filter.items()}

        printv(
            f"Done! Took: {time_keeper.lap():.2f}s. Elapsed: {time_keeper.differential():.2f}s. Writing to DB",
            args.verbose,
        )

        db.put_bytes(
            "getall:target_variants",
            json_encoder.encode(variant_filter),
        )

        del variant_filter
        nt_db.put("getall:valid_refs", ",".join(list(top_refs)))

        head_to_seq = get_head_to_seq(nt_db, recipe)

        if is_assembly_or_genome:
            arguments = []
            for gene, hits in output:
                arguments.append(
                    ConvertArgs(
                        hits,
                        gene,
                    ),
                )

            if post_threads > 1:
                # with Pool(post_threads) as pool:
                output = pool.map(convert_and_cull, arguments)
            else:
                output = [convert_and_cull(arg) for arg in arguments]
        passes = 0
        encoder = json.Encoder()
        for result in output:
            if is_assembly_or_genome:
                hits, gene = result.hits, result.gene
            else:
                gene, hits = result

            out = []
            for hit in hits:
                if not is_assembly_or_genome:
                    hit = ReporterHit(
                        hit.node,
                        hit.frame,
                        hit.qstart,
                        hit.qend,
                        hit.gene,
                        hit.ref,
                        hit.uid,
                        hit.refs,
                    )
                hit.seq = head_to_seq[hit.node][hit.qstart - 1 : hit.qend]
                if hit.frame < 0:
                    hit.seq = bio_revcomp(hit.seq)
                out.append(hit)
                dupe_divy_headers[gene].add(hit.node)

            passes += len(out)
            db.put_bytes(f"gethits:{gene}", encoder.encode(out))

        del head_to_seq
        if global_log:
            with open(path.join(input_path, "multi.log"), "w") as fp:
                fp.write("\n".join(global_log))
        printv(
            f"{multi_kicks} multi kicks",
            args.verbose,
        )
        printv(
            f"{internal_kicks} internal kicks",
            args.verbose,
        )
        printv(
            f"Wrote {passes} results after {multi_kicks+internal_kicks} kicks.",
            args.verbose,
        )
        orthoset_raw_path = path.join(orthosets_dir, orthoset, "raw")
        orthoset_aln_path = path.join(orthosets_dir, orthoset, "aln")
        top_path = path.join(input_path, "top")

        printv(
            f"Writing top reference alignments",
            args.verbose,
        )

        if path.exists(top_path):
            rmtree(top_path)
        makedirs(top_path, exist_ok=True)

        arguments = []
        for gene in present_genes:
            arguments.append(
                (orthoset_raw_path, orthoset_aln_path, top_refs, gene_target_to_taxa[gene], top_path, gene)
            )

        if post_threads > 1:
            pool.starmap(top_reference_realign, arguments)
        else:
            for arg in arguments:
                top_reference_realign(*arg)
        if post_threads > 1:
            pool.close()
            pool.terminate()

        gene_dupe_count = defaultdict(dict)
        for gene, headers in dupe_divy_headers.items():
            for base_header in headers:
                if base_header in dupe_counts:
                    gene_dupe_count[gene][base_header] = dupe_counts[base_header]

        db.put("getall:presentgenes", ",".join(present_genes))

        key = "getall:gene_dupes"
        data = json_encoder.encode(gene_dupe_count)  # type=dict[str, dict[str, int]]
        nt_db.put_bytes(key, data)

    del top_refs
    del db
    del nt_db

    printv(f"Done! Took {time_keeper.differential():.2f}s overall.", args.verbose)

    return True


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    results = []
    for input_path in args.INPUT:
        results.append(run_process(args, input_path))
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return all(results)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Diamond"
    raise Exception(
        msg,
    )
