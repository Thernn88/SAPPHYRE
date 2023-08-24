import decimal
import itertools
import os
from argparse import Namespace
from collections import Counter, defaultdict
from math import ceil
from multiprocessing.pool import Pool
from shutil import rmtree
from tempfile import NamedTemporaryFile, TemporaryDirectory
from time import time
from typing import Union
import phymmr_tools
import numpy as np
import pandas as pd
import wrap_rocks
from msgspec import Struct, json

from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, printv

# call the phymmr_tools.Hit class using an unpacked list of values found in the row input.
# def create_hit(row) -> Hit:


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
    grouped_data: pd.DataFrame
    target_to_taxon: dict[str, tuple[str, str, int]]
    debug: bool
    result_fp: str
    pairwise_refs: set[str]
    is_assembly: bool
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


class ContainmentArgs(Struct, frozen=True):
    hits: list[Hit]
    target_to_taxon: dict[str, tuple[str, str, int]]
    debug: bool
    gene: str


class ContainmentReturn(Struct, frozen=True):
    kicks: set
    kick_count: int
    log: list[str]
    gene: str


class ConvertArgs(Struct, frozen=True):
    hits: list[Hit]
    gene: str


class ConvertReturn(Struct, frozen=True):
    gene: str
    hits: list[ReporterHit]


def get_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Get the overlap between two ranges.

    Args:
    ----
        a_start (int): The starting position of range A.
        a_end (int): The ending position of range A.
        b_start (int): The starting position of range B.
        b_end (int): The ending position of range B.

    Returns:
    -------
        int: The number of elements in the overlap between the two ranges.
    """
    overlap_coords = phymmr_tools.get_overlap(
        a_start,
        a_end,
        b_start,
        b_end,
        0,
    )
    if overlap_coords == None:
        return 0
    
    return (overlap_coords[1] - overlap_coords[0]) + 1


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

    # Assign the highest scoring hit as the master
    master = this_args.hits[0]

    # Assign all the lower scoring hits as candidates
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
        amount_of_overlap = get_overlap(
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
                # If the score difference is not greater than 5% trigger miniscule score
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
    kicks = itertools.count()
    this_kicks = set()

    for header_hits in this_args.this_hits:
        for hit_a, hit_b in itertools.combinations(header_hits, 2):
            if hit_a.uid in this_kicks or hit_b.uid in this_kicks:
                continue

            # Calculate overlap percent over the internal overlap's length
            overlap_amount = get_overlap(
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
                                f"{this_args.gene}, {hit_b.uid}, {round(hit_b.score, 2)}, {hit_b.qstart}, ",
                                f"{hit_b.qend} Internal kicked out by {this_args.gene}, {hit_a.uid}, ",
                                f"{round(hit_a.score, 2)}, {hit_a.qstart}, {hit_a.qend}, {round(percent, 3)}",
                            ],
                        )
                    )

    return InternalReturn(this_args.gene, next(kicks), log, this_kicks)


def containments(this_args: ContainmentArgs):
    PARALOG_SCORE_DIFF = 0.1
    KICK_PERCENT = 0.8
    this_args.hits.sort(key=lambda x: x.score, reverse=True)

    log = []
    kicks = set()

    for i, top_hit in enumerate(this_args.hits):
        if top_hit.uid in kicks:
            continue
        for top_ref_hit in top_hit.refs:
            top_hit_reference = top_ref_hit.ref

            for hit in this_args.hits[i + 1 :]:
                if hit.uid in kicks:
                    continue

                for ref_hit in hit.refs:
                    if this_args.target_to_taxon[ref_hit.query][1] != top_hit_reference:
                        continue

                    if (
                        get_score_difference(top_hit.score, hit.score)
                        < PARALOG_SCORE_DIFF
                    ):
                        continue

                    overlap_percent = get_overlap(
                        top_ref_hit.start,
                        top_ref_hit.end,
                        ref_hit.start,
                        ref_hit.end,
                    ) / min(
                        (top_ref_hit.end - top_ref_hit.start),
                        (ref_hit.end - ref_hit.start),
                    )

                    if overlap_percent > KICK_PERCENT:
                        kicks.add(hit.uid)
                        if this_args.debug:
                            log.append(
                                hit.gene
                                + " "
                                + hit.node
                                + "|"
                                + str(hit.frame)
                                + " kicked out by "
                                + top_hit.node
                                + "|"
                                + str(top_hit.frame)
                                + f" overlap: {overlap_percent:.2f}",
                            )
                        break

    return ContainmentReturn(kicks, len(kicks), log, this_args.gene)


def convert_and_cull(this_args: ConvertArgs) -> ConvertReturn:
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

    for _, header_df in pargs.grouped_data.groupby("header"):
        frame_to_hits = defaultdict(list)
        hits = []
        for row in header_df.values:
            target = row[1]
            sstart = row[7]
            send = row[8]
            frame = row[2]
            qstart = row[5]
            qend = row[6]
            evalue = row[3]

            if evalue > pargs.evalue_threshold:
                continue

            if frame < 0:
                qend, qstart = qstart, qend
            gene, ref, _ = pargs.target_to_taxon[target]

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
                hash(time()),
                ref,
                refs,
            )

            frame_to_hits[this_hit.frame].append(this_hit)

        for hits in frame_to_hits.values():
            hits.sort(key=lambda x: x.score, reverse=True)
            genes_present = {hit.gene for hit in hits}

            kicks = set()
            if len(genes_present) > 1:
                result = multi_filter(MultiArgs(hits, pargs.debug))
                kicks = result.this_kicks
                multi_kicks += result.kick_count
                if pargs.debug:
                    this_log.extend(result.log)

            if any(hits):
                if len(genes_present) > 1:
                    gene_hits = defaultdict(list)
                    for hit in hits:
                        if hit.uid not in kicks:
                            gene_hits[hit.gene].append(hit)
                else:
                    gene_hits = {
                        hits[0].gene: [hit for hit in hits if hit.uid not in kicks]
                    }

                for _, hits in gene_hits.items():
                    top_hit = hits[0]

                    if pargs.is_assembly:
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

    with open(pargs.result_fp, "wb") as result_file:
        result_file.write(json.encode(result))


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
    OVERSHOOT_AMOUNT = 1.75
    # Minimum headers to try to estimate thread distrubution
    MINIMUM_HEADERS = 32000
    # Minimum amount of hits to delegate to a process
    MINIMUM_CHUNKSIZE = 50
    # Hard coded values for assembly datasets
    ASSEMBLY_EVALUE = 5
    ASSEMBLY_MIN_ORF = 30
    # Default min orf value
    NORMAL_MIN_ORF = 20

    json_encoder = json.Encoder()

    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    sensitivity = args.sensitivity
    top_amount = args.top

    taxa = os.path.basename(input_path)
    printv(f"Processing: {taxa}", args.verbose, 0)
    printv("Grabbing reference data from Orthoset DB.", args.verbose)
    # make dirs
    diamond_path = os.path.join(input_path, "diamond")
    if args.overwrite and os.path.exists(diamond_path):
        rmtree(diamond_path)
    os.makedirs(diamond_path, exist_ok=True)

    num_threads = args.processes
    post_threads = args.processes if args.processes < THREAD_CAP else THREAD_CAP

    orthoset_db_path = os.path.join(orthosets_dir, orthoset, "rocksdb")
    diamond_db_path = os.path.join(
        orthosets_dir,
        orthoset,
        "diamond",
        orthoset + ".dmnd",
    )
    # Legacy db check
    if not os.path.exists(orthoset_db_path) or not os.path.exists(diamond_db_path):
        gen = False
        print("Could not find orthoset directory or diamond database.")
        if os.path.exists(f"{orthoset}.sqlite"):
            if (
                input(
                    f"Could not find orthoset DB at {orthoset_db_path}.\nDetected legacy orthoset: {orthoset}.sqlite\nWould you like to generate new Orthoset DB? Y/N: ",
                ).lower()
                == "y"
            ):
                print("Attempting to generate DB")
                os.system(
                    f"python3 -m sapphyre -p {num_threads} Makeref {orthoset}.sqlite -s {orthoset} -all",
                )
                gen = True
        if not gen:
            print("Aborting")
            return False
    orthoset_db = wrap_rocks.RocksDB(orthoset_db_path)

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
    time_keeper.lap()  # Reset timer

    printv("Done! Grabbing NT sequences.", args.verbose)

    nt_db_path = os.path.join(input_path, "rocksdb", "sequences", "nt")

    nt_db = wrap_rocks.RocksDB(nt_db_path)

    # Check if sequence is assembly
    dbis_assembly = nt_db.get("get:isassembly")
    evalue = args.evalue
    min_orf = NORMAL_MIN_ORF
    is_assembly = False
    if dbis_assembly and dbis_assembly == "True":
        is_assembly = True
        evalue = ASSEMBLY_EVALUE
        min_orf = ASSEMBLY_MIN_ORF

    precision = decimal.Decimal("5") * decimal.Decimal("0.1") ** decimal.Decimal(
        evalue,
    )

    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
    else:
        msg = (
            "Nucleotide sequence not found in database. Did Prepare succesfully finish?"
        )
        raise ValueError(
            msg,
        )
    out = [nt_db.get_bytes(f"ntbatch:{i}") for i in recipe]

    dupe_counts = json.decode(nt_db.get_bytes("getall:dupes"), type=dict[str, int])

    out_path = os.path.join(diamond_path, f"{sensitivity}")
    extension_found = False
    for extension in [".tsv", ".tsv.tar.gz", ".gz", ".tsv.gz"]:
        if os.path.exists(out_path + extension):
            out_path += extension
            extension_found = True
            break

    if not os.path.exists(out_path) or os.stat(out_path).st_size == 0:
        with TemporaryDirectory(dir=gettempdir()) as dir, NamedTemporaryFile(
            dir=dir,
        ) as input_file:
            input_file.write(b"".join(out))
            input_file.flush()

            quiet = "--quiet" if args.verbose == 0 else ""

            printv(
                f"Done! Running Diamond. Elapsed time {time_keeper.differential():.2f}s.",
                args.verbose,
            )
            time_keeper.lap()  # Reset timer
            if not extension_found:
                out_path += ".tsv"
            os.system(
                f"diamond blastx -d {diamond_db_path} -q {input_file.name} -o {out_path} --{sensitivity}-sensitive --masking 0 -e {precision} --compress 1 --outfmt 6 qseqid sseqid qframe evalue bitscore qstart qend sstart send {quiet} --top {top_amount} --min-orf {min_orf} --max-hsps 0 -p {num_threads}",
            )
            out_path += ".gz"
            input_file.seek(0)

        printv(
            f"Diamond done. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Reading file to memory",
            args.verbose,
        )
    else:
        printv(
            f"Found existing Diamond output. Elapsed time {time_keeper.differential():.2f}s. Reading file to memory",
            args.verbose,
        )

    db = wrap_rocks.RocksDB(os.path.join(input_path, "rocksdb", "hits"))
    output = defaultdict(list)
    multi_kicks = 0

    global_log = []
    dupe_divy_headers = defaultdict(set)

    if os.stat(out_path).st_size == 0:
        printv("Diamond returned zero hits.", args.verbose, 0)
        return True

    df = pd.read_csv(
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
            "frame": np.int8,
            "evalue": np.float64,
            "score": np.float32,
            "qstart": np.uint16,
            "qend": np.uint16,
            "sstart": np.uint16,
            "send": np.uint16,
        },
    )

    target_counts = df["target"].value_counts()
    combined_count = Counter()
    taxon_to_targets = defaultdict(list)
    for target, count in target_counts.to_dict().items():
        _, ref_taxa, _ = target_to_taxon[target]
        taxon_to_targets[ref_taxa].append(target)
        combined_count[ref_taxa] += count

    top_refs = set()
    pairwise_refs = set()
    top_targets = set()
    most_common = combined_count.most_common()
    min_count = min(most_common[0:5], key=lambda x: x[1])[1]
    target_count = min_count - (min_count * args.top_ref)
    for taxa, count in most_common:
        if count >= target_count:
            top_refs.add(taxa)
            top_targets.update(taxon_to_targets[taxa])
        if count >= min_count:
            pairwise_refs.add(taxa)
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
                        np.where(
                            df[start_index:]["header"].values
                            == last_header,
                        )[0][-1]
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
                    end_index = (
                        np.where(
                            df[start_index : (start_index + estimated_end)]["header"].values
                            == last_header,
                        )[0][-1]
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
                is_assembly,
                precision,
            )
            for i, index in enumerate(indices)
        )

        printv(
            f"Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Processing data.",
            args.verbose,
        )
        if post_threads > 1:
            with Pool(post_threads) as pool:
                pool.map(process_lines, arguments)
        else:
            for arg in arguments:
                process_lines(arg)

        printv(
            f"Processed. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Reading thread outputs",
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
            f"Done reading outputs. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Doing internal filters",
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
            with Pool(post_threads) as pool:
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

        next_step = (
            "Writing to db" if not is_assembly else "Doing Assembly Containments"
        )
        printv(
            f"Filtering done. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. {next_step}",
            args.verbose,
        )

        # containment_kicks = 0
        # containment_log = []
        # arguments = []
        # present_genes = []
        # if is_assembly:
        #     for gene, hits in output.items():
        #         arguments.append(
        #             ContainmentArgs(
        #                 hits,
        #                 target_to_taxon,
        #                 args.debug,
        #                 gene,
        #             ),
        #         )
        #
        #     if post_threads > 1:
        #         with Pool(post_threads) as pool:
        #             results = pool.map(containments, arguments)
        #     else:
        #         results = [containments(arg) for arg in arguments]
        #
        #     next_output = []
        #     for result in results:
        #         containment_kicks += result.kick_count
        #         if args.debug:
        #             containment_log.extend(result.log)
        #
        #         next_output.append(
        #             (
        #                 result.gene,
        #                 [i for i in output[result.gene] if i.uid not in result.kicks],
        #             ),
        #         )
        #         present_genes.append(result.gene)
        #     output = next_output
        # else:
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

        db.put_bytes(
            "getall:target_variants",
            json_encoder.encode(variant_filter),
        )  # type=dict[str, list[str]]

        del variant_filter
        nt_db.put("getall:valid_refs", ",".join(list(top_refs)))

        head_to_seq = {}
        for content in out:
            lines = content.decode().split("\n")
            head_to_seq.update(
                {
                    lines[i][1:]: lines[i + 1]
                    for i in range(0, len(lines), 2)
                    if lines[i] != ""
                },
            )
        del out

        if is_assembly:
            arguments = []
            for gene, hits in output:
                arguments.append(
                    ConvertArgs(
                        hits,
                        gene,
                    ),
                )

            if post_threads > 1:
                with Pool(post_threads) as pool:
                    output = pool.map(convert_and_cull, arguments)
            else:
                output = [convert_and_cull(arg) for arg in arguments]

        passes = 0
        encoder = json.Encoder()
        for result in output:
            if is_assembly:
                hits, gene = result.hits, result.gene
            else:
                gene, hits = result

            out = []
            for hit in hits:
                if not is_assembly:
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
                    hit.seq = phymmr_tools.bio_revcomp(hit.seq)
                out.append(hit)
                dupe_divy_headers[gene].add(hit.node)

            passes += len(out)
            db.put_bytes(f"gethits:{gene}", encoder.encode(out))

        del head_to_seq
        if global_log:
            with open(os.path.join(input_path, "multi.log"), "w") as fp:
                fp.write("\n".join(global_log))
        # if containment_log:
        #     with open(os.path.join(input_path, "containment.log"), "w") as fp:
        #         fp.write("\n".join(containment_log))

        printv(
            f"{multi_kicks} multi kicks",
            args.verbose,
        )
        printv(
            f"{internal_kicks} internal kicks",
            args.verbose,
        )
        # printv(
        #     f"{containment_kicks} containment kicks",
        #     args.verbose,
        # )
        printv(
            # f"Took {time_keeper.lap():.2f}s for {multi_kicks+internal_kicks+containment_kicks} kicks leaving {passes} results. Writing to DB",
            f"Took {time_keeper.lap():.2f}s for {multi_kicks + internal_kicks} kicks leaving {passes} results. Writing to DB",
            args.verbose,
        )

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

    printv(f"Took {time_keeper.differential():.2f}s overall.", args.verbose)

    return True


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
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
