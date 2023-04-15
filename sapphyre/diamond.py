from argparse import Namespace
from collections import Counter, defaultdict, namedtuple
import itertools
import json
from math import ceil
import os
from shutil import rmtree
import sys
from tempfile import TemporaryDirectory, NamedTemporaryFile
from multiprocessing.pool import Pool
from time import time
import numpy as np
import pandas as pd
import orjson
import wrap_rocks

# from phymmr_tools import Hit, ReferenceHit
from .utils import printv, gettempdir
from .timekeeper import TimeKeeper, KeeperMode

# namedtuple for the processing arguments.
ProcessingArgs = namedtuple(
    "ProcessingArgs",
    [
        "i",
        "grouped_data",
        "target_to_taxon",
        "debug",
        "result_fp"
    ],
)


# Define a function lst that returns an empty list when called.
# This function will be used as the default_factory argument for a defaultdict object.
def lst() -> list:
    return []


# call the phymmr_tools.Hit class using an unpacked list of values found in the row input.
# def create_hit(row) -> Hit:
#     return Hit(*row)


class ReferenceHit:
    __slots__ = (
        "target",
        "sstart",
        "send",
    )

    def __init__(self, target, sstart, send):
        self.target = target
        self.sstart = int(sstart)
        self.send = int(send)

    def to_json(self):
        return {
            "target": self.target,
            "sstart": self.sstart,
            "send": self.send,
        }


class Hit:
    __slots__ = (
        "header",
        "qstart",
        "qend",
        "gene",
        "score",
        "evalue",
        "length",
        "kick",
        "seq",
        "frame",
        "uid",
        "reftaxon",
        "target",
        "sstart",
        "send",
        "reference_hits",
        "json",
    )

    def __init__(self, df):
        self.uid = hash(time())
        (
            self.header,
            self.target,
            self.frame,
            self.evalue,
            self.score,
            self.qstart,
            self.qend,
            self.sstart,
            self.send,
        ) = df
        self.gene = None
        self.reftaxon = None
        self.kick = False
        self.seq = None
        self.json = None
        self.reference_hits = [ReferenceHit(self.target, self.sstart, self.send)]
        if self.frame < 0:
            self.qend, self.qstart = self.qstart, self.qend
        self.length = self.qend - self.qstart + 1

    def convert_reference_hits(self):
        self.reference_hits = [i.to_json() for i in self.reference_hits]

    def to_json(self):
        return {
            "header": self.header,
            "frame": self.frame,
            "seq": self.seq,
            "taxon": self.reftaxon,
            "ali_start": self.qstart,
            "ali_end": self.qend,
            "ref_hits": self.reference_hits,
            "kick": False,
            "uid":self.uid,
            "score":self.score,
        }


def get_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """
    Get the overlap between two ranges.

    Args:
        a_start (int): The starting position of range A.
        a_end (int): The ending position of range A.
        b_start (int): The starting position of range B.
        b_end (int): The ending position of range B.

    Returns:
        int: The number of elements in the overlap between the two ranges.
    """
    # Calculate the left most position out of each range's end.
    overlap_end = min(a_end, b_end)
    # Calculate the right most position out of each range's start.
    overlap_start = max(a_start, b_start)

    # If the ranges do not overlap, return 0.
    if overlap_end < overlap_start:
        return 0

    # Calculate the number of elements in the overlap.
    amount = (overlap_end - overlap_start) + 1
    return amount


def get_score_difference(score_a: float, score_b: float) -> float:
    """
    Get the decimal difference between two scores.

    Args:
        score_a (float): The first score.
        score_b (float): The second score.

    Returns:
        float: The decimal difference between the two scores.
    """
    # If either score is zero return zero.
    if score_a == 0.0 or score_b == 0.0:
        return 0.0

    # Return the decimal difference between the largest score and the smallest score.
    return max(score_a, score_b) / min(score_a, score_b)


def multi_filter(hits: list, debug: bool) -> tuple[list, int, list]:
    """
    Filter out hits that imprecisely map to multiple genes.

    Args:
        hits (list): A list of Hit objects.
        debug (bool): A boolean indicating whether or not to print debug messages.
    Returns:
        tuple[list, int, list]: A tuple containing a list of filtered hits, the number of hits filtered, and a list of debug rows.
    """
    # Global variables for the multi_filter function.
    MULTI_PERCENTAGE_OF_OVERLAP = 0.3
    MULTI_SCORE_DIFFERENCE = 1.05

    log = []

    # Assign the highest scoring hit as the master
    master = hits[0]

    # Assign all the lower scoring hits as candidates
    candidates = hits[1:]

    for i, candidate in enumerate(candidates, 1):
        # Skip if the candidate has been kicked already or if master and candidate are internal hits:
        if not candidate or master.gene == candidate.gene:
            continue

        # Get the distance between the master and candidate
        distance = master.qend - master.qstart

        # Get the amount and percentage overlap between the master and candidate
        amount_of_overlap = get_overlap(
            master.qstart, master.qend, candidate.qstart, candidate.qend
        )
        percentage_of_overlap = amount_of_overlap / distance

        # If the overlap is greater than 30% and the score difference is greater than 5%
        if percentage_of_overlap >= MULTI_PERCENTAGE_OF_OVERLAP:
            score_difference = get_score_difference(master.score, candidate.score)
            if score_difference >= MULTI_SCORE_DIFFERENCE:
                # Kick the candidate
                hits[i].kick = True
                candidates[i - 1].kick = True
                if debug:
                    log.append(
                        (
                            candidate.gene,
                            candidate.header,
                            candidate.reftaxon,
                            candidate.score,
                            candidate.qstart,
                            candidate.qend,
                            "Kicked out by",
                            master.gene,
                            master.header,
                            master.reftaxon,
                            master.score,
                            master.qstart,
                            master.qend,
                        )
                    )
            else:
                # If the score difference is not greater than 5% trigger miniscule score
                if debug:
                    log.extend(
                        [
                            (
                                hit.gene,
                                hit.header,
                                hit.reftaxon,
                                hit.score,
                                hit.qstart,
                                hit.qend,
                                "Kicked due to miniscule score",
                                master.gene,
                                master.header,
                                master.reftaxon,
                                master.score,
                                master.qstart,
                                master.qend,
                            )
                            for hit in hits
                            if hit
                        ]
                    )
                for hit in hits:
                    hit.kick = True
                return len(hits), log

    # Filter out the None values
    passes = [hit for hit in hits if not hit.kick]
    return len(hits) - len(passes), log


def internal_filter(
    gene: str, header_based: dict, debug: bool, internal_percent: float
) -> tuple[set, list, int]:
    """
    Filters out overlapping hits who map to the same gene.

    Args:
        gene (str): The gene to filter.
        header_based (dict): A dictionary of hits grouped by header.
        debug (bool): Whether to print debug information.
        internal_percent (float): The percentage of overlap required to constitute a kick.
    Returns:
        tuple[set, list, int]:
            A tuple containing the set of kicked hit headers,
            the log of kicked hits, and the number of hits kicked.
    """
    log = []
    kicks = itertools.count()
    this_kicks = set()

    for hits in header_based.values():
        for hit_a, hit_b in itertools.combinations(hits, 2):
            if hit_a["kick"] or hit_b["kick"]:
                continue

            # Calculate overlap percent over the internal overlap's length
            overlap_amount = get_overlap(
                hit_a["ali_start"], hit_a["ali_end"], hit_b["ali_start"], hit_b["ali_end"]
            )
            internal_start = min(hit_a["ali_start"], hit_b["ali_start"])
            internal_end = max(hit_a["ali_end"], hit_b["ali_end"])
            internal_length = internal_end - internal_start
            percent = overlap_amount / internal_length

            # If bp overlap and the percent of overlap is greater than or equal to the internal percent arg
            if overlap_amount > 0 and percent >= internal_percent:
                # Iterate the kicks counter and set the kick attribute to True
                next(kicks)
                hit_b["kick"] = True

                # Add the hit to the kicked set
                this_kicks.add(hit_b["uid"])

                if debug:
                    log.append(
                        (
                            gene,
                            hit_b["uid"],
                            hit_b["taxon"],
                            round(hit_b["score"], 2),
                            hit_b["ali_start"],
                            hit_b["ali_end"],
                            "Internal kicked out by",
                            gene,
                            hit_a["uid"],
                            hit_a["taxon"],
                            round(hit_a["score"], 2),
                            hit_a["ali_start"],
                            hit_a["ali_end"],
                            round(percent, 3),
                        )
                    )

    return this_kicks, log, next(kicks)


def internal_filtering(
    gene: str, hits: list, debug: bool, internal_percent: float
) -> tuple[str, int, list, list]:
    """
    Performs the internal filter on a list of hits and returns a count of the number of hits
    kicked, a log of the hits kicked, and the list of hits that passed the filter.

    Args:
        gene (str): The gene that the hits belong to.
        hits (list): The list of hits to filter.
        debug (bool): Whether to log debug information or not.

    Returns:
        tuple[str, int, list, list]:
            A tuple containing the gene, the number of hits kicked, a log of the hits kicked,
            and the list of hits that passed the filter.
    """
    # Perform the internal filter.
    kicked_hits, this_log, this_kicks = internal_filter(gene, hits, debug, internal_percent)

    # Return the gene, the number of hits kicked, a log of the hits kicked, and the list of hits
    return (gene, this_kicks, this_log, kicked_hits)


def process_lines(pargs: ProcessingArgs) -> tuple[dict[str, Hit], int, list[str]]:
    """
    Process a subset of lines from the Diamond tsv result.

    Args:
        pargs (ProcessingArgs): The arguments for processing the lines.
    Returns:
        tuple[dict[str, Hit], int, list[str]]:
            A tuple containing the dictionary of hits grouped by header,
            the number of hits kicked, and the list of kicked hits.
    """
    output = defaultdict(list)
    multi_kicks = 0
    this_log = []

    for _, header_df in pargs.grouped_data.groupby("header"):
        frame_to_hits = defaultdict(list)
        # Convert each row in the dataframe to a Hit() object
        hits = [Hit(row) for row in header_df.values]
        hits.sort(key=lambda hit: hit.score, reverse=True)
        # Add reference data to each
        for hit in hits:
            hit.gene, hit.reftaxon, _ = pargs.target_to_taxon[hit.target]
            frame_to_hits[hit.frame].append(hit)

        for hits in frame_to_hits.values():
            genes_present = {hit.gene for hit in hits}

            if len(genes_present) > 1:
                this_kicks, log = multi_filter(hits, pargs.debug)
                multi_kicks += this_kicks
                if pargs.debug:
                    this_log.extend(log)

            if any(hits):
                if len(genes_present) > 1:
                    gene_hits = defaultdict(lst)
                    for hit in hits:
                        if not hit.kick:
                            gene_hits[hit.gene].append(hit)
                else:
                    gene_hits = {hit.gene: [hit for hit in hits if not hit.kick]}

                for _, hits in gene_hits.items():
                    top_hit = hits[0]

                    ref_seqs = [
                        ReferenceHit(hit.target, hit.sstart, hit.send)
                        for hit in hits
                        if hit.gene == top_hit.gene and hit != top_hit
                    ]
                    top_hit.reference_hits.extend(ref_seqs)
                    top_hit.convert_reference_hits()

                    output[top_hit.gene].append(top_hit.to_json())

    with open(pargs.result_fp, "w") as result_file:
        result_file.write(json.dumps(output))

    return pargs.i, multi_kicks, this_log


def run_process(args: Namespace, input_path: str) -> bool:
    """
    Run the main process on the input path.

    Args:
        args (Namespace): The arguments for the program.
        input_path (str): The path to the input directory.
    Returns:
        bool: Whether the process was successful or not.
    """
    # Estimated size of header chunks in a sub-section of the dataframe.
    CHUNK_ESTIMATED_SIZE = 6

    # Due to the thread bottleneck of chunking a ceiling is set on the threads post reporter
    THREAD_CAP = 32

    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    strict_search_mode = args.strict_search_mode
    sensitivity = args.sensitivity
    top_amount = args.top

    taxa = os.path.basename(input_path)
    printv(f"Processing: {taxa}", args.verbose, 0)
    printv("Grabbing reference data from Orthoset DB.", args.verbose)
    # make dirs
    diamond_path = os.path.join(input_path, "diamond")
    if args.overwrite:
        if os.path.exists(diamond_path):
            rmtree(diamond_path)
    os.makedirs(diamond_path, exist_ok=True)

    num_threads = args.processes
    post_threads = args.processes if args.processes < THREAD_CAP else THREAD_CAP

    orthoset_db_path = os.path.join(orthosets_dir, orthoset, "rocksdb")
    diamond_db_path = os.path.join(
        orthosets_dir, orthoset, "diamond", orthoset + ".dmnd"
    )
    # Potato brain check
    if not os.path.exists(orthoset_db_path) or not os.path.exists(diamond_db_path):
        if (
            input(
                f"Could not find orthoset DB at {orthoset_db_path}. Would you like to generate it? Y/N: "
            ).lower()
            == "y"
        ):
            print("Attempting to generate DB")
            os.system(
                f"python3 -m sapphyre -p {num_threads} Makeref {orthoset}.sqlite -s {orthoset}"
            )
        else:
            print("Aborting")
            sys.exit(1)
    orthoset_db = wrap_rocks.RocksDB(orthoset_db_path)

    target_to_taxon = orjson.loads(orthoset_db.get("getall:targetreference"))

    del orthoset_db
    time_keeper.lap()  # Reset timer

    printv("Done! Grabbing NT sequences.", args.verbose)

    nt_db_path = os.path.join(input_path, "rocksdb", "sequences", "nt")

    nt_db = wrap_rocks.RocksDB(nt_db_path)

    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
    else:
        raise ValueError(
            "Nucleotide sequence not found in database. Did Prepare succesfully finish?"
        )
    out = [nt_db.get(f"ntbatch:{i}") for i in recipe]

    dupe_counts = orjson.loads(nt_db.get("getall:dupes"))

    out_path = os.path.join(diamond_path, f"{sensitivity}.tsv")
    if not os.path.exists(out_path) or os.stat(out_path).st_size == 0:
        with TemporaryDirectory(dir=gettempdir()) as dir, NamedTemporaryFile(
            dir=dir
        ) as input_file:
            input_file.write("".join(out).encode())
            input_file.flush()

            quiet = "--quiet" if args.verbose == 0 else ""

            printv(
                f"Done! Running Diamond. Elapsed time {time_keeper.differential():.2f}s.",
                args.verbose,
            )
            time_keeper.lap()  # Reset timer
            os.system(
                f"diamond blastx -d {diamond_db_path} -q {input_file.name} -o {out_path} --{sensitivity}-sensitive --masking 0 -e {args.evalue} --outfmt 6 qseqid sseqid qframe evalue bitscore qstart qend sstart send {quiet} --top {top_amount} --min-orf 20 --max-hsps 0 -p {num_threads}"
            )
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
    output = {}
    multi_kicks = 0

    global_log = []
    dupe_divy_headers = {}
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
            "frame": "int8",
            "evalue": np.float64,
            "score": np.float32,
            "qstart": "int16",
            "qend": "int16",
            "sstart": "int32",
            "send": "int32",
        },
    )
    target_counts = df["target"].value_counts()
    combined_count = Counter()
    taxon_to_targets = {}
    for target, count in target_counts.to_dict().items():
        _, ref_taxa, _ = target_to_taxon[target]
        taxon_to_targets.setdefault(ref_taxa, []).append(target)
        combined_count[ref_taxa] += count

    top_refs = set()
    top_targets = set()
    most_common = combined_count.most_common()
    min_count = min(most_common[0:5], key=lambda x: x[1])[1]
    target_count = min_count - (min_count * args.top_ref)
    for taxa, count in most_common:
        if count >= target_count:
            top_refs.add(taxa)
            top_targets.update(taxon_to_targets[taxa])
    target_has_hit = set(df["target"].unique())
    df = df[(df["target"].isin(top_targets))]
    headers = df["header"].unique()
    if len(headers) > 0:
        per_thread = ceil(len(headers) / args.processes)
        estimated_end = ceil(CHUNK_ESTIMATED_SIZE * per_thread)
        arguments = []
        indices = []
        for x, i in enumerate(np.arange(0, len(headers), per_thread), 1):
            if x == 1:
                start_index = 0
            else:
                start_index = end_index + 1

            if x != post_threads:
                last_header = headers[i + per_thread - 1]
                end_index = (
                    np.where(df[start_index:(start_index+estimated_end)]["header"].values == last_header)[0][-1]
                    + start_index
                )
            else:
                end_index = len(df) - 1

            indices.append((start_index, end_index))

        temp_files = [NamedTemporaryFile(dir=gettempdir()) for _ in range(post_threads)]

        arguments = (
            ProcessingArgs(
                i, 
                df.iloc[index[0] : index[1] + 1],
                target_to_taxon,
                args.debug,
                temp_files[i].name
            )
            for i, index in enumerate(indices)
        )

        printv(
            f"Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Processing data.",
            args.verbose,
        )

        with Pool(post_threads) as p:
            result = p.map(process_lines, arguments)

        del arguments
        for process_index, mkicks, this_log in result:
            this_result = temp_files[process_index]
            this_output = json.load(open(this_result.name, "r"))
            this_result.close()
            
            for gene, hits in this_output.items():
                if gene not in output:
                    output[gene] = hits
                else:
                    output[gene].extend(hits)

            multi_kicks += mkicks

            if args.debug:
                global_log.extend(this_log)
        del this_output, mkicks, this_log, result
        printv(
            f"Processed. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Doing internal filters",
            args.verbose,
        )

        requires_internal = {}
        internal_order = []
        for gene, hits in output.items():
            this_counter = Counter([i["header"] for i in hits]).most_common()
            requires_internal[gene] = {}
            if this_counter[0][1] > 1:
                this_hits = np.sum(i[1] for i in this_counter if i[1] > 1)
                this_common = {i[0] for i in this_counter if i[1] > 1}
                for hit in [i for i in hits if i["header"] in this_common]:
                    requires_internal[gene].setdefault(hit["header"], []).append(hit)

                internal_order.append((gene, this_hits))

        internal_order.sort(key=lambda x: x[1], reverse=True)

        with Pool(post_threads) as pool:
            internal_results = pool.starmap(
                internal_filtering,
                [
                    (gene, requires_internal[gene], args.debug, args.internal_percent)
                    for gene, _ in internal_order
                ],
            )

        internal_result = defaultdict(lst)
        for gene, kick_count, this_log, kicked_hits in internal_results:
            internal_result[gene] = (kick_count, this_log, kicked_hits)

        printv(
            f"Filtering done. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Writing to db",
            args.verbose,
        )

        head_to_seq = {}
        for content in out:
            lines = content.split("\n")
            head_to_seq.update(
                {
                    lines[i][1:]: lines[i + 1]
                    for i in np.arange(0, len(lines), 2)
                    if lines[i] != ""
                }
            )
        del out

        passes = 0
        internal_kicks = 0
        for gene, hits in output.items():
            dupe_divy_headers[gene] = {}
            kicks = set()
            out = []
            if gene in internal_result:
                this_kicks, this_log, kicks = internal_result[gene]

                internal_kicks += this_kicks
                if args.debug:
                    global_log.extend(this_log)

            for hit in hits:
                if hit["uid"] in kicks:
                    continue
                hit["seq"] = head_to_seq[hit["header"]]
                out.append(hit)
                dupe_divy_headers[gene][hit["header"]] = 1

            passes += len(out)
            db.put_bytes(f"gethits:{gene}", orjson.dumps(out))
        del head_to_seq
        if global_log:
            with open(os.path.join(input_path, "multi.log"), "w") as fp:
                fp.write(
                    "\n".join([",".join([str(i) for i in line]) for line in global_log])
                )

        printv(
            f"{multi_kicks} multi kicks",
            args.verbose,
        )
        printv(
            f"{internal_kicks} internal kicks",
            args.verbose,
        )
        printv(
            f"Took {time_keeper.lap():.2f}s for {multi_kicks+internal_kicks} kicks leaving {passes} results. Writing dupes and present gene data",
            args.verbose,
        )

        gene_dupe_count = {}
        for gene, headers in dupe_divy_headers.items():
            for base_header in headers.keys():
                if base_header in dupe_counts:
                    gene_dupe_count.setdefault(gene, {})[base_header] = dupe_counts[
                        base_header
                    ]

        db.put("getall:presentgenes", ",".join(list(output.keys())))

        key = "getall:gene_dupes"
        data = orjson.dumps(gene_dupe_count)
        nt_db.put_bytes(key, data)

    # DOING VARIANT FILTER
    variant_filter = {}

    for target, ref_tuple in target_to_taxon.items():
        gene, ref_taxon, data_length = ref_tuple
        if ref_taxon in top_refs:
            variant_filter.setdefault(gene, []).append((ref_taxon, target, data_length))

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
                variants_with_hits = np.sum(
                    i[1] in target_has_hit for i in this_targets
                )
                all_variants_kicked = variants_with_hits == 0
                if all_variants_kicked:
                    reintroduce = max(this_targets, key=lambda x: x[2])
                    out_targets.append(reintroduce[1])
                    continue

                out_targets.extend(
                    [i[1] for i in this_targets if i[1] in target_has_hit]
                )

            variant_filter[gene] = out_targets
        else:
            variant_filter.pop(gene, -1)

    variant_filter = {k: list(v) for k, v in variant_filter.items()}
    db.put_bytes("getall:target_variants", orjson.dumps(variant_filter))

    del variant_filter
    nt_db.put("getall:valid_refs", ",".join(list(top_refs)))
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
    raise Exception(
        "Cannot be called directly, please use the module:\nsapphyre Diamond"
    )
