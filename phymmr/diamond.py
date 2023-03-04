from collections import Counter, namedtuple
import itertools
from math import ceil
import os
from shutil import rmtree
import sys
from tempfile import TemporaryDirectory, NamedTemporaryFile
import json
from multiprocessing.pool import Pool
import wrap_rocks

from .utils import printv, gettempdir
from .timekeeper import TimeKeeper, KeeperMode


def reciprocal_check(hits, strict, taxa_present):
    first = None
    hit_on_taxas = {i: 0 for i in taxa_present}

    for hit in hits:
        if not strict:
            return hit
        first = hit

        hit_on_taxas[hit.reftaxon] = 1  # Just need 1 hit

        if strict:
            if all(hit_on_taxas.values()):
                return first

    if any(hit_on_taxas.values()):
        return first

    return None


class Hit:
    __slots__ = (
        "header",
        "qstart",
        "qend",
        "gene",
        "score",
        "reftaxon",
        "seq",
        "evalue",
        "length",
        "kick",
        "frame",
        "full_header",
        "target"
    )

    def __init__(
        self, header, ref_header, frame, evalue, score, qstart, qend, gene, reftaxon
    ):
        self.header = header
        self.target = ref_header
        self.gene = gene
        self.reftaxon = reftaxon
        self.seq = None
        self.score = float(score)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.evalue = float(evalue)
        self.kick = False
        self.frame = int(frame)
        if self.frame < 0:
            self.qend, self.qstart = self.qstart, self.qend
        self.length = self.qend - self.qstart + 1

        if self.frame < 0:
            self.full_header = (
                self.header + f"|[revcomp]:[translate({abs(self.frame)})]"
            )
        else:
            self.full_header = self.header + f"|[translate({self.frame})]"

    def to_json(self):
        return {
            "header": self.full_header,
            "seq": self.seq,
            "ref_taxon": self.reftaxon,
            "ali_start": self.qstart,
            "ali_end": self.qend,
        }


def get_overlap(a_start, a_end, b_start, b_end):
    overlap_end = min(a_end, b_end)
    overlap_start = max(a_start, b_start)
    amount = (overlap_end - overlap_start) + 1  # inclusive

    return 0 if amount < 0 else amount


def get_difference(scoreA, scoreB):
    """
    Returns decimal difference of two scores
    """
    if scoreA == 0 or scoreB == 0:
        return 0
    return max(scoreA, scoreB) / min(scoreA, scoreB)


def multi_filter(hits, debug):
    kick_happend = True

    log = []

    while kick_happend:
        kick_happend = False
        master = hits[0]
        candidates = hits[1:]

        master_env_start = master.qstart
        master_env_end = master.qend

        miniscule_score = False
        for i, candidate in enumerate(candidates, 1):
            if candidate:
                if master.gene != candidate.gene:
                    distance = (master_env_end - master_env_start) + 1  # Inclusive
                    amount_of_overlap = get_overlap(
                        master_env_start,
                        master_env_end,
                        candidate.qstart,
                        candidate.qend,
                    )
                    percentage_of_overlap = amount_of_overlap / distance

                    if percentage_of_overlap >= 0.3:  # min_overlap_multi:
                        score_difference = get_difference(master.score, candidate.score)
                        if score_difference >= 1.05:
                            kick_happend = True
                            hits[i] = None
                            candidates[i - 1] = None
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
                            miniscule_score = True
                            break
        if miniscule_score:
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
            return [], len(hits), log

    passes = [i for i in hits if i]
    return passes, len(hits) - len(passes), log

def internal_filter(header_based: dict, debug: bool, internal_percent: float) -> list:
    log = []
    kicks = 0
    this_kicks = set()

    for hits in header_based.values():
        for hit_a, hit_b in itertools.combinations(hits, 2):
            if hit_a.kick or hit_b.kick:
                continue

            overlap_amount = get_overlap(
                hit_a.qstart, hit_a.qend, hit_b.qstart, hit_b.qend
            )
            percent = overlap_amount / hit_a.length if overlap_amount > 0 else 0

            if percent >= internal_percent:
                kicks += 1
                hit_b.kick = True

                this_kicks.add(hit_b.full_header)

                if debug:
                    log.append(
                        (
                            hit_b.gene,
                            hit_b.header,
                            hit_b.reftaxon,
                            hit_b.score,
                            hit_b.qstart,
                            hit_b.qend,
                            "Internal kicked out by",
                            hit_a.gene,
                            hit_a.header,
                            hit_a.reftaxon,
                            hit_a.score,
                            hit_a.qstart,
                            hit_a.qend,
                        )
                    )

    return this_kicks, log, kicks


def internal_filtering(gene, hits, debug, internal_percent):
    kicked_hits, this_log, this_kicks = internal_filter(hits, debug, internal_percent)

    return {gene: (this_kicks, this_log, kicked_hits)}


def count_reftaxon(file_pointer, taxon_lookup: dict, percent: float) -> list:
    """
    Counts reference taxon in very.tsv file. Returns a list
    containg the 5 most popular names.
    Possible speed if we cut the string at the : and only
    lookup names after we know the counts?
    """
    rextaxon_count = {}
    header_lines = {}
    current_header = None
    target_has_hit = set()

    for line in file_pointer:
        header, ref_header_pair, frame = line.split("\t")[:3]
        if hash(header) != current_header:
            current_header = hash(header)
            ref_variation_filter = set()
            
        ref_gene, ref_taxon, _ = taxon_lookup[ref_header_pair]

        target_has_hit.add(ref_header_pair)

        ref_variation_key = header+ref_taxon
        if not ref_variation_key in ref_variation_filter:
            ref_variation_filter.add(ref_variation_key)
            rextaxon_count[ref_taxon] = rextaxon_count.get(ref_taxon, 0) + 1
        header_lines.setdefault(header + frame, []).append(
            line.strip() + f"\t{ref_gene}\t{ref_taxon}"
        )
    sorted_counts = [x for x in rextaxon_count.items()]
    top_names = []
    total_references = 0
    if sorted_counts:
        sorted_counts.sort(key=lambda x: x[1], reverse=True)
        target_count = min([i[1] for i in sorted_counts[0:5]])
        target_count = target_count - (target_count * percent)
        total_references = len(sorted_counts)
        top_names = set([x[0] for x in sorted_counts if x[1] >= target_count])

    return top_names, total_references, list(header_lines.values()), target_has_hit

ProcessingArgs = namedtuple(
    "ProcessingArgs",
    [
        "lines",
        "debug",
        "top_refs",
        "total_references",
        "strict_search_mode",
        "reference_taxa",
        "evalue",

    ],
)

def process_lines(
    pargs: ProcessingArgs
):
    output = {}
    multi_kicks = 0
    this_log = []
    for this_lines in pargs.lines:
        hits = [Hit(*hit.split("\t")) for hit in this_lines]  # convert to Hit object
        hits = list(filter(lambda x: x.evalue <= pargs.evalue, hits))
        hits.sort(key=lambda x: x.score, reverse=True)
        genes_present = {hit.gene for hit in hits}

        if len(genes_present) > 1:
            hits, this_kicks, log = multi_filter(hits, pargs.debug)
            multi_kicks += this_kicks
            if pargs.debug:
                this_log.extend(log)

        best_hit = reciprocal_check(hits, pargs.strict_search_mode, pargs.reference_taxa)

        if best_hit:
            output.setdefault(best_hit.gene, []).append(best_hit)

    return output, multi_kicks, this_log


def run_process(args, input_path) -> None:
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
                f"python3 -m phymmr -p {num_threads} Makeref {orthoset}.sqlite -s {orthoset}"
            )
        else:
            print("Aborting")
            sys.exit(1)
    orthoset_db = wrap_rocks.RocksDB(orthoset_db_path)

    reference_taxa = json.loads(orthoset_db.get("getall:taxainset"))
    target_to_taxon = json.loads(orthoset_db.get("getall:targetreference"))

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

    dupe_counts = json.loads(nt_db.get("getall:dupes"))

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
                f"diamond blastx -d {diamond_db_path} -q {input_file.name} -o {out_path} --{sensitivity}-sensitive --masking 0 -e {args.evalue} --outfmt 6 qseqid sseqid qframe evalue bitscore qstart qend {quiet} --top {top_amount} --max-hsps 0 -p {num_threads}"
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

    chunks = args.chunks
    chunk_count = itertools.count(1)

    global_log = []
    dupe_divy_headers = {}
    with open(out_path) as fp:
        top_refs, total_references, lines, target_has_hit = count_reftaxon(
            fp, target_to_taxon, args.top_ref
        )
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
                variants_with_hits = sum([i[1] in target_has_hit for i in this_targets])
                all_variants_kicked = variants_with_hits == 0
                if all_variants_kicked:
                    reintroduce = max(this_targets, key = lambda x : x[2])
                    out_targets.append(reintroduce[1])
                    continue
                
                out_targets.extend([i[1] for i in this_targets if i[1] in target_has_hit])
            
            variant_filter[gene] = out_targets
        else:
            variant_filter.pop(gene, -1)

    variant_filter = {k: list(v) for k, v in variant_filter.items()}
    db.put("getall:target_variants", json.dumps(variant_filter))

    del variant_filter


    printv(
        f"Processing {chunks} chunk(s). Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s",
        args.verbose,
    )
    if lines:
        per_thread = ceil(ceil(len(lines) / chunks) / num_threads)
        for i in range(0, chunks*num_threads, num_threads):
            printv(
                f"Processing chunk {next(chunk_count)}. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s",
                args.verbose,
            )
            
            arguments = [
                (
                    ProcessingArgs(
                        lines[(per_thread * j) : (per_thread * (j + 1))],
                        args.debug,
                        top_refs,
                        total_references,
                        strict_search_mode,
                        reference_taxa,
                    ),
                )
                for j in range(i, i + num_threads)
            ]

            with Pool(num_threads) as p:
                result = p.starmap(process_lines, arguments)
            del arguments
            for this_output, mkicks, this_log in result:
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

        nt_db.put("getall:valid_refs", ",".join(list(top_refs)))
        del top_refs

        requires_internal = {}
        internal_order = []
        for gene, hits in output.items():
            this_counter = Counter([i.header for i in hits]).most_common()
            requires_internal[gene] = {}
            if this_counter[0][1] > 1:
                this_hits = sum([i[1] for i in this_counter if i[1] > 1])
                this_common = {i[0] for i in this_counter if i[1] > 1}
                for hit in [i for i in hits if i.header in this_common]:
                    requires_internal[gene].setdefault(hit.header, []).append(hit)

                internal_order.append((gene, this_hits))

        internal_order.sort(key=lambda x: x[1], reverse=True)

        with Pool(args.processes) as pool:
            internal_results = pool.starmap(
                internal_filtering,
                [
                    (gene, requires_internal[gene], args.debug, args.internal_percent)
                    for gene, _ in internal_order
                ],
            )

        internal_result = {}
        for result in internal_results:
            internal_result.update(result)

        printv(
            f"Filtering done. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Writing to db",
            args.verbose,
        )

        head_to_seq = {}
        for content in out:
            lines = content.split("\n")
            for i in range(0, len(lines), 2):
                if lines[i] == "":
                    continue

                header = lines[i][1:]
                head_to_seq[header] = lines[i + 1]
        del out

        passes = 0
        internal_kicks = 0
        for gene, hits in output.items():
            dupe_divy_headers[gene] = {}
            kicks = {}
            out = []
            if gene in internal_result:
                this_kicks, this_log, kicks = internal_result[gene]

                internal_kicks += this_kicks
                if args.debug:
                    global_log.extend(this_log)

                out = []
                if kicks:
                    for hit in hits:
                        if hit.full_header not in kicks:
                            hit.seq = head_to_seq[hit.header.split("|")[0]]
                            out.append(hit.to_json())
                            dupe_divy_headers[gene][hit.header] = 1
            if not kicks:
                for hit in hits:
                    hit.seq = head_to_seq[hit.header.split("|")[0]]
                    out.append(hit.to_json())
                    dupe_divy_headers[gene][hit.header] = 1

            passes += len(out)
            db.put(f"gethits:{gene}", json.dumps(out))
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
        data = json.dumps(gene_dupe_count)
        nt_db.put(key, data)

        del db
        del nt_db

        printv(f"Took {time_keeper.differential():.2f}s overall.", args.verbose)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for input_path in args.INPUT:

        run_process(args, input_path)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr diamondPal"
    )
