import itertools
import os
from shutil import rmtree
import sys
from tempfile import TemporaryDirectory, NamedTemporaryFile
import json
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

    return None, None


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
        "kick"
    )

    def __init__(self, header, gene, reftaxon, score, qstart, qend, length, evalue=None):
        self.header = header
        self.gene = gene
        self.reftaxon = reftaxon
        self.seq = None
        self.score = float(score)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.evalue = float(evalue)
        self.length = int(length)
        self.kick = False

    def to_json(self):
        return {
            "header": self.header,
            "seq": self.seq,
            "ref_taxon": self.reftaxon,
            "ali_start": self.qstart,
            "ali_end": self.qend,
        }

    def __lt__(self, other):
        return self.score < other.score

    def __le__(self, other):
        return self.score <= other.score

    def __eq__(self, other):
        return self.score == other.score

    def __ne__(self, other):
        return self.score != other.score

    def __gt__(self, other):
        return self.score > other.score

    def __ge__(self, other):
        return self.score >= other.score


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


def get_sequence_results(fp, target_to_taxon, min_length, min_score, debug, e_min_size, e_min_evalue, internal_percent):
    header_hits = {}
    header_maps_to = {}
    this_header = None
    this_header_hits = []

    this_log = []
    intermediate_multi_kicks = 0
    intermediate_internal_kicks = 0
    intermediate_evalue_kicks = 0
    for line in fp:
        (
            raw_header,
            ref_header,
            frame,
            evalue,
            score,
            qstart,
            qend,
        ) = line.split("\t")
        qstart = int(qstart)
        qend = int(qend)
        gene, reftaxon = target_to_taxon[ref_header]

        if int(frame) < 0:
            length = qstart - qend + 1
            qend, qstart = qstart, qend
            header = raw_header + f"|[revcomp]:[translate({frame[1:]})]"
        else:
            length = qend - qstart + 1
            header = raw_header + f"|[translate({frame})]"

        if length < min_length:
            continue

        if float(score) < min_score:
            continue

        if not this_header:
            this_header = raw_header
        else:
            if this_header != raw_header:

                fail_e_check, log = hits_are_bad(this_header_hits, debug, e_min_size, e_min_evalue)
                if fail_e_check:
                    intermediate_evalue_kicks += len(this_header_hits)
                    if debug:
                        this_log.extend(log)
                    continue

                if sum(header_maps_to.values()) > 1:
                    this_header_hits, this_kicks, log = multi_filter(this_header_hits, debug)
                    if debug:
                        this_log.extend(log)
                    
                    intermediate_multi_kicks += this_kicks

                this_kicks, this_header_hits, log = internal_filter(this_header_hits, debug, internal_percent)
                if debug:
                    this_log.extend(log)
                
                intermediate_internal_kicks += this_kicks

                header_hits = {}
                for hit in this_header_hits:
                    header_hits.setdefault(hit.header, []).append(hit)

                yield header_hits, intermediate_multi_kicks, intermediate_internal_kicks, intermediate_evalue_kicks, this_log

                header_maps_to = {}
                this_header = raw_header
                this_log = []
                intermediate_multi_kicks = 0
                intermediate_internal_kicks = 0
                intermediate_evalue_kicks = 0

        header_maps_to[gene] = 1
        this_header_hits.append(Hit(header, gene, reftaxon, score, qstart, qend, length, evalue=evalue))
        


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

                    if percentage_of_overlap >= 0.5:  # min_overlap_multi:
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
                                        "Multi kicked out by",
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
                        if hit is not None
                    ]
                )
            return [], len(hits), log

    passes = [i for i in hits if i]
    return passes, len(hits) - len(passes), log

def internal_filter(hits, debug, interal_percent):
    log = []
    kicks = 0

    for hit_a, hit_b in itertools.combinations(hits, 2):
        if hit_a.kick or hit_b.kick:
            continue

        if hit_a.header == hit_b.header:
            continue

        if hit_a.gene != hit_b.gene:
            continue

        overlap_amount = get_overlap(hit_a.qstart, hit_a.qend, hit_b.qstart, hit_b.qend)
        percent = overlap_amount / hit_a.length
        if percent >= interal_percent:
            kicks += 1
            if debug:
                log.append(
                    (
                        hit_a.gene,
                        hit_a.header,
                        hit_a.reftaxon,
                        hit_a.score,
                        hit_a.qstart,
                        hit_a.qend,
                        "Internal kicked out by",
                        hit_b.gene,
                        hit_b.header,
                        hit_b.reftaxon,
                        hit_b.score,
                        hit_b.qstart,
                        hit_b.qend,
                    )
                )
            hit_b.kick = True

    return kicks, [i for i in hits if not i.kick], log

def hits_are_bad(
    hits: list, debug: bool, min_size, min_evalue
) -> bool:
    """
    Checks a list of hits for minimum size and score. If below both, kick.
    Returns True is checks fail, otherwise returns False.
    """
    evalue_log = []
    if len(hits) > min_size:
        return False, evalue_log
    for hit in hits:
        if hit.evalue < min_evalue:
            return False, evalue_log
    if debug:
        evalue_log = [
            (
                candidate.gene,
                candidate.header,
                candidate.reftaxon,
                candidate.score,
                candidate.qstart,
                candidate.qend,
                str(candidate.evalue),
                str(len(hits)),
                "Group kicked due to length and evalue",
            )
            for candidate in hits
        ]
    return True, evalue_log


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

    head_to_seq = {}
    for content in out:
        lines = content.split("\n")
        for i in range(0, len(lines), 2):
            if lines[i] == "":
                continue

            header = lines[i][1:]
            head_to_seq[header] = lines[i + 1]

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
            f"Diamond done. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s",
            args.verbose,
        )
    else:
        printv(
            f"Found existing Diamond output. Elapsed time {time_keeper.differential():.2f}s",
            args.verbose,
        )
    del out

    db = wrap_rocks.RocksDB(os.path.join(input_path, "rocksdb", "hits"))
    output = {}
    kicks = 0
    passes = 0
    evalue_kicks = 0
    global_log = []
    dupe_divy_headers = {}
    multi_kicks = 0
    internal_kicks = 0
    if args.skip_multi:
        printv("Skipping multi-filtering", args.verbose)
    with open(out_path) as fp:
        for header_dict, intermediate_multi_kicks, intermediate_internal_kicks, intermediate_evalue_kicks, this_log in get_sequence_results(
                fp,
                target_to_taxon,
                args.min_length,
                args.min_score,
                args.debug, 
                args.min_amount, 
                args.min_evalue, 
                args.internal_percent
            ):

            kicks += (intermediate_multi_kicks + intermediate_internal_kicks + intermediate_evalue_kicks)
            internal_kicks += intermediate_internal_kicks
            evalue_kicks += intermediate_evalue_kicks
            multi_kicks += intermediate_multi_kicks
            if args.debug:
                global_log.extend(this_log)

            for hits in header_dict.values():
                hits.sort(key=lambda x: x.score, reverse=True)
                passes += len(hits)
                best_hit = reciprocal_check(
                    hits, strict_search_mode, reference_taxa
                )

                if best_hit:
                    base_header = best_hit.header.split("|")[0]
                    dupe_divy_headers.setdefault(best_hit.gene, {})[base_header] = 1

                    best_hit.seq = head_to_seq[base_header]

                    output.setdefault(best_hit.gene, []).append(
                        best_hit.to_json()
                    )

    for gene, hits in output.items():
        db.put(f"gethits:{gene}", json.dumps(hits))

    if global_log:
        with open(os.path.join(input_path, "multi.log"), "w") as fp:
            fp.write(
                "\n".join([",".join([str(i) for i in line]) for line in global_log])
            )
    print(
        f"Took {time_keeper.lap():.2f}s for {kicks} kicks leaving {passes} results. Writing to DB"
    )
    print(f"{evalue_kicks} kicks were due to Evalue Filtering")
    print(f"{internal_kicks} kicks were due to Internal Filtering")
    print(f"{multi_kicks} kicks were due to Multi Filtering")

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
