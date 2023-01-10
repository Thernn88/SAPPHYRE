
import os
from math import ceil
from shutil import rmtree
import sys
from tempfile import TemporaryDirectory, NamedTemporaryFile
import json
from multiprocessing.pool import Pool
import wrap_rocks
import phymmr_tools

from . import rocky
from .utils import printv, gettempdir, parseFasta
from .timekeeper import TimeKeeper, KeeperMode

def reciprocal_check(hits, strict, taxa_present):
    first = None
    second = None
    hit_on_taxas = {i:0 for i in taxa_present}

    for hit in hits:
        if not first:
            first = hit
        elif not second:
            if first.reftaxon != hit.reftaxon and first.gene == hit.gene:
                second = hit

        hit_on_taxas[hit.reftaxon] = 1 # Just need 1 hit
        
        if strict:
            if all(hit_on_taxas.values()):
                return first, second
        else:
            if first and second:
                return first, second #Return early if we have 2 hits
    
    if any(hit_on_taxas.values()):
        return first, second

    return None, None

class Hit:
    __slots__ = (
        "header",
        "qstart",
        "qend",
        "gene",
        "score",
        "reftaxon",
        "full_seq",
        "trim_seq"
    )
    def __init__(self, header, gene, reftaxon, score, qstart, qend):
        self.header = header
        self.gene = gene
        self.reftaxon = reftaxon
        self.full_seq = None
        self.trim_seq = None
        self.score =score
        self.qstart = int(qstart)
        self.qend = int(qend)

    def to_first_json(self):
        return  {
                    "header": self.header,
                    "gene": self.gene,
                    "full_seq": self.full_seq,
                    "trim_seq": self.trim_seq,
                    "ref_taxon": self.reftaxon,
                    "ali_start": self.qstart,
                    "ali_end": self.qend,
                }

    def to_second_json(self):
        return  {
                    "trim_seq": self.trim_seq,
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

def get_sequence_results(fp, target_to_taxon, head_to_seq):
    header_hits = {}
    header_maps_to = {}
    this_header = None
    for line in fp:
        raw_header, ref_header, frame, score, qstart, qend, qlen = line.strip().split("\t")
        qstart = int(qstart)
        qend = int(qend)
        gene, reftaxon = target_to_taxon[ref_header]
        

        if int(frame) < 0:
            qstart = int(qlen) - qstart
            qend = int(qlen) - qend

            header = raw_header + f"|[revcomp]:[translate({abs(int(frame))})]"
        else:
            header = raw_header + f"|[translate({abs(int(frame))})]"

        if not this_header: 
            this_header = raw_header
        else:
            if this_header != raw_header:
                for hits in header_hits.values():
                    yield sorted(hits, key = lambda x: x.score, reverse=True), sum(list(header_maps_to.values()))    

                header_hits = {}
                header_maps_to = {}
                this_header = raw_header

        header_maps_to[gene] = 1
        
        header_hits.setdefault(header, []).append(Hit(header, gene, reftaxon, score, qstart, qend))
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
                    amount_of_overlap = get_overlap(master_env_start, master_env_end, candidate.qstart,
                                                    candidate.qend)
                    percentage_of_overlap = amount_of_overlap / distance

                    if percentage_of_overlap >= 0.5:#min_overlap_multi:
                        score_difference = get_difference(master.score, candidate.score)
                        if score_difference >= 1.05:
                            kick_happend = True
                            hits[i] = None
                            candidates[i-1] = None
                            if debug:
                                log.append((candidate.gene, candidate.header, candidate.reftaxon, candidate.score, candidate.qstart, candidate.qend, "Kicked out by", master.gene, master.header, master.reftaxon, master.score, master.qstart, master.qend))
                        else:
                            miniscule_score = True
                            break
        if miniscule_score:
            if debug:
                log.extend([(candidate.gene, candidate.header, candidate.reftaxon, candidate.score, candidate.qstart, candidate.qend, "Kicked due to miniscule score", master.gene, master.header, master.reftaxon, master.score, master.qstart, master.qend) for hit in hits if hit])
            return [], len(hits), log

    passes = [i for i in hits if i]
    return passes, len(hits) - len(passes), log

def run_process(args, input_path) -> None:
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input
    
    strict_search_mode = args.strict_search_mode
    sensitivity = args.sensitivity
    top_amount = args.top

    taxa = os.path.basename(input_path)
    printv(f'Processing: {taxa}', args.verbose, 0)
    printv("Grabbing reference data from Orthoset DB.", args.verbose)
    # make dirs
    diamond_path = os.path.join(input_path, "diamond")
    if args.overwrite:
        if os.path.exists(diamond_path):
            rmtree(diamond_path)
    os.makedirs(diamond_path, exist_ok=True)

    num_threads = args.processes

    orthoset_db_path = os.path.join(orthosets_dir, orthoset, "rocksdb")
    diamond_db_path = os.path.join(orthosets_dir, orthoset, "diamond", orthoset+'.dmnd')
    # Potato brain check
    if not os.path.exists(orthoset_db_path) or not os.path.exists(diamond_db_path):
        if input(f"Could not find orthoset DB at {orthoset_db_path}. Would you like to generate it? Y/N: ").lower() == "y":
            print("Attempting to generate DB")
            os.system(f"python3 -m phymmr -p {num_threads} Makeref {orthoset}.sqlite -s {orthoset}")
        else:
            print("Aborting")
            sys.exit(1)
    orthoset_db = wrap_rocks.RocksDB(orthoset_db_path)

    reference_taxa = json.loads(orthoset_db.get("getall:taxainset"))
    target_to_taxon = json.loads(orthoset_db.get("getall:targetreference"))

    del orthoset_db
    time_keeper.lap() #Reset timer

    printv("Done! Grabbing NT sequences.", args.verbose)

    nt_db_path = os.path.join(input_path, "rocksdb", "sequences", "nt")

    nt_db = wrap_rocks.RocksDB(nt_db_path)

    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
    else:
        raise ValueError("Nucleotide sequence not found in database. Did Prepare succesfully finish?")
    out = [nt_db.get(f"ntbatch:{i}") for i in recipe]

    dupe_counts = json.loads(nt_db.get("getall:dupes"))

    head_to_seq = {}
    for content in out:
        lines = content.split("\n")
        for i in range(0, len(lines), 2):
            if lines[i] == "":
                continue

            header = lines[i][1:]
            head_to_seq[header] = lines[i+1]

    out_path = os.path.join(diamond_path, f"{sensitivity}.tsv")
    if not os.path.exists(out_path) or os.stat(out_path).st_size == 0:
        with TemporaryDirectory(dir=gettempdir()) as dir, NamedTemporaryFile(dir=dir) as input_file:
            input_file.write("".join(out).encode())
            input_file.flush()
            
            quiet = "--quiet" if args.verbose == 0 else ""

            printv(f"Done! Running Diamond. Elapsed time {time_keeper.differential():.2f}s.", args.verbose)
            time_keeper.lap() #Reset timer
            os.system(f"diamond blastx -d {diamond_db_path} -q {input_file.name} -o {out_path} --{sensitivity}-sensitive --masking 0 -e {args.evalue} --outfmt 6 qseqid sseqid qframe bitscore qstart qend qlen {quiet} --top {top_amount} --max-hsps 0 -p {num_threads}")
            input_file.seek(0)

        printv(f"Diamond done. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s", args.verbose)
    else:
        printv(f"Found existing Diamond output. Elapsed time {time_keeper.differential():.2f}s", args.verbose)
    del out

    db = wrap_rocks.RocksDB(os.path.join(input_path, "rocksdb", "hits"))
    output = {}
    kicks = 0 
    passes = 0
    global_log = []
    dupe_divy_headers = {}
    if not args.multi:
        printv("Skipping multi-filtering", args.verbose)
    with open(out_path) as fp:
        for hits, requires_multi in get_sequence_results(fp, target_to_taxon, head_to_seq):
            if requires_multi and args.multi:
                hits, this_kicks, log = multi_filter(hits, args.debug)
                kicks += this_kicks
                if args.debug:
                    global_log.extend(log)

            passes += len(hits)
            first_hit, rerun_hit = reciprocal_check(hits, strict_search_mode, reference_taxa)

            if first_hit:
                base_header = first_hit.header.split("|")[0]
                dupe_divy_headers.setdefault(first_hit.gene, {})[base_header] = 1

                nuc_seq = head_to_seq[base_header] #Todo move later
                if "revcomp" in first_hit.header:
                    nuc_seq = phymmr_tools.bio_revcomp(nuc_seq)
                
                first_hit.full_seq = nuc_seq

                first_hit.trim_seq = nuc_seq[first_hit.qstart : first_hit.qend]
                if rerun_hit:
                    rerun_hit.trim_seq = nuc_seq[rerun_hit.qstart : rerun_hit.qend]

                rerun_hit = rerun_hit.to_second_json() if rerun_hit else None
                output.setdefault(first_hit.gene, []).append({"f":first_hit.to_first_json(), "s":rerun_hit})

    for gene, hits in output.items():
        db.put(f"gethits:{gene}", json.dumps(hits))

    if global_log:
        with open(os.path.join(input_path, "multi.log"), "w") as fp:
            fp.write("\n".join([",".join([str(i) for i in line]) for line in global_log]))

    print(f"Took {time_keeper.lap():.2f}s for {kicks} kicks leaving {passes} results. Writing to DB")
            
    gene_dupe_count = {}
    for gene, headers in dupe_divy_headers.items():
        for base_header in headers.keys():
            if base_header in dupe_counts:
                gene_dupe_count.setdefault(gene, {})[base_header] = dupe_counts[base_header]

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
