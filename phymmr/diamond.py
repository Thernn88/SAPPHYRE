
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
        elif not second and first.reftaxon != hit.reftaxon:
            second = hit

        hit_on_taxas[hit.reftaxon] = 1 # Just need 1 hit
        
        if strict:
            if all(hit_on_taxas.values()):
                return [first, second]
        else:
            if any(hit_on_taxas.values()):
                return [first, second]
    
    return None

class Hit:
    __slots__ = (
        "header",
        "ref_header",
        "evalue",
        "score",
        "qstart",
        "qend",
        "sstart",
        "send",
        "length",
        "gene",
        "reftaxon",
        "full_seq",
        "trim_seq",
        "sstart",
        "send"
    )
    def __init__(self, header, ref_header, evalue, score, full_seq, trim_seq, gene, reftaxon, qstart, qend, sstart, send):
        self.header = header
        self.ref_header = ref_header
        self.evalue = evalue
        self.score = float(score)
        self.gene = gene
        self.reftaxon = reftaxon
        self.full_seq = full_seq
        self.trim_seq = trim_seq
        self.qstart = int(qstart)
        self.qend = int(qend)

        self.sstart = int(sstart)
        self.send = int(send)


    def to_json(self):
        return  {
                    "header": self.header,
                    "gene": self.gene,
                    "full_seq": self.full_seq,
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
    try:
        if scoreA / scoreB > 1:
            return scoreA / scoreB
        if scoreB / scoreA > 1:
            return scoreB / scoreA
        if scoreA == scoreB:
            return 1
    except ZeroDivisionError:
        return 0

def multi_filter(hits, debug):
    kick_happend = True

    log = []

    while kick_happend:
        kick_happend = False
        master = hits[0]
        candidates = hits[1:]

        master_env_start = master.sstart
        master_env_end = master.send

        miniscule_score = False
        for i, candidate in enumerate(candidates, 1):
            if candidate:
                if master.gene != candidate.gene:
                    distance = (master_env_end - master_env_start) + 1  # Inclusive
                    amount_of_overlap = get_overlap(master_env_start, master_env_end, candidate.sstart,
                                                    candidate.send)

                    percentage_of_overlap = amount_of_overlap / distance

                    if percentage_of_overlap >= 0.5:#min_overlap_multi:
                        score_difference = get_difference(master.score, candidate.score)
                        if score_difference >= 1.05:
                            kick_happend = True
                            hits[i] = None
                            candidates[i-1] = None
                            if debug:
                                log.append((candidate.gene, candidate.header, candidate.reftaxon, candidate.score, "Kicked out by", master.gene, master.header, master.reftaxon, master.score))
                        else:
                            miniscule_score = True
                            break
        if miniscule_score:
            if debug:
                log.extend([(candidate.gene, candidate.header, candidate.reftaxon, candidate.score, "Kicked due to miniscule score", master.gene, master.header, master.reftaxon, master.score) for hit in hits if hit])
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
            os.system(f"diamond blastx -d {diamond_db_path} -q {input_file.name} -o {out_path} --{sensitivity}-sensitive --masking 0 --outfmt 6 qseqid sseqid qframe evalue bitscore qstart qend sstart send length qlen {quiet} --top {top_amount} --max-hsps 0 -p {num_threads}")
            input_file.seek(0)

        printv(f"Diamond done. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s", args.verbose)
    else:
        printv(f"Found existing Diamond output. Elapsed time {time_keeper.differential():.2f}s", args.verbose)
    del out

    this_header_out = {}
    header_maps_to = {}
    with open(out_path) as fp:
        for line in fp:
            raw_header, ref_header, frame, evalue, score, qstart, qend, sstart, send, length, qlen = line.strip().split("\t")
            qstart = int(qstart)
            qend = int(qend)
            gene, reftaxon = target_to_taxon[ref_header]
            nuc_seq = head_to_seq[raw_header] #Todo move later
            

            if int(frame) < 0:
                qstart = int(qlen) - qstart
                qend = int(qlen) - qend

                nuc_seq = phymmr_tools.bio_revcomp(nuc_seq)
                
                header = raw_header + f"|[revcomp]:[translate({abs(int(frame))})]"
            else:
                header = raw_header + f"|[translate({abs(int(frame))})]"

            header_maps_to.setdefault(header, []).append(gene)

            trimmed_sequence = nuc_seq[qstart : qend]

            this_header_out.setdefault(header, []).append(Hit(header, ref_header, evalue, score, nuc_seq, trimmed_sequence, gene, reftaxon, qstart, qend, sstart, send))

    for header in this_header_out:
        this_header_out[header] = sorted(this_header_out[header], key = lambda x: x.score, reverse=True)

    #Multi filter
    time_keeper.lap()
    if not args.multi:
        print("Skipping multi")
        requires_multi = []
    else:
        requires_multi = [header for header, genes in header_maps_to.items() if len(genes) > 1]
    kicks = 0 
    passes = 0
    global_log = []
    for header in requires_multi:
        results, this_kicks, log = multi_filter(this_header_out[header], args.debug)
        if args.debug:
            global_log.extend(log)
        this_header_out[header] = results
        passes += len(results)
        kicks += this_kicks
    
    if global_log:
        with open(os.path.join(input_path, "multi.log"), "w") as fp:
            fp.write("\n".join([",".join([str(i) for i in line]) for line in global_log]))

    print(f"Took {time_keeper.lap():.2f} for {kicks} kicks leaving {passes} results")
    #Internal filter
    #Reciprocal filter
    printv(f"Read diamond output. Took {time_keeper.lap():.2f}s. Elapsed time {time_keeper.differential():.2f}s. Doing Reciprocal Check", args.verbose)
    arguments = [(data, strict_search_mode, reference_taxa,) for data in this_header_out.values()]

    with Pool(num_threads) as pool:
        results = pool.starmap(reciprocal_check, arguments)

    final_gene_out = {}
    for result in results:
        if result:
            final_gene_out.setdefault(result[0].gene, []).append(result)

    printv(f"Reciprocal done. Took {time_keeper.lap():.2f}s. Doing output.", args.verbose)
            
    db = wrap_rocks.RocksDB(os.path.join(input_path, "rocksdb", "hits"))

    present_genes = []
    gene_dupe_count = {}
    for gene, hits in final_gene_out.items():
        output = []
        for first_hit, rerun_hit in hits:
            base_header = first_hit.header.split("|")[0]
            if base_header in dupe_counts:
                gene_dupe_count.setdefault(gene, {})[base_header] = dupe_counts[base_header]
            

            rerun_hit = rerun_hit.to_json() if rerun_hit else None
            output.append({"f":first_hit.to_json(), "s":rerun_hit})
        db.put(f"gethits:{gene}", json.dumps(output))
        present_genes.append(gene)

    db.put("getall:presentgenes", ",".join(present_genes))

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
