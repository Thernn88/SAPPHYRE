from decimal import Decimal
import math
from time import time
import warnings
from collections import defaultdict
from shutil import rmtree
from tempfile import NamedTemporaryFile
from os import path, system, stat
from multiprocessing import Pool
from wrap_rocks import RocksDB
from msgspec import Struct, json
from sapphyre_tools import bio_revcomp, get_overlap
from Bio import BiopythonWarning
from Bio.Seq import Seq
from .utils import parseFasta, printv, gettempdir, writeFasta
from .diamond import ReferenceHit, ReporterHit as Hit
from .timekeeper import TimeKeeper, KeeperMode


class HmmHit(Struct):
    node: str
    score: float
    frame: int
    evalue: float|None
    qstart: int
    qend: int
    gene: str
    query: str
    uid: int|None
    refs: list[ReferenceHit]
    seq: str = None
    

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


def get_diamondhits(
    rocks_hits_db: RocksDB,
) -> dict[str, list[Hit]]:
    """Returns a dictionary of gene to corresponding hits.

    Args:
    ----
        rocks_hits_db (RocksDB): RocksDB instance
    Returns:
        dict[str, list[Hit]]: Dictionary of gene to corresponding hits
    """
    present_genes = rocks_hits_db.get("getall:presentgenes")
    if not present_genes:
        printv("ERROR: No genes found in hits database", 0)
        printv("Please make sure Diamond completed successfully", 0)
        return None
    genes_to_process = present_genes.split(",")

    gene_based_results = []
    for gene in genes_to_process:
        gene_result = rocks_hits_db.get_bytes(f"gethits:{gene}")
        if not gene_result:
            printv(
                f"WARNING: No hits found for {gene}. If you are using a gene list file this may be a non-issue",
                0,
            )
            continue
        gene_based_results.append((gene, gene_result))

    return genes_to_process, gene_based_results


def shift(frame, by):
    coeff = 1
    if frame <= 0:
        coeff = -1

    frame = abs(frame) + by

    if frame > 3:
        frame = frame - 3
    if frame < 1:
        frame = frame + 3

    return frame * coeff


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
    if overlap_coords is None:
        return 0

    return overlap_coords[1] - overlap_coords[0]


def internal_filter_gene(this_gene_hits, debug, min_overlap_internal=0.9, score_diff_internal=1.5):
    this_gene_hits.sort(key=lambda hit: hit[0].score, reverse=True)
    filtered_sequences_log = []

    internal_log_template = "{},{},{},{},{},{},Internal Overlapped with Highest Score,{},{},{},{},{},{}"

    for i, (hit_a, start_a, end_a) in enumerate(this_gene_hits):
        if not hit_a:
            continue
        for j in range(len(this_gene_hits) - 1, i, -1):
            hit_b, start_b, end_b = this_gene_hits[j]
            if hit_b:
                if hit_a.node != hit_b.node:
                    if ((hit_a.score / hit_b.score) if hit_b.score != 0 else 0) < score_diff_internal:
                        break

                    amount_of_overlap = get_overlap_amount(start_a, end_a, start_b, end_b)
                    distance = (end_b - start_b) + 1  # Inclusive
                    percentage_of_overlap = amount_of_overlap / distance

                    if percentage_of_overlap >= min_overlap_internal:
                        this_gene_hits[j] = None, None, None
                        if debug:
                            filtered_sequences_log.append(
                                internal_log_template.format(
                                    hit_b.gene,
                                    hit_b.node,
                                    hit_b.frame,
                                    hit_b.score,
                                    start_b,
                                    end_b,
                                    hit_a.gene,
                                    hit_a.node,
                                    hit_a.frame,
                                    hit_a.score,
                                    start_a,
                                    end_a,
                                )
                            )


    return [i[0] for i in this_gene_hits if i[0] is not None], filtered_sequences_log


def hmm_search(batches, this_seqs, is_full, hmm_output_folder, top_location, overwrite, map_mode, debug, verbose, evalue_threshold, chomp_max_distance):
    batch_result = []
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    for gene, diamond_hits in batches:
        diamond_hits = json.decode(diamond_hits, type=list[Hit])
        printv(f"Processing: {gene}", verbose, 2)
        aligned_sequences = []
        this_hmm_output = path.join(hmm_output_folder, f"{gene}.hmmout")

        get_id = lambda header: int(header.split("_")[1])
        
        hits_have_frames_already = defaultdict(set)
        diamond_ids = list()
        for hit in diamond_hits:
            diamond_ids.append(get_id(hit.node))
            hits_have_frames_already[hit.node].add(hit.frame)


        diamond_ids.sort()
        clusters = []
        current_cluster = []
        for child_index in diamond_ids:
            if not current_cluster:
                current_cluster.append(child_index)
                current_index = child_index
            else:
                if child_index - current_index <= chomp_max_distance:
                    current_cluster.append(child_index)
                    current_index = child_index
                else:
                    if len(current_cluster) > 2:
                        clusters.append((current_cluster[0], current_cluster[-1]))
                    current_cluster = [child_index]
                    current_index = child_index
        
        if current_cluster:
            if len(current_cluster) > 2:
                clusters.append((current_cluster[0], current_cluster[-1]))
                
        cluster_dict = {}
        for cluster in clusters:
            for i in range(cluster[0]-chomp_max_distance, cluster[1] + chomp_max_distance + 1):
                cluster_dict[i] = cluster
        
        nt_sequences = {}
        parents = {}

        cluster_full = set()
        cluster_queries = defaultdict(list)

        children = {}
        unaligned_sequences = []
        if is_full:
            fallback = {}
            nodes_in_gene = set()
            for hit in diamond_hits:
                id = get_id(hit.node)


                this_crange = cluster_dict.get(id)
                if this_crange:
                    cluster_queries[this_crange].append(hit.query)

                nodes_in_gene.add(hit.node)
                query = f"{hit.node}|{hit.frame}"
                parents[query] = hit

                if hit.node not in fallback:
                    fallback[hit.node] = hit

            #grab most occuring query
            if clusters:
                cluster_queries = {k: max(set(v), key=v.count) for k, v in cluster_queries.items()}
                smallest_cluster_in_range = min([i[0] for i in clusters]) - chomp_max_distance
                largest_cluster_in_range = max([i[1] for i in clusters]) + chomp_max_distance


                source_clusters = {}
                for header, seq in this_seqs.items():
                    if header in nodes_in_gene:
                        continue
                    id = get_id(header)
                    if id < smallest_cluster_in_range or id > largest_cluster_in_range:
                        continue
                    
                    this_crange = cluster_dict.get(id)
                    
                    if this_crange:
                        source_clusters[header] = this_crange
                        cluster_full.add(header)
                        nodes_in_gene.add(header)
                
            for node in nodes_in_gene:
                parent_seq = this_seqs[node]

                # Forward frame 1
                query = f"{node}|1"
                nt_sequences[query] = parent_seq
                unaligned_sequences.append((query, parent_seq))

                # Forward 2 & 3
                for shift_by in [1, 2]:
                    shifted = shift(1, shift_by)
                    new_query = f"{node}|{shifted}"
                    nt_sequences[new_query] = parent_seq[shift_by:]
                    unaligned_sequences.append((new_query, parent_seq[shift_by:]))
                    if new_query in parents:
                        continue
                    if node in cluster_full:
                        continue
                    children[new_query] = fallback[node]

                # Reversed frame 1
                bio_revcomp_seq = bio_revcomp(parent_seq)
                query = f"{node}|-1"
                nt_sequences[query] = bio_revcomp_seq
                unaligned_sequences.append((query, bio_revcomp_seq))

                # Reversed 2 & 3
                for shift_by in [1, 2]:
                    shifted = shift(-1, shift_by)
                    new_query = f"{node}|{shifted}"
                    nt_sequences[new_query] = bio_revcomp_seq[shift_by:]
                    unaligned_sequences.append((new_query, bio_revcomp_seq[shift_by:]))
                    if new_query in parents:
                        continue
                    if node in cluster_full:
                        continue
                    children[new_query] = fallback[node]


        else:
            for hit in diamond_hits:
                raw_sequence = hit.seq
                frame = hit.frame
                query = f"{hit.node}|{frame}"
                unaligned_sequences.append((query, raw_sequence))
                parents[query] = hit

                for shift_by in [1, 2]:
                    shifted = shift(frame, shift_by)
                    if not shifted in hits_have_frames_already[hit.node]:
                        new_query = f"{hit.node}|{shifted}"
                        unaligned_sequences.append((new_query, raw_sequence[shift_by:]))
                        hits_have_frames_already[hit.node].add(shifted)
                        children[new_query] = hit

        for header, seq in unaligned_sequences:
            nt_sequences[header] = seq
            

        top_file = path.join(top_location, f"{gene}.aln.fa")
        output = []
        new_outs = []
        parents_done = set()
        passed_ids = set()
        hmm_log = []
        hmm_log_template = "{},{},{},{}"
        seq_template = ">{}\n{}\n"

        if debug > 2 or not path.exists(this_hmm_output) or stat(this_hmm_output).st_size == 0 or overwrite:
            for header, seq in unaligned_sequences:
                translate = str(Seq(seq).translate())
                aligned_sequences.append((header, translate))
                
            with NamedTemporaryFile(dir=gettempdir()) as hmm_temp_file:
                system(f"hmmbuild '{hmm_temp_file.name}' '{top_file}' > /dev/null")

                if debug > 2:
                    this_hmm_in = path.join(hmm_output_folder, f"{gene}_input.fa")
                    writeFasta(this_hmm_in, aligned_sequences)
                    system(
                    f"hmmsearch -o {this_hmm_output} --nobias --domT 10.0 {hmm_temp_file.name} {this_hmm_in} > /dev/null",
                    )
                else:
                    with NamedTemporaryFile(dir=gettempdir()) as aligned_files:
                        writeFasta(aligned_files.name, aligned_sequences)
                        aligned_files.flush()
                        system(
                        f"hmmsearch --nobias --domtblout {this_hmm_output} --domT 10.0 {hmm_temp_file.name} {aligned_files.name} > /dev/null",
                        )

        if debug > 2:
            continue#return "", [], [], [], []

        data = defaultdict(list)
        high_score = 0
        with open(this_hmm_output) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                while "  " in line:
                    line = line.replace("  ", " ")
                line = line.strip().split()

                query = line[0]

                start, end, ali_start, ali_end, score = int(line[17]), int(line[18]), int(line[15]), int(line[16]), float(line[13])

                high_score = max(high_score, score)

                data[query].append((start - 1, end, score, ali_start, ali_end))

        if map_mode:
            score_thresh = high_score * 0.9

            queries = []
            for query, results in data.items():
                queries.append((query, [i for i in results if i[2] >= score_thresh]))
        else:
            queries = data.items()

        for query, results in queries:
            node, frame = query.split("|")
            if node in cluster_full:
                if not query in parents_done:
                    for result in results:
                        start, end, score, ali_start, ali_end = result
                        start = start * 3
                        end = end * 3

                        sequence = nt_sequences[query][start: end]

                        new_qstart = start

                        passed_ids.add(get_id(node))
                        parents_done.add(query)

                        cquery = cluster_queries[source_clusters[node]]

                        new_hit = HmmHit(node=node, score=score, frame=int(frame), evalue=0, qstart=new_qstart, qend=new_qstart + len(sequence), gene=gene, query=cquery, uid=None, refs=[], seq=sequence)
                        hmm_log.append(hmm_log_template.format(new_hit.gene, new_hit.node, new_hit.frame, "Found in Cluster Full"))
                        output.append(new_hit)
                continue
            
            if query in parents:
                hit = parents[query]
                if not f"{hit.node}|{hit.frame}" in parents_done:
                    for result in results:
                        start, end, score, ali_start, ali_end = result
                        start = start * 3
                        end = end * 3

                        if map_mode:
                            sequence = nt_sequences[query]
                            if len(sequence) % 3 != 0:
                                sequence += ("N" * (3 - len(sequence) % 3))
                            start = 0
                        else:
                            sequence = nt_sequences[query][start: end]

                        if is_full:
                            new_qstart = start
                        else:
                            new_qstart = hit.qstart + start

                        new_uid = hit.uid + start + end

                        passed_ids.add(get_id(hit.node))
                        parents_done.add(f"{hit.node}|{hit.frame}")
                        new_hit = HmmHit(node=hit.node, score=score, frame=hit.frame, evalue=hit.evalue, qstart=new_qstart, qend=new_qstart + len(sequence), gene=hit.gene, query=hit.query, uid=new_uid, refs=hit.refs, seq=sequence)
                        output.append(new_hit)

            if query in children:
                _, frame = query.split("|")
                parent = children[query]
                for result in results:
                    start, end, score, ali_start, ali_end = result
                    start = start * 3
                    end = end * 3

                    if map_mode:
                        sequence = nt_sequences[query]
                        if len(sequence) % 3 != 0:
                            sequence += ("N" * (3 - len(sequence) % 3))
                        start = 0
                    else:
                        sequence = nt_sequences[query][start: end]

                    if is_full:
                        new_qstart = start
                    else:
                        new_qstart = hit.qstart + start

                    new_uid = parent.uid + start + end

                    passed_ids.add(get_id(hit.node))
                    clone = HmmHit(node=parent.node, score=score, frame=int(frame), evalue=parent.evalue, qstart=new_qstart, qend=new_qstart + len(sequence), gene=parent.gene, query=parent.query, uid=new_uid, refs=parent.refs, seq=sequence)
                    new_outs.append((f"{clone.gene},{clone.node},{clone.frame}"))
                    
                    output.append(clone)

        diamond_kicks = []
        for hit in diamond_hits:
            if not f"{hit.node}|{hit.frame}" in parents_done:
                if math.floor(math.log10(abs(hit.evalue))) <= -evalue_threshold:
                    this_id = get_id(hit.node)
                    has_neighbour = False
                    neighbour = None
                    for id in passed_ids:
                        if abs(this_id - id) <= chomp_max_distance:
                            neighbour = str(id)
                            has_neighbour = True
                            break

                    if has_neighbour:
                        printv(f"Rescued {hit.node}", verbose, 2)
                        new_hit = HmmHit(node=hit.node, score=0, frame=hit.frame, evalue=hit.evalue, qstart=hit.qstart, qend=hit.qend, gene=hit.gene, query=hit.query, uid=hit.uid, refs=hit.refs, seq=hit.seq)
                        output.append(new_hit)
                        parents_done.add(f"{hit.node}|{hit.frame}")
                        hmm_log.append(hmm_log_template.format(hit.gene, hit.node, hit.frame, f"Rescued by NODE_{neighbour}. Evalue: {hit.evalue}"))
                        continue

                diamond_kicks.append(hmm_log_template.format(hit.gene, hit.node, hit.frame, f"Kicked. Evalue: {hit.evalue}"))
                
        batch_result.append((gene, output, new_outs, hmm_log, diamond_kicks))
    return batch_result

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


def miniscule_multi_filter(hits, debug):
    hits.sort(key=lambda hit: hit.score, reverse=True)
    
    MULTI_PERCENTAGE_OF_OVERLAP = 0.3
    MULTI_SCORE_DIFFERENCE = 1.1

    edge_log = []
    log = []
    kicked_indices = set()

    # Assign the highest scoring hit as the 'master'
    master = hits[0]

    # Assign all the lower scoring hits as 'candidates'
    candidates = hits[1:]

    internal_template = "{},{},{},{},{},Kicked due to miniscule score,{},{},{},{},{}"
    internal_template_kick = "{},{},{},{},{},Kicked due to not being in highest scoring gene: {},{},{},{},{},{}"

    for i, candidate in enumerate(candidates, 1):
        # Skip if the candidate has been kicked already or if master and candidate are internal hits:
        if i in kicked_indices or master.gene == candidate.gene:
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
            if master.score and candidate.score:
                score_difference = get_score_difference(master.score, candidate.score)
                if score_difference < MULTI_SCORE_DIFFERENCE:
                    # If the score difference is not greater than 10% trigger miniscule score.
                    # Miniscule score means hits map to loosely to multiple genes thus
                    # we can't determine a viable hit on a single gene and must kick all.
                    if debug:
                        log = [internal_template.format(hit.gene, hit.node, hit.score, hit.qstart, hit.qend, master.gene, master.node, master.score, master.qstart, master.qend) for hit in hits if hit]

                    return [], len(hits), log, edge_log
            else:
                edge_log.append(f"{candidate.gene},{candidate.node},{candidate.frame},{candidate.evalue},{master.gene}")

        # If the score difference is great enough
        # kick the candidate and log the kick.
        log.append(internal_template_kick.format(candidate.gene, candidate.node, candidate.score, candidate.qstart, candidate.qend, master.gene, master.gene, master.node, master.score, master.qstart, master.qend))
        kicked_indices.add(i)

    return [hit for i, hit in enumerate(hits) if i not in kicked_indices], len(list(kicked_indices)), log, edge_log


def do_folder(input_folder, args):
    printv(f"Processing {input_folder}", args.verbose, 1)
    hits_db = RocksDB(path.join(input_folder, "rocksdb", "hits"))
    tk = TimeKeeper(KeeperMode.DIRECT)
    diamond_genes, transcripts_mapped_to = get_diamondhits(
        hits_db
    )

    head_to_seq = {}
    seq_db = RocksDB(path.join(input_folder, "rocksdb", "sequences", "nt"))
    is_genome = seq_db.get("get:isgenome")
    is_genome = is_genome == "True"
    is_assembly = seq_db.get("get:isassembly")
    is_assembly = is_assembly == "True"
    is_full = is_genome or is_assembly or args.full
    if is_full:
        recipe = seq_db.get("getall:batches").split(",")
        head_to_seq = get_head_to_seq(seq_db, recipe)
    del seq_db

    hmm_output_folder = path.join(input_folder, "hmmsearch")
    
    if (args.debug > 1 or args.overwrite) and path.exists(hmm_output_folder):
        rmtree(hmm_output_folder, ignore_errors=True)
    if not path.exists(hmm_output_folder):
        system(f"mkdir {hmm_output_folder}")

    top_location = path.join(input_folder, "top")

    per_batch = math.ceil(len(transcripts_mapped_to) / args.processes)
    batches = [(transcripts_mapped_to[i:i + per_batch], head_to_seq, is_full, hmm_output_folder, top_location, args.overwrite, args.map, args.debug, args.verbose, args.evalue_threshold, args.chomp_max_distance) for i in range(0, len(transcripts_mapped_to), per_batch)]

    if args.processes <= 1:
        for this_arg in batches:
            all_hits.append(hmm_search(*this_arg))
    else:
        with Pool(args.processes) as p:
            all_hits = p.starmap(hmm_search, batches)

    printv(f"Running miniscule score filter", args.verbose, 1)

    log = ["Gene,Node,Frame"]
    klog = ["Gene,Node,Frame"]
    kick_log = ["Gene,Node,Frame"]
    multi_causing_log = ["Gene,Node,Frame,Evalue,Master Gene"]
    header_based_results = defaultdict(list)
    printv(f"Processing results", args.verbose, 1)
    for batch in all_hits:
        for gene, hits, logs, klogs, dkicks in batch:
            if not gene:
                continue

            for hit in hits:
                header_based_results[hit.node].append(hit)
            
            log.extend(logs)
            klog.extend(klogs)
            kick_log.extend(dkicks)

    mlog = ["Gene,Node,Score,Start,End,Reason,Master Gene,Header,Score,Start,End"]
    mkicks = 0

    multi_filter_args = []

    gene_based_results = {gene: [] for gene in diamond_genes}
    for header, hits in header_based_results.items():
        genes_present = {hit.gene for hit in hits}
        if len(genes_present) > 1:
            multi_filter_args.append((hits, args.debug))
        else:
            for hit in hits:
                gene_based_results[hit.gene].append(hit)
            
    if args.processes <= 1:
        result = []
        for hits, debug in multi_filter_args:
            result.append(miniscule_multi_filter(hits, debug))
    else:
        with Pool(args.processes) as p:
            result = p.starmap(miniscule_multi_filter, multi_filter_args)
    
    for hits, kicked_count, this_log, edge_log in result:
        mlog.extend(this_log)
        multi_causing_log.extend(edge_log)
        if kicked_count:
            mkicks += kicked_count

        for hit in hits:
            gene_based_results[hit.gene].append(hit)
        
    printv(f"Kicked {mkicks} hits due to miniscule score", args.verbose, 1)
    printv(f"Writing results to db", args.verbose, 1)
    get_id = lambda header: int(header.split("_")[1])
    hmmsearch_cluster_log = []
    for gene, hits in gene_based_results.items():
        ids = [get_id(hit.node) for hit in hits]
        ids.sort()

        clusters = []
        current_cluster = []
        for child_index in ids:
            if not current_cluster:
                current_cluster.append(child_index)
                current_index = child_index
            else:
                if child_index - current_index <= args.chomp_max_distance:
                    current_cluster.append(child_index)
                    current_index = child_index
                else:
                    if len(current_cluster) > 2:
                        clusters.append((current_cluster[0], current_cluster[-1]))
                    current_cluster = [child_index]
                    current_index = child_index

        if current_cluster:
            if len(current_cluster) > 2:
                clusters.append((current_cluster[0], current_cluster[-1]))

        clusters.sort(key=lambda x: x[0])

        cluster_string = ", ".join([f"{cluster[0]}-{cluster[1]}" for cluster in clusters])

        hmmsearch_cluster_log.append((len(clusters), f"{gene},{len(hits)},{len(clusters)},{cluster_string}"))

        hits_db.put_bytes(f"gethmmhits:{gene}", json.encode(hits))

    del hits_db

    hmmsearch_cluster_log.sort(key=lambda x: x[0], reverse=True)
    with open(path.join(input_folder, "hmm_cluster.csv"), "w") as f:
        f.write("Gene,Sequence Count,Cluster Count,Clusters\n")
        f.write("\n".join([i[1] for i in hmmsearch_cluster_log]))
    with open(path.join(input_folder, "hmmsearch_log.log"), "w") as f:
        f.write("\n".join(klog))
    with open(path.join(input_folder, "hmmsearch_kick.log"), "w") as f:
        f.write("\n".join(kick_log))
    with open(path.join(input_folder, "hmmsearch_multi.log"), "w") as f:
        f.write("\n".join(mlog))
    with open(path.join(input_folder, "hmmsearch_zero_multis.log"), "w") as f:
        f.write("\n".join(multi_causing_log))

    printv(f"Done with {input_folder}. Took {tk.lap():.2f}s", args.verbose, 1)
    return True

def main(args):
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exist.", args.verbose, 0)
        return False
    result = []
    for input_path in args.INPUT:
        result.append(
            do_folder(
                input_path,
                args,
            ),
        )
    
    return all(result)
