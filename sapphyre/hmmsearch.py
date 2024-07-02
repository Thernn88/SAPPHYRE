import math
import warnings
import subprocess
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
    node: int
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


def merge_clusters(clusters):
    # Sort the clusters by their starting point
    sorted_clusters = sorted(clusters, key=lambda x: x[0])
    merged_clusters = []
    
    current_start, current_end = sorted_clusters[0]
    
    for start, end in sorted_clusters[1:]:
        if get_overlap(current_start, current_end, start, end, 1) is not None:
            current_start = min(current_start, start)  # Handle left overlap
            current_end = max(current_end, end)        # Handle right overlap and containment
        else:
            merged_clusters.append((current_start, current_end))
            current_start, current_end = start, end
    
    merged_clusters.append((current_start, current_end))  # Add the last cluster
    
    return merged_clusters

def load_Sequences(source_seqs):
    if type(source_seqs) is dict:
        this_seqs = source_seqs
    else:
        this_seqs = {int(header): seq for header, seq in parseFasta(source_seqs)}

    return this_seqs

def shift_targets(is_full, query_template, nodes_in_gene, diamond_hits, cluster_full, fallback, hits_have_frames_already, unaligned_sequences, nt_sequences, parents, children, required_frames, this_seqs, bio_revcomp):
    if is_full:
        for node in nodes_in_gene:
            parent_seq = this_seqs.get(node)
            if parent_seq is None:
                # Full cluster search introduced header not found in db.
                # Quicker to filter after the fact 
                continue

            strands_present = cluster_full.get(node, {"+", "-"})

            unaligned_sequences.append((node, parent_seq))
            if "+" in strands_present:
                # Forward frame 1
                query = query_template.format(node, 1)
                nt_sequences[query] = parent_seq
                required_frames[node].update({1,2,3})
                if query not in parents and node not in cluster_full:
                    children[query] = fallback[node]

                # Forward 2 & 3
                for shift_by in [1, 2]:
                    shifted = shift(1, shift_by)
                    new_query = query_template.format(node, shifted)
                    nt_sequences[new_query] = parent_seq[shift_by:]
                    if new_query in parents:
                        continue
                    if node in cluster_full:
                        continue
                    children[new_query] = fallback[node]

            if "-" in strands_present:
                # Reversed frame 1
                bio_revcomp_seq = bio_revcomp(parent_seq)
                query = query_template.format(node, -1)
                nt_sequences[query] = bio_revcomp_seq
                required_frames[node].update({-1,-2,-3})
                if query not in parents and node not in cluster_full:
                    children[query] = fallback[node]

                # Reversed 2 & 3
                for shift_by in [1, 2]:
                    shifted = shift(-1, shift_by)
                    new_query = query_template.format(node, shifted)
                    nt_sequences[new_query] = bio_revcomp_seq[shift_by:]
                    if new_query in parents:
                        continue
                    if node in cluster_full:
                        continue
                    children[new_query] = fallback[node]
    else:
        for hit in diamond_hits:
            raw_sequence = hit.seq
            frame = hit.frame
            query = query_template.format(hit.node, frame)
            unaligned_sequences.append((query, raw_sequence))
            parents[query] = hit

            for shift_by in [1, 2]:
                shifted = shift(frame, shift_by)
                if not shifted in hits_have_frames_already[hit.node]:
                    new_query = query_template.format(hit.node, shifted)
                    unaligned_sequences.append((new_query, raw_sequence[shift_by:]))
                    hits_have_frames_already[hit.node].add(shifted)
                    children[new_query] = hit


def add_full_cluster_search(clusters, chomp_max_distance, nodes_in_genes, source_clusters, cluster_full, nodes_in_gene, cluster_dict):
    for cluster in clusters:
        for i in range(cluster[0][0] - chomp_max_distance, cluster[0][1] + chomp_max_distance + 1):
            header = i
            
            if header in nodes_in_gene:
                continue
            
            this_crange, strands_present = cluster_dict.get(header, (None, None))
            
            if this_crange:
                source_clusters[header] = this_crange
                cluster_full[header] = strands_present
                nodes_in_gene.add(header)


def get_results(hmm_output, map_mode):
    data = defaultdict(list)
    high_score = 0
    with open(hmm_output) as f:
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

    return queries


def add_new_result(map_mode, gene, query, results, is_full, cluster_full, cluster_queries, source_clusters, nt_sequences, parents, children, parents_done, passed_ids, output, hmm_log, hmm_log_template):
    node, frame = query.split("|")
    id = int(node)
    if id in cluster_full:
        if not query in parents_done:
            frame = int(frame)
            for result in results:
                start, end, score, ali_start, ali_end = result
                start = start * 3
                end = end * 3

                sequence = nt_sequences[query][start: end]

                new_qstart = start
                if frame < 0:
                    new_qstart -= (3 - abs(frame))
                else:
                    new_qstart += frame

                passed_ids.add(id)
                parents_done.add(query)

                cluster_range = source_clusters[id]
                where_from = "Cluster Full"
                cquery = cluster_queries[cluster_range]

                new_hit = HmmHit(node=id, score=score, frame=int(frame), evalue=0, qstart=new_qstart, qend=new_qstart + len(sequence), gene=gene, query=cquery, uid=None, refs=[], seq=sequence)
                hmm_log.append(hmm_log_template.format(new_hit.gene, new_hit.node, new_hit.frame, where_from))
                output.append(new_hit)
        return
    
    if query in parents or query in children:
        if query in parents:
            hit = parents[query]
            if f"{hit.node}|{hit.frame}" in parents_done:
                return
            
            frame = hit.frame
        else:
            hit = children[query]
            frame = int(query.split("|")[1])
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
                if frame < 0:
                    new_qstart -= (3 - abs(frame))

                else:
                    new_qstart += frame

            else:
                new_qstart = hit.qstart + start

            new_uid = hit.uid + start + end

            passed_ids.add(hit.node)
            parents_done.add(f"{hit.node}|{frame}")
            new_hit = HmmHit(node=hit.node, score=score, frame=frame, evalue=hit.evalue, qstart=new_qstart, qend=new_qstart + len(sequence), gene=hit.gene, query=hit.query, uid=new_uid, refs=hit.refs, seq=sequence)
            output.append(new_hit)


def hmm_search(batches, source_seqs, is_full, is_genome, hmm_output_folder, aln_ref_location, overwrite, map_mode, debug, verbose, evalue_threshold, chomp_max_distance):
    batch_result = []
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    this_seqs = load_Sequences(source_seqs)#, nodes_in_gene)
    for gene, diamond_hits in batches:
        diamond_hits = json.decode(diamond_hits, type=list[Hit])
        printv(f"Processing: {gene}", verbose, 2)
        aligned_sequences = []
        this_hmm_output = path.join(hmm_output_folder, f"{gene}.hmmout")

        
        hits_have_frames_already = defaultdict(set)
        diamond_ids = []
        for hit in diamond_hits:
            diamond_ids.append((hit.node, hit.frame))
            hits_have_frames_already[hit.node].add(hit.frame)

        if is_genome:
            diamond_ids.sort()
            clusters = []
            current_cluster = []

            strands_present = set()

            for child_index, frame in diamond_ids:
                if not current_cluster:
                    current_cluster.append(child_index)
                    current_index = child_index
                    strand = "+" if frame > 0 else "-"
                    strands_present.add(strand)
                else:
                    if child_index - current_index <= chomp_max_distance:
                        current_cluster.append(child_index)
                        current_index = child_index
                        strand = "+" if frame > 0 else "-"
                        strands_present.add(strand)
                    else:
                        if len(current_cluster) > 2:
                            clusters.append(((current_cluster[0], current_cluster[-1]), strands_present))
                        current_cluster = [child_index]
                        strand = "+" if frame > 0 else "-"
                        strands_present = {strand}
                        current_index = child_index
            
            if current_cluster:
                if len(current_cluster) > 2:
                    clusters.append(((current_cluster[0], current_cluster[-1]), strands_present))
                    
            cluster_dict = {}
            for cluster, strands_present in clusters:
                for i in range(cluster[0]-chomp_max_distance, cluster[1] + chomp_max_distance + 1):
                    cluster_dict[i] = (cluster, strands_present)
        
        nt_sequences = {}
        parents = {}

        cluster_full = {}
        cluster_queries = defaultdict(list)
        source_clusters = {}
        
        children = {}
        unaligned_sequences = []
        required_frames = defaultdict(set)
        query_template = "{}|{}"
        if is_full:
            fallback = {}
            nodes_in_gene = set()
            for hit in diamond_hits:
                nodes_in_gene.add(hit.node)
                query = query_template.format(hit.node, hit.frame)
                parents[query] = hit

                if hit.node not in fallback:
                    fallback[hit.node] = hit
                
                if not is_genome:
                    continue

                this_crange, _ = cluster_dict.get(hit.node, (None, None))
                if this_crange:
                    cluster_queries[this_crange].append(hit.query)
                
            if is_genome and clusters:
                cluster_queries = {k: max(set(v), key=v.count) for k, v in cluster_queries.items()}
                add_full_cluster_search(clusters, chomp_max_distance, nodes_in_gene, source_clusters, cluster_full, nodes_in_gene, cluster_dict)
                      
        shift_targets(is_full, query_template, nodes_in_gene, diamond_hits, cluster_full, fallback, hits_have_frames_already, unaligned_sequences, nt_sequences, parents, children, required_frames, this_seqs, bio_revcomp)

        aln_file = path.join(aln_ref_location, f"{gene}.aln.fa")
        output = []
        new_outs = []
        parents_done = set()
        passed_ids = set()
        hmm_log = []
        hmm_log_template = "{},{},{},{}"
        
        if debug > 2 or not path.exists(this_hmm_output) or stat(this_hmm_output).st_size == 0 or overwrite:
            if is_full:
                with NamedTemporaryFile(dir=gettempdir()) as unaligned_tmp, NamedTemporaryFile(dir=gettempdir()) as aln_tmp:
                    writeFasta(unaligned_tmp.name, unaligned_sequences)
                    ft = ["fastatranslate", unaligned_tmp.name]
                    subprocess.run(ft, stdout=aln_tmp)
                    #system(f"fastatranslate {unaligned_tmp.name} > {aln_tmp.name}")

                    for header, seq in parseFasta(aln_tmp.name, True):
                        frame = int(header[-3])
                        if "rev" in header:
                            frame = -frame
                        header = int(header.split(" ")[0])
                        if frame in required_frames[header]:
                            query = query_template.format(header, frame)
                            aligned_sequences.append((query, seq))
            else:
                for header, seq in unaligned_sequences:
                    nt_sequences[header] = seq
                    aligned_sequences.append((header, str(Seq(seq).translate())))
                    
                        
            with NamedTemporaryFile(dir=gettempdir()) as hmm_temp_file:
                hmmb = ["hmmbuild", hmm_temp_file.name, aln_file]
                subprocess.run(hmmb, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                #system(f"hmmbuild '{hmm_temp_file.name}' '{aln_file}' > /dev/null")

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
                        hmms = ["hmmsearch", "--nobias", "--domtblout", this_hmm_output, "--domT", "10.0", hmm_temp_file.name, aligned_files.name]
                        subprocess.run(hmms, stdout=subprocess.DEVNULL)
                        #system(
                        #f"hmmsearch --nobias --domtblout {this_hmm_output} --domT 10.0 {hmm_temp_file.name} {aligned_files.name} > /dev/null",
                        #)

        if debug > 2:
            continue#return "", [], [], [], []

        queries = get_results(this_hmm_output, map_mode)

        for query, results in queries:
            add_new_result(map_mode, gene, query, results, is_full, cluster_full, cluster_queries, source_clusters, nt_sequences, parents, children, parents_done, passed_ids, output, hmm_log, hmm_log_template)

        diamond_kicks = []
        for hit in diamond_hits:
            if not query_template.format(hit.node, hit.frame) in parents_done:
                if is_genome and math.floor(math.log10(abs(hit.evalue))) <= -evalue_threshold:
                    this_id = hit.node
                    has_neighbour = False
                    neighbour = None
                    for id in passed_ids:
                        if abs(this_id - id) <= chomp_max_distance:
                            neighbour = str(id)
                            has_neighbour = True
                            break

                    if has_neighbour:
                        printv(f"Rescued {hit.node}", verbose, 2)
                                     
                        if hit.frame < 0:      
                            new_start = hit.qend
                            new_start = len(this_seqs[hit.node]) - 2 - new_start # -1 for len inclusive and -1 for 0 based
                            new_end = new_start + len(hit.seq)
                        else:
                            new_start = hit.qstart
                            new_end = new_start + len(hit.seq)
                        
                        new_hit = HmmHit(node=hit.node, score=0, frame=hit.frame, evalue=hit.evalue, qstart=new_start, qend=new_end, gene=hit.gene, query=hit.query, uid=hit.uid, refs=hit.refs, seq=hit.seq)
                        output.append(new_hit)
                        parents_done.add(query_template.format(hit.node, hit.frame))
                        hmm_log.append(hmm_log_template.format(hit.gene, hit.node, hit.frame, f"Rescued by {neighbour}. Evalue: {hit.evalue}"))
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
        lines = nt_db.get_bytes(f"ntbatch:{i}").decode().splitlines()
        head_to_seq.update(
            {
                int(lines[i][1:]): lines[i + 1]
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

    aln_ref_location = path.join(input_folder, "top")
    if not path.exists(aln_ref_location):
        orthoset_path = path.join(args.orthoset_input, args.orthoset)
        aln_ref_location = None
        for ortho_path in ["final", "cleaned", "trimmed", "aln"]:
            if path.exists(path.join(orthoset_path, ortho_path)):
                aln_ref_location = ortho_path
                break
    if aln_ref_location is None:
        printv("ERROR: Could not find alignment reference", args.verbose, 0)
        return False

    per_batch = math.ceil(len(transcripts_mapped_to) / args.processes)
    temp_source_file = None
    if args.processes > 1:
        temp_source_file = NamedTemporaryFile(dir=gettempdir(), prefix="seqs_", suffix=".fa")
        writeFasta(temp_source_file.name, head_to_seq.items())
        seq_source = temp_source_file.name
    else:
        seq_source = head_to_seq
    batches = [(transcripts_mapped_to[i:i + per_batch], seq_source, is_full, is_genome, hmm_output_folder, aln_ref_location, args.overwrite, args.map, args.debug, args.verbose, args.evalue_threshold, args.chomp_max_distance) for i in range(0, len(transcripts_mapped_to), per_batch)]

    if args.processes <= 1:
        all_hits = []
        for this_arg in batches:
            all_hits.append(hmm_search(*this_arg))
    else:
        with Pool(args.processes) as p:
            all_hits = p.starmap(hmm_search, batches)

    printv("Running miniscule score filter", args.verbose, 1)

    log = ["Gene,Node,Frame"]
    klog = ["Gene,Node,Frame"]
    kick_log = ["Gene,Node,Frame"]
    multi_causing_log = ["Gene,Node,Frame,Evalue,Master Gene"]
    header_based_results = defaultdict(list)
    printv("Processing results", args.verbose, 1)
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
    for hits in header_based_results.values():
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
    printv("Writing results to db", args.verbose, 1)
    for gene, hits in gene_based_results.items():
        hits_db.put_bytes(f"gethmmhits:{gene}", json.encode(hits))

    del temp_source_file
    del hits_db

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
