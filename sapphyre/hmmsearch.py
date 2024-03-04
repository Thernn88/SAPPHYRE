import warnings
from Bio import BiopythonWarning
from collections import defaultdict
from shutil import rmtree
from tempfile import NamedTemporaryFile
from .diamond import ReferenceHit, ReporterHit as Hit
from wrap_rocks import RocksDB
from .utils import printv, gettempdir, parseFasta, writeFasta
from os import path, system, stat
from msgspec import Struct, json
from multiprocessing import Pool
from sapphyre_tools import translate, bio_revcomp, get_overlap
from Bio.Seq import Seq
from .timekeeper import TimeKeeper, KeeperMode


class HmmHit(Struct):
    node: str
    score: float
    frame: int
    qstart: int
    qend: int
    gene: str
    query: str
    uid: int
    refs: list[ReferenceHit]
    seq: str = None
    


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
        this_hits = json.decode(gene_result, type=list[Hit])
        if not gene_result:
            printv(
                f"WARNING: No hits found for {gene}. If you are using a gene list file this may be a non-issue",
                0,
            )
            continue
        gene_based_results.append((gene, this_hits))

    return gene_based_results


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
    if overlap_coords == None:
        return 0

    return overlap_coords[1] - overlap_coords[0]


def internal_filter_gene(this_gene_hits, debug, min_overlap_internal=0.9, score_diff_internal=1.5):
    this_gene_hits.sort(key=lambda hit: hit[0].score, reverse=True)
    filtered_sequences_log = []

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
                                f"{hit_b.gene},{hit_b.node},{hit_b.frame},{hit_b.score},{start_b},{end_b},Internal Overlapped with Highest Score,{hit_a.gene},{hit_a.node},{hit_a.frame},{hit_a.score},{start_a},{end_a}"
                            )


    return [i[0] for i in this_gene_hits if i[0] is not None], filtered_sequences_log


def hmm_search(gene, diamond_hits, this_seqs, is_full, hmm_output_folder, top_location, overwrite, debug, verbose):
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    printv(f"Processing: {gene}", verbose, 2)
    aligned_sequences = []
    this_hmm_output = path.join(hmm_output_folder, f"{gene}.hmmout")
    
    hits_have_frames_already = defaultdict(set)
    for hit in diamond_hits:
        hits_have_frames_already[hit.node].add(hit.frame)

    nt_sequences = {}
    parents = {}
    children = {}
    unaligned_sequences = []
    if is_full:
        fallback = {}
        nodes_in_gene = set()
        for hit in diamond_hits:
            nodes_in_gene.add(hit.node)
            query = f"{hit.node}|{hit.frame}"
            parents[query] = hit

            if hit.node not in fallback:
                fallback[hit.node] = hit
        
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
        aligned_sequences.append((header, str(Seq(seq).translate())))
        
    top_file = path.join(top_location, f"{gene}.aln.fa")
    if debug > 1 or not path.exists(this_hmm_output) or stat(this_hmm_output).st_size == 0 or overwrite:
        with NamedTemporaryFile(dir=gettempdir()) as hmm_temp_file:
            system(f"hmmbuild '{hmm_temp_file.name}' '{top_file}' > /dev/null")

            if debug > 1:
                this_hmm_in = path.join(hmm_output_folder, f"{gene}_input.fa")
                writeFasta(this_hmm_in, aligned_sequences)
                system(
                f"hmmsearch -o {this_hmm_output} --domT 10.0 {hmm_temp_file.name} {this_hmm_in} > /dev/null",
                )
            else:
                with NamedTemporaryFile(dir=gettempdir()) as aligned_files:
                    writeFasta(aligned_files.name, aligned_sequences)
                    aligned_files.flush()
                    system(
                    f"hmmsearch --domtblout {this_hmm_output} --domT 10.0 {hmm_temp_file.name} {aligned_files.name} > /dev/null",
                    )

    if debug > 1:
        return "", [], [], [], []
    data = defaultdict(list)
    with open(this_hmm_output) as f:
        for line in f:
            if line.startswith("#"):
                continue
            while "  " in line:
                line = line.replace("  ", " ")
            line = line.strip().split()

            query = line[0]

            start, end, ali_start, ali_end, score = int(line[17]), int(line[18]), int(line[15]), int(line[16]), float(line[7])

            data[query].append((start - 1, end, score, ali_start, ali_end))

    output = []
    new_outs = []
    parents_done = set()
    for query, results in data.items():
        
        if query in parents:
            hit = parents[query]
            if not f"{hit.node}|{hit.frame}" in parents_done:
                for result in results:
                    start, end, score, ali_start, ali_end = result
                    start = start * 3
                    end = end * 3

                    sequence = nt_sequences[query][start: end]

                    if is_full:
                        new_qstart = start
                    else:
                        new_qstart = hit.qstart + start


                    parents_done.add(f"{hit.node}|{hit.frame}")
                    new_hit = HmmHit(node=hit.node, score=score, frame=hit.frame, qstart=new_qstart, qend=new_qstart + len(sequence), gene=hit.gene, query=hit.query, uid=hit.uid, refs=hit.refs, seq=sequence)
                    output.append((new_hit, ali_start, ali_end))

        if query in children:
            _, frame = query.split("|")
            parent = children[query]
            for result in results:
                start, end, score, ali_start, ali_end = result
                start = start * 3
                end = end * 3

                sequence = nt_sequences[query][start: end]

                if is_full:
                    new_qstart = start
                else:
                    new_qstart = hit.qstart + start


                clone = HmmHit(node=parent.node, score=score, frame=int(frame), qstart=new_qstart, qend=new_qstart + len(sequence), gene=parent.gene, query=parent.query, uid=parent.uid, refs=parent.refs, seq=sequence)
                new_outs.append((f"{clone.gene},{clone.node},{clone.frame}"))
                
                output.append((clone, ali_start, ali_end))

    kick_log = []
    for hit in diamond_hits:
        if not f"{hit.node}|{hit.frame}" in parents_done:
            kick_log.append(f"{hit.gene},{hit.node},{hit.frame}")

    # output, filtered_sequences_log = internal_filter_gene(output, debug)
    filtered_sequences_log = []

    return gene, [i[0] for i in output], new_outs, kick_log, filtered_sequences_log

def get_arg(transcripts_mapped_to, head_to_seq, is_full, hmm_output_folder, top_location, overwrite, debug, verbose):
    for gene, transcript_hits in transcripts_mapped_to:
        this_seqs = {}
        if is_full:
            for hit in transcript_hits:
                this_seqs[hit.node] = head_to_seq[hit.node]
        yield gene, transcript_hits, this_seqs, is_full, hmm_output_folder, top_location, overwrite, debug, verbose


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

def do_folder(input_folder, args):
    hits_db = RocksDB(path.join(input_folder, "rocksdb", "hits"))
    tk = TimeKeeper(KeeperMode.DIRECT)
    transcripts_mapped_to = get_diamondhits(
        hits_db
    )

    head_to_seq = {}
    if args.full:
        seq_db = RocksDB(path.join(input_folder, "rocksdb", "sequences", "nt"))
        recipe = seq_db.get("getall:batches").split(",")
        head_to_seq = get_head_to_seq(seq_db, recipe)

    hmm_output_folder = path.join(input_folder, "hmmsearch")
    
    if (args.debug or args.overwrite) and path.exists(hmm_output_folder):
        rmtree(hmm_output_folder, ignore_errors=True)
    if not path.exists(hmm_output_folder):
        system(f"mkdir {hmm_output_folder}")

    hmm_location = path.join(args.orthoset_input, args.orthoset, "hmms")
    top_location = path.join(input_folder, "top")

    arguments = get_arg(transcripts_mapped_to, head_to_seq, args.full, hmm_output_folder, top_location, args.overwrite, args.debug, args.verbose)

    all_hits = []

    printv(f"Processing {input_folder}", args.verbose, 1)

    if args.processes <= 1:
        for this_arg in arguments:
            all_hits.append(hmm_search(*this_arg))
    else:
        with Pool(args.processes) as p:
            all_hits = p.starmap(hmm_search, arguments)

    log = ["Gene,Node,Frame"]
    klog = ["Gene,Node,Frame"]
    ilog = ["Kicked Gene,Header,Frame,Score,Start,End,Reason,Master Gene,Header,Frame,Score,Start,End"]
    for gene, hits, logs, klogs, ilogs in all_hits:
        if not gene:
            continue
        
        hits_db.put_bytes(f"gethmmhits:{gene}", json.encode([i for i in hits]))
        log.extend(logs)
        klog.extend(klogs)
        ilog.extend(ilogs)

    del hits_db

    with open(path.join(input_folder, "hmmsearch_new.log"), "w") as f:
        f.write("\n".join(log))
    with open(path.join(input_folder, "hmmsearch_kick.log"), "w") as f:
        f.write("\n".join(klog))
    with open(path.join(input_folder, "hmmsearch_internal.log"), "w") as f:
        f.write("\n".join(ilog))

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