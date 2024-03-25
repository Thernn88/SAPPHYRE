from time import time
import warnings
from Bio import BiopythonWarning
from collections import defaultdict
from shutil import rmtree
from tempfile import NamedTemporaryFile, TemporaryDirectory
from .diamond import ReferenceHit, ReporterHit as Hit
from wrap_rocks import RocksDB
from .utils import printv, gettempdir, parseFasta, writeFasta
from os import path, system, stat
from msgspec import Struct, json
from multiprocessing import Pool
from sapphyre_tools import translate, bio_revcomp, get_overlap
from pyfastx import Fastq
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

    gene_based_results = {}
    for gene in genes_to_process:
        gene_result = rocks_hits_db.get_bytes(f"gethits:{gene}")
        this_hits = json.decode(gene_result, type=list[Hit])
        if not gene_result:
            printv(
                f"WARNING: No hits found for {gene}. If you are using a gene list file this may be a non-issue",
                0,
            )
            continue
        gene_based_results[gene] = this_hits

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


def hmm_search(gene, diamond_hits, this_seqs, top_location, verbose):
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    printv(f"Processing: {gene}", verbose, 2)

    unaligned_sequences = []
    query_to_hit = {}
    hits_have_frames_already = defaultdict(set)
    for hit in diamond_hits:
        raw_sequence = hit.seq
        frame = hit.frame
        query = f"{hit.node}|{frame}"
        unaligned_sequences.append((query, raw_sequence))
        query_to_hit[query] = hit

        for shift_by in [1, 2]:
            shifted = shift(frame, shift_by)
            if not shifted in hits_have_frames_already[hit.node]:
                new_query = f"{hit.node}|{shifted}"
                unaligned_sequences.append((new_query, raw_sequence[shift_by:]))
                hits_have_frames_already[hit.node].add(shifted)
                query_to_hit[new_query] = hit
            
    nt_sequences = {}
    aligned_sequences = []
    for header, seq in unaligned_sequences:
        nt_sequences[header] = seq
        aligned_sequences.append((header, str(Seq(seq).translate())))

    output = []
    
    with TemporaryDirectory(dir=gettempdir()) as dir:
        fastq_gene_file = path.join(dir, gene+".fastq")
        with open(fastq_gene_file, "w") as fp:
            for line in this_seqs:
                fp.write(line)
        system(f"spades --only-assembler --sc -k 33 -s {fastq_gene_file} -o {path.join(dir,gene)} > {path.join(dir, 'spades_log.txt')}")

        this_hmm_fasta = path.join(dir, gene,"six_fold_and_aln.fasta")
        contig_path = path.join(dir, gene,"contigs.fasta")

        system(f"fastatranslate {contig_path} > {this_hmm_fasta}")

        six_fold = list(parseFasta(this_hmm_fasta, True)) + aligned_sequences
        writeFasta(this_hmm_fasta, six_fold)
        

        this_hmm_output = path.join(dir, gene,"hmmresult.hmm")
        hmm_temp_file = path.join(dir, gene,"temp_hmm.hmmdb")

        top_file = path.join(top_location ,gene+".aln.fa")

        system(f"hmmbuild '{hmm_temp_file}' '{top_file}' > /dev/null")

        system(
            f"hmmsearch --domtblout {this_hmm_output} --domT 10.0 {hmm_temp_file} {this_hmm_fasta} > {path.join(dir, 'hmm_log.txt')}",
            )
        
        if not path.exists(this_hmm_output):
            print(gene,"failed!")
            return None, []

        highest_score_headers = {}
        data = defaultdict(list)
        with open(this_hmm_output) as fp:
            for line in fp:
                if line[0] == "#":
                    continue
                while "  " in line:
                    line = line.replace("  ", " ")
                line = line.strip().split()

                if line[-1] != "-":
                    if line[0] in highest_score_headers:
                        _, _, _, score = highest_score_headers[line[0]]
                        if float(line[7]) > score:
                            highest_score_headers[line[0]] = line[-1], int(line[17]), int(line[18]), float(line[7])
                    else:
                        highest_score_headers[line[0]] = line[-1], int(line[17]), int(line[18]), float(line[7])
                else:
                    query = line[0]

                    start, end, ali_start, ali_end, score = int(line[17]), int(line[18]), int(line[15]), int(line[16]), float(line[7])

                    data[query].append((start - 1, end, score, ali_start, ali_end))

        headers = {header: (frame, start, end, score) for header, (frame, start, end, score) in highest_score_headers.items()}

        queries = data.items()
        hit_queries = []

        for query, results in queries:
            if query in query_to_hit:
                hit = query_to_hit[query]
                for result in results:
                    start, end, score, ali_start, ali_end = result
                    start = start * 3
                    end = end * 3

                    sequence = nt_sequences[query][start: end]
                    new_qstart = hit.qstart + start

                    new_hit = HmmHit(node=hit.node, score=score, frame=hit.frame, qstart=new_qstart, qend=new_qstart + len(sequence), gene=hit.gene, query=hit.query, uid=hit.uid, refs=hit.refs, seq=sequence)
                    hit_queries.append(hit.query)
                    output.append(new_hit)

        for header, seq in parseFasta(contig_path, True):
            if header in headers:
                frame,start,end,score = headers[header]
                is_revcomp = False
                if "rev" in frame:
                    seq = bio_revcomp(seq)

                    is_revcomp = True

                frame = int(frame.split("translate(")[1].split(")")[0])
                if is_revcomp:
                    frame = frame // - 1

                start = start * 3
                end = end * 3
                
                seq = seq[(start):(end)]

                output.append(
                    HmmHit(
                        None,
                        score,
                        frame,
                        start,
                        end,
                        gene,
                        max(set(hit_queries), key=hit_queries.count),
                        abs(hash(time())),
                        [],
                        seq
                    )
                )

    return gene, output

def get_arg(transcripts_mapped_to, gene_fastq_paths, top_location, verbose):
    for gene, transcript_hits in transcripts_mapped_to.items():
        # if "99" in gene:
        yield gene, transcript_hits, gene_fastq_paths[gene], top_location, verbose


def get_head_to_seq(nt_db, recipe, is_full):
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

    last_id = int(list(head_to_seq.keys())[-1].split("_")[1])

    if not is_full:
        return {}, last_id

    return head_to_seq, last_id

def do_folder(input_folder, args):
    hits_db = RocksDB(path.join(input_folder, "rocksdb", "hits"))
    tk = TimeKeeper(KeeperMode.DIRECT)
    transcripts_mapped_to = get_diamondhits(
        hits_db
    )

    positions_to_keep = defaultdict(dict)
    keep_set = defaultdict(set)
    nt_db_path = path.join(input_folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(nt_db_path)

    recipe = nt_db.get("getall:batches").split(",")
    head_to_seq, last_id = get_head_to_seq(nt_db, recipe, args.full)
    print(last_id)
    try:
        original_positons = json.decode(nt_db.get("getall:original_positions"), type=dict[str, int])
    except:
        try:
            original_positons = json.decode(nt_db.get("getall:original_positions"), type=dict[str, tuple[int, int, tuple[int, int]|int|None]])
        except:
            original_positons = json.decode(nt_db.get("getall:original_positions"), type=dict[str, tuple[int, int]])

    original_inputs = json.decode(nt_db.get("getall:original_inputs"), type=list[str])
    present_genes = hits_db.get("getall:presentgenes").split(',')
    file_index = 0
    print("Generating Fastq dict")
    for gene in present_genes:
        this_hits = transcripts_mapped_to[gene]
        for i, hit in enumerate(this_hits):
            og_pos_entry = original_positons[hit.node]
            if type(og_pos_entry) != tuple:
                line_index = og_pos_entry
            elif len(og_pos_entry) == 2:
                file_index, line_index = og_pos_entry
            else:
                file_index, line_index, _ = og_pos_entry
                    
            keep_set[file_index].add(line_index)


            positions_to_keep[file_index][line_index] = hit.gene, hit.node

    out_dict = defaultdict(list)
    for file_index, file in enumerate(original_inputs):
        for i, (header, seq, qual) in enumerate(Fastq(file, build_index = False)):
            if i in keep_set[file_index]:
                gene, node = positions_to_keep[file_index][i]
                out_dict[gene].append((f"@{node}\n{seq}\n+{node}\n{qual}\n"))

    hmm_output_folder = path.join(input_folder, "hmmsearch")
    
    if (args.debug or args.overwrite) and path.exists(hmm_output_folder):
        rmtree(hmm_output_folder, ignore_errors=True)
    if not path.exists(hmm_output_folder):
        system(f"mkdir {hmm_output_folder}")

    top_location = path.join(input_folder, "top")

    arguments = get_arg(transcripts_mapped_to, out_dict, top_location, args.verbose)

    all_hits = []

    printv(f"Processing {input_folder}", args.verbose, 1)

    if args.processes <= 1:
        for this_arg in arguments:
            all_hits.append(hmm_search(*this_arg))
    else:
        with Pool(args.processes) as p:
            all_hits = p.starmap(hmm_search, arguments)

    template = "NODE_{}"
    for gene, hits in all_hits:
        if not gene:
            continue

        for hit in hits:
            if hit.node is None:
                hit.node = template.format(last_id + 1)
                last_id += 1

        input(hits)
        
        hits_db.put_bytes(f"gethmmhits:{gene}", json.encode([i for i in hits]))

    del hits_db

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