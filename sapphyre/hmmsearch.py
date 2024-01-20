import warnings
from Bio import BiopythonWarning
from collections import defaultdict
from shutil import rmtree
from tempfile import NamedTemporaryFile
from .diamond import ReporterHit as Hit
from wrap_rocks import RocksDB
from .utils import printv, gettempdir, parseFasta, writeFasta
from os import path, system, stat
from msgspec import json
from multiprocessing import Pool
from phymmr_tools import translate, bio_revcomp
from Bio.Seq import Seq
from .timekeeper import TimeKeeper, KeeperMode

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


def hmm_search(gene, diamond_hits, hmm_output_folder, hmm_location, overwrite, debug, verbose):
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
    node_has_frames_already = hits_have_frames_already.copy()
    for hit in diamond_hits:
        unaligned_sequences = []
        raw_sequence = hit.seq
        frame = hit.frame
        query = f"{hit.node}|{frame}"
        unaligned_sequences.append((query, raw_sequence))
        parents[query] = hit

        for shift_by in [1, 2]:
            shifted = shift(frame, shift_by)
            if not shifted in node_has_frames_already[hit.node]:
                new_query = f"{hit.node}|{shifted}"
                unaligned_sequences.append((new_query, raw_sequence[shift_by:]))
                node_has_frames_already[hit.node].add(shifted)
                children[new_query] = hit
            
        
        new_sequences = []
        for header, seq in unaligned_sequences:
            nt_sequences[header] = seq
            new_sequences.append((header, str(Seq(seq).translate())))

        aligned_sequences.extend(new_sequences)  #
        
    hmm_file = path.join(hmm_location, f"{gene}.hmm")
    if debug or not path.exists(this_hmm_output) or stat(this_hmm_output).st_size == 0 or overwrite:
            if debug:
                this_hmm_in = path.join(hmm_output_folder, f"{gene}_input.fa")
                writeFasta(this_hmm_in, aligned_sequences)
                system(
                f"hmmsearch -o {this_hmm_output} --domT 10.0 {hmm_file} {this_hmm_in} > /dev/null",
                )
            else:
                with NamedTemporaryFile(dir=gettempdir()) as aligned_files:
                    writeFasta(aligned_files.name, aligned_sequences)
                    aligned_files.flush()
                    system(
                    f"hmmsearch --domtblout {this_hmm_output} --domT 10.0 {hmm_file} {aligned_files.name} > /dev/null",
                    )

    if debug:
        return "", [], []
    data = defaultdict(list)
    with open(this_hmm_output) as f:
        for line in f:
            if line.startswith("#"):
                continue
            while "  " in line:
                line = line.replace("  ", " ")
            line = line.strip().split()

            query = line[0]

            start, end = int(line[17]), int(line[18])

            data[query].append((start - 1, end))

    output = []
    new_outs = []
    parents_done = set()
    for query, results in data.items():
        
        if query in parents:
            hit = parents[query]
            if not f"{hit.node}|{hit.frame}" in parents_done:
                for result in results:
                    start, end = result
                    start = start * 3
                    end = end * 3

                    sequence = nt_sequences[query][start: end]

                    new_qstart = hit.qstart + start


                    parents_done.add(f"{hit.node}|{hit.frame}")
                    new_hit = Hit(node=hit.node, frame=int(frame), qstart=new_qstart, qend=new_qstart + len(sequence), gene=hit.gene, query=hit.query, uid=hit.uid, refs=hit.refs, seq=sequence)
                    output.append(new_hit)

        if query in children:
            _, frame = query.split("|")
            parent = children[query]
            for result in results:
                start, end = result
                start = start * 3
                end = end * 3

                sequence = nt_sequences[query][start: end]

                new_qstart = parent.qstart + start


                clone = Hit(node=parent.node, frame=int(frame), qstart=new_qstart, qend=new_qstart + len(sequence), gene=parent.gene, query=parent.query, uid=parent.uid, refs=parent.refs, seq=sequence)
                new_outs.append((f"{clone.gene},{clone.node},{clone.frame}"))
                output.append(clone)

    kick_log = []
    for hit in diamond_hits:
        if not f"{hit.node}|{hit.frame}" in parents_done:
            kick_log.append(f"{hit.gene},{hit.node},{hit.frame}")

    return gene, output, new_outs, kick_log

def get_arg(transcripts_mapped_to, hmm_output_folder, hmm_location, overwrite, debug, verbose):
    for gene, transcript_hits in transcripts_mapped_to:
        yield gene, transcript_hits, hmm_output_folder, hmm_location, overwrite, debug, verbose


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

    hmm_output_folder = path.join(input_folder, "hmmsearch")
    
    if (args.debug or args.overwrite) and path.exists(hmm_output_folder):
        rmtree(hmm_output_folder, ignore_errors=True)
    if not path.exists(hmm_output_folder):
        system(f"mkdir {hmm_output_folder}")

    hmm_location = path.join(args.orthoset_input, args.orthoset, "hmms")

    arguments = get_arg(transcripts_mapped_to, hmm_output_folder, hmm_location, args.overwrite, args.debug, args.verbose)

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
    for gene, hits, logs, klogs in all_hits:
        if not gene:
            continue
        hits_db.put_bytes(f"gethmmhits:{gene}", json.encode(hits))
        log.extend(logs)
        klog.extend(klogs)

    del hits_db

    with open(path.join(input_folder, "hmmsearch_new.log"), "w") as f:
        f.write("\n".join(log))
    with open(path.join(input_folder, "hmmsearch_kick.log"), "w") as f:
        f.write("\n".join(klog))

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