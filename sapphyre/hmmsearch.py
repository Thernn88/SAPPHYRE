from collections import defaultdict
from shutil import rmtree
from tempfile import NamedTemporaryFile
from .diamond import ReporterHit as Hit
from wrap_rocks import RocksDB
from .utils import printv, gettempdir, parseFasta, writeFasta
from . import rocky
from os import path, system, stat
from msgspec import json
from multiprocessing import Pool
from phymmr_tools import translate, bio_revcomp
from Bio.Seq import Seq

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


def hmm_search(gene, diamond_hits, parent_sequences, hmm_output_folder, hmm_location, overwrite, debug):
    printv(f"Processing: {gene}", 1)
    aligned_sequences = []
    this_hmm_output = path.join(hmm_output_folder, f"{gene}.hmmout")
    if debug or not path.exists(this_hmm_output) or stat(this_hmm_output).st_size == 0 or overwrite:
        hits_have_frames_already = defaultdict(set)
        for hit in diamond_hits:
            hits_have_frames_already[hit.node].add(hit.frame)

        nt_sequences = {}
        parents = {}
        children = {}
        node_has_frames_already = hits_have_frames_already.copy()
        for hit in diamond_hits:
            unaligned_sequences = []
            sequence = parent_sequences[hit.node]
            raw_sequence = sequence[hit.qstart - 1 : hit.qend]
            frame = hit.frame
            query = f"{hit.node}|{frame}"
            unaligned_sequences.append((query, raw_sequence))
            parents[query] = hit
            if abs(frame) == 1:
                next_shift = [1,2]
            elif abs(frame) == 2:
                next_shift = [-1,1]
            elif abs(frame) == 3:
                next_shift = [-1,-2]

            for shift_by in next_shift:
                shifted = shift(frame, shift_by)
                if not shifted in node_has_frames_already[hit.node]:
                    new_query = f"{hit.node}|{shifted}"
                    if shifted < 0:
                        unaligned_sequences.append((new_query, sequence[hit.qstart - 1 : hit.qend - shift_by]))
                    else:
                        unaligned_sequences.append((new_query, bio_revcomp(sequence[hit.qstart + shift_by - 1 : hit.qend])))
                    node_has_frames_already[hit.node].add(shifted)
                    children[new_query] = hit
            
            
        
            new_sequences = []
            for header, seq in unaligned_sequences:
                if frame < 0:
                    seq = bio_revcomp(seq)
                nt_sequences[header] = seq
                new_sequences.append((header, str(Seq(seq).translate())))

            aligned_sequences.extend(new_sequences)  #
        hmm_file = path.join(hmm_location, f"{gene}.hmm")
        with NamedTemporaryFile(dir=gettempdir()) as aligned_files:
            writeFasta(aligned_files.name, aligned_sequences)
            aligned_files.flush()
            if debug:
                system(
                f"hmmsearch -o {this_hmm_output} --domT 10.0 {hmm_file} {aligned_files.name} > /dev/null",
                )
            else:
                system(
                f"hmmsearch --domtblout {this_hmm_output} --domT 10.0 {hmm_file} {aligned_files.name} > /dev/null",
                )

    if debug:
        return "", []
    data = {}
    with open(this_hmm_output) as f:
        for line in f:
            if line.startswith("#"):
                continue
            while "  " in line:
                line = line.replace("  ", " ")
            line = line.strip().split()

            query = line[0]

            if query not in data:
                start, end = int(line[17]), int(line[18])

                data[query] = (start - 1, end)

    output = []
    parents_done = set()
    for query, result in data.items():
        
        if query in parents:
            hit = parents[query]
            if not f"{hit.node}|{hit.frame}" in parents_done:
                parents_done.add(f"{hit.node}|{hit.frame}")
                output.append(hit)

        if query in children:
            _, frame = query.split("|")
            parent = children[query]
            
            sequence = nt_sequences[query]

            start, end = result
            start = start * 3
            end = end * 3

            clone = Hit(node=parent.node, frame=int(frame), qstart=hit.qstart, qend=hit.qend, gene=parent.gene, query=parent.query, uid=parent.uid, refs=parent.refs, seq=sequence[start: end])

            output.append(clone)

    return gene, output

def get_arg(transcripts_mapped_to, raw_db_sequences, hmm_output_folder, hmm_location, overwrite, debug):
    for gene, transcript_hits in transcripts_mapped_to:
        this_seqs = {}
        for hit in transcript_hits:
            this_seqs[hit.node] = raw_db_sequences[hit.node]
    
        yield gene, transcript_hits, this_seqs, hmm_output_folder, hmm_location, overwrite, debug


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
    transcripts_mapped_to = get_diamondhits(
        rocky.get_rock("rocks_hits_db"),
    )

    nt_db = rocky.get_rock("rocks_nt_db")
    recipe = nt_db.get("getall:batches").split(",")
    raw_db_sequences = get_head_to_seq(nt_db, recipe)

    hmm_output_folder = path.join(input_folder, "hmmsearch")
    
    if (args.debug or args.overwrite) and path.exists(hmm_output_folder):
        rmtree(hmm_output_folder, ignore_errors=True)
    if not path.exists(hmm_output_folder):
        system(f"mkdir {hmm_output_folder}")

    hmm_location = path.join(args.orthoset_input, args.orthoset, "hmms")

    arguments = get_arg(transcripts_mapped_to, raw_db_sequences, hmm_output_folder, hmm_location, args.overwrite, args.debug)

    all_hits = []

    if args.processes <= 1:
        for this_arg in arguments:
            all_hits.append(hmm_search(*this_arg))
    else:
        with Pool(args.processes) as p:
            all_hits = p.starmap(hmm_search, arguments)


    hits_db = rocky.get_rock("rocks_hits_db")
    for gene, hits in all_hits:
        if not gene:
            continue
        hits_db.put_bytes(f"gethmmhits:{gene}", json.encode(hits))

    return True

def main(args):
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exist.", args.verbose, 0)
        return False
    rocky.create_pointer(
        "rocks_orthoset_db",
        path.join(args.orthoset_input, args.orthoset, "rocksdb"),
    )
    result = []
    for input_path in args.INPUT:
        rocks_db_path = path.join(input_path, "rocksdb")
        rocky.create_pointer(
            "rocks_nt_db",
            path.join(rocks_db_path, "sequences", "nt"),
        )
        rocky.create_pointer("rocks_hits_db", path.join(rocks_db_path, "hits"))
        result.append(
            do_folder(
                input_path,
                args,
            ),
        )
        rocky.close_pointer("rocks_nt_db")
        rocky.close_pointer("rocks_hits_db")
    rocky.close_pointer("rocks_orthoset_db")
    
    return all(result)