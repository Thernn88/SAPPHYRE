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

def hmm_search(gene, diamond_hits, parent_sequences, hmm_output_folder, hmm_location, overwrite, debug):
    # printv(f"Processing: {gene}", 1)
    aligned_sequences = []
    this_hmm_output = path.join(hmm_output_folder, f"{gene}.hmmout")
    if debug or not path.exists(this_hmm_output) or stat(this_hmm_output).st_size == 0 or overwrite:
        for hit in diamond_hits:
            # raw_sequence = parent_sequences[hit.node]
            # frame = hit.frame

            # if frame < 0:
            #     raw_sequence = bio_revcomp(raw_sequence)

            # frame_offset = abs(int(frame))-1
            # raw_sequence = raw_sequence[frame_offset:]

            # this_aa = str(Seq(raw_sequence).translate())
            this_aa = str(Seq(hit.seq).translate())

            aligned_sequences.append((hit.node+"_"+str(hit.frame), this_aa))

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
    data = defaultdict(list)
    with open(this_hmm_output) as f:
        for line in f:
            if line.startswith("#"):
                continue
            while "  " in line:
                line = line.replace("  ", " ")
            line = line.strip().split()

            query, start, end, index = line[0], int(line[17]), int(line[18]), line[9]

            data[query].append((start - 1, end, index))

    output = []
    for hit in diamond_hits:
        query = hit.node+"_"+str(hit.frame)
        if query not in data:
            continue

        # raw_sequence = parent_sequences[hit.node]
        # frame = hit.frame

        # if frame < 0:
        #     raw_sequence = bio_revcomp(raw_sequence)

        # for start, end, index in data[query]:
            
        #     node = hit.node+'-'+index
            
        #     start = start * 3
        #     end = end * 3

        #     hstart = start
        #     hend = end
            
        #     seq = hit.seq[start:end]
        #     this_hit = Hit(node, hit.frame, hit.qstart, hit.qend, hit.gene, hit.query, hit.uid, hit.refs, seq, hstart, hend)
        output.append(hit)

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