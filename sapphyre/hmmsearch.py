from collections import defaultdict
from shutil import rmtree
from tempfile import NamedTemporaryFile, TemporaryDirectory
from .diamond import ReporterHit as Hit
from wrap_rocks import RocksDB
from .utils import printv, gettempdir, parseFasta, writeFasta
from . import rocky
import os
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

def hmm_search(gene, path, flagged_sequences, hmm_output_folder, hmm_location, overwrite, debug):
    # printv(f"Processing: {gene}", 1)
    aligned_sequences = []
    sequence_data = {}
    this_aa_output = os.path.join(hmm_output_folder, "aa", f"{gene}.aa.fa")
    this_nt_output = os.path.join(hmm_output_folder, "nt", f"{gene}.nt.fa")

    this_hmm_output = os.path.join(hmm_output_folder, f"{gene}.txt")
    if debug or not os.path.exists(this_aa_output) or os.stat(this_aa_output).st_size == 0 or overwrite:
        for header, sequence in flagged_sequences:
            this_aa = str(Seq(sequence).translate())

            sequence_data[header] = sequence

            aligned_sequences.append((header, this_aa))

        hmm_file = os.path.join(hmm_location, f"{gene}.hmm")

        with TemporaryDirectory(dir=gettempdir()) as temp_dir:
            with NamedTemporaryFile(dir=temp_dir) as aligned_files, NamedTemporaryFile(dir=temp_dir) as result_file:
                writeFasta(aligned_files.name, aligned_sequences)
                aligned_files.flush()
                if debug == 2:
                    os.system(
                    f"hmmsearch -o {this_hmm_output} --domT 10.0 {hmm_file} {aligned_files.name} > /dev/null",
                    )
                    return "", []
                else:
                    os.system(
                    f"hmmsearch --domtblout {result_file.name} --domT 10.0 {hmm_file} {aligned_files.name} > /dev/null",
                    )
                    
            
                result = []

                with open(result_file.name) as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        while "  " in line:
                            line = line.replace("  ", " ")
                        line = line.strip().split()

                        query, start, end, score = line[0], int(line[17]), int(line[18]), float(line[13])
                        start -= 1
                        
                        result.append((query, start, end, score))

    result.sort(key=lambda x: x[3], reverse=True)
    done = set()
    aa_output = []
    nt_output = []
    for query, start, end, score in result:
        if query not in done:
            done.add(query)
        else:
            continue

        seq = sequence_data[query][start*3:end*3]

        aa_output.append((query, str(Seq(seq).translate())))
        nt_output.append((query, seq))

    if aa_output:
        writeFasta(this_nt_output, nt_output)
        if debug:
            writeFasta(this_aa_output, aa_output)
            return True
        
        else:
            with TemporaryDirectory(dir=gettempdir()) as temp_dir:
                with NamedTemporaryFile(dir=temp_dir) as aligned_files, NamedTemporaryFile(dir=temp_dir) as temp_result:
                    writeFasta(aligned_files.name, aa_output)
                    aligned_files.flush()
                    os.system(
                        f"mafft --anysymbol --quiet --jtt 1 --addfragments {aligned_files.name} --thread 1 {path} > {temp_result.name}"
                    )

                    writeFasta(this_aa_output, parseFasta(temp_result.name, True))

    return True

def get_arg(gene_paths, hmm_output_folder, hmm_location, overwrite, debug, og_aa):
    for this_path in gene_paths:
        flagged_sequences = list(parseFasta(this_path))

        gene_raw = os.path.basename(this_path)
        
        yield gene_raw.split(".")[0], os.path.join(og_aa, gene_raw), flagged_sequences, hmm_output_folder, hmm_location, overwrite, debug


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
    flagged_aa = os.path.join(input_folder, "trimmed", "flagged")

    og_aa = os.path.join(input_folder, "trimmed", "aa")

    hmm_output_folder = os.path.join(input_folder, "hmmsearch")
    hmm_output_folder_aa = os.path.join(input_folder, "hmmsearch", "aa")
    hmm_output_folder_nt = os.path.join(input_folder, "hmmsearch", "nt")
    
    if (args.debug or args.overwrite) and os.path.exists(hmm_output_folder):
        rmtree(hmm_output_folder, ignore_errors=True)
    if not os.path.exists(hmm_output_folder):
        os.mkdir(hmm_output_folder)
        os.mkdir(hmm_output_folder_aa)
        os.mkdir(hmm_output_folder_nt)

    hmm_location = os.path.join(args.orthoset_input, args.orthoset, "hmms")

    gene_paths = [os.path.join(flagged_aa, i) for i in os.listdir(flagged_aa) if i.endswith(".fa")]
    arguments = get_arg(gene_paths, hmm_output_folder, hmm_location, args.overwrite, args.debug, og_aa)

    all_hits = []

    if args.processes <= 1:
        for this_arg in arguments:
            all_hits.append(hmm_search(*this_arg))
    else:
        with Pool(args.processes) as p:
            all_hits = p.starmap(hmm_search, arguments)

    return True

def main(args):
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exist.", args.verbose, 0)
        return False
    rocky.create_pointer(
        "rocks_orthoset_db",
        os.path.join(args.orthoset_input, args.orthoset, "rocksdb"),
    )
    result = []
    for input_path in args.INPUT:
        rocks_db_path = os.path.join(input_path, "rocksdb")
        rocky.create_pointer(
            "rocks_nt_db",
            os.path.join(rocks_db_path, "sequences", "nt"),
        )
        rocky.create_pointer("rocks_hits_db", os.path.join(rocks_db_path, "hits"))
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