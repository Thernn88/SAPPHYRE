from __future__ import annotations

import argparse
import gzip
import json
import math
import mmap
import os
import re
from itertools import count
from multiprocessing.pool import Pool
from pathlib import Path
from queue import Queue
from shutil import rmtree
from subprocess import call
from sys import argv
from threading import Thread
from time import time

import phymmr_tools
import wrap_rocks
import xxhash
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm

from .timekeeper import KeeperMode, TimeKeeper
from .utils import ConcurrentLogger




def truncate_taxa(header: str, extension=None) -> str:
    """
    Given a fasta header, checks the end for problematic tails.
    If found, truncates the string.
    Returns the string + suffix to check for name matches.
    """
    # search for _# and _R#, where # is digits
    result = header
    m = re.search(r"_R?\d+$", header)
    if m:
        tail_length = m.end() - m.start()
        result = result[0:-tail_length]
    if extension:
        result = result + extension
    return result

def get_seq_count(filename: Path, is_gz: bool):
    count = 0
    if is_gz:
        f = gzip.open(filename, "r+")
        for line in f:
            if line[0] == ord(">") or line[0] == ord("+"):
                count += 1
    else: # We can grab it cheaper
        by = filename.read_bytes()
        ascii_le = ord('>')
        ascii_pl = ord('+')

        char_count = {ascii_le: 0, ascii_pl: 0}

        def map_char_count(l: bytearray):
            char_count.setdefault(l[0], 0)
            char_count[l[0]] += 1

        by_splt_lines = by.splitlines()
        list(map(map_char_count, by_splt_lines))

        if char_count[ascii_le] > 0:
            count = char_count[ascii_le]
        elif char_count[ascii_pl] > 0:
            count = char_count[ascii_pl]

    return count

def N_trim(parent_sequence, MINIMUM_SEQUENCE_LENGTH):
    t1 = time()
    if "N" in parent_sequence:
        # Get N indices and start and end of sequence
        indices = (
            [0]
            + [i for i, ltr in enumerate(parent_sequence) if ltr == "N"]
            + [len(parent_sequence)]
        )

        for i in range(0, len(indices) - 1):
            start = indices[i]
            end = indices[i + 1]

            length = end - start
            if length >= MINIMUM_SEQUENCE_LENGTH:
                raw_seq = parent_sequence[start + 1: end]
                yield raw_seq, time() - t1
                t1 = time()
    else:
        yield parent_sequence, time() - t1


def add_pc_to_db(db, key: int, data: list) -> str:
    """Join data into a string and add it to the database.

    Args:
        db: instance of rocks db
        key: integer representing the index
        data: list of stuff

    Returns:
        Formatted key string."""
    kstr = f"prepared:{key}"
    db.put(kstr, "".join(data))
    return kstr


def translate(in_path: Path, out_path: Path, translate_program="fastatranslate", genetic_code=1):
    call(
        " ".join([
            translate_program,
            '--geneticcode',
            str(genetic_code),
            str(in_path), 
            '>',
            str(out_path)
        ]),
        shell=True
    )
    in_path.unlink()


def main(args):
    msgq = Queue()
    printv = ConcurrentLogger(msgq)
    printv.start()
    input_path = Path(args.INPUT)

    if not input_path.exists():
        printv("ERROR: An existing directory must be provided.", args.verbose, 0)
        return False

    trim_times = []  # Append computed time for each loop.
    dedup_time = 0
    global_time_keeper = TimeKeeper(KeeperMode.DIRECT)

    PROT_MAX_SEQS_PER_LEVEL = args.sequences_per_level
    MINIMUM_SEQUENCE_LENGTH = args.minimum_sequence_length
    num_threads = args.processes
    folder = args.INPUT
    project_root = Path(__file__).parent.parent

    allowed_filetypes = ["fa", "fas", "fasta", "fq", "fastq", "fq", "fastq.gz", "fq.gz"]

    rocksdb_folder_name = "rocksdb"
    sequences_db_name = "sequences"
    core_directory = Path("PhyMMR").joinpath(input_path.parts[-1])
    secondary_directory = project_root.joinpath(core_directory)

    # Create necessary directories
    printv("Creating directories", args.verbose)
    secondary_directory.mkdir(parents=True, exist_ok=True)

    taxa_runs = {}

    # Scan all the files in the input. Remove _R# and merge taxa
    
    allowed_files_glob = list(input_path.glob("**/*.f[aq]*"))
    for f in allowed_files_glob:
        file_is_file = f.is_file()
        
        file_suffix = f.suffix
        file_suffix_allowed = file_suffix[1:] in allowed_filetypes
        
        if (file_is_file and file_suffix_allowed):
            taxa = f.stem.split(".")[0]

            formatted_taxa = truncate_taxa(taxa, extension=".fa")

            taxa_runs.setdefault(formatted_taxa, [])
            taxa_runs[formatted_taxa].append(f)
    
    for formatted_taxa_out, components in taxa_runs.items():
        global_time_keeper.lap()
        taxa_time_keeper = TimeKeeper(KeeperMode.DIRECT)
        printv(f"Preparing: {formatted_taxa_out}", args.verbose, 0)

        taxa_destination_directory = secondary_directory.joinpath(formatted_taxa_out)

        rocksdb_path = taxa_destination_directory.joinpath(rocksdb_folder_name)
        sequences_db_path = rocksdb_path.joinpath(sequences_db_name)

        if args.clear_database and rocksdb_path.exists():
            printv("Clearing old database", args.verbose)
            rmtree(rocksdb_path)

        printv("Creating rocksdb database", args.verbose)

        rocksdb_path.mkdir(parents=True, exist_ok=True)

        db = wrap_rocks.RocksDB(str(sequences_db_path))

        taxa_destination_directory.mkdir(parents=True, exist_ok=True)

        prepared_file_destination = taxa_destination_directory.joinpath(formatted_taxa_out)

        prot_path = taxa_destination_directory.joinpath(formatted_taxa_out.replace(".fa", "_prot.fa"))

        duplicates = {}
        rev_comp_save = {}
        transcript_mapped_to = {}
        dupes = count()
        this_index = 1

        fa_file_out = []
        printv("Formatting input sequences and inserting into database", args.verbose)
        for fa_file_path in components:           
            if any([fa_file_path.suffix in (".fq", ".fastq")]):
                read_method = FastqGeneralIterator
            else:
                read_method = SimpleFastaParser

            is_compressed = fa_file_path.suffix == '.gz'
            if is_compressed:
                fasta_file = read_method(gzip.open(fa_file_path, "rt"))
            else:
                fasta_file = read_method(open(fa_file_path, "r"))
                                
            if args.verbose: 
                seq_grab_count = get_seq_count(fa_file_path, is_compressed)
            for_loop = tqdm(fasta_file, total=seq_grab_count) if args.verbose else fasta_file
            for row in for_loop:
                header, parent_seq = row[:2]

                if not len(parent_seq) >= MINIMUM_SEQUENCE_LENGTH:
                    continue
                parent_seq = parent_seq.upper()
                # for seq in N_trim(parent_seq, MINIMUM_SEQUENCE_LENGTH):
                for seq, tt in N_trim(parent_seq, MINIMUM_SEQUENCE_LENGTH):
                    trim_times.append(tt)
                    length = len(seq)
                    header = f"NODE_{this_index}_length_{length}"

                    # Check for dupe, if so save how many times that sequence occured
                    seq_start = time()

                    if seq in transcript_mapped_to:
                        duplicates.setdefault(transcript_mapped_to[seq], 1)
                        duplicates[transcript_mapped_to[seq]] += 1
                        next(dupes)
                        continue
                    else:
                        transcript_mapped_to[seq] = header

                    # Rev-comp sequence. Save the reverse compliment in a hashmap with the original
                    # sequence so we don't have to rev-comp this unique sequence again
                    if seq in rev_comp_save:
                        rev_seq = rev_comp_save[seq]
                    else:
                        rev_seq = phymmr_tools.bio_revcomp(seq)
                        rev_comp_save[seq] = rev_seq

                    # Check for revcomp dupe, if so save how many times that sequence occured
                    if rev_seq in transcript_mapped_to:
                        duplicates.setdefault(transcript_mapped_to[rev_seq], 1)
                        duplicates[transcript_mapped_to[rev_seq]] += 1
                        next(dupes)
                        continue
                    else:
                        transcript_mapped_to[rev_seq] = header

                    seq_end = time()
                    dedup_time += seq_end - seq_start
                    this_index += 1

                    # If no dupe, write to prepared file and db
                    line = ">" + header + "\n" + seq + "\n"

                    fa_file_out.append(line)

                    # Get rid of space and > in header (blast/hmmer doesn't like it) Need to push modified external to remove this. ToDo.
                    preheader = header.replace(" ", "|")  # pre-hash header
                    # Key is a hash of the header
                    db.put(xxhash.xxh64_hexdigest(preheader), f"{preheader}\n{seq}")

        printv("Translating prepared file", args.verbose)

        aa_dupe_count = count()

        path_run_tmp = Path("/run/shm")
        path_dev_tmp = Path("/dev/shm")
        path_tmp_tmp = Path("/tmp")

        if path_dev_tmp.exists() and path_dev_tmp.is_dir():
            tmp_path = path_dev_tmp
        elif path_run_tmp.exists() and path_run_tmp.is_dir():
            tmp_path = path_run_tmp
        elif path_tmp_tmp.exists() and path_tmp_tmp.is_dir():
            tmp_path = path_tmp_tmp
        else:
            tmp_path = project_root.joinpath("tmp")
            tmp_path.mkdir(parents=True, exist_ok=True)

        sequences_per_thread = math.ceil((len(fa_file_out) / 2) / num_threads) * 2
        translate_files = []
        for i in range(0, len(fa_file_out), sequences_per_thread):
            in_path = tmp_path.joinpath(
                formatted_taxa_out + f'_{len(translate_files)}.fa'
            )
            out_path = tmp_path.joinpath(
                formatted_taxa_out + f'_prot_{len(translate_files)}.fa'
            )
            translate_files.append((in_path, out_path))
            in_path.write_text("\n".join(fa_file_out[i:i + sequences_per_thread]))

        with Pool(num_threads) as translate_pool:
            translate_pool.starmap(translate, translate_files)

        if args.keep_prepared:
            prepared_file_destination.write_text("\n".join(fa_file_out))
        del fa_file_out

        prot_components = []

        printv(f"Storing translated file in DB. Translate took {taxa_time_keeper.lap():.2f}s", args.verbose)

        out_lines = []
        for prepare_file, translate_file in translate_files:
            with open(translate_file, "r+", encoding="utf-8") as fp:
                for header, seq in SimpleFastaParser(fp):
                    header = header.replace(" ", "|")
                    if seq in transcript_mapped_to:
                        duplicates.setdefault(transcript_mapped_to[seq], 1)
                        duplicates[transcript_mapped_to[seq]] += 1
                        next(aa_dupe_count)
                        continue
                    else:
                        transcript_mapped_to[seq] = header
                    out_lines.append(">" + header + "\n" + seq + "\n")
            os.remove(translate_file)

        if args.keep_prepared:
            open(prot_path, 'w').writelines(out_lines)

        aa_dupes = next(aa_dupe_count)
        printv(f"AA dedupe took {taxa_time_keeper.lap():.2f}s. Kicked {aa_dupes} dupes", args.verbose, 2)

        levels = math.ceil(len(out_lines) / PROT_MAX_SEQS_PER_LEVEL)
        per_level = math.ceil(len(out_lines) / levels)

        component = 0
        for i in range(0, len(out_lines), per_level):
            component += 1

            data = out_lines[i: i + per_level]
            prot_components.append(str(component))
            db.put(f"getprot:{component}", "".join(data))

        db.put("getall:prot", ",".join(prot_components))

        printv(f"Translation and storing done! Took {taxa_time_keeper.lap():.2f}s", args.verbose)

        sequence_count = this_index - 1

        printv(
            "Inserted {} sequences over {} batches. Found {} NT and {} AA dupes.".format(
                sequence_count, levels, next(dupes), aa_dupes
            ),
            args.verbose,
        )

        del transcript_mapped_to  # Clear mem

        # Store the count of dupes in the database
        db.put("getall:dupes", json.dumps(duplicates))

        printv(f"{formatted_taxa_out} took {global_time_keeper.lap():.2f}s overall\n", args.verbose)

    printv(f"Finished! Took {global_time_keeper.differential():.2f}s overall.", args.verbose, 0)
    printv("N_trim time: {} seconds".format(sum(trim_times)), args.verbose, 2)
    printv(f"Dedupe time: {dedup_time}", args.verbose, 2)
    
    msgq.join()
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr Prepare"
    )
