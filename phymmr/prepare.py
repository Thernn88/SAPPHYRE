from __future__ import annotations

import gzip
import json
import math
import mmap
import os
import re
from itertools import count
from multiprocessing.pool import Pool
from shutil import rmtree
from time import time

import phymmr_tools
import wrap_rocks
import xxhash
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm

from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv


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


def get_seq_count(filename, is_gz):
    counter = count()
    if is_gz:
        f = gzip.open(filename, "r+")
        for line in f:
            if line[0] == ord(">") or line[0] == ord("+"):
                next(counter)
    else: # We can grab it cheaper
        f = open(filename, "r+")
        buf = mmap.mmap(f.fileno(), 0)
        readline = buf.readline
        line = readline()
        while line:
            if line[0] == ord(">") or line[0] == ord("+"):
                next(counter)
            line = readline()
    return next(counter)


def N_trim(parent_sequence, minimum_sequence_length, tike):
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
            if length >= minimum_sequence_length:
                raw_seq = parent_sequence[start + 1: end]
                tike.timer1("now")
                yield raw_seq
    else:
        yield parent_sequence


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


def translate(in_path, out_path, translate_program="fastatranslate", genetic_code=1):
    os.system(
        f"{translate_program} --geneticcode {genetic_code} '{in_path}' > '{out_path}'"
    )
    os.remove(in_path)


def main(args):
    if not os.path.exists(args.INPUT):
        printv("ERROR: An existing directory must be provided.", args.verbose, 0)
        return False

    # Append computed time for each loop.
    trim_times = TimeKeeper(KeeperMode.SUM)
    dedup_time = 0
    global_time_keeper = TimeKeeper(KeeperMode.DIRECT)

    PROT_MAX_SEQS_PER_LEVEL = args.sequences_per_level
    MINIMUM_SEQUENCE_LENGTH = args.minimum_sequence_length
    num_threads = args.processes
    folder = args.INPUT

    allowed_filetypes = ["fa", "fas", "fasta", "fq", "fastq", "fq", "fastq.gz", "fq.gz"]

    rocksdb_folder_name = "rocksdb"
    sequences_db_name = "sequences"
    core_directory = "PhyMMR"
    secondary_directory = os.path.join(core_directory, os.path.basename(os.path.normpath(folder)))

    # Create necessary directories
    printv("Creating directories", args.verbose)

    os.makedirs(secondary_directory, exist_ok=True)

    taxa_runs = {}

    # Scan all the files in the input. Remove _R# and merge taxa
    for file in os.listdir(folder):
        file_type = '.'.join(file.split(".")[1:])
        if (
            os.path.isfile(os.path.join(folder, file))
            and file_type in allowed_filetypes
        ):
            taxa = file.split(".")[0]

            formatted_taxa = truncate_taxa(taxa, extension=".fa")

            taxa_runs.setdefault(formatted_taxa, [])
            taxa_runs[formatted_taxa].append(file)
    
    for formatted_taxa_out, components in taxa_runs.items():
        global_time_keeper.lap()
        taxa_time_keeper = TimeKeeper(KeeperMode.DIRECT)
        printv(f"Preparing: {formatted_taxa_out}", args.verbose, 0)

        taxa_destination_directory = os.path.join(
            secondary_directory, formatted_taxa_out
        )
        rocksdb_path = os.path.join(taxa_destination_directory, rocksdb_folder_name)
        sequences_db_path = os.path.join(rocksdb_path, sequences_db_name)

        if args.clear_database and os.path.exists(rocksdb_path):
            printv("Clearing old database", args.verbose)
            rmtree(rocksdb_path)

        printv("Creating rocksdb database", args.verbose)

        os.makedirs(rocksdb_path, exist_ok=True)

        db = wrap_rocks.RocksDB(sequences_db_path)

        os.makedirs(taxa_destination_directory, exist_ok=True)

        prepared_file_destination = os.path.join(
            taxa_destination_directory, formatted_taxa_out
        )

        prot_path = os.path.join(
            taxa_destination_directory, formatted_taxa_out.replace(".fa", "_prot.fa")
        )

        duplicates = {}
        rev_comp_save = {}
        transcript_mapped_to = {}
        dupes = count()
        this_index = 1

        fa_file_out = []
        printv("Formatting input sequences and inserting into database", args.verbose)
        for file in components:
            fa_file_path = os.path.join(folder, file)
            
            if ".fq" in fa_file_path or ".fastq" in fa_file_path:
                read_method = FastqGeneralIterator
            else:
                read_method = SimpleFastaParser

            is_compressed = file.split('.')[-1] == 'gz'
            if is_compressed:
                fasta_file = read_method(gzip.open(fa_file_path, "rt"))
            else:
                fasta_file = read_method(open(fa_file_path, encoding="utf-8"))
            
            if args.verbose: 
                seq_grab_count = get_seq_count(fa_file_path, is_compressed)
            for_loop = tqdm(fasta_file, total=seq_grab_count) if args.verbose else fasta_file
            for row in for_loop:
                header, parent_seq = row[:2]

                if not len(parent_seq) >= MINIMUM_SEQUENCE_LENGTH:
                    continue
                parent_seq = parent_seq.upper()
                # for seq in N_trim(parent_seq, MINIMUM_SEQUENCE_LENGTH):
                for seq, tt in N_trim(parent_seq, MINIMUM_SEQUENCE_LENGTH, trim_times):
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

        if os.path.exists("/run/shm"):
            tmp_path = "/run/shm"
        elif os.path.exists("/dev/shm"):
            tmp_path = "/dev/shm"
        else:
            tmp_path = os.path.join(folder, "tmp")
            os.makedirs(tmp_path, exist_ok=True)

        sequences_per_thread = math.ceil((len(fa_file_out) / 2) / num_threads) * 2
        translate_files = []
        for i in range(0, len(fa_file_out), sequences_per_thread):
            in_path = os.path.join(tmp_path, formatted_taxa_out + f'_{len(translate_files)}.fa')
            out_path = os.path.join(tmp_path, formatted_taxa_out + f'_prot_{len(translate_files)}.fa')
            translate_files.append((in_path, out_path))
            open(in_path, 'w').writelines(fa_file_out[i:i + sequences_per_thread])

        with Pool(num_threads) as translate_pool:
            translate_pool.starmap(translate, translate_files)

        if args.keep_prepared:
            open(prepared_file_destination, "w", encoding="UTF-8").writelines(fa_file_out)
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
    printv("N_trim time: {} seconds".format(trim_times.time1), args.verbose, 2)
    printv(f"Dedupe time: {dedup_time}", args.verbose, 2)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr Prepare"
    )