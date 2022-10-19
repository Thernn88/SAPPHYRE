import argparse
import json
import math
import os
import re
from itertools import count
from shutil import rmtree
from sys import argv
from time import time
from multiprocessing.pool import Pool
import blosum_distance
import wrap_rocks
import xxhash
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm


def printv(msg, verbosity):
    if verbosity:
        print(msg)


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
                raw_seq = parent_sequence[start + 1 : end]
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

def translate(in_path, out_path, translate_program = "fastatranslate", genetic_code = 1):
    os.system(
            f"{translate_program} --geneticcode {genetic_code} '{in_path}' > '{out_path}'"
        )

def main(argv):
    trim_time = 0
    dedup_time = 0
    global_start = time()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", default="Snails", type=str, help="Path to input directory."
    )
    parser.add_argument(
        "-c",
        "--clear_database",
        action="store_true",
        help="Overwrite existing rocksdb database.",
    )
    parser.add_argument(
        "-ml",
        "--minimum_sequence_length",
        default=90,
        type=int,
        help="Minimum input sequence length.",
    )
    parser.add_argument(
        "-sl",
        "--sequences_per_level",
        default=1000000,
        type=int,
        help="Amount of sequences to store per database entry.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Number of threads used to call processes.",
    )
    parser.add_argument(
        "-k",
        "--keep_prepared",
        action="store_true",
        help="Writes the prepared input fasta into the output taxa directory.",
    )
    parser.add_argument("-v", "--verbose", default=0, type=int, help="Verbose debug.")
    args = parser.parse_args()

    PROT_MAX_SEQS_PER_LEVEL = args.sequences_per_level
    MINIMUM_SEQUENCE_LENGTH = args.minimum_sequence_length
    num_threads = args.processes

    allowed_filetypes = ["fa", "fas", "fasta"]

    rocksdb_folder_name = "rocksdb"
    sequences_db_name = "sequences"
    core_directory = "PhyMMR"
    secondary_directory = os.path.join(core_directory, os.path.basename(args.input))

    # Create necessary directories
    printv("Creating directories", args.verbose)

    os.makedirs(secondary_directory, exist_ok=True)

    taxa_runs = {}
    rev_comp_save = {}

    # Scan all the files in the input. Remove _R# and merge taxa
    for file in os.listdir(args.input):
        if (
            os.path.isfile(os.path.join(args.input, file))
            and file.split(".")[-1] in allowed_filetypes
        ):
            taxa = file.split(".")[0]

            formatted_taxa = truncate_taxa(taxa, extension=".fa")

            taxa_runs.setdefault(formatted_taxa, [])
            taxa_runs[formatted_taxa].append(file)

    for formatted_taxa_out, components in taxa_runs.items():
        taxa_start = time()
        printv(f"Preparing {formatted_taxa_out}", args.verbose)

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
            taxa_destination_directory, formatted_taxa_out.replace(".fa","_prot.fa")
        )

        duplicates = {}
        transcript_mapped_to = {}
        dupes = count()
        this_index = 1
        trim_times = []  # Append computed time for each loop.

        fa_file_out = []
        printv("Formatting input sequences and inserting into database", args.verbose)
        for file in components:
            if ".fa" not in file:
                continue

            fa_file_directory = os.path.join(args.input, file)
            fasta_file = SimpleFastaParser(open(fa_file_directory, encoding="UTF-8"))

            for header, parent_seq in tqdm(fasta_file) if args.verbose else fasta_file:
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
                        rev_seq = blosum_distance.bio_revcomp(seq)
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
        prot_start = time()

        if os.path.exists("/run/shm"):
            tmp_path = "/run/shm"
        elif os.path.exists("/dev/shm"):
            tmp_path = "/dev/shm"
        else:
            tmp_path = os.path.join(args.input, "tmp")
            os.makedirs(tmp_path, exist_ok=True)

        sequences_per_thread = math.ceil((len(fa_file_out) / 2) / num_threads) * 2
        translate_files = []
        for i in range(0,len(fa_file_out),sequences_per_thread):
            in_path = os.path.join(tmp_path,formatted_taxa_out+f'_{len(translate_files)}.fa')
            out_path = os.path.join(tmp_path,formatted_taxa_out+f'_prot_{len(translate_files)}.fa')
            translate_files.append((in_path, out_path))
            open(in_path, 'w').writelines(fa_file_out[i:i + sequences_per_thread])

        with Pool(num_threads) as translate_pool:
            translate_pool.starmap(translate, translate_files)  
            
        if args.keep_prepared:
            open(prepared_file_destination, "w", encoding="UTF-8").writelines(fa_file_out)
        del fa_file_out

        prot_components = []

        printv("Storing translated file in DB. Translate took {:.2f}s".format(time()-prot_start), args.verbose)

        out_lines = []
        aa_dedupe_time = time()
        for prepare_file, translate_file in translate_files:
            with open(translate_file, "r+", encoding="utf-8") as fp:
                for header, seq in SimpleFastaParser(fp):
                    header = header.replace(" ","|")
                    if seq in transcript_mapped_to:
                        duplicates.setdefault(transcript_mapped_to[seq], 1)
                        duplicates[transcript_mapped_to[seq]] += 1
                        next(aa_dupe_count)
                        continue
                    else:
                        transcript_mapped_to[seq] = header
                    out_lines.append(">"+header+"\n"+seq+"\n")
            os.remove(prepare_file)
            os.remove(translate_file)

        if args.keep_prepared:
            open(prot_path,'w').writelines(out_lines)
        aa_dupes = next(aa_dupe_count)
        printv("AA dedupe took {:.2f}s. Kicked {} dupes".format(time()-aa_dedupe_time, aa_dupes), args.verbose)

        levels = math.ceil(len(out_lines) / PROT_MAX_SEQS_PER_LEVEL)

        component = 0
        for i in range(0,len(out_lines),levels):
            component += 1

            data = out_lines[i: i+levels]
            prot_components.append(str(component))
            db.put(f"getprot:{component}", "".join(data))

        db.put("getall:prot", ",".join(prot_components))

        printv("Translation and storing done! Took {:.2f}s".format(time()-prot_start), args.verbose)

        sequence_count = this_index - 1

        printv(
            "Inserted {} sequences for {} found {} NT & {} AA dupes.\n".format(
                sequence_count, formatted_taxa, next(dupes), aa_dupes
            ),
            args.verbose,
        )

        del transcript_mapped_to  # Clear mem


        # Store the count of dupes in the database
        db.put("getall:dupes", json.dumps(duplicates))

        printv("Took {:.2f}s for {}".format(time() - taxa_start, file), args.verbose)

    print("Finished took {:.2f}s overall.".format(time() - global_start))
    print("N_trim time: {} seconds".format(sum(trim_times)))
    print(f"Dedupe time: {dedup_time}")


if __name__ == "__main__":
    main(argv)
