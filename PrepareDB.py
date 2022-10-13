import argparse
import json
import math
import os
import re
from itertools import count
from shutil import rmtree
from sys import argv
from time import time

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
        "-m",
        "--max_prepared_batch_size",
        default=100000,
        type=int,
        help="Max sequences per prepared batch in db. Default: 100 thousand.",
    )
    parser.add_argument(
        "-k",
        "--keep_prepared",
        action="store_true",
        help="Writes the prepared input fasta into the output taxa directory.",
    )
    parser.add_argument("-v", "--verbose", default=0, type=int, help="Verbose debug.")
    args = parser.parse_args()

    MAX_PREPARED_LEVEL_SIZE = args.max_prepared_batch_size
    MINIMUM_SEQUENCE_LENGTH = args.minimum_sequence_length

    allowed_filetypes = ["fa", "fas", "fasta"]

    rocksdb_folder_name = "rocksdb"
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

        if args.clear_database and os.path.exists(rocksdb_path):
            printv("Clearing old database", args.verbose)
            rmtree(rocksdb_path)

        printv("Creating rocksdb database", args.verbose)

        db = wrap_rocks.RocksDB(rocksdb_path)

        os.makedirs(taxa_destination_directory, exist_ok=True)

        prepared_file_destination = os.path.join(
            taxa_destination_directory, formatted_taxa_out
        )

        duplicates = {}
        transcript_mapped_to = {}
        dupes = (count())
        prepared_component_all = []
        this_index = 1
        component_i = (count())
        prepared_recipe = []
        trim_times = []  # Append computed time for each loop.

        if args.keep_prepared:
            fa_file_out = open(prepared_file_destination, "w", encoding="UTF-8")
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

                    if args.keep_prepared:
                        fa_file_out.write(line)

                    # Save prepared file lines in a list to ration into the db
                    prepared_component_all.append(line)

                    if len(prepared_component_all) >= MAX_PREPARED_LEVEL_SIZE:
                        key = add_pc_to_db(
                            db, next(component_i), prepared_component_all
                        )
                        prepared_recipe.append(key)
                        prepared_component_all = []

                    # Get rid of space and > in header (blast/hmmer doesn't like it) Need to push modified external to remove this. ToDo.
                    preheader = header.replace(" ", "|")  # pre-hash header
                    # Key is a hash of the header
                    db.put(xxhash.xxh64_hexdigest(preheader), f"{preheader}\n{seq}")

        if args.keep_prepared:
            fa_file_out.close()
        sequence_count = this_index - 1

        printv(
            "Inserted {} sequences for {} found {} dupes.\n".format(
                sequence_count, formatted_taxa, next(dupes)
            ),
            args.verbose,
        )

        del transcript_mapped_to  # Clear mem

        printv("Storing prepared file in database", args.verbose)

        if len(prepared_component_all) != 0:
            key = add_pc_to_db(
                db, next(component_i), prepared_component_all
            )
            prepared_component_all = []
            prepared_recipe.append(key)            

        # Store the keys to each component of the prepared file
        db.put("getall:prepared", ",".join(prepared_recipe))

        # Store the count of dupes in the database
        db.put("getall:dupes", json.dumps(duplicates))

        printv("Took {:.2f}s for {}".format(time() - taxa_start, file), args.verbose)

    print("Finished took {:.2f}s overall.".format(time() - global_start))
    print("N_trim time: {} seconds".format(sum(trim_times)))
    print(f"Dedupe time: {dedup_time}")


if __name__ == "__main__":
    main(argv)
