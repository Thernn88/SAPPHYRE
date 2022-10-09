import argparse
from sys import argv
import os
import wrap_rocks
import json
import re
from time import time
from shutil import rmtree
import hashlib
from tqdm import tqdm
import math
#from Bio.Seq import Seq
import blosum_distance


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
        tail_length = m.end()-m.start()
        result = result[0:-tail_length]
    if extension:
        result = result + extension
    return result


def rev_comp(seq):
    return blosum_distance.bio_revcomp(seq)

def N_trim(parent_sequence, MINIMUM_SEQUENCE_LENGTH):
    if 'N' in parent_sequence:
        # Get N indices and start and end of sequence
        indices = [0] + [i for i, ltr in enumerate(parent_sequence) if ltr == "N"] + [len(parent_sequence)]

        for i in range(0, len(indices)-1):
            start = indices[i]
            end = indices[i+1]

            length = end - start
            if length >= MINIMUM_SEQUENCE_LENGTH:
                # raw_seq = parent_sequence[start:end].strip('N')
                # yield raw_seq
                raw_seq = parent_sequence[start+1:end]
                yield raw_seq
    else:
        yield parent_sequence

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
    parser.add_argument(
        "-v", 
        "--verbose", 
        default=0, 
        type=int, 
        help="Verbose debug.")
    args = parser.parse_args()

    MAX_PREPARED_LEVEL_SIZE = args.max_prepared_batch_size
    MINIMUM_SEQUENCE_LENGTH = args.minimum_sequence_length

    allowed_filetypes = ['fa', 'fas', 'fasta']

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
        if os.path.isfile(os.path.join(args.input, file)) and file.split('.')[-1] in allowed_filetypes:
            taxa = file.split('.')[0]

            formatted_taxa = truncate_taxa(taxa, extension='.fa')

            taxa_runs.setdefault(formatted_taxa, [])
            taxa_runs[formatted_taxa].append(file)
    
    for formatted_taxa_out, components in taxa_runs.items():
        taxa_start = time()
        printv(f"Preparing {formatted_taxa_out}", args.verbose)

        taxa_destination_directory = os.path.join(secondary_directory, formatted_taxa_out)
        rocksdb_path = os.path.join(taxa_destination_directory, rocksdb_folder_name)

        if args.clear_database and os.path.exists(rocksdb_path):
            printv("Clearing old database", args.verbose)
            rmtree(rocksdb_path)

        printv("Creating rocksdb database", args.verbose)

        db = wrap_rocks.RocksDB(rocksdb_path)

        os.makedirs(taxa_destination_directory, exist_ok=True)

        prepared_file_destination = os.path.join(taxa_destination_directory, formatted_taxa_out)
        
        duplicates = {}
        dupe_set = set()
        dupes = 0
        prepared_component_all = []
        this_index = 1
        component_i = 0
        prepared_recipe = []

        if args.keep_prepared:
            fa_file_out = open(prepared_file_destination, "w", encoding="UTF-8")
        printv("Formatting input sequences and inserting into database", args.verbose)
        for file in components:
            fa_file_directory = os.path.join(args.input, file)
            if ".fa" in file:
                with open(fa_file_directory, encoding="UTF-8") as fa_file_in:
                    lines = fa_file_in.readlines()

                    if lines[-1] == "\n": lines = lines[:-1]

                    for_loop_range = (
                        tqdm(range(0, len(lines), 2))
                        if args.verbose
                        else range(0, len(lines), 2)
                    )
                    
                    for i in for_loop_range:
                        parent_seq = lines[i+1].strip()

                        if len(parent_seq) >= MINIMUM_SEQUENCE_LENGTH:
                            trim_start = time()
                            gen_object = N_trim(parent_seq, MINIMUM_SEQUENCE_LENGTH)
                            trim_end = time()
                            trim_time = trim_end - trim_start
                            # for seq in N_trim(parent_seq, MINIMUM_SEQUENCE_LENGTH):
                            for seq in gen_object:
                                # Check for dupe, if so save how many times that sequence occured
                                seq_start = time()
                                seq_hash = hash(seq)
                                if seq_hash in dupe_set:
                                    if seq_hash in duplicates:
                                        duplicates[seq_hash] += 1
                                    else:
                                        duplicates[seq_hash] = 2
                                    dupes += 1
                                    continue
                                else:
                                    dupe_set.add(seq_hash)

                                # Rev-comp sequence. Save the reverse compliment in a hashmap with the original
                                # sequence so we don't have to rev-comp this unique sequence again
                                if seq in rev_comp_save:
                                    rev_seq = rev_comp_save[seq]
                                else:
                                    rev_seq = rev_comp(seq)
                                    rev_comp_save[seq] = rev_seq

                                # Check for revcomp dupe, if so save how many times that sequence occured
                                seq_hash = hash(rev_seq)
                                if seq_hash in dupe_set:
                                    if seq_hash in duplicates:
                                        duplicates[seq_hash] += 1
                                    else:
                                        duplicates[seq_hash] = 2
                                    dupes += 1
                                    continue
                                else:
                                    dupe_set.add(seq_hash)
                                

                                length = len(seq)
                                header = f">NODE_{this_index}_length_{length}"
                                seq_end = time()
                                dedup_time += seq_end - seq_start
                                this_index += 1
                                
                                # If no dupe, write to prepared file and db
                                line = header+'\n'+seq+'\n'

                                if args.keep_prepared is True:
                                    fa_file_out.write(line)

                                # Save prepared file lines in a list to ration into the db
                                prepared_component_all.append(line)

                                if len(prepared_component_all) >= MAX_PREPARED_LEVEL_SIZE:
                                    this_batch = "".join(prepared_component_all)

                                    del prepared_component_all
                                    prepared_component_all = []
                                    component_i += 1

                                    key = f"prepared:{component_i}"

                                    db.put(key, this_batch)

                                    prepared_recipe.append(key)

                                # Get rid of space and > in header (blast/hmmer doesn't like it) Need to push modified external to remove this. ToDo.
                                preheader = header.replace(" ", "|").replace(">", "") # pre-hash header
                                
                                # Data that will be stored in the database
                                data = f"{preheader}\n{seq}"

                                # Hash the header
                                key = hashlib.sha256(preheader.encode()).hexdigest()

                                # Write to rocksdb
                                db.put(key, data)

        if args.keep_prepared:
            fa_file_out.close()
        sequence_count = (this_index - 1)

        printv(
            "Inserted {} sequences for {} found {} dupes.\n".format(sequence_count, formatted_taxa, dupes),
            args.verbose
        )

        del dupe_set # Clear mem

        printv("Storing prepared file in database", args.verbose)

        if len(prepared_component_all) != 0:
            this_batch = "".join(prepared_component_all)
            prepared_component_all = []
            key = f"prepared:{component_i}"

            db.put(key, this_batch)

            prepared_recipe.append(key)

        # Store the keys to each component of the prepared file
        key = "getall:prepared"
        data = ",".join(prepared_recipe)

        db.put(key, data)

        # Store the count of dupes in the database
        key = "getall:dupes"
        data = json.dumps(duplicates)

        db.put(key, data)

        printv("Took {:.2f}s for {}".format(time() - taxa_start, file), args.verbose)

    print("Finished took {:.2f}s overall.".format(time() - global_start))
    print(f"N_trim time: {trim_time}")
    print(f"Dedupe time: {dedup_time}")


if __name__ == "__main__":
    main(argv)
