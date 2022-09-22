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


def truncate_taxa(header: str, extension=None) -> str:
    """
    Given a fasta header, checks the end for problematic tails.
    If found, truncates the string.
    Returns the string + suffix to check for name matches.
    """
    # search for _# and _R#, where # is digits
    result = header
    m = re.search(r"_R?\d+$", header)
    if m is not None:
        tail_length = m.end()-m.start()
        result = result[0:-tail_length]
    if extension is not None:
        result = result + extension
    return result


def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def main(argv):
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
        "-m",
        "--max_prepared_batch_size",
        default=1000000, 
        type=int, 
        help="Max sequences per prepared batch in db. Default: 1 million.",
    )
    parser.add_argument(
        "-v", 
        "--verbose", 
        default=0, 
        type=int, 
        help="Verbose debug.")
    args = parser.parse_args()

    MAX_PREPARED_LEVEL_SIZE = args.max_prepared_batch_size

    allowed_filetypes = ['fa', 'fas', 'fasta']

    rocksdb_folder_name = "rocksdb"
    core_directory = "PhyMMR"
    secondary_directory = os.path.join(core_directory, os.path.basename(args.input))

    # Create necessary directories
    if args.verbose != 0:
        print("Creating directories")

    if not os.path.exists(core_directory):
        os.mkdir(core_directory)

    if not os.path.exists(secondary_directory):
        os.mkdir(secondary_directory)

    taxa_runs = {}
    rev_comp_save = {}

    # Scan all the files in the input. Remove _R# and merge taxa
    for file in os.listdir(args.input):
        if os.path.isfile(os.path.join(args.input, file)) and file.split('.')[-1] in allowed_filetypes:
            taxa = file.split('.')[0]
            # if taxa[-1].isnumeric() and taxa[-2] == 'R' and taxa[-3] == '_': # Contains "_R#"
            #     formatted_taxa = taxa[:-3]+'.fa'
            # else:
            #     formatted_taxa = file
            formatted_taxa = truncate_taxa(taxa, extension='.fa')

            if formatted_taxa in taxa_runs:
                taxa_runs[formatted_taxa].append(file)
            else:
                taxa_runs[formatted_taxa] = [file]
    
    for formatted_taxa_out, components in taxa_runs.items():
        taxa_start = time()
        if args.verbose != 0:
            print(f"Preparing {formatted_taxa_out}")

        taxa_destination_directory = os.path.join(secondary_directory, formatted_taxa_out)
        rocksdb_path = os.path.join(taxa_destination_directory, rocksdb_folder_name)

        if args.clear_database:
            if os.path.exists(rocksdb_path):
                if args.verbose != 0:
                    print("Clearing old database")
                rmtree(rocksdb_path)

        if args.verbose >= 1:
            print("Creating rocksdb database")

        db = wrap_rocks.RocksDB(rocksdb_path)

        if not os.path.exists(taxa_destination_directory):
                os.mkdir(taxa_destination_directory)

        prepared_file_destination = os.path.join(taxa_destination_directory, formatted_taxa_out)
        
        duplicates = {}
        dupe_set = set()
        dupes = 0
        prepared_component_all = []
        sequence_count = 0

        if args.verbose >= 1:
            print("Formatting input sequences and inserting into database")
        with open(
            prepared_file_destination, "w", encoding="UTF-8"
        ) as fa_file_out:
            for file in components:
                fa_file_directory = os.path.join(args.input, file)
                if ".fa" in file:
                    with open(fa_file_directory, encoding="UTF-8") as fa_file_in:
                        lines = fa_file_in.readlines()

                        if lines[-1] == "\n": lines = lines[:-1]

                        sequence_count += int(len(lines) / 2)

                        for_loop_range = (
                            tqdm(range(0, len(lines), 2))
                            if args.verbose != 0
                            else range(0, len(lines), 2)
                        )
                        
                        for i in for_loop_range:
                            header = lines[i].strip()
                            seq = lines[i+1].strip()

                            # Check for dupe, if so save how many times that sequence occured
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
                            
                            # If no dupe, write to prepared file and db
                            line = header+'\n'+seq+'\n'

                            fa_file_out.write(line)

                            # Save prepared file lines in a list to ration into the db
                            prepared_component_all.append(line)

                            # Get rid of space and > in header (blast/hmmer doesn't like it)
                            preheader = header.replace(" ", "|").replace(">", "") # pre-hash header
                            
                            # Data that will be stored in the database
                            data = f"{preheader}\n{seq}"

                            # Hash the header
                            header = hashlib.sha256(preheader.encode()).hexdigest()

                            # Write to rocksdb
                            db.put(header, data)

        if args.verbose != 0:
            print(
                "Inserted {} sequences for {} found {} dupes.\n".format(sequence_count-dupes, formatted_taxa, dupes)
            )

        del dupe_set # Clear mem

        if args.verbose != 0:
            print("Storing prepared file in database")
        
        # Figure out how many levels are needed to split sequence_count into MAX_PREPARED_LEVEL_SIZE sized chunks 
        prepared_file_levels = math.ceil(sequence_count / MAX_PREPARED_LEVEL_SIZE)
        sequences_per_level = math.ceil(sequence_count / prepared_file_levels)
        prepared_recipe = []

        # Grab the chunks of the prepared file lines and store in db
        for i in range(0, prepared_file_levels):
            start_index = i * sequences_per_level
            end_index = (i + 1) * sequences_per_level

            this_batch = "".join(prepared_component_all[start_index:end_index])
            key = f"prepared:{i}"

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

        if args.verbose != 0:
            print(
                "Took {:.2f}s for {}".format(time() - taxa_start, file)
            )

    print("Finished took {:.2f}s overall.".format(time() - global_start))


if __name__ == "__main__":
    main(argv)
