import argparse
from sys import argv
import os
import wrap_rocks
import json
from time import time
from shutil import rmtree
import hashlib
from tqdm import tqdm
import math

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
        "-k",
        "--keep_prepared",
        action="store_true",
        help="Keep prepared input file in new directory.",
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

    rocksdb_folder_name = "rocksdb"
    core_directory = "PhyMMR"
    secondary_directory = os.path.join(core_directory, os.path.basename(args.input))

    if args.verbose != 0:
        print("Creating directories")

    if not os.path.exists(core_directory):
        os.mkdir(core_directory)

    if not os.path.exists(secondary_directory):
        os.mkdir(secondary_directory)

    for file in os.listdir(args.input):
        if ".fa" in file:
            taxa_start = time()
            if args.verbose != 0:
                print(f"Preparing {file}")

            fa_file_directory = os.path.join(args.input, file)
            taxa_destination_directory = os.path.join(secondary_directory, file)

            if not os.path.exists(taxa_destination_directory):
                os.mkdir(taxa_destination_directory)

            prepared_file_destination = os.path.join(taxa_destination_directory, file)

            if not os.path.exists(prepared_file_destination):
                with open(
                    prepared_file_destination, "w", encoding="UTF-8"
                ) as fa_file_out:
                    with open(fa_file_directory, encoding="UTF-8") as fa_file_in:
                        for line in fa_file_in:
                            line = line.strip()
                            if line != "":
                                fa_file_out.write(line + "\n")

            rocksdb_path = os.path.join(taxa_destination_directory, rocksdb_folder_name)

            if args.clear_database:
                if os.path.exists(rocksdb_path):
                    if args.verbose != 0:
                        print("Clearing old database")
                    rmtree(rocksdb_path)

            if args.verbose >= 1:
                print("Creating rocksdb database")

            db = wrap_rocks.RocksDB(rocksdb_path)

            prepared_component_all = []
            duplicates = {}
            dupe_set = set()

            if args.verbose != 0:
                print("Inserting sequences into database")

            with open(prepared_file_destination) as file_in:
                lines = file_in.read().split("\n")

            if lines[-1] == "":
                lines = lines[:-1]

            sequence_count = int(len(lines)/2)

            for_loop_range = (
                tqdm(range(0, len(lines), 2))
                if args.verbose != 0
                else range(0, len(lines), 2)
            )

            for i in for_loop_range:
                header = lines[i]
                sequence = lines[i + 1]

                prepared_component_all.append(header+'\n'+sequence+'\n')

                this_hash = hash(sequence)
                if this_hash not in dupe_set:
                    dupe_set.add(this_hash)
                else:
                    if sequence in duplicates:
                        duplicates[sequence] += 1
                    else:
                        duplicates[sequence] = 2

                # get rid of space and > in header (blast/hmmer doesn't like it)
                preheader = header.replace(" ", "|").replace(">", "") # pre-hash header
                
                data = f"{preheader}\n{sequence}"

                header = hashlib.sha256(preheader.encode()).hexdigest()
                db.put(header, data)

            if args.verbose != 0:

                print(
                    "Inserted {} sequences for {}.\n".format(sequence_count, file)
                )

            del dupe_set # Clear mem

            if args.verbose != 0:
                print("Storing prepared file in database")

            # store the keys to each component of the prepared file
            
            prepared_file_levels = math.ceil(sequence_count / MAX_PREPARED_LEVEL_SIZE)
            sequences_per_level = math.ceil(sequence_count / prepared_file_levels)
            prepared_recipe = []

            for i in range(0, prepared_file_levels):
                start_index = i * sequences_per_level
                end_index = (i + 1) * sequences_per_level

                this_batch = "".join(prepared_component_all[start_index:end_index])
                key = f"prepared:{i}"

                db.put(key, this_batch)

                prepared_recipe.append(key)

            key = "getall:prepared"
            data = ",".join(prepared_recipe)

            db.put(key, data)

            # store the count of dupes in the database
            key = "getall:dupes"
            data = json.dumps(duplicates)

            db.put(key, data)

            if args.verbose != 0:
                print(
                    "Cleaning up. Took {:.2f}s for {}".format(time() - taxa_start, file)
                )

            if args.keep_prepared is not True:
                # We can delete these as it's stored in memory for future use if needed
                os.remove(prepared_file_destination)

    print("Finished took {:.2f}s overall.".format(time() - global_start))


if __name__ == "__main__":
    main(argv)
