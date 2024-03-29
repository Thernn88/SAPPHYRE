from collections import Counter, defaultdict
import os
import re
from collections.abc import Callable, Generator
from itertools import count
from pathlib import Path
from queue import Queue
from shutil import rmtree
from tempfile import NamedTemporaryFile
from time import time
from typing import Any
from tqdm import tqdm

import sapphyre_tools
import wrap_rocks
import xxhash
from msgspec import json

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta


class IndexIter:
    """
    Utilize the itertools count function but with the ability to call the index without incrementing it
    """

    def __init__(self) -> None:
        self.counter = count(1)
        self.x = next(self.counter)

    def __next__(self):
        self.x = next(self.counter)


def truncate_taxa(taxa: str, extension=None) -> str:
    """Given a fasta header, checks the end for problematic tails.
    If found, truncates the string.
    Returns the string + suffix to check for name matches.
    """

    result = taxa.replace("_001", "")
    m = re.search(r"(_\d.fa)|(_R\d.fa)|(_part\d.fa)", result + extension)

    if m:
        tail_length = m.end() - m.start() - len(extension)
        result = result[0:-tail_length]
    if extension:
        result += extension

    return result


def group_taxa_in_glob(
    globbed: Generator[Path, Any, Any],
) -> dict[str, list[Path]]:
    """Filters and truncates fasta files for processing. Returns a dict of taxa to files."""
    ALLOWED_FILETYPES_NORMAL = [
        ".fa",
        ".fas",
        ".fasta",
        ".fq",
        ".fastq",
        ".fq",
        ".fsa_nt",
    ]
    ALLOWED_FILETYPES_GZ = [f"{ft}.gz" for ft in ALLOWED_FILETYPES_NORMAL]

    ALLOWED_FILETYPES = ALLOWED_FILETYPES_NORMAL + ALLOWED_FILETYPES_GZ
    taxa_runs = {}

    for f in globbed:
        file_is_file = f.is_file()

        file_name = f.name
        file_suffix_allowed = any(suff in file_name for suff in ALLOWED_FILETYPES)

        if file_is_file and file_suffix_allowed:
            taxa = f.stem.split(".")[0]
            formatted_taxa = truncate_taxa(taxa, extension=".fa")

            taxa_runs.setdefault(formatted_taxa, [])
            taxa_runs[formatted_taxa].append(f)

    return taxa_runs

class SeqDeduplicator:
    def __init__(
        self,
        minimum_sequence_length: int,
        verbose: int,
        overlap_length,
        rename,
    ) -> None:
        self.original_positions = {}
        self.original_inputs = []
        self.lines = []
        self.file_index = 0
        self.minimum_sequence_length = minimum_sequence_length
        self.verbose = verbose
        self.this_assembly = False
        self.this_genome = False
        self.overlap_length = overlap_length
        self.rename = rename
        self.transcript_mapped_to = defaultdict(dict)

    def __call__(
        self,
        fa_file_path: str,
        duplicates: dict[str, int],
        rev_comp_save: dict[str, int],
        dupes: count,
        this_index: IndexIter,
    ):
        CHOMP_LEN = 750
        CHOMP_CUTOFF = 50000
        ASSEMBLY_LEN = 750

        self.original_inputs.append(str(fa_file_path))

        if not self.this_genome:
            for _, seq in parseFasta(fa_file_path, True):
                if len(seq) > CHOMP_CUTOFF:
                    self.this_genome = True
                    self.this_assembly = False
                    break

                    

        for_loop = enumerate(parseFasta(fa_file_path, True))
        if self.verbose > 1:
            for_loop = list(for_loop)

        if self.verbose>1:
            for_loop = tqdm(for_loop)

        requires = False
        for line_index, (header, parent_seq) in for_loop:
            if line_index == 0:
                requires = any(l.islower() for l in parent_seq)
            if len(parent_seq) < self.minimum_sequence_length:
                continue

            if requires:
                parent_seq = parent_seq.upper()

            if not self.rename:
                n_sequences = [chunk for chunk in parent_seq.split("N") if len(chunk) >= self.minimum_sequence_length]
            else:
                n_sequences = (chunk for chunk in parent_seq.split("N") if len(chunk) >= self.minimum_sequence_length)
            individual_index = IndexIter()
            header_template = "NODE_{}"
            append_index_template = "{}_{}"
            sequence_template = ">{}\n{}\n"

            for seq in n_sequences:
                seq_len = len(seq)
                if self.rename:
                    header = header_template.format(this_index.x)
                else:
                    header = header.split(" ")[0]
                seq_hash = xxhash.xxh3_64(seq).hexdigest()

                this_hash_subset = self.transcript_mapped_to[seq_len]

                # Check for dupe, if so save how many times that sequence occured
                if seq_hash in this_hash_subset:
                    duplicates[this_hash_subset[seq_hash]] += 1
                    next(dupes)
                    continue
                this_hash_subset[seq_hash] = header

                # Rev-comp sequence. Save the reverse compliment in a hashmap with the original
                # sequence so we don't have to rev-comp this unique sequence again
                if seq_hash in rev_comp_save:
                    rev_seq_hash = rev_comp_save[seq_hash]
                else:
                    rev_seq_hash = xxhash.xxh3_64(
                        sapphyre_tools.bio_revcomp(seq),
                    ).hexdigest()
                    rev_comp_save[seq_hash] = rev_seq_hash

                # Check for revcomp dupe, if so save how many times that sequence occured
                if rev_seq_hash in this_hash_subset:
                    duplicates[this_hash_subset[rev_seq_hash]] += 1
                    next(dupes)
                    continue
                this_hash_subset[rev_seq_hash] = header

                if not self.this_genome:
                    if not self.this_assembly and len(parent_seq) >= ASSEMBLY_LEN:
                        self.this_assembly = True
                    this_header = header
                    if not self.rename and len(n_sequences) > 1:
                        this_header = append_index_template.format(header, individual_index)
                        next(individual_index)

                    self.original_positions[this_header] = (self.file_index, line_index)

                    self.lines.append(sequence_template.format(this_header, seq))
                    next(this_index)
                else:
                    individual_index = IndexIter()
                    for i in range(0, len(seq), CHOMP_LEN - self.overlap_length):
                        if self.rename:
                            this_header = header_template.format(this_index.x)
                        else:
                            this_header = append_index_template.format(header, individual_index.x)
                            next(individual_index)
                        if len(seq[i:i+CHOMP_LEN]) < self.minimum_sequence_length:
                            continue
                        self.lines.append(sequence_template.format(this_header, seq[i:i+CHOMP_LEN]))
                        next(this_index)

        self.file_index += 1

def map_taxa_runs(
    formatted_taxa_out,
    secondary_directory,
    verbose,
    overwrite,
    minimum_sequence_length,
    overlap_length,
    rename,
    components,
    keep_prepared,
    chunk_size,
    skip_entropy,
):
    """
    Removes duplicate sequences, renames them and inserts them into the taxa database.
    """
    ROCKSDB_FOLDER_NAME = "rocksdb"
    SEQUENCES_FOLDER_NAME = "sequences"
    NT_DB_NAME = "nt"

    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    printv(f"Preparing: {formatted_taxa_out}", verbose)

    taxa_destination_directory = secondary_directory.joinpath(formatted_taxa_out)

    rocksdb_path = taxa_destination_directory.joinpath(ROCKSDB_FOLDER_NAME)
    sequences_folder = rocksdb_path.joinpath(SEQUENCES_FOLDER_NAME)
    nt_db_path = sequences_folder.joinpath(NT_DB_NAME)

    if overwrite and rocksdb_path.exists():
        printv("Clearing old database", verbose)
        rmtree(rocksdb_path)

    printv("Creating rocksdb database", verbose)
    rocksdb_path.mkdir(parents=True, exist_ok=True)
    taxa_destination_directory.mkdir(parents=True, exist_ok=True)

    nt_db = wrap_rocks.RocksDB(str(nt_db_path))
    prepared_file_destination = taxa_destination_directory.joinpath(formatted_taxa_out)

    duplicates = Counter()
    rev_comp_save = {}
    dupes = count()
    this_index = IndexIter()

    printv(
        f"Got rocksdb database. Took {time_keeper.lap():.2f}s. Processing sequences",
        verbose,
    )

    deduper = SeqDeduplicator(
        minimum_sequence_length,
        verbose,
        overlap_length,
        rename,
    )
    for fa_file_path in components:
        deduper(
            fa_file_path,
            duplicates,
            rev_comp_save,
            dupes,
            this_index,
        )

    fa_file_out = deduper.lines
    this_is_assembly = deduper.this_assembly
    this_is_genome = deduper.this_genome
    original_positions = deduper.original_positions
    original_inputs = deduper.original_inputs
    del deduper
    if not skip_entropy:
        fa_file_out = sapphyre_tools.entropy_filter(fa_file_out, 0.7)
    recipe = []
    recipe_index = IndexIter()
    final = IndexIter()

    current_count = IndexIter()
    current_batch = []
    prepared_file = None
    if keep_prepared:
        prepared_file = prepared_file_destination.open("w")

    printv(
        f"Done! Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Inserting sequences into database.",
        verbose,
    )
    
    for line in fa_file_out:
        if not skip_entropy:
            header, seq = line
            line = f"{header}\n{seq}\n"
        
        next(final)
        current_batch.append(line)

        if prepared_file:
            prepared_file.write(line)

        next(current_count)

        if current_count.x > chunk_size:
            current_count = IndexIter()
            next(recipe_index)
            recipe.append(str(recipe_index.x))

            nt_db.put(f"ntbatch:{recipe_index.x}", "".join(current_batch))
            current_batch = []

    if current_batch:
        next(recipe_index)
        recipe.append(str(recipe_index.x))

        nt_db.put(f"ntbatch:{recipe_index.x}", "".join(current_batch))

    if prepared_file:
        prepared_file.close()

    recipe_data = ",".join(recipe)
    nt_db.put("getall:batches", recipe_data)

    nt_db.put("get:isassembly", str(this_is_assembly))
    nt_db.put("get:isgenome", str(this_is_genome))

    printv(
        f"Inserted {this_index.x} sequences. Found {next(dupes)} duplicates. Entropy removed {this_index.x - final.x}. Took {time_keeper.lap():.2f}s.",
        verbose,
    )

    # Store the original positions
    if not (this_is_assembly or this_is_genome):
        nt_db.put_bytes("getall:original_positions", json.encode(original_positions))
        nt_db.put_bytes("getall:original_inputs", json.encode(original_inputs))

    # Store the count of dupes in the database
    nt_db.put_bytes("getall:dupes", json.encode(duplicates))

    printv(
        f"Done! {formatted_taxa_out} took {time_keeper.differential():.2f}s overall.",
        verbose,
    )


def main(args):
    input_path = Path(args.INPUT)

    if not input_path.exists():
        printv("ERROR: An existing directory must be provided.", args.verbose, 0)
        return False

    minimum_sequence_length = args.minimum_sequence_length

    core_directory = Path(args.out).joinpath(input_path.parts[-1])
    secondary_directory = Path(os.getcwd()).joinpath(core_directory)

    # Create necessary directories
    printv("Creating directories", args.verbose)
    secondary_directory.mkdir(parents=True, exist_ok=True)

    global_time_keeper = TimeKeeper(KeeperMode.DIRECT)

    globbed = list(input_path.glob("*"))
    globbed.sort()
    taxa_runs = group_taxa_in_glob(globbed)

    for formatted_taxa_out, components in taxa_runs.items():
        map_taxa_runs(
            formatted_taxa_out,
            secondary_directory,
            args.verbose,
            args.overwrite,
            minimum_sequence_length,
            args.overlap_length,
            not args.no_rename,
            components,
            args.keep_prepared,
            args.chunk_size,
            args.skip_entropy,
        )

    printv(
        f"Finished! Took {global_time_keeper.differential():.2f}s overall.",
        args.verbose,
        0,
    )

    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Prepare"
    raise Exception(
        msg,
    )
