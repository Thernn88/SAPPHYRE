from __future__ import annotations

import json
import os
import re

from itertools import count
from pathlib import Path
from queue import Queue
from shutil import rmtree
import subprocess
from tempfile import NamedTemporaryFile
from time import time
from typing import Any, Callable, Dict, Generator, List, Tuple

import phymmr_tools
import wrap_rocks
import xxhash
from tqdm import tqdm

from .timekeeper import TimeKeeper, KeeperMode
from .utils import ConcurrentLogger, gettempdir, parseFasta, writeFasta

ROCKSDB_FOLDER_NAME = "rocksdb"
SEQUENCES_FOLDER_NAME = "sequences"
NT_DB_NAME = "nt"
CORE_FOLDER = "PhyMMR"
ALLOWED_FILETYPES_NORMAL = [
    ".fa",
    ".fas",
    ".fasta",
    ".fq",
    ".fastq",
    ".fq",
]
ALLOWED_FILETYPES_GZ = [f"{ft}.gz" for ft in ALLOWED_FILETYPES_NORMAL]

ALLOWED_FILETYPES = ALLOWED_FILETYPES_NORMAL + ALLOWED_FILETYPES_GZ


class IndexIter:
    def __init__(self):
        self.counter = count(1)
        self.x = next(self.counter)

    def __str__(self):
        return str(self.x)

    def __next__(self):
        self.x = next(self.counter)


def truncate_taxa(header: str, extension=None) -> str:
    """
    Given a fasta header, checks the end for problematic tails.
    If found, truncates the string.
    Returns the string + suffix to check for name matches.
    """
    # search for _# and _R#, where # is digits
    result = header
    m = re.search(r"(_L\d\d\d){0,1}(_R\d+)(_\d\d\d){0,1}$", header)
    if m:
        tail_length = m.end() - m.start()
        result = result[0:-tail_length]
    if extension:
        result = result + extension
    return result


def N_trim(parent_sequence: str, minimum_sequence_length: int, tike: TimeKeeper):
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
                raw_seq = parent_sequence[start + 1 : end]
                tike[0].timer1("now")
                yield raw_seq
    else:
        yield parent_sequence


def glob_for_fasta_and_save_for_runs(
    globbed: Generator[Path, Any, Any]
) -> Dict[str, List[Path]]:
    """
    Scan all the files in the input. Remove _R# and merge taxa
    """

    taxa_runs = {}

    for f in globbed:
        file_is_file = f.is_file()

        file_name = f.name
        file_suffix_allowed = any([suff in file_name for suff in ALLOWED_FILETYPES])

        if file_is_file and file_suffix_allowed:
            taxa = f.stem.split(".")[0]

            formatted_taxa = truncate_taxa(taxa, extension=".fa")

            taxa_runs.setdefault(formatted_taxa, [])
            taxa_runs[formatted_taxa].append(f)

    return taxa_runs


class SeqDeduplicator:
    def __init__(
        self, fa_file_path: Path, temp_path: str, db: Any, minimum_sequence_length: int, verbose: int
    ):
        self.fa_file_path = fa_file_path
        self.minimum_sequence_length = minimum_sequence_length
        self.verbose = verbose
        self.nt_db = db
        self.recipe = []
        self.temp_path = temp_path

    def __call__(
        self,
        trim_times: TimeKeeper,
        duplicates: Dict[str, int],
        rev_comp_save: Dict[str, int],
        transcript_mapped_to: Dict[str, str],
        dupes: count[int],
        this_index: IndexIter,
        fa_file_out: List[str],
        dedup_time: List[int],
    ):
        fasta_file = parseFasta(self.fa_file_path)

        for_loop = tqdm(fasta_file) if self.verbose else fasta_file

        with open(self.temp_path, "w") as fp:
            for row in for_loop:
                header, parent_seq = row[:2]

                if not len(parent_seq) >= self.minimum_sequence_length:
                    continue
                parent_seq = parent_seq.upper()

                for seq in N_trim(parent_seq, self.minimum_sequence_length, trim_times):
                    header = f"NODE_{this_index}"
                    seq_hash = xxhash.xxh64(seq).hexdigest()

                    # Check for dupe, if so save how many times that sequence occured
                    seq_start = time()

                    if seq_hash in transcript_mapped_to:
                        duplicates.setdefault(transcript_mapped_to[seq_hash], 1)
                        duplicates[transcript_mapped_to[seq_hash]] += 1
                        next(dupes)
                        continue
                    transcript_mapped_to[seq_hash] = header

                    # Rev-comp sequence. Save the reverse compliment in a hashmap with the original
                    # sequence so we don't have to rev-comp this unique sequence again
                    if seq_hash in rev_comp_save:
                        rev_seq_hash = rev_comp_save[seq_hash]
                    else:
                        rev_seq_hash = xxhash.xxh64(
                            phymmr_tools.bio_revcomp(seq)
                        ).hexdigest()
                        rev_comp_save[seq_hash] = rev_seq_hash

                    # Check for revcomp dupe, if so save how many times that sequence occured
                    if rev_seq_hash in transcript_mapped_to:
                        duplicates.setdefault(transcript_mapped_to[rev_seq_hash], 1)
                        duplicates[transcript_mapped_to[rev_seq_hash]] += 1
                        next(dupes)
                        continue
                    transcript_mapped_to[rev_seq_hash] = header

                    seq_end = time()
                    dedup_time[0] += seq_end - seq_start

                    # If no dupe, write to prepared file and db
                    next(this_index)

                    fa_file_out.append((header, seq))

                    fp.write(f">{header}\n{seq}\n")

def return_clean_lines(clean_file_path, batch_size, nt_db, keep_prepared, path):
    recipe_index = 0
    current_count = IndexIter()
    recipe = []
    current_batch = []

    if keep_prepared:
        fp = open(path, "w")

    for header, sequence in parseFasta(clean_file_path):
        line = f">{header}\n{sequence}\n"

        if keep_prepared:
            fp.write(line)

        current_batch.append(line)
        next(current_count)
        if current_count.x >= batch_size:
            current_count = IndexIter()
            recipe.append(str(recipe_index))
            nt_db.put(f"ntbatch:{recipe_index}", "".join(current_batch))
            recipe_index += 1
            current_batch = []
    
    if keep_prepared:
        fp.close()

    if current_batch:
        recipe.append(str(recipe_index))
        nt_db.put(f"ntbatch:{recipe_index}", "".join(current_batch))
        recipe_index += 1

    recipe_data = ",".join(recipe)
    nt_db.put("getall:batches", recipe_data)

def return_dupe_table(dupe_file_path, nt_db):
    dupe_table = {}
    dupe_count = IndexIter()
    for header, sequence in parseFasta(dupe_file_path):
        if header.endswith('.'):
            master = header[:-1]
        else:
            next(dupe_count)
            dupe_table.setdefault(master, []).append((header, sequence)) 

    nt_db.put("getall:hamm_dupes", json.dumps(dupe_table))
    return dupe_count.x

def hamming_dedupe(fa_file_dest, batch_size, threads, nt_db, keep_prepared, prepared_file_destination):
    with NamedTemporaryFile(dir = gettempdir()) as clean_file, NamedTemporaryFile(dir = gettempdir()) as dupe_file:
        os.system(f"minirmd/minirmd -i {fa_file_dest} -o {clean_file.name} -l {dupe_file.name} -t {threads} -d 1")

        return_clean_lines(clean_file.name, batch_size, nt_db, keep_prepared, prepared_file_destination)
        return return_dupe_table(dupe_file.name, nt_db)

class DatabasePreparer:
    def __init__(
        self,
        formatted_taxa_out: str,
        components: List[Path],
        verbose: bool,
        printv: Callable,
        clear_database: bool,
        keep_prepared: bool,
        minimum_sequence_length: int,
        dedup_time: List[int],
        trim_times: TimeKeeper,
        chunk_size: int,
        threads,
    ):
        self.fto = formatted_taxa_out
        self.comp = components
        self.trim_times = trim_times
        self.fa_file_out = []
        self.verbose = verbose
        self.printv = printv
        self.clear_database = clear_database
        self.minimum_sequence_length = minimum_sequence_length
        self.keep_prepared = keep_prepared
        self.duplicates = {}
        self.rev_comp_save = {}
        self.transcript_mapped_to = {}
        self.dupes = count()
        self.this_index = IndexIter()
        self.dedup_time = dedup_time
        self.chunk_size = chunk_size
        self.threads = threads

    def init_db(
        self,
        secondary_directory: Path,
    ) -> Any:
        self.taxa_time_keeper = TimeKeeper(KeeperMode.DIRECT)
        self.printv(f"Preparing: {self.fto}", self.verbose, 0)
        self.printv(
            "Formatting input sequences and inserting into database", self.verbose
        )

        taxa_destination_directory = secondary_directory.joinpath(self.fto)

        rocksdb_path = taxa_destination_directory.joinpath(ROCKSDB_FOLDER_NAME)
        sequences_folder = rocksdb_path.joinpath(SEQUENCES_FOLDER_NAME)
        nt_db_path = sequences_folder.joinpath(NT_DB_NAME)

        if self.clear_database and rocksdb_path.exists():
            self.printv("Clearing old database", self.verbose)
            rmtree(rocksdb_path)

        self.printv("Creating rocksdb database", self.verbose)
        rocksdb_path.mkdir(parents=True, exist_ok=True)
        self.nt_db = wrap_rocks.RocksDB(str(nt_db_path))
        taxa_destination_directory.mkdir(parents=True, exist_ok=True)
        self.prepared_file_destination = taxa_destination_directory.joinpath(self.fto)

    def dedup(self):
        with NamedTemporaryFile(dir = gettempdir()) as temp:
            for fa_file_path in self.comp:
                SeqDeduplicator(
                    fa_file_path, temp.name, self.nt_db, self.minimum_sequence_length, self.verbose
                )(
                    self.trim_times,
                    self.duplicates,
                    self.rev_comp_save,
                    self.transcript_mapped_to,
                    self.dupes,
                    self.this_index,
                    self.fa_file_out,
                    self.dedup_time,
                )
            temp.flush()

            hamm_dupes = hamming_dedupe(temp.name, self.chunk_size, self.threads, self.nt_db, self.keep_prepared, self.prepared_file_destination)

        if self.keep_prepared:
            writeFasta(self.prepared_file_destination, self.fa_file_out)

        sequence_count = str(self.this_index)

        self.printv(
            f"Inserted {sequence_count} sequences. Found {next(self.dupes)} exact duplicates. Found {hamm_dupes} hamm duplicates.",
            self.verbose,
        )

        # Store the count of dupes in the database
        self.nt_db.put("getall:dupes", json.dumps(self.duplicates))

        self.printv(
            f"{self.fto} took {self.taxa_time_keeper.differential():.2f}s overall\n",
            self.verbose,
        )

    def __call__(
        self,
        secondary_directory: Path,
    ):
        self.init_db(secondary_directory)
        self.dedup()


def map_taxa_runs(
    tuple_in: Tuple[str, str],
    verbose: int,
    printv: Callable,
    clear_database: int,
    keep_prepared: int,
    minimum_sequence_length: int,
    secondary_directory: Path,
    dedup_time: List[int],
    trim_times: TimeKeeper,
    chunk_size,
    threads,
):
    formatted_taxa_out, components = tuple_in

    DatabasePreparer(
        formatted_taxa_out,
        components,
        verbose,
        printv,
        clear_database,
        keep_prepared,
        minimum_sequence_length,
        dedup_time,
        trim_times,
        chunk_size,
        threads,
    )(
        secondary_directory,
    )


def main(args):
    msgq = Queue()
    printv = ConcurrentLogger(msgq)
    printv.start()

    input_path = Path(args.INPUT)

    if not input_path.exists():
        printv("ERROR: An existing directory must be provided.", args.verbose, 0)
        return False

    minimum_sequence_length = args.minimum_sequence_length

    project_root = Path(__file__).parent.parent
    core_directory = Path(CORE_FOLDER).joinpath(input_path.parts[-1])
    secondary_directory = project_root.joinpath(core_directory)

    # Create necessary directories
    printv("Creating directories", args.verbose)
    secondary_directory.mkdir(parents=True, exist_ok=True)

    trim_times = [TimeKeeper(KeeperMode.SUM)]  # Append computed time for each loop.
    dedup_time = [0]
    global_time_keeper = TimeKeeper(KeeperMode.DIRECT)

    globbed = list(input_path.glob("*.f[aq]*"))
    globbed.sort()
    taxa_runs = glob_for_fasta_and_save_for_runs(globbed)

    ls_args = [
        (
            tuple_in,
            args.verbose,
            printv,
            args.clear_database,
            args.keep_prepared,
            minimum_sequence_length,
            secondary_directory,
            dedup_time,
            trim_times,
            args.chunk_size,
            args.processes
        )
        for tuple_in in taxa_runs.items()
    ]

    list(map(lambda x: map_taxa_runs(x[0], *x[1:]), ls_args))

    printv(
        f"Finished! Took {global_time_keeper.differential():.2f}s overall.",
        args.verbose,
        0,
    )
    printv(f"N_trim time: {trim_times[0].time1} seconds", args.verbose, 2)
    printv(f"Dedupe time: {dedup_time[0]}", args.verbose, 2)

    msgq.join()
    return True


if __name__ == "__main__":
    raise Exception("Cannot be called directly, please use the module:\nphymmr Prepare")