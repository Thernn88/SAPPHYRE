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
from subprocess import call
from time import time
from typing import Any, Callable, Dict, Generator, List, Tuple
from shutil import rmtree

import phymmr_tools
import wrap_rocks
import xxhash
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm

from .timekeeper import KeeperMode, TimeKeeper
from .utils import ConcurrentLogger

ROCKSDB_FOLDER_NAME = "rocksdb"
SEQUENCES_DB_NAME = "sequences"
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

def N_trim(
    parent_sequence: str, 
    minimum_sequence_length: int, 
    tike: TimeKeeper
):
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
                tike[0].timer1("now")
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


def translate(in_path: Path, out_path: Path, translate_program="fastatranslate", genetic_code=1):
    """
    Call FASTA translate utility using subprocess.call
    """

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


def glob_for_fasta_and_save_for_runs(globbed: Generator[Path, Any, Any]) -> Dict[str, List[Path]]:
    """
    Scan all the files in the input. Remove _R# and merge taxa
    """

    taxa_runs = {}

    for f in globbed:
        file_is_file = f.is_file()

        file_name = f.name
        file_suffix_allowed = any([suff in file_name for suff in ALLOWED_FILETYPES])

        if (file_is_file and file_suffix_allowed):
            taxa = f.stem.split(".")[0]

            formatted_taxa = truncate_taxa(taxa, extension=".fa")

            taxa_runs.setdefault(formatted_taxa, [])
            taxa_runs[formatted_taxa].append(f)

    return taxa_runs


class SeqDeduplicator:
    def __init__(
        self,
        fa_file_path: Path,
        db: Any,
        minimum_sequence_length: int,
        verbose: int
    ):
        self.fa_file_path = fa_file_path
        self.minimum_sequence_length = minimum_sequence_length
        self.verbose = verbose
        self.db = db       

    def __call__(
        self,
        trim_times: TimeKeeper,
        duplicates: Dict[str, int],
        rev_comp_save: Dict[str, int],
        transcript_mapped_to: Dict[str, str],
        dupes: List[count[int]],
        this_index: List[int],
        fa_file_out: List[str],
        dedup_time: List[int],
    ):

        suffix = self.fa_file_path.suffix
        if any([s in self.fa_file_path.stem for s in (".fq", ".fastq")]):
            read_method = FastqGeneralIterator
        else:
            read_method = SimpleFastaParser

        is_compressed = suffix == '.gz'
        if is_compressed:
            fasta_file = list(read_method(gzip.open(str(self.fa_file_path), "rt")))
        else:
            fasta_file = read_method(open(str(self.fa_file_path), "r"))

        if self.verbose:
            seq_grab_count = get_seq_count(
                self.fa_file_path, is_compressed)

        for_loop = tqdm(
            fasta_file, total=seq_grab_count) if self.verbose else fasta_file

        for row in for_loop:
            header, parent_seq = row[:2]

            if not len(parent_seq) >= self.minimum_sequence_length:
                continue
            parent_seq = parent_seq.upper()

            for seq in N_trim(parent_seq, self.minimum_sequence_length, trim_times):
                length = len(seq)
                header = f"NODE_{this_index[0]}_length_{length}"

                # Check for dupe, if so save how many times that sequence occured
                seq_start = time()

                if seq in transcript_mapped_to:
                    duplicates.setdefault(
                        transcript_mapped_to[seq], 1)
                    duplicates[transcript_mapped_to[seq]] += 1
                    next(dupes[0])
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
                    duplicates.setdefault(
                        transcript_mapped_to[rev_seq], 1)
                    duplicates[transcript_mapped_to[rev_seq]] += 1
                    next(dupes[0])
                    continue
                else:
                    transcript_mapped_to[rev_seq] = header

                seq_end = time()
                dedup_time[0] += seq_end - seq_start
                this_index[0] += 1

                # If no dupe, write to prepared file and db
                line = f">{header}\n{seq}\n"

                fa_file_out.append(line)

                # Get rid of space and > in header (blast/hmmer doesn't like it) Need to push modified external to remove this. ToDo.
                preheader = header.replace(" ", "|")  # pre-hash header
                # Key is a hash of the header
                self.db.put(xxhash.xxh64_hexdigest(
                    preheader), f"{preheader}\n{seq}")


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
        trim_times: TimeKeeper
    ):
        self.fto = formatted_taxa_out
        self.comp = components
        self.trim_times = trim_times
        self.db = None
        self.prepared_file_destination = None
        self.prot_path = None
        self.fa_file_out = []
        self.verbose = verbose
        self.printv = printv
        self.clear_database = clear_database
        self.seq_grab_count = 0
        self.minimum_sequence_length = minimum_sequence_length
        self.taxa_time_keeper = None
        self.prot_components = []
        self.translate_files = []
        self.keep_prepared = keep_prepared
        self.duplicates = {}
        self.rev_comp_save = {}
        self.transcript_mapped_to = {}
        self.dupes = [count()]
        self.this_index = [1]
        self.fa_file_out = []
        self.dedup_time = dedup_time

    def init_db(
            self,
            secondary_directory: Path,
            global_time_keeper: List[TimeKeeper]
    ) -> Any:
        global_time_keeper[0].lap()
        self.taxa_time_keeper = TimeKeeper(KeeperMode.DIRECT)
        self.printv(f"Preparing: {self.fto}", self.verbose, 0)
        self.printv(
            "Formatting input sequences and inserting into database", self.verbose)

        taxa_destination_directory = secondary_directory.joinpath(self.fto)

        rocksdb_path = taxa_destination_directory.joinpath(ROCKSDB_FOLDER_NAME)
        sequences_db_path = rocksdb_path.joinpath(SEQUENCES_DB_NAME)

        if self.clear_database and rocksdb_path.exists():
            self.printv("Clearing old database", self.verbose)
            rmtree(rocksdb_path)

        self.printv("Creating rocksdb database", self.verbose)
        rocksdb_path.mkdir(parents=True, exist_ok=True)
        self.db = wrap_rocks.RocksDB(str(sequences_db_path))
        taxa_destination_directory.mkdir(parents=True, exist_ok=True)
        self.prepared_file_destination = taxa_destination_directory.joinpath(
            self.fto)
        self.prot_path = taxa_destination_directory.joinpath(
            self.fto.replace(".fa", "_prot.fa"))

    def dedup(self):
        for fa_file_path in self.comp:
            SeqDeduplicator(
                fa_file_path,
                self.db,
                self.minimum_sequence_length,
                self.verbose
            )(
                self.trim_times,
                self.duplicates,
                self.rev_comp_save,
                self.transcript_mapped_to,
                self.dupes,
                self.this_index,
                self.fa_file_out,
                self.dedup_time
            )

    def translate_fasta_and_write_to_disk(
        self,
        num_threads: int,
        tmp_path: Path,
    ):
        self.printv("Translating prepared file", self.verbose)
        
        self.taxa_time_keeper.lap()
        
        sequences_per_thread = math.ceil(
            (len(self.fa_file_out) // 2) // num_threads) * 2
        self.translate_files = []

        for i in range(0, len(self.fa_file_out), sequences_per_thread):
            in_path = tmp_path.joinpath(
                self.fto + f'_{len(self.translate_files)}.fa'
            )
            out_path = tmp_path.joinpath(
                self.fto + f'_prot_{len(self.translate_files)}.fa'
            )
            self.translate_files.append((in_path, out_path))
            in_path.write_text(
                "".join(self.fa_file_out[i:i + sequences_per_thread]))

        with Pool(num_threads) as translate_pool:
            translate_pool.starmap(translate, self.translate_files)

        if self.keep_prepared:
            self.prepared_file_destination.write_text(
                "".join(self.fa_file_out))

        self.printv(
            f"Storing translated file in DB. Translate took {self.taxa_time_keeper.lap():.2f}s", self.verbose)

    def store_in_db(
        self,
        prot_max_seqs_per_level: int,
        global_time_keeper: List[TimeKeeper]

    ):        

        prot_components = []
        aa_dupe_count = count()

        out_lines = []
        for _, translate_file in self.translate_files:
            with open(translate_file, "r+", encoding="utf-8") as fp:
                for header, seq in SimpleFastaParser(fp):
                    header = header.replace(" ", "|")
                    if seq in self.transcript_mapped_to:
                        self.duplicates.setdefault(
                            self.transcript_mapped_to[seq], 1)
                        self.duplicates[self.transcript_mapped_to[seq]] += 1
                        next(aa_dupe_count)
                        continue
                    else:
                        self.transcript_mapped_to[seq] = header
                    out_lines.append(f">{header}\n{seq}\n")
            os.remove(translate_file)

        if self.keep_prepared:
            open(self.prot_path, 'w').writelines(out_lines)

        aa_dupes = next(aa_dupe_count)
        self.printv(
            f"AA dedupe took {self.taxa_time_keeper.lap():.2f}s. Kicked {aa_dupes} dupes", self.verbose, 2)

        levels = math.ceil(len(out_lines) / prot_max_seqs_per_level)
        per_level = math.ceil(len(out_lines) / levels)

        component = 0
        for i in range(0, len(out_lines), per_level):
            component += 1

            data = out_lines[i: i + per_level]
            prot_components.append(str(component))
            self.db.put(f"getprot:{component}", "".join(data))

        self.db.put("getall:prot", ",".join(prot_components))

        self.printv(
            f"Translation and storing done! Took {self.taxa_time_keeper.lap():.2f}s", self.verbose)

        sequence_count = self.this_index[0] - 1

        self.printv(
            "Inserted {} sequences over {} batches. Found {} NT and {} AA dupes.".format(
                sequence_count, levels, next(self.dupes[0]), aa_dupes
            ),
            self.verbose,
        )

        # Store the count of dupes in the database
        self.db.put("getall:dupes", json.dumps(self.duplicates))

        self.printv(
            f"{self.fto} took {global_time_keeper[0].lap():.2f}s overall\n", self.verbose)

    def __call__(
        self,
        secondary_directory: Path,
        global_time_keeper: List[TimeKeeper],
        num_threads: int,
        tmp_path: Path,
        prot_max_seqs_per_level: int
    ):
        self.init_db(secondary_directory, global_time_keeper)
        self.dedup()
        self.translate_fasta_and_write_to_disk(num_threads, tmp_path)
        self.store_in_db(prot_max_seqs_per_level, global_time_keeper)


def map_taxa_runs(
    tuple_in: Tuple[str, str],
    verbose: int,
    printv: Callable,
    clear_database: int,
    keep_prepared: int,
    minimum_sequence_length: int,
    secondary_directory: Path,
    global_time_keeper: List[TimeKeeper],
    num_threads: int,
    tmp_path: Path,
    prot_max_seqs_per_level: Path,
    dedup_time: List[int],
    trim_times: TimeKeeper
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
        trim_times
    )(
        secondary_directory,
        global_time_keeper,
        num_threads,
        tmp_path,
        prot_max_seqs_per_level
    )


def main(args):
    msgq = Queue()
    printv = ConcurrentLogger(msgq)
    printv.start()

    input_path = Path(args.INPUT)

    if not input_path.exists():
        printv("ERROR: An existing directory must be provided.", args.verbose, 0)
        return False

    prot_max_seqs_per_level = args.sequences_per_level
    minimum_sequence_length = args.minimum_sequence_length
    num_threads = args.processes

    project_root = Path(__file__).parent.parent
    core_directory = Path(CORE_FOLDER).joinpath(input_path.parts[-1])
    secondary_directory = project_root.joinpath(core_directory)

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

    # Create necessary directories
    printv("Creating directories", args.verbose)
    secondary_directory.mkdir(parents=True, exist_ok=True)
   

    trim_times = [TimeKeeper(KeeperMode.SUM)]  # Append computed time for each loop.
    dedup_time = [0]
    global_time_keeper = [TimeKeeper(KeeperMode.DIRECT)]

    globbed = input_path.glob("*.f[aq]*")
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
            global_time_keeper,
            num_threads,
            tmp_path,
            prot_max_seqs_per_level,
            dedup_time,
            trim_times
        )
        for tuple_in in taxa_runs.items()
    ]

    list(map(lambda x: map_taxa_runs(x[0], *x[1:]), ls_args))

    printv(
        f"Finished! Took {global_time_keeper[0].differential():.2f}s overall.", args.verbose, 0)
    printv(f"N_trim time: {trim_times[0].time1} seconds", args.verbose, 2)
    printv(f"Dedupe time: {dedup_time[0]}", args.verbose, 2)

    msgq.join()
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr Prepare"
    )
