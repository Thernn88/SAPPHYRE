from __future__ import annotations

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
from shutil import rmtree
from tempfile import TemporaryDirectory, NamedTemporaryFile
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
from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv, gettempdir,ConcurrentLogger

ROCKSDB_FOLDER_NAME = "rocksdb"
SEQUENCES_FOLDER_NAME  = "sequences"
AA_DB_NAME = "aa"
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


def translate(sequences: list, translate_program="fastatranslate", genetic_code=1):
    """
    Call FASTA translate utility using subprocess.call
    """

    with TemporaryDirectory(dir=gettempdir()) as tmp_dir, NamedTemporaryFile(mode = "r" , dir = tmp_dir) as out_path:
        with NamedTemporaryFile(mode = "w+", dir = tmp_dir) as in_path:
            in_path.writelines(sequences)
            in_path.flush()

            call(
                " ".join([
                    translate_program,
                    '--geneticcode',
                    str(genetic_code),
                    str(in_path.name),
                    '>',
                    str(out_path.name)
                ]),
                shell=True
            )

        for header, seq in SimpleFastaParser(out_path):
            yield header, seq

class SeqBuilder:
    def __init__(self, tmt: dict, dup: dict, prot_file: Path):
        self.tmt = tmt
        self.dup = dup
        self.aa_dupe_count = count()
        self.num_sequences = count()
        self.out_lines = []
        self.prot_file = prot_file

        self.prot_file.write_text("")
        self.prot_handle = self.prot_file.open(mode="a")


    def add_sequence(self, header, seq):
        header = header.replace(" ", "|")
        seq_hash = xxhash.xxh64(seq).hexdigest()
        if seq_hash in self.tmt:
            self.dup.setdefault(self.tmt[seq_hash], 1)
            self.dup[self.tmt[seq_hash]] += 1
            next(self.aa_dupe_count)
            return
        else:
            self.tmt[seq_hash] = header
        self.prot_handle.write(f">{header}\n{seq}\n")
        next(self.num_sequences)
    
    def finalise(self):
        self.prot_handle.close()
        return next(self.num_sequences), self.tmt, self.dup, next(self.aa_dupe_count)


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
        self.nt_db = db       

    def __call__(
        self,
        trim_times: TimeKeeper,
        duplicates: Dict[str, int],
        rev_comp_save: Dict[str, int],
        transcript_mapped_to: Dict[str, str],
        dupes: count[int],
        this_index: int,
        fa_file_out: List[str],
        dedup_time: List[int],
    ):

        suffix = self.fa_file_path.suffix

        read_method = SimpleFastaParser
        if any([s in self.fa_file_path.stem for s in (".fq", ".fastq")]):
            read_method = FastqGeneralIterator

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
                header_index = next(this_index)
                header = f"NODE_{header_index}_length_{length}"
                seq_hash = xxhash.xxh64(seq).hexdigest()

                # Check for dupe, if so save how many times that sequence occured
                seq_start = time()

                if seq_hash in transcript_mapped_to:
                    duplicates.setdefault(
                        transcript_mapped_to[seq_hash], 1)
                    duplicates[transcript_mapped_to[seq_hash]] += 1
                    next(dupes)
                    continue
                else:
                    transcript_mapped_to[seq_hash] = header

                # Rev-comp sequence. Save the reverse compliment in a hashmap with the original
                    # sequence so we don't have to rev-comp this unique sequence again
                if seq_hash in rev_comp_save:
                    rev_seq_hash = rev_comp_save[seq_hash]
                else:
                    rev_seq_hash = xxhash.xxh64(phymmr_tools.bio_revcomp(seq)).hexdigest()
                    rev_comp_save[seq_hash] = rev_seq_hash

                # Check for revcomp dupe, if so save how many times that sequence occured
                if rev_seq_hash in transcript_mapped_to:
                    duplicates.setdefault(
                        transcript_mapped_to[rev_seq_hash], 1)
                    duplicates[transcript_mapped_to[rev_seq_hash]] += 1
                    next(dupes)
                    continue
                else:
                    transcript_mapped_to[rev_seq_hash] = header

                seq_end = time()
                dedup_time[0] += seq_end - seq_start

                # If no dupe, write to prepared file and db
                line = f">{header}\n{seq}\n"

                fa_file_out.append(line)

                # Get rid of space and > in header (blast/hmmer doesn't like it) Need to push modified external to remove this.
                preheader = header.replace(" ", "|")  # pre-hash header
                # Key is a hash of the header
                self.nt_db.put(xxhash.xxh64_hexdigest(
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
        self.this_index = count(1)
        self.dedup_time = dedup_time

    def init_db(
            self,
            secondary_directory: Path,
    ) -> Any:
        self.taxa_time_keeper = TimeKeeper(KeeperMode.DIRECT)
        self.printv(f"Preparing: {self.fto}", self.verbose, 0)
        self.printv(
            "Formatting input sequences and inserting into database", self.verbose)

        taxa_destination_directory = secondary_directory.joinpath(self.fto)

        rocksdb_path = taxa_destination_directory.joinpath(ROCKSDB_FOLDER_NAME)
        sequences_folder = rocksdb_path.joinpath(SEQUENCES_FOLDER_NAME)
        aa_db_path = sequences_folder.joinpath(AA_DB_NAME)
        nt_db_path = sequences_folder.joinpath(NT_DB_NAME)

        if self.clear_database and rocksdb_path.exists():
            self.printv("Clearing old database", self.verbose)
            rmtree(rocksdb_path)

        self.printv("Creating rocksdb database", self.verbose)
        rocksdb_path.mkdir(parents=True, exist_ok=True)
        self.aa_db = wrap_rocks.RocksDB(str(aa_db_path))
        self.nt_db = wrap_rocks.RocksDB(str(nt_db_path))
        taxa_destination_directory.mkdir(parents=True, exist_ok=True)
        self.prepared_file_destination = taxa_destination_directory.joinpath(
            self.fto)
        self.prot_path = taxa_destination_directory.joinpath(
            self.fto.replace(".fa", "_prot.fa"))

    def dedup(self):
        for fa_file_path in self.comp:
            SeqDeduplicator(
                fa_file_path,
                self.nt_db,
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
    ):
        self.printv("Translating and Deduping prepared files", self.verbose)
        
        self.taxa_time_keeper.lap()

        self.tsequences = SeqBuilder(self.transcript_mapped_to, self.duplicates, self.prot_path)

        for header, sequence in translate(self.fa_file_out):
            self.tsequences.add_sequence(header, sequence)

        if self.keep_prepared:
            self.prepared_file_destination.write_text(
                "".join(self.fa_file_out))

        del self.fa_file_out

    def store_in_db(
        self,
        prot_max_seqs_per_level: int,
    ):        

        self.printv(
            f"Storing translated file in DB. Translate and Dedupe took {self.taxa_time_keeper.lap():.2f}s", self.verbose)

        prot_components = []
        num_sequences, self.transcript_mapped_to, self.duplicates, aa_dupes = self.tsequences.finalise()
        del self.tsequences

        num_lines = num_sequences * 2

        self.printv(
            f"AA dedupe kicked {aa_dupes} dupes", self.verbose, 2)

        levels = math.ceil(num_lines / prot_max_seqs_per_level)
        per_level = math.ceil(num_lines / levels)

        prot_components = []
        component = 1

        with open(self.prot_path) as fp:
            current_out = []
            current_count = count()
            for line in fp:
                current_out.append(line)
                if next(current_count) == per_level:
                    component += 1
                    prot_components.append(str(component))
                    self.aa_db.put(f"getprot:{component}", "".join(current_out))

                    #Reset counter and output
                    current_out = []
                    current_count = count()
            
            if current_out:
                component += 1
                prot_components.append(str(component))
                self.aa_db.put(f"getprot:{component}", "".join(current_out))
            
        if not self.keep_prepared:
            self.prot_path.unlink()

        self.aa_db.put("getall:prot", ",".join(prot_components))

        self.printv(
            f"Translation and storing done! Took {self.taxa_time_keeper.lap():.2f}s", self.verbose)

        sequence_count = next(self.this_index)

        self.printv(
            "Inserted {} sequences over {} batches. Found {} NT and {} AA dupes.".format(
                sequence_count, levels, next(self.dupes), aa_dupes
            ),
            self.verbose,
        )

        # Store the count of dupes in the database
        self.aa_db.put("getall:dupes", json.dumps(self.duplicates))

        self.printv(
            f"{self.fto} took {self.taxa_time_keeper.differential():.2f}s overall\n", self.verbose)

    def __call__(
        self,
        secondary_directory: Path,
        prot_max_seqs_per_level: int
    ):
        self.init_db(secondary_directory)
        self.dedup()
        self.translate_fasta_and_write_to_disk()
        self.store_in_db(prot_max_seqs_per_level)


def map_taxa_runs(
    tuple_in: Tuple[str, str],
    verbose: int,
    printv: Callable,
    clear_database: int,
    keep_prepared: int,
    minimum_sequence_length: int,
    secondary_directory: Path,
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

    project_root = Path(__file__).parent.parent
    core_directory = Path(CORE_FOLDER).joinpath(input_path.parts[-1])
    secondary_directory = project_root.joinpath(core_directory)

    # Create necessary directories
    printv("Creating directories", args.verbose)
    secondary_directory.mkdir(parents=True, exist_ok=True)
   

    trim_times = [TimeKeeper(KeeperMode.SUM)]  # Append computed time for each loop.
    dedup_time = [0]
    global_time_keeper = TimeKeeper(KeeperMode.DIRECT)

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
            prot_max_seqs_per_level,
            dedup_time,
            trim_times
        )
        for tuple_in in taxa_runs.items()
    ]

    list(map(lambda x: map_taxa_runs(x[0], *x[1:]), ls_args))

    printv(
        f"Finished! Took {global_time_keeper.differential():.2f}s overall.", args.verbose, 0)
    printv(f"N_trim time: {trim_times[0].time1} seconds", args.verbose, 2)
    printv(f"Dedupe time: {dedup_time[0]}", args.verbose, 2)

    msgq.join()
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr Prepare"
    )
