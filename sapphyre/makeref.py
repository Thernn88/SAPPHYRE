from collections import Counter, defaultdict
from math import ceil
import os
import sqlite3
from itertools import count
from multiprocessing.pool import Pool
from pathlib import Path
from tempfile import NamedTemporaryFile

import wrap_rocks
from Bio import SeqIO
from msgspec import json
from .utils import printv
from .timekeeper import TimeKeeper, KeeperMode


class Sequence:
    __slots__ = ("header", "aa_sequence", "nt_sequence", "taxon", "gene", "id")

    def __init__(self, header, aa_sequence, nt_sequence, taxon, gene, id) -> None:
        self.header = header
        self.aa_sequence = aa_sequence
        self.nt_sequence = nt_sequence
        self.taxon = taxon
        self.gene = gene
        self.id = id

    def to_tuple(self):
        return self.header, self.aa_sequence

    def __str__(self) -> str:
        return f">{self.header}\n{self.aa_sequence}\n"


class Sequence_Set:
    def __init__(self, name) -> None:
        self.name = name
        self.sequences = []
        self.aligned_sequences = {}

    def add_sequence(self, seq: Sequence) -> None:
        self.sequences.append(seq)

    def add_aligned_sequences(self, gene: str, aligned_sequences: str) -> None:
        self.aligned_sequences[gene] = aligned_sequences

    def get_aligned_sequences(self) -> str:
        return self.aligned_sequences
    
    def get_last_id(self) -> int:
        if not self.sequences:
            return 0
        return self.sequences[-1].id

    def get_gene_taxon_count(self) -> Counter:
        data = defaultdict(set)

        for seq in self.sequences:
            data[seq.taxon].add(seq.gene)

        return Counter({taxon: len(genes) for taxon, genes in data.items()})

    def get_taxon_in_set(self) -> list:
        """Returns every taxon contained in the set."""
        data = set()
        for sequence in self.sequences:
            data.add(sequence.taxon)

        return list(data)

    def get_genes(self) -> list:
        """Returns every gene contained in the set."""
        data = set()
        for sequence in self.sequences:
            data.add(sequence.gene)

        return list(data)

    def get_gene_dict(self, raw=False) -> dict:
        """Returns a dictionary with every gene and its corresponding sequences."""
        data = {}
        for sequence in self.sequences:
            if raw:
                data.setdefault(sequence.gene, []).append(sequence)
            else:
                data.setdefault(sequence.gene, []).append(str(sequence))

        return data

    def get_core_sequences(self) -> dict:
        """Returns a dictionary of every gene and its corresponding taxon, headers and sequences."""
        core_sequences = {}
        for sequence in self.sequences:
            core_sequences.setdefault(sequence.gene, {}).setdefault("aa", []).append(
                (sequence.taxon, sequence.header, sequence.aa_sequence),
            )
            core_sequences.setdefault(sequence.gene, {}).setdefault("nt", [])
            if sequence.nt_sequence:
                core_sequences[sequence.gene]["nt"].append(
                    (sequence.taxon, sequence.header, sequence.nt_sequence),
                )

        return core_sequences

    def get_diamond_data(self):
        diamond_data = []
        target_to_taxon = {}
        taxon_to_sequences = {}

        for seq in self.sequences:
            diamond_data.append(f">{seq.header}\n{seq.aa_sequence}\n")
            target_to_taxon[seq.header] = seq.gene, seq.taxon, len(seq.aa_sequence)

            taxon_to_sequences.setdefault(seq.gene, {})[seq.taxon] = seq.aa_sequence

        return "".join(diamond_data), target_to_taxon, taxon_to_sequences

    def __str__(self) -> str:
        return "".join(map(str, self.sequences))
    
    def absorb(self, other):
        current_cursor = self.get_last_id() + 1
        for i, seq in enumerate(other.sequences):
            seq.id = current_cursor + i
            self.sequences.append(seq)
        for i, seq in enumerate(other.aligned_sequences):
            seq.id = current_cursor + i
            self.aligned_sequences[i] = seq


def generate_aln(set: Sequence_Set, align_method, overwrite, pool, verbosity, set_path):
    sequences = set.get_gene_dict()

    aln_path = set_path.joinpath("aln")
    aln_path.mkdir(exist_ok=True)

    arguments = []
    for gene, fasta in sequences.items():
        arguments.append(
            (
                gene,
                fasta,
                aln_path,
                align_method,
                overwrite,
                verbosity,
            ),
        )

    aligned_sequences_components = pool.starmap(aln_function, arguments)

    for gene, aligned_sequences in aligned_sequences_components:
        set.add_aligned_sequences(gene, aligned_sequences)


def aln_function(gene, content, aln_path, align_method, overwrite, verbosity):
    fa_file = aln_path.joinpath(gene + ".fa")
    aln_file = fa_file.with_suffix(".aln.fa")
    if not aln_file.exists() or overwrite:
        printv(f"Generating: {gene}", verbosity, 2)
        with fa_file.open(mode="w") as fp:
            fp.writelines(content)

        if align_method == "clustal":
            os.system(
                f"clustalo -i '{fa_file}' -o '{aln_file}'  --full --iter=5 --full-iter --threads=4 --force",
            )  # --verbose
        else:
            os.system(f"mafft-linsi '{fa_file}' > '{aln_file}'")

    aligned_result = []
    file = aln_file if aln_file.exists() else fa_file  # No alignment required
    for seq_record in SeqIO.parse(file, "fasta"):
        header = seq_record.description
        seq = str(seq_record.seq)
        aligned_result.append((header, seq))

    return gene, aligned_result

def make_diamonddb(set: Sequence_Set, overwrite, threads):
    diamond_dir = Path(SETS_DIR, set.name, "diamond")
    diamond_dir.mkdir(exist_ok=True)

    db_file = diamond_dir.joinpath(set.name + ".dmnd")

    diamond_db_data, target_to_taxon, taxon_to_sequences = set.get_diamond_data()

    if db_file.exists() and not overwrite:
        return target_to_taxon, taxon_to_sequences

    with NamedTemporaryFile(mode="w") as fp:
        fp.write(diamond_db_data)
        fp.flush()

        os.system(
            f"diamond makedb --in '{fp.name}' --db '{db_file}' --threads {threads}",
        )

    return target_to_taxon, taxon_to_sequences


SETS_DIR = None

def generate_subset(file_paths, taxon_to_kick:set):
    subset = Sequence_Set("subset")
    index = count()
    for file in file_paths:
        for seq_record in SeqIO.parse(file, "fasta"):
            header, data = seq_record.description.split(" ", 1)
            data = json.decode(data)
            seq = str(seq_record.seq)
            taxon = data["organism_name"].replace(" ", "_")
            if taxon not in taxon_to_kick:
                gene = data["pub_og_id"]

                subset.add_sequence(Sequence(header, seq, "", taxon, gene, next(index)))

    return subset

def main(args):
    tk = TimeKeeper(KeeperMode.DIRECT)
    global SETS_DIR
    SETS_DIR = Path(args.orthoset_dir)

    verbosity = args.verbose  # 0

    set_name = args.set  # eg "Ortholog_set_Mecopterida_v4"
    if not set_name:
        printv("Fatal: Missing set name (-s)", verbosity, 0)
        assume = os.path.basename(args.INPUT).split(".")[0].replace(" ", "_")
        if input(f"Would you like to use the input file name as the set name? ({assume}) (y/n): ").lower() == "y":
            set_name = assume
        else:
            return False
        
    kick = args.kick
    if kick:
        if not os.path.exists(kick):
            printv(f"Warning: Kick file {kick} does not exist, Ignoring", verbosity, 0)
            kick = None
        else:
            with open(kick) as fp:
                kick = set(fp.read().split("\n"))
            printv(f"Found {len(kick)} taxon to kick.", verbosity)
        
    input_file = args.INPUT  # "Ortholog_set_Mecopterida_v4.sqlite"
    if not input_file or not os.path.exists(input_file):
        printv("Fatal: Input file not defined or does not exist (-i)", verbosity, 0)
        return False
    align_method = args.align_method  # "clustal"
    threads = args.processes  # 2
    overwrite = args.overwrite  # False
    do_align = args.align or args.all
    do_count = args.count or args.all
    do_diamond = args.diamond or args.all
    this_set = Sequence_Set(set_name)

    index = count()

    set_path = SETS_DIR.joinpath(set_name)
    set_path.mkdir(exist_ok=True)

    if input_file.split(".")[-1] == "fa":
        printv("Input Detected: Single fasta", verbosity)
        for seq_record in SeqIO.parse(input_file, "fasta"):
            header, data = seq_record.description.split(" ", 1)
            data = json.decode(data)
            seq = str(seq_record.seq)
            taxon = data["organism_name"].replace(" ", "_")
            if taxon not in kick:
                gene = data["pub_og_id"]

                this_set.add_sequence(Sequence(header, seq, "", taxon, gene, next(index)))
    elif input_file.split(".")[-1] in {"sql", "sqlite", "sqlite3", "db", "db3", "s3db", "sl3"}:
        printv("Input Detected: Legacy SQL database", verbosity)
        input_path = Path(input_file)
        if not input_path.exists():
            input_path = SETS_DIR.joinpath(input_file)

        orthoset_db_con = sqlite3.connect(input_path)
        cursor = orthoset_db_con.cursor()

        nt_data = {}
        query = """SELECT n.id, n.sequence FROM orthograph_ntseqs AS n"""

        rows = cursor.execute(query)

        nt_data = dict(rows)

        cursor = orthoset_db_con.cursor()

        query = """SELECT p.nt_seq, t.name, o.ortholog_gene_id, a.header, a.sequence,  a.id
                FROM orthograph_orthologs         AS o
            INNER JOIN orthograph_sequence_pairs    AS p
            ON o.sequence_pair = p.id
            INNER JOIN orthograph_aaseqs      AS a
            ON a.id = p.aa_seq
            INNER JOIN orthograph_taxa AS t
            ON a.taxid = t.id
           """

        rows = cursor.execute(query)

        for row in rows:
            nt_seq = None
            nt_id, taxon, gene, header, aa_seq, id = row
            if taxon not in kick:
                if nt_id in nt_data:
                    nt_seq = nt_data[nt_id]
                this_set.add_sequence(Sequence(header, aa_seq, nt_seq, taxon, gene, id))
    else:
        printv("Input Detected: Folder containing Fasta", verbosity)
        file_paths = []
        for file in os.listdir(input_file):
            if file.endswith(".fa"):
                file_paths.append(os.path.join(input_file, file))

        per_thread = ceil(len(file_paths) / threads)
        distributed_files = [(file_paths[i: i + per_thread], kick,) for i in range(0, len(file_paths), per_thread)]

        printv(f"Generating {threads} subsets containing {per_thread} fasta files each", verbosity)
        with Pool(threads) as pool:
            subsets = pool.starmap(generate_subset, distributed_files)

        printv("Merging subsets", verbosity)
        for subset in subsets:
            this_set.absorb(subset)

    if do_count:
        printv("Generating taxon stats", verbosity)
        with open(os.path.join(set_path, "taxon_stats.csv"), "w") as fp:
            fp.write("taxon,count\n")
            counter = this_set.get_gene_taxon_count()

            for taxon, tcount in counter.most_common():
                fp.write(f"{taxon},{tcount}\n")

    if do_align:
        with Pool(threads) as pool:
            printv("Generating aln", verbosity)
            generate_aln(this_set, align_method, overwrite, pool, verbosity, set_path)

    if do_diamond:
        printv("Making Diamond DB", verbosity)
        target_to_taxon, taxon_to_sequences = make_diamonddb(this_set, overwrite, threads)

    printv("Writing to RocksDB", verbosity, 1)
    rocks_db_path = set_path.joinpath("rocksdb")
    rocksdb_db = wrap_rocks.RocksDB(str(rocks_db_path))

    encoder = json.Encoder()
    rocksdb_db.put_bytes("getall:taxoninset", encoder.encode(this_set.get_taxon_in_set()))

    if do_diamond:
        rocksdb_db.put_bytes("getall:refseqs", encoder.encode(taxon_to_sequences))
        rocksdb_db.put_bytes("getall:targetreference", encoder.encode(target_to_taxon))

    for gene, data in this_set.get_core_sequences().items():
        rocksdb_db.put_bytes(f"getcore:{gene}", encoder.encode(data))

    printv(f"Done! Took {tk.differential():.2f}s", verbosity, 1)
    return True
