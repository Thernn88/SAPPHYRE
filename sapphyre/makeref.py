import gzip
import os
import sqlite3
from collections import Counter, defaultdict
from itertools import count
from math import ceil
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import copyfileobj, rmtree
from tempfile import NamedTemporaryFile

import wrap_rocks
import xxhash
from Bio import SeqIO
from msgspec import json
from sapphyre_tools import constrained_distance, find_index_pair, get_overlap

from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta


class Sequence:
    __slots__ = (
        "raw_head",
        "header",
        "aa_sequence",
        "nt_sequence",
        "taxon",
        "gene",
        "id",
    )

    def __init__(
        self, raw_head, header, aa_sequence, nt_sequence, taxon, gene, id
    ) -> None:
        self.raw_head = raw_head
        self.header = header
        self.aa_sequence = aa_sequence.replace("-","")
        self.nt_sequence = nt_sequence if nt_sequence is None else nt_sequence.replace("-","")
        self.taxon = taxon
        self.gene = gene
        self.id = id

    def raw_seq(self):
        """
        Returns the unaligned aa sequence
        """
        return self.aa_sequence.replace("-", "")

    def to_tuple(self):
        """
        Returns the record tuple (header, aa_sequence)
        """
        return self.header, self.aa_sequence

    def seq_with_regen_data(self, nt = False):
        """
        Returns the sequence with the organism name and pub_og_id present for regeneration
        """

        sequence = self.nt_sequence if nt else self.aa_sequence

        if self.raw_head:
            return f">{self.raw_head}\n{sequence}\n"

        return (
            f">{self.header}"
            + ' {"organism_name": "'
            + self.taxon
            + '", "pub_og_id": "'
            + self.gene
            + '"}'
            + f"\n{sequence}\n"
        )

    def __str__(self) -> str:
        return f">{self.header}\n{self.aa_sequence}\n"


class Sequence_Set:
    def __init__(self, name) -> None:
        self.name = name
        self.sequences = []
        self.aligned_sequences = {}
        self.has_aligned = False
        self.has_nt = False

    def add_sequence(self, seq: Sequence) -> None:
        """Adds a sequence to the set."""
        if not self.has_nt and seq.nt_sequence:            
            self.has_nt = True
        self.sequences.append(seq)

    def add_aligned_sequences(
        self, gene: str, aligned_sequences: list[Sequence]
    ) -> None:
        """
        Adds a list of aligned sequences to the set.
        """
        if not self.has_aligned:
            self.has_aligned = True
        self.aligned_sequences[gene] = aligned_sequences

    def get_aligned_sequences(self) -> dict:
        """
        Returns the aligned sequences dict
        """
        return self.aligned_sequences

    def get_last_id(self) -> int:
        """
        Gets the last id of the current sequences in the set
        """
        if not self.sequences:
            return 0
        return self.sequences[-1].id

    def get_gene_taxon_count(self) -> Counter:
        """
        Returns the count of genes present for each taxon
        """
        data = defaultdict(set)

        for seq in self.sequences:
            data[seq.taxon].add(seq.gene)

        return Counter({taxon: len(genes) for taxon, genes in data.items()})

    def kick_dupes(self, headers):
        """
        Kicks sequences with headers in the headers list
        """
        self.sequences = [i for i in self.sequences if i.header not in headers]

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

    def get_gene_dict(self, raw=False, nt=False) -> dict:
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
        from_sequence_list = (
            [item for sublist in self.aligned_sequences.values() for item in sublist]
            if self.has_aligned
            else self.sequences
        )
        for sequence in from_sequence_list:
            core_sequences.setdefault(sequence.gene, {}).setdefault("aa", []).append(
                (
                    sequence.taxon,
                    sequence.header,
                    sequence.aa_sequence
                    if not self.has_aligned
                    else sequence.raw_seq(),
                ),
            )
            core_sequences.setdefault(sequence.gene, {}).setdefault("nt", [])
            if sequence.nt_sequence:
                core_sequences[sequence.gene]["nt"].append(
                    (
                        sequence.taxon,
                        sequence.header,
                        sequence.aa_sequence
                        if not self.has_aligned
                        else sequence.raw_seq(),
                    ),
                )

        return core_sequences

    def get_diamond_data(self):
        """Returns the data required to generate a diamond database."""
        diamond_data = []
        target_to_taxon = {}
        taxon_to_sequences = {}

        from_sequence_list = (
            [item for sublist in self.aligned_sequences.values() for item in sublist]
            if self.has_aligned
            else self.sequences
        )
        for seq in from_sequence_list:
            aaseq = seq.raw_seq() if self.has_aligned else seq.aa_sequence
            key = f"{seq.gene}|{seq.header}"
            diamond_data.append(f">{key}\n{aaseq}\n")
            target_to_taxon[key] = seq.gene, seq.taxon, len(aaseq)

            taxon_to_sequences.setdefault(seq.gene, {})[seq.taxon] = aaseq

        return "".join(diamond_data), target_to_taxon, taxon_to_sequences

    def __str__(self) -> str:
        return "".join(map(str, self.sequences))

    def absorb(self, other):
        """Merges two sets together."""
        current_cursor = self.get_last_id() + 1
        if not self.has_aligned and other.has_aligned:
            self.has_aligned = other.has_aligned
        if not self.has_nt and other.has_nt:
            self.has_nt = other.has_nt
        for i, seq in enumerate(other.sequences):
            seq.id = current_cursor + i
            self.sequences.append(seq)
        for i, seq in enumerate(other.aligned_sequences):
            seq.id = current_cursor + i
            self.aligned_sequences[i] = seq


def cull(sequences, cull_result, percent):
    """
    Culls each edge of the sequences to a column where the percentage of non-gap characters is greater than or equal to the percent argument.
    """
    msa_length = len(sequences[0][1])

    for i in range(msa_length):
        cull_start = i
        if sum(1 for seq in sequences if seq[1][i] != "-") / len(sequences) >= percent:
            break

    for i in range(msa_length - 1, 0, -1):
        cull_end = i
        if sum(1 for seq in sequences if seq[1][i] != "-") / len(sequences) >= percent:
            break
    
    cull_end += 1 # include last bp

    out = []
    for header, seq in sequences:
        left_flank = seq[:cull_start]
        right_flank = seq[cull_end:]
        seq = seq[cull_start: cull_end]

        left_bp = len(left_flank) - left_flank.count("-")
        right_bp = len(right_flank) - right_flank.count("-")

        cull_result[header] = left_bp, right_bp

        out.append((header, seq))

    return out


def generate_hmm(set: Sequence_Set, overwrite, threads, verbosity, set_path):
    """
    Generates the .hmm files for each gene in the set.
    """
    aligned_sequences = set.get_aligned_sequences()
    hmm_path = set_path.joinpath("hmms")
    hmm_path.mkdir(exist_ok=True)

    arguments = []
    for gene, sequences in aligned_sequences.items():
        arguments.append((gene, sequences, hmm_path, overwrite, verbosity))

    with Pool(threads) as pool:
        pool.starmap(hmm_function, arguments)


def hmm_function(gene, sequences, hmm_path, overwrite, verbosity):
    """
    Calls the hmm build function and returns the result.
    """
    hmm_file = hmm_path.joinpath(gene + ".hmm")
    if not hmm_file.exists() or overwrite:
        printv(f"Generating: {gene}", verbosity, 2)
        with NamedTemporaryFile(mode="w") as fp:
            fp.write("".join([i.seq_with_regen_data() for i in sequences]))
            fp.flush()

            os.system(f"hmmbuild '{hmm_file}' '{fp.name}'")


def generate_raw(
    set: Sequence_Set,
    overwrite,
    pool,
    verbosity,
    raw_path,
    raw_nt_path,
):
    """
    Generates the .fa files for each gene in the set.
    """
    sequences = set.get_gene_dict(True)

    arguments = []
    for gene, fasta in sequences.items():
        arguments.append((gene, fasta, raw_path, raw_nt_path, overwrite, verbosity))

    pool.starmap(raw_function, arguments)


def raw_function(gene, sequences, raw_path, raw_nt_path, overwrite, verbosity):
    """
    Writes the sequences to a .fa file.
    """
    raw_fa_file = raw_path.joinpath(gene + ".fa")
    if not raw_fa_file.exists() or overwrite:
        printv(f"Generating: {gene}", verbosity, 2)
        with raw_fa_file.open(mode="w") as fp:
            fp.write("".join([i.seq_with_regen_data() for i in sequences]))
    if not raw_nt_path:
        return
    raw_nt_fa_file = raw_nt_path.joinpath(gene + ".nt.fa")
    if not raw_nt_fa_file.exists() or overwrite:
        with raw_nt_fa_file.open(mode="w") as fp:
            fp.write("".join([i.seq_with_regen_data(nt=True) for i in sequences]))


def generate_aln(
    set: Sequence_Set,
    align_method,
    overwrite,
    pool,
    verbosity,
    set_path,
    raw_path,
    do_cull,
    cull_percent,
):
    """
    Generates the .aln.fa files for each gene in the set.
    """
    sequences = set.get_gene_dict(True).copy()

    aln_path = set_path.joinpath("aln")
    aln_path.mkdir(exist_ok=True)

    trimmed_path = set_path.joinpath("trimmed")
    if os.path.exists(trimmed_path):
        rmtree(trimmed_path)

    trimmed_path.mkdir(exist_ok=True)

    arguments = []
    for gene, fasta in sequences.items():
        arguments.append(
            (
                gene,
                fasta,
                raw_path,
                aln_path,
                trimmed_path,
                align_method,
                overwrite,
                verbosity,
                do_cull,
                cull_percent,
            ),
        )

    aligned_sequences_components = pool.starmap(aln_function, arguments)

    for gene, aligned_sequences, dupe_headers in aligned_sequences_components:
        set.kick_dupes(dupe_headers)
        set.add_aligned_sequences(gene, aligned_sequences)


def do_merge(sequences):
    """
    Merges perfectly overlapping sequences.
    """
    merged_indices = set()
    merge_occured = True
    while merge_occured:
        merge_occured = False
        for i, (header_a, seq_a) in enumerate(sequences):
            for j, (header_b, seq_b) in enumerate(sequences):
                if i == j or i in merged_indices or j in merged_indices:
                    continue
                if header_a.split(":")[0] != header_b.split(":")[0]:
                    continue

                start_a, end_a = find_index_pair(seq_a, "-")
                start_b, end_b = find_index_pair(seq_b, "-")

                overlap_coords = get_overlap(start_a, end_a, start_b, end_b, 1)
                if overlap_coords is None:
                    if end_a < start_b:
                        new_seq = seq_a[:end_a] + seq_b[end_a:]
                    else:
                        new_seq = seq_b[:end_b] + seq_a[end_b:]

                    seq_a = new_seq

                    header_a = header_a + "&&" + header_b
                    merged_indices.add(j)
                    sequences[i] = (header_a, seq_a)
                    merge_occured = True
                else:
                    kmer_a = seq_a[overlap_coords[0] : overlap_coords[1]]
                    kmer_b = seq_b[overlap_coords[0] : overlap_coords[1]]

                    if constrained_distance(kmer_a, kmer_b) == 0:
                        overlap_coord = overlap_coords[0]
                        if start_b >= start_a and end_b <= end_a:
                            new_seq = (
                                seq_a[:overlap_coord]
                                + seq_b[overlap_coord:end_b]
                                + seq_a[end_b:]
                            )
                        elif start_a >= start_b and end_a <= end_b:
                            new_seq = (
                                seq_b[:overlap_coord]
                                + seq_a[overlap_coord:end_a]
                                + seq_b[end_a:]
                            )
                        elif start_b >= start_a:
                            new_seq = seq_a[:overlap_coord] + seq_b[overlap_coord:]
                        else:
                            new_seq = seq_b[:overlap_coord] + seq_a[overlap_coord:]

                        seq_a = new_seq

                        header_a = header_a + "&&" + header_b
                        merged_indices.add(j)
                        sequences[i] = (header_a, seq_a)
                        merge_occured = True

    return sequences


def aln_function(
    gene,
    sequences,
    raw_path,
    aln_path,
    trimmed_path,
    align_method,
    overwrite,
    verbosity,
    do_cull,
    cull_percent,
):
    """
    Calls the alignment program, runs some additional logic on the result and returns the aligned sequences
    """
    raw_fa_file = raw_path.joinpath(gene + ".fa")
    aln_file = aln_path.joinpath(gene + ".aln.fa")
    trimmed_path = trimmed_path.joinpath(gene + ".aln.fa")
    if not aln_file.exists() or overwrite:
        printv(f"Generating: {gene}", verbosity, 2)

        if not raw_fa_file.exists():
            msg = f"Raw file {raw_fa_file} does not exist"
            raise FileNotFoundError(
                msg
            )
        
        if len(list(parseFasta(raw_fa_file))) == 1:
            writeFasta(aln_file, parseFasta(raw_fa_file))
        elif align_method == "clustal":
            os.system(
                f"clustalo -i '{raw_fa_file}' -o '{aln_file}'  --full --iter=5 --threads=1 --full-iter --force",
            )  # --verbose
        else:
            os.system(f"mafft-linsi --thread 1 '{raw_fa_file}' > '{aln_file}'")

    aligned_result = []
    aligned_dict = {}
    file = aln_file if aln_file.exists() else raw_fa_file  # No alignment required
    for header, seq in parseFasta(file, True):
        header = header.split(" ")[0]
        aligned_result.append((header, seq))

    cull_result = {}
    if do_cull:
        aligned_result = cull(aligned_result, cull_result, cull_percent)

    duped_headers = set()
    seq_hashes = set()
    for header, seq in aligned_result:
        component = header.split(":")[0]
        seq_hash = xxhash.xxh3_64(component + seq).hexdigest()
        if seq_hash in seq_hashes:
            duped_headers.add(header)
        else:
            seq_hashes.add(seq_hash)

    if duped_headers:
        aligned_result = [i for i in aligned_result if i[0] not in duped_headers]

    aligned_result = [i for i in aligned_result if len(i[1]) != i[1].count("-")]

    aligned_result = do_merge(aligned_result)

    writeFasta(trimmed_path, aligned_result, False)
    for header, seq in aligned_result:
        aligned_dict[header] = seq

    output = []
    for seq in sequences:
        if seq.header in aligned_dict:
            seq.aa_sequence = aligned_dict[seq.header]
            
            if seq.nt_sequence and cull_result != {}:
                left_bp_remove, right_bp_remove = cull_result[seq.header]
                seq.nt_sequence = seq.nt_sequence[left_bp_remove:][:-right_bp_remove] # Ugly but works

            output.append(seq)

    return gene, output, duped_headers


def make_diamonddb(set: Sequence_Set, overwrite, threads):
    """
    Calls the diamond makedb function and returns the data to insert into the rocksdb.
    """
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


def generate_subset(file_paths, taxon_to_kick: set, nt_input: str = None):
    """
    Grabs alls sequences in a fasta file and inserts them into a subset.
    """
    subset = Sequence_Set("subset")
    index = count()
    for raw_file in file_paths:

        nt_seqs = {}
        if nt_input:
            gene_as_nt = os.path.join(nt_input, os.path.basename(raw_file).replace(".aa.",".nt."))
            if os.path.exists(gene_as_nt):
                nt_seqs = dict(parseFasta(gene_as_nt))

        temp_file = None
        if raw_file.endswith(".gz"):
            temp_file = NamedTemporaryFile(dir = gettempdir())
            with gzip.open(raw_file, 'rb') as f_in, open(temp_file.name, 'wb') as f_out:
                copyfileobj(f_in, f_out)
            file = temp_file.name
        else:
            file = raw_file
            
        this_counter = Counter()

        #Grab headers
        headers = set()
        multiheaders = set()
        with open(file) as fp:
            for line in fp:
                if line[0] != ">" or line[-1] == ".":
                    continue
                if "{" in line:
                    break
                
                header = line.split("|")[2]
                if header in headers:
                    multiheaders.add(header)
                    continue
                headers.add(header)

        for seq_record in SeqIO.parse(file, "fasta"):
            if "|" in seq_record.description and " " not in seq_record.description:
                #Assuming legacy GENE|TAXON|ID
                if seq_record.description.endswith("."):
                    continue

                gene, taxon, header = seq_record.description.split("|")[:3]
                if header in multiheaders:
                    this_counter[header] += 1
                    header = f"{header}_{this_counter[header]}"
                    
                # header = gene+"_"+header
                if taxon.lower() not in taxon_to_kick:
                    subset.add_sequence(
                        Sequence(
                            None, # Want to generate new style
                            header,
                            str(seq_record.seq),
                            nt_seqs.get(seq_record.description),
                            taxon,
                            gene,
                            next(index),
                        )
                    )
            else:
                header, data = seq_record.description.split(" ", 1)
                #Assuming ID {JSON}
                data = json.decode(data)
                seq = str(seq_record.seq)
                taxon = data["organism_name"].replace(" ", "_")
                if (
                    taxon.lower() not in taxon_to_kick
                    and data["organism_name"].lower() not in taxon_to_kick
                ):
                    gene = data["pub_og_id"]

                    subset.add_sequence(
                        Sequence(
                            seq_record.description,
                            header,
                            seq,
                            nt_seqs.get(seq_record.description),
                            taxon,
                            gene,
                            next(index),
                        )
                    )
        del temp_file
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
        if (
            input(
                f"Would you like to use the input file name as the set name? ({assume}) (y/n): "
            ).lower()
            == "y"
        ):
            set_name = assume
        else:
            return False

    kick = args.kick
    if kick:
        if not os.path.exists(kick):
            printv(f"Warning: Kick file {kick} does not exist, Ignoring", verbosity, 0)
            kick = set()
        else:
            with open(kick) as fp:
                kick = set([i.lower() for i in fp.read().split("\n") if i.strip()])
            printv(f"Found {len(kick)} taxon to kick.", verbosity)
    else:
        kick = set()

    nc_genes = args.non_coding_genes
    if nc_genes:
        if not os.path.exists(nc_genes):
            printv(
                f"Warning: NCD Gene file {nc_genes} does not exist, Ignoring",
                verbosity,
                0,
            )
            nc_genes = set()
        else:
            with open(nc_genes) as fp:
                nc_genes = set(fp.read().split("\n"))
            printv(f"Found {len(nc_genes)} non-coding genes.", verbosity)
    else:
        nc_genes = set()

    input_file = args.INPUT  # "Ortholog_set_Mecopterida_v4.sqlite"
    if not input_file or not os.path.exists(input_file):
        printv("Fatal: Input file not defined or does not exist (-i)", verbosity, 0)
        return False
    align_method = args.align_method  # "clustal"
    nt_input = args.nt_input
    threads = args.processes  # 2
    overwrite = args.overwrite  # False
    do_align = args.align or args.all
    do_count = args.count or args.all
    do_diamond = args.diamond or args.all
    do_hmm = args.hmmer or args.all
    cull_percent = args.cull_percent
    do_cull = cull_percent != 0
    this_set = Sequence_Set(set_name)
    set_path = SETS_DIR.joinpath(set_name)
    set_path.mkdir(exist_ok=True)

    #special case for reconcile inputs
    inputs = None
    if os.path.exists(os.path.join(input_file, "aa_merged")):
        nt_input = os.path.join(input_file, "nt_merged")
        input_file = os.path.join(input_file, "aa_merged")
    else:
        #detect super dir
        if any(os.path.exists(os.path.join(input_file, i, "aa_merged")) for i in os.listdir(input_file)):
            inputs = [(os.path.join(
                            input_file, i, "aa_merged"
                        ), os.path.join(
                            input_file, i, "nt_merged"
                        )) 
                        for i in os.listdir(input_file) 
                        if os.path.exists(
                            os.path.join(input_file, i, "aa_merged")
                        )]

    if not inputs:
        inputs = [(input_file, nt_input)]

    for input_file, nt_input in inputs:
        printv(f"Reading files from {input_file}", verbosity)
        if input_file.split(".")[-1] == "fa" and os.path.isfile(input_file):
            printv("Input Detected: Single fasta", verbosity)
            subset = generate_subset([input_file], kick)
            this_set.absorb(subset)

        elif input_file.split(".")[-1] in {
            "sql",
            "sqlite",
            "sqlite3",
            "db",
            "db3",
            "s3db",
            "sl3",
        }:
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
                    this_set.add_sequence(
                        Sequence(None, header, aa_seq, nt_seq, taxon, gene, id)
                    )
        else:
            printv("Input Detected: Folder containing Fasta", verbosity)
            file_paths = []
            for file in os.listdir(input_file):
                if file.endswith(".fa") or file.endswith(".fa.gz"):
                    file_paths.append(os.path.join(input_file, file))

            per_thread = ceil(len(file_paths) / threads)
            distributed_files = [
                (
                    file_paths[i : i + per_thread],
                    kick,
                    nt_input,
                )
                for i in range(0, len(file_paths), per_thread)
            ]

            printv(
                f"Generating {threads} subsets containing {per_thread} fasta files each",
                verbosity,
            )
            with Pool(threads) as pool:
                subsets = pool.starmap(generate_subset, distributed_files)

            printv("Merging subsets", verbosity)
            for subset in subsets:
                this_set.absorb(subset)

    printv("Got input!", verbosity)

    if do_count:
        printv("Generating taxon stats", verbosity)
        with open(os.path.join(set_path, "taxon_stats.csv"), "w") as fp:
            fp.write("taxon,count\n")
            counter = this_set.get_gene_taxon_count()

            for taxon, tcount in counter.most_common():
                fp.write(f"{taxon},{tcount}\n")

    raw_path = set_path.joinpath("raw")
    raw_path.mkdir(exist_ok=True)

    raw_nt_path = None

    if this_set.has_nt:
        raw_nt_path = set_path.joinpath("raw_nt")
        raw_nt_path.mkdir(exist_ok=True)

    with Pool(threads) as pool:
        printv("Generating raw", verbosity)
        generate_raw(
            this_set,
            overwrite,
            pool,
            verbosity,
            raw_path,
            raw_nt_path,
        )

        if do_align or do_cull or do_hmm:
            printv("Generating aln", verbosity)
            generate_aln(
                this_set,
                align_method,
                overwrite,
                pool,
                verbosity,
                set_path,
                raw_path,
                do_cull,
                cull_percent,
            )

    if do_hmm:
        generate_hmm(this_set, overwrite, threads, verbosity, set_path)

    if do_diamond:
        printv("Making Diamond DB", verbosity)
        target_to_taxon, taxon_to_sequences = make_diamonddb(
            this_set, overwrite, threads
        )

    printv("Writing to RocksDB", verbosity, 1)
    rocks_db_path = set_path.joinpath("rocksdb")
    rocksdb_db = wrap_rocks.RocksDB(str(rocks_db_path))

    encoder = json.Encoder()
    rocksdb_db.put_bytes(
        "getall:taxoninset", encoder.encode(this_set.get_taxon_in_set())
    )

    rocksdb_db.put_bytes("getall:nc_genes", ",".join(list(nc_genes)).encode())

    if do_diamond:
        rocksdb_db.put_bytes("getall:refseqs", encoder.encode(taxon_to_sequences))
        rocksdb_db.put_bytes("getall:targetreference", encoder.encode(target_to_taxon))

    for gene, data in this_set.get_core_sequences().items():
        rocksdb_db.put_bytes(f"getcore:{gene}", encoder.encode(data))

    printv(f"Done! Took {tk.differential():.2f}s", verbosity, 1)
    return True
