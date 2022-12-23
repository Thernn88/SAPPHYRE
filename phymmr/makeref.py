from itertools import count
import json
from pathlib import Path
import os
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from multiprocessing.pool import Pool
import sqlite3
import argparse
import wrap_rocks

class Sequence:
    __slots__ = (
        "header",
        "aa_sequence",
        "nt_sequence",
        "taxa",
        "gene",
        "id"
    )
    def __init__(self, header, aa_sequence, nt_sequence, taxa, gene, id):
        self.header = header
        self.aa_sequence = aa_sequence
        self.nt_sequence = nt_sequence
        self.taxa = taxa
        self.gene = gene
        self.id = id

    def to_tuple(self):
        return self.header, self.aa_sequence

    def __str__(self) -> str:
        return f">{self.header}\n{self.aa_sequence}\n"

class Sequence_Set:
    def __init__(self, name):
        self.name = name
        self.sequences = []
        self.aligned_sequences = {}

    def add_sequence(self, seq: Sequence) -> None:
        self.sequences.append(seq)

    def add_aligned_sequences(self, gene:str, aligned_sequences: str) -> None:
        self.aligned_sequences[gene] = aligned_sequences

    def get_aligned_sequences(self) -> str:
        return self.aligned_sequences

    def get_taxa_in_set(self) -> list:
        """
        Returns every taxa contained in the set
        """
        data = {}
        for sequence in self.sequences:
            data[sequence.taxa] = 1
        
        return list(data.keys())

    def get_genes(self) -> list:
        """
        Returns every gene contained in the set
        """
        data = {}
        for sequence in self.sequences:
            data[sequence.gene] = 1
    
        return list(data.keys())

    def get_gene_dict(self, raw = False) -> dict:
        """
        Returns a dictionary with every gene and its corresponding sequences
        """
        data = {}
        for sequence in self.sequences:
            if raw:
                data.setdefault(sequence.gene, []).append(sequence)
            else:
                data.setdefault(sequence.gene, []).append(str(sequence))

        return data

    def get_core_sequences(self) -> dict:
        """
        Returns a dictionary of every gene and its corresponding taxa, headers and sequences
        """
        core_sequences = {}
        for sequence in self.sequences:
            core_sequences.setdefault(sequence.gene, {}).setdefault("aa", []).append((sequence.taxa, sequence.header, sequence.aa_sequence))
            core_sequences.setdefault(sequence.gene, {}).setdefault("nt", []).append((sequence.taxa, sequence.header, sequence.aa_sequence))

        return core_sequences

    def get_blast_data(self):
        blast_data = []
        target_to_taxon = {}
        taxon_to_sequences = {}

        for seq in self.sequences:
            blast_data.append(f">{seq.id}\n{seq.aa_sequence}\n")
            target_to_taxon.setdefault(seq.gene,{})[seq.id] = seq.taxa

            taxon_to_sequences.setdefault(seq.gene, {})[seq.taxa] = seq.aa_sequence

        return "".join(blast_data), target_to_taxon, taxon_to_sequences

    def __str__(self):
        return "".join(map(str, self.sequences))

def get_set_path(set_name):
    set_path = SETS_DIR.joinpath(set_name)
    set_path.mkdir(exist_ok=True)
    return set_path

def generate_aln(set: Sequence_Set, align_method, threads, overwrite):
    sequences = set.get_gene_dict()

    set_path = get_set_path(set.name)

    aln_path = set_path.joinpath("aln")
    aln_path.mkdir(exist_ok=True)

    arguments = []
    for gene, fasta in sequences.items():
        arguments.append((gene, fasta, aln_path, align_method, overwrite,))

    with Pool(threads) as pool:
        aligned_sequences_components = pool.starmap(aln_function, arguments)

    for gene, aligned_sequences in aligned_sequences_components:
        set.add_aligned_sequences(gene, aligned_sequences)
        
def aln_function(gene, content, aln_path, align_method, overwrite):
    fa_file = aln_path.joinpath(gene+".fa")
    aln_file = fa_file.with_suffix('.aln.fa')
    if not aln_file.exists() or overwrite:
        print(f"Generating: {gene}")
        with fa_file.open(mode="w") as fp:
            fp.writelines(content)

        if align_method == "clustal":
            os.system(f"clustalo -i '{fa_file}' -o '{aln_file}'  --full --iter=5 --full-iter --threads=4 --force")#--verbose
        else:
            os.system(f"mafft-linsi '{fa_file}' > '{aln_file}'")
    
    aligned_result = []
    file = aln_file if aln_file.exists() else fa_file # No alignment required
    for seq_record in SeqIO.parse(file, "fasta"):
        header = seq_record.description
        seq = str(seq_record.seq)
        aligned_result.append((header, seq))

    return gene, aligned_result


def generate_stockholm(set: Sequence_Set, threads, overwrite):
    aligned_sequences = set.get_aligned_sequences()

    set_path = get_set_path(set.name)

    aln_path = set_path.joinpath("aln")
    aln_path.mkdir(exist_ok=True)

    arguments = []
    for gene, raw_fasta in aligned_sequences.items():
        arguments.append((aln_path, gene, raw_fasta, overwrite, ))
        
    with Pool(threads) as pool:
        stockh_files = pool.starmap(f2s, arguments)

    return stockh_files

def f2s(aln_path, gene, raw_fasta, overwrite):
    stockh_file = aln_path.joinpath(gene+".stockh")
    if not stockh_file.exists() or overwrite:
        with stockh_file.open(mode="w") as fp:
            fp.write("# STOCKHOLM 1.0\n")
            for seq_tuple in raw_fasta:
                header,aa_sequence = seq_tuple
                fp.write(f"{header: <50}{aa_sequence}\n")
            fp.write("//")
    return (gene,stockh_file)

def generate_hmms(set:Sequence_Set, stockh_files, threads, overwrite):
    set_path = get_set_path(set.name)
    hmm_path = set_path.joinpath("hmms")
    hmm_path.mkdir(exist_ok=True)
    arguments = []
    for gene, stockh_file in stockh_files:
        hmm_file = hmm_path.joinpath(gene).with_suffix('.hmm')

        arguments.append((stockh_file, hmm_file, overwrite, ))

    with Pool(threads) as pool:
        pool.starmap(hmmbuild, arguments)

def hmmbuild(stockhfile: Path, hmm_file: Path, overwrite):
    hmmname = hmm_file.with_suffix("").name

    if not hmm_file.exists() or overwrite:
        os.system(f"hmmbuild -n '{hmmname}' '{hmm_file}' '{stockhfile}'")

def make_blastdb(set: Sequence_Set, overwrite):
    blast_dir = Path(SETS_DIR, set.name, 'blast')
    blast_dir.mkdir(exist_ok=True)

    db_file = blast_dir.joinpath(set.name)

    blast_db_data, target_to_taxon, taxon_to_sequences = set.get_blast_data()

    if db_file.with_suffix('.psq').exists():
        if not overwrite:
            return target_to_taxon, taxon_to_sequences
    
    with NamedTemporaryFile(mode="w") as fp:
        fp.write(blast_db_data)
        fp.flush()

        os.system(f"makeblastdb -in '{fp.name}' -out '{db_file}' -input_type fasta -dbtype prot -title '{set.name}'")#-parse_seqids

    return target_to_taxon, taxon_to_sequences

def main(args):
    global SETS_DIR
    SETS_DIR = Path(args.orthoset_dir)

    if not args.INPUT:
        print("Fatal: Missing input file")
    if not args.set:
        print("Fatal: Missing set name (-s)")

    set_name = args.set#"Ortholog_set_Mecopterida_v4"
    input_file = args.INPUT#"Ortholog_set_Mecopterida_v4.sqlite"
    align_method = args.align_method#"clustal"
    threads = args.processes#2 
    overwrite = args.overwrite#False
    this_set = Sequence_Set(set_name)

    index = count()

    if input_file.split('.')[-1] == "fa":
        for seq_record in SeqIO.parse(input_file, "fasta"):
            header, data = seq_record.description.split(' ', 1)
            data = json.loads(data)
            seq = str(seq_record.seq)
            taxa = data["organism_name"].replace(" ", "_")
            gene = data["pub_og_id"]
            
            this_set.add_sequence(Sequence(header, seq, "", taxa, gene, next(index)))
    else:
        input_path = Path(input_file)
        if not input_path.exists():
            input_path = SETS_DIR.joinpath(input_file)

        orthoset_db_con = sqlite3.connect(input_path)
        cursor = orthoset_db_con.cursor()

        query = f'''SELECT t.name, o.ortholog_gene_id, a.header, a.sequence,  n.sequence, a.id
                FROM orthograph_orthologs         AS o
            INNER JOIN orthograph_sequence_pairs    AS p
            ON o.sequence_pair = p.id
            INNER JOIN orthograph_aaseqs      AS a
            ON a.id = p.aa_seq
            INNER JOIN orthograph_ntseqs      AS n
            ON n.id = p.nt_seq
            INNER JOIN orthograph_taxa AS t
            ON a.taxid = t.id
           '''
            
        rows = cursor.execute(query)

        for row in rows:
            taxa, gene, header, aa_seq, nt_seq, id = row
            this_set.add_sequence(Sequence(header, aa_seq, nt_seq, taxa, gene, id))
    #

    print('Generating aln')
    generate_aln(this_set, align_method, threads, overwrite)

    print("Creating HMM's Stockholm Alignment Files")
    stockh_files = generate_stockholm(this_set, threads, overwrite)

    print("Generating Hmms")
    generate_hmms(this_set, stockh_files, threads, overwrite)

    print("Making blast DB")
    target_to_taxon, taxon_to_sequences = make_blastdb(this_set, overwrite)

    print("Writing to DB")
    rocks_db_path = get_set_path(set_name).joinpath("rocksdb")
    rocksdb_db = wrap_rocks.RocksDB(str(rocks_db_path))

    rocksdb_db.put("getall:taxainset", json.dumps(this_set.get_taxa_in_set()))

    rocksdb_db.put("getall:refseqs", json.dumps(taxon_to_sequences))
    rocksdb_db.put("getall:targetreference", json.dumps(target_to_taxon))

    for gene, data in this_set.get_core_sequences().items():
        rocksdb_db.put(f"getcore:{gene}", json.dumps(data))

    print("Done!")
    return True