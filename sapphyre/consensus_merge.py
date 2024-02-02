"""Merges all sequences per taxa into single sequence in each gene.

PyLint 9.61/10
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import islice
from multiprocessing.pool import Pool
from os import makedirs, mkdir, path, remove, listdir
from pathlib import Path
from shutil import rmtree
from typing import Any, Literal

from Bio.Seq import Seq
from msgspec import Struct, json
from numpy import uint8
from phymmr_tools import find_index_pair, get_overlap, score_splits
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

def get_IUPAC(combination_set):
    combinations = {
        "A":"A",
        "C":"C",
        "G":"G",
        "T":"T",
        "AG": "R",
        "CT": "Y",
        "CG": "S",
        "AT": "W",
        "GT": "K",
        "AC": "M",
        "CGT": "B",
        "AGT": "D",
        "ACT": "H",
        "ACG": "V",
    }
    
    return combinations.get("".join(sorted(list(combination_set))), "N")


def make_nt(aa_gene):
    return aa_gene.replace(".aa.",".nt.")


def expand_region(original: tuple, expansion: tuple) -> tuple:
    """Expands two (start, end) tuples to cover the entire region."""
    start = min(original[0], expansion[0])
    end = max(original[1], expansion[1])
    return start, end


def disperse_into_overlap_groups(taxa_pair: list) -> list[tuple]:
    """Splits list of (header,sequence) into overlap based groups.

    Returns (overlap region, sequences in overlap region)
    """
    result = []
    current_group = []
    current_region = None

    for sequence in taxa_pair:
        if (
            current_region is None
            or get_overlap(
                sequence.start, sequence.end, current_region[0], current_region[1], 0
            )
            is None
        ):
            if current_group:
                result.append((current_region, current_group))
            current_region = (sequence.start, sequence.end)
            current_group = [sequence]
        else:
            current_group.append(sequence)
            current_region = expand_region(
                current_region, (sequence.start, sequence.end)
            )

    if current_group:
        result.append((current_region, current_group))

    return result


class Node(Struct):
    header: str
    sequence: str
    start: int
    end: int
    

def get_node(header):
    return header.split("|")[3]


class do_gene():
    def __init__(self, aa_gene_input, nt_gene_input, aa_gene_output, nt_gene_output, compress) -> None:
        self.aa_gene_input = aa_gene_input
        self.nt_gene_input = nt_gene_input

        self.aa_gene_output = aa_gene_output
        self.nt_gene_output = nt_gene_output

        self.compress = compress

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        self.do_gene(*args, **kwds)

    def do_gene(self, aa_gene, nt_gene):
        raw_gene = aa_gene.split('.')[0]
        print("Doing:",raw_gene)

        aa_gene_path = path.join(self.aa_gene_input, aa_gene)
        nt_gene_path = path.join(self.nt_gene_input, nt_gene)
        

        nt_seqs = parseFasta(nt_gene_path)

        nt_out = []
        candidates = []

        for header, seq in nt_seqs:
            if header.endswith("."):
                nt_out.append((header, seq))
            else:
                candidates.append(Node(header, seq, *find_index_pair(seq, "-")))
        
        overlap_groups = disperse_into_overlap_groups(candidates)
        header_to_ambig = {}
        for region, group in overlap_groups:
            columns = defaultdict(set)
            for node in group:
                for i, let in enumerate(node.sequence):
                    columns[i].add(let)

            for i, column in columns.items():
                if column == {"-"}:
                    columns[i] = "-"
                    continue
                
                column = {char for char in column if char != '-'}

                if column == {"A"}:
                    columns[i] = "A"
                    continue

                if column == {"C"}:
                    columns[i] = "C"
                    continue

                if column == {"G"}:
                    columns[i] = "G"
                    continue
                
                
                if column == {"T"}:
                    columns[i] = "T"
                    continue

                columns[i] = get_IUPAC(column)
            

            new_seq = "".join(columns[i] for i in range(len(columns)))
            new_header = "&&".join(get_node(i.header) for i in group)
            
            codons = [new_seq[i:i+3] for i in range(0, len(new_seq), 3)]
            acgt_gap = {"A", "C", "G", "T", "-"}
            ambig_cols = set()
            for i, codon in enumerate(codons):
                if {i for i in set(codon) if i not in acgt_gap}:
                    ambig_cols.add(i)

            header_to_ambig[new_header] = ambig_cols

            nt_out.append((new_header, new_seq))

        writeFasta(path.join(self.nt_gene_output, nt_gene), nt_out, self.compress)

        aa_seqs = parseFasta(aa_gene_path)

        aa_out = []
        candidates = []
        for header, seq in aa_seqs:
            if header.endswith("."):
                aa_out.append((header, seq))
            else:
                candidates.append(Node(header, seq, *find_index_pair(seq, "-")))

        overlap_groups = disperse_into_overlap_groups(candidates)

        for region, group in overlap_groups:
            new_header = "&&".join(get_node(i.header) for i in group)
            ambig_cols = header_to_ambig[new_header]
            new_seq = []

            for i in range(len(group[0].sequence)):
                if i in ambig_cols:
                    new_seq.append("X")
                else:
                    amino = list({x.sequence[i] for x in group if x.start <= i <= x.end and x.sequence[i] != "-"})
                    if len(amino) == 0:
                        new_seq.append("-")
                        continue

                    new_seq.append(amino[0])
            
            new_seq = "".join(new_seq)
            aa_out.append((new_header, new_seq))

        writeFasta(path.join(self.aa_gene_output, aa_gene), aa_out, self.compress)


def do_folder(input_folder, args):
    print("Processing:", input_folder)
    gene_input_folder = None
    for folder in ["hmmfilter", "excise", "blosum"]:
        if path.exists(path.join(input_folder, "outlier", folder)):
            gene_input_folder = path.join(input_folder, "outlier", folder)
            break
    
    aa_gene_input = path.join(gene_input_folder, "aa")
    nt_gene_input = path.join(gene_input_folder, "nt")

    aa_gene_output = path.join(input_folder, "aa_merged")
    nt_gene_output = path.join(input_folder, "nt_merged")

    if path.exists(aa_gene_output):
        rmtree(aa_gene_output)
    if path.exists(nt_gene_output):
        rmtree(nt_gene_output)

    mkdir(aa_gene_output)
    mkdir(nt_gene_output)

    gene_func = do_gene(aa_gene_input, nt_gene_input, aa_gene_output, nt_gene_output, args.compress)

    arguments = []
    for aa_gene in listdir(aa_gene_input):
        nt_gene = make_nt(aa_gene)
        
        arguments.append((aa_gene, nt_gene))
    
    if args.processes > 1:
        with Pool(args.processes) as pool:
            pool.starmap(gene_func, arguments)
    else:
        for argument in arguments:
            gene_func(*argument)



def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for folder in args.INPUT:
        do_folder(folder, args)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre MergeOverlap"
    raise Exception(
        msg,
    )
