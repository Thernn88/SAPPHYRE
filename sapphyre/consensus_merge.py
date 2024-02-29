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
from sapphyre_tools import find_index_pair, get_overlap, score_splits
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
    
    return combinations.get("".join(sorted(list(i for i in combination_set if i != "-"))), "N")


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
    count: int
    start: int
    end: int

def get_header_parts(headers):
    headers.sort()
    node_part = "&&".join(header.split("|")[3] for header in headers)
    ref_part = Counter(header.split("|")[1] for header in headers).most_common(1)[0][0]
    return node_part, ref_part, headers[0].split("|")[2]


def do_consensus(nodes, threshold):
    if not nodes:
        return ""

    length = len(nodes[0].sequence)
    consensus_sequence = ""
    cand_coverage = {}

    for i in range(length):
        counts = {}

        for node in nodes:
            if i >= node.start and i < node.end:
                counts.setdefault(node.sequence[i], 0)
                counts[node.sequence[i]] += node.count

        if not counts:
            consensus_sequence += "-"
            cand_coverage[i] = 0
            continue


        max_count = max(counts.values())
        total_count = sum(counts.values())

        cand_coverage[i] = total_count

        if max_count / total_count > threshold:
            consensus_sequence += max(counts, key=counts.get)
        elif len(set(counts.keys()) - {"-"}) == 1:
            consensus_sequence += "-"
        else:
            consensus_sequence += 'X'

    return consensus_sequence, cand_coverage


class do_gene():
    def __init__(self, aa_gene_input, nt_gene_input, aa_gene_output, nt_gene_output, compress, debug, prepare_dupe_counts, reporter_dupe_counts) -> None:
        self.aa_gene_input = aa_gene_input
        self.nt_gene_input = nt_gene_input

        self.aa_gene_output = aa_gene_output
        self.nt_gene_output = nt_gene_output

        self.compress = compress
        self.debug = debug

        self.threshold = 0.5

        self.prepare_dupes = prepare_dupe_counts
        self.reporter_dupes = reporter_dupe_counts

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        self.do_gene(*args, **kwds)

    def do_gene(self, aa_gene, nt_gene):
        raw_gene = aa_gene.split('.')[0]
        # print("Doing:",raw_gene)

        

        aa_gene_path = path.join(self.aa_gene_input, aa_gene)
        nt_gene_path = path.join(self.nt_gene_input, nt_gene)
        

        nt_seqs = parseFasta(nt_gene_path)

        nt_out = []
        candidates = []

        for header, seq in nt_seqs:
            if header.endswith("."):
                nt_out.append((header, seq))
            else:
                node = header.split("|")[3]
                count = self.prepare_dupes.get(node, 1) + sum(
                    self.prepare_dupes.get(node, 1)
                    for node in self.reporter_dupes.get(node, [])
                )
                candidates.append(Node(header, seq, count, *find_index_pair(seq, "-")))
        
        candidates.sort(key=lambda x: x.start)

        overlap_groups = disperse_into_overlap_groups(candidates)
        
        for region, group in overlap_groups:
            triplets = defaultdict(list)
            
            min_start = min(node.start for node in group)
            msa_end = len(group[0].sequence)

            for node in group:

                for i in range(node.start, node.end, 3):
                    triplets[i].extend([node.sequence[i: i+3]] * node.count)

            out_seq = []
            for i, column in triplets.items():
                column_counts = Counter(column)

                most_common = column_counts.most_common(1)[0]
                if most_common[1] / sum(column_counts.values()) > self.threshold:
                    out_seq.append(most_common[0])
                    continue

                column = set(column)
                if column == set():
                    out_seq.append("-")
                    continue

                for i in range(0,3):
                    out_seq.append(get_IUPAC({char[i] for char in column}))
            
            new_seq = ("-" * min_start) + "".join(out_seq)
            new_seq = new_seq + ("-" * (len(new_seq) - msa_end))
            

            new_node, new_ref, old_taxa = get_header_parts([i.header for i in group])

            new_header = f"{raw_gene}|{new_ref}|{old_taxa}|{new_node}"
            nt_out.append((new_header, new_seq))
            if self.debug:
                for node in group:
                    nt_out.append((node.header.split("|")[3], node.sequence))

            

        writeFasta(path.join(self.nt_gene_output, nt_gene), nt_out, self.compress)

        aa_seqs = parseFasta(aa_gene_path)

        aa_out = []
        candidates = []
        for header, seq in aa_seqs:
            if header.endswith("."):
                aa_out.append((header, seq))
            else:
                node = header.split("|")[3]
                count = self.prepare_dupes.get(node, 1) + sum(
                    self.prepare_dupes.get(node, 1)
                    for node in self.reporter_dupes.get(node, [])
                )
                candidates.append(Node(header, seq, count, *find_index_pair(seq, "-")))

        candidates.sort(key=lambda x: x.start)

        overlap_groups = disperse_into_overlap_groups(candidates)

        for region, group in overlap_groups:
            new_node, new_ref, old_taxa = get_header_parts([i.header for i in group])

            new_header = f"{raw_gene}|{new_ref}|{old_taxa}|{new_node}"

            new_seq = []
            

            cand_seq, x = do_consensus(group, self.threshold)

            aa_out.append((new_header, cand_seq))

            if self.debug:
                for node in group:
                    aa_out.append((node.header.split("|")[3], node.sequence))

        writeFasta(path.join(self.aa_gene_output, aa_gene), aa_out, self.compress)


def do_folder(input_folder, args):
    print("Processing:", input_folder)
    gene_input_folder = None
    for folder in ["excise","internal","hmmfilter", "blosum"]:
        if path.exists(path.join(input_folder, "outlier", folder)):
            gene_input_folder = path.join(input_folder, "outlier", folder)
            break

    if gene_input_folder is None:
        printv(f"No gene input folder found in {input_folder}", args.verbose, 2)
    else:
        printv(f"Current gene: {gene_input_folder}", args.verbose, 2)

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

    nt_db_path = path.join(input_folder, "rocksdb", "sequences", "nt")
    prepare_dupe_counts, reporter_dupe_counts = {}, {}
    if path.exists(nt_db_path):
        nt_db = RocksDB(nt_db_path)
        prepare_dupe_counts = json.decode(
            nt_db.get("getall:gene_dupes"), type=dict[str, dict[str, int]]
        )
        reporter_dupe_counts = json.decode(
            nt_db.get("getall:reporter_dupes"), type=dict[str, dict[str, list]]
        )
        del nt_db

    gene_func = do_gene(aa_gene_input, nt_gene_input, aa_gene_output, nt_gene_output, args.compress, args.debug, prepare_dupe_counts, reporter_dupe_counts)

    arguments = []
    for aa_gene in listdir(aa_gene_input):
        # if "EOG5CZ8WT" in aa_gene:
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
