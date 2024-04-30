"""Merges all sequences per taxa into single sequence in each gene.

PyLint 9.61/10
"""
from __future__ import annotations

from collections import Counter, defaultdict
from multiprocessing.pool import Pool
from os import mkdir, path, listdir
from shutil import rmtree
from typing import Any

from msgspec import Struct, json
from sapphyre_tools import find_index_pair, get_overlap
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


        else:
            key_set = set(counts.keys()) - {"-"}
            if len(key_set) == 1:
                consensus_sequence += key_set.pop()
            else:
                consensus_sequence += 'X'

    return consensus_sequence


def detect_ambig_with_gaps(nodes, gap_percentage=0.25):
    if not nodes:
        return []

    length = len(nodes[0].sequence)
    indices = []
    prev_indx_target_or_internal_gap = False

    for i in range(length):
        
        counts = {}

        for node in nodes:
            if i >= node.start and i < node.end:
                counts.setdefault(node.sequence[i], 0)
                counts[node.sequence[i]] += node.count

        if not counts:
            continue

        if prev_indx_target_or_internal_gap and "-" in counts and len(counts.keys()) == 1:
            indices.append((i,False))
            prev_indx_target_or_internal_gap = True
            continue

        prev_indx_target_or_internal_gap = False

        total_count = sum(counts.values())

        if "-" in counts and len(counts.keys()) > 1:
            if counts["-"] / total_count >= gap_percentage:
                indices.append((i,True))
                prev_indx_target_or_internal_gap = True
        
        
    
    if not indices:
        return []
    
    flag_start = None
    flag_end = None
    flag_regions = []

    start = None
    end = None

    regions = []

    # print(indices)

    for index, flag in indices:
        if flag:
            if flag_start is None:
                flag_start = index
                flag_end = index
            elif flag_end + 1 == index:
                flag_end = index
            else:
                flag_regions.append((flag_start, flag_end))
                flag_start = index
                flag_end = index

        elif flag_end is not None:
            flag_regions.append((flag_start, flag_end))
            flag_start = None
            flag_end = None
        
        if start is None:
            start = index
            end = index
        elif end + 1 == index:
            end = index
        else:
            regions.append(((start, end), flag_regions))
            start = None
            end = None
            flag_regions = []
    
    if start is not None:
        if flag_start is not None:
            flag_regions.append((flag_start, flag_end))
        regions.append(((start, end), flag_regions))

    return regions

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
        return self.do_gene(*args, **kwds)

    def do_gene(self, aa_gene, nt_gene):
        raw_gene = aa_gene.split('.')[0]
        # print("Doing:",raw_gene)

        

        aa_gene_path = path.join(self.aa_gene_input, aa_gene)
        nt_gene_path = path.join(self.nt_gene_input, nt_gene)
        

        nt_seqs = parseFasta(nt_gene_path)
        aa_seqs = parseFasta(aa_gene_path)

        aa_out = []
        aa_candidates = []
        for header, seq in aa_seqs:
            if header.endswith("."):
                aa_out.append((header, seq))
            else:
                node = header.split("|")[3]
                count = self.prepare_dupes.get(node, 1) + sum(
                    self.prepare_dupes.get(node, 1)
                    for node in self.reporter_dupes.get(node, [])
                )
                aa_candidates.append(Node(header, seq, count, *find_index_pair(seq, "-")))

        aa_candidates.sort(key=lambda x: x.start)

        aa_overlap_groups = disperse_into_overlap_groups(aa_candidates)

        move_record = defaultdict(list)
        for region, group in aa_overlap_groups:
            gap_regions = detect_ambig_with_gaps(group)
            
            for gap_region, kmers in gap_regions:
                for i in range(len(kmers)-1, -1, -1):
                    if i - 1 < 0:
                        continue

                    kmer_right_start, kmer_right_end = kmers[i]
                    kmer_left_start, kmer_left_end = kmers[i-1]

                    kmer_left_distance = kmer_left_end - kmer_left_start
                    kmer_right_distance = kmer_right_end - kmer_right_start

                    if kmer_left_distance != kmer_right_distance:
                        continue

                    # If sequence has overlap with the gap_region generate a flex consensus using left kmer and check if right kmers merge
                    left_consensus = set()
                    check_again_for_right = []
                    for i, node in enumerate(group):
                        if get_overlap(node.start, node.end, gap_region[0], gap_region[1], 0) is not None:
                            kmer = node.sequence[kmer_left_start:kmer_left_end + 1].replace("-", "")
                            if kmer:
                                left_consensus.add(kmer)
                            else:
                                check_again_for_right.append(i)

                    for i in check_again_for_right:
                        node = group[i]
                        
                        right_kmer = node.sequence[kmer_right_start:kmer_right_end + 1].replace("-", "")
                        if right_kmer in left_consensus:

                            this_sequence = list(node.sequence)

                            for x, i in zip(range(kmer_right_start, kmer_right_end+1), range(kmer_left_start, kmer_left_end+1)):
                                this_sequence[i] = this_sequence[x]
                                this_sequence[x] = "-"
                                move_record[node.header].append((i,x))
                        
                            node.sequence = "".join(this_sequence)
                            node.start, node.end = find_index_pair(node.sequence, "-")
                        

        nt_out = []
        candidates = []
        move_log = []

        for header, seq in nt_seqs:
            if header.endswith("."):
                nt_out.append((header, seq))
            else:
                node = header.split("|")[3]
                count = self.prepare_dupes.get(node, 1) + sum(
                    self.prepare_dupes.get(node, 1)
                    for node in self.reporter_dupes.get(node, [])
                )

                moves = move_record[header]
                if moves:
                    seq_as_triplets = [seq[i:i+3] for i in range(0, len(seq), 3)]

                    for i, x in moves:
                        seq_as_triplets[i] = seq_as_triplets[x]
                        seq_as_triplets[x] = "---"

                    if self.debug:
                        move_log.append((header+"_Before", seq))
                        move_log.append((header+"_After", "".join(seq_as_triplets)))

                    seq = "".join(seq_as_triplets)
                

                candidates.append(Node(header, seq, count, *find_index_pair(seq, "-")))
        
        candidates.sort(key=lambda x: x.start)

        overlap_groups = disperse_into_overlap_groups(candidates)

        for _, group in overlap_groups:
            new_node, new_ref, old_taxa = get_header_parts([i.header for i in group])

            new_header = f"{raw_gene}|{new_ref}|{old_taxa}|{new_node}"
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

                key_set = set(column) - {"---"}

                if len(key_set) == 1:
                    out_seq.append(key_set.pop())
                    continue

                for x in range(0,3):
                    triplet_column = {char[x] for char in column}
                    if len(triplet_column) == 1:
                        out_seq.append(triplet_column.pop())
                    else:
                        out_seq.append(get_IUPAC(triplet_column))

            
            new_seq = ("-" * min_start) + "".join(out_seq)
            new_seq = new_seq + ("-" * (msa_end - len(new_seq)))
            

            
            nt_out.append((new_header, new_seq))
            if self.debug > 1:
                for node in group:
                    nt_out.append((node.header.split("|")[3], node.sequence))

            

        writeFasta(path.join(self.nt_gene_output, nt_gene), nt_out, self.compress)

        for _, group in aa_overlap_groups:
            new_node, new_ref, old_taxa = get_header_parts([i.header for i in group])

            new_header = f"{raw_gene}|{new_ref}|{old_taxa}|{new_node}"

            new_seq = []
            

            cand_seq = do_consensus(group, self.threshold)

            aa_out.append((new_header, cand_seq))

            if self.debug > 1:
                for node in group:
                    aa_out.append((node.header.split("|")[3], node.sequence))

        writeFasta(path.join(self.aa_gene_output, aa_gene), aa_out, self.compress)

        return move_log

def do_folder(input_folder, args):
    print("Processing:", input_folder)
    gene_input_folder = path.join(input_folder, "outlier", "blosum") if args.map else None
    if gene_input_folder is None:
        for folder in ["internal", "excise", "clusters", "hmmfilter", "collapsed", "blosum"]:
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
        nt_gene = make_nt(aa_gene)
        
        arguments.append((aa_gene, nt_gene))
        
    if args.processes > 1:
        with Pool(args.processes) as pool:
            logs = pool.starmap(gene_func, arguments)
    else:
        logs = []
        for argument in arguments:
            logs.append(gene_func(*argument))

    if args.debug:
        log_output = []
        for log in logs:
            if log:
                log_output.extend(log)
        writeFasta(path.join(input_folder, "move_log.fa"), log_output)



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
