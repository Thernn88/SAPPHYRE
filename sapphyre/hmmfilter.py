from collections import defaultdict
from itertools import combinations
from math import ceil
from multiprocessing import Pool
from os import listdir, mkdir, path
import re

from msgspec import Struct
from phymmr_tools import (
    constrained_distance,
    dumb_consensus,
    find_index_pair,
    get_overlap,
    is_same_kmer,
)
from wrap_rocks import RocksDB
from msgspec import json

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta


class HmmfilterArgs(Struct):
    compress: bool
    uncompress_intermediates: bool
    processes: int

    verbose: int
    debug: int
    consensus: float

    from_folder: str
    min_overlap_internal: float


class BatchArgs(Struct):
    args: HmmfilterArgs
    genes: list
    nt_input_path: str
    nt_out_path: str
    aa_input_path: str
    aa_out_path: str
    compress: bool
    gene_scores: dict


class NODE(Struct):
    header: str
    base_header: str
    score: float
    sequence: str | list
    start: int
    end: int


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
                counts[node.sequence[i]] += 1

        if not counts:
            consensus_sequence += "-"
            cand_coverage[i] = 0
            continue

        max_count = max(counts.values())
        total_count = sum(counts.values())

        cand_coverage[i] = total_count

        if max_count / total_count >= threshold:
            consensus_sequence += max(counts, key=counts.get)
        else:
            consensus_sequence += 'X'

    return consensus_sequence, cand_coverage


def del_cols(sequence, columns, nt=False):
    if nt:
        seq = [sequence[i: i+3] for i in range(0, len(sequence), 3)]
        for i in columns:
            seq[i] = "---"
        return "".join(seq)
    else:
        seq = list(sequence)
        for i in columns:
            seq[i] = "-"
        return "".join(seq)

def internal_filter_gene(nodes, debug, gene, min_overlap_internal, score_diff_internal=1.5):
    nodes.sort(key=lambda hit: hit.score, reverse=True)
    filtered_sequences_log = []
    kicks = set()

    for i, hit_a in enumerate(nodes):
        if not hit_a:
            continue
        for j in range(len(nodes) - 1, i, -1):
            hit_b = nodes[j]
            if hit_b:
                if hit_a.base_header != hit_b.base_header:
                    if ((hit_a.score / hit_b.score) if hit_b.score != 0 else 0) < score_diff_internal:
                        break

                    overlap_coords = get_overlap(
                        hit_a.start,
                        hit_a.end,
                        hit_b.start,
                        hit_b.end,
                        1,
                    )
                    amount_of_overlap = 0 if overlap_coords is None else overlap_coords[1] - overlap_coords[0]

                    distance = (hit_b.end - hit_b.start) + 1  # Inclusive
                    percentage_of_overlap = amount_of_overlap / distance

                    if percentage_of_overlap >= min_overlap_internal:

                        kmer_a = hit_a.sequence[overlap_coords[0]: overlap_coords[1]]
                        kmer_b = hit_b.sequence[overlap_coords[0]: overlap_coords[1]]

                        if not is_same_kmer(kmer_a, kmer_b):
                            kicks.add(hit_b.header)
                            nodes[j] = None
                            if debug:
                                filtered_sequences_log.append(
                                    f"{gene},{hit_b.header},{hit_b.score},{hit_b.start},{hit_b.end},Internal Overlapped with Highest Score,{gene},{hit_a.header},{hit_a.score},{hit_a.start},{hit_a.end}\nHit A Kmer: {kmer_a}\nHit B Kmer: {kmer_b}\n"
                                )


    return [i for i in nodes if i is not None], filtered_sequences_log, kicks




def process_batch(
    batch_args: BatchArgs,
):
    """
    Processes a batch of genes

    Args:
    ----
        batch_args (BatchArgs): The arguments for the batch
    Returns:
    -------
        bool: True if the batch was successful, False otherwise
        list: The list of kicks for debug
        int: The total number of kicks
        list: The list of genes that were kicked
        list: The list of consensus kicks for debug
        int: The total number of sequences that passed
    """
    args = batch_args.args

    consensus_kicks = []
    reported_regions = []
    internal_kicks = []

    total_x_before = 0
    total_x_after = 0

    for input_gene in batch_args.genes:
        gene = input_gene.split('.')[0]
        kicked_headers = set()
        printv(f"Doing: {gene}", args.verbose, 2)

        nt_out = path.join(batch_args.nt_out_path, input_gene)
        aa_out = path.join(batch_args.aa_out_path, input_gene.replace(".nt.", ".aa."))

        nt_in = path.join(batch_args.nt_input_path, input_gene)
        aa_in = path.join(batch_args.aa_input_path, input_gene.replace(".nt.", ".aa."))

        aa_sequences = parseFasta(aa_in)

        aa_output = []
        
        nodes = []

        this_gene_scores = batch_args.gene_scores.get(gene, {})

        # make nodes out of the non reference sequences for processing
        aa_count = 0
        for header, sequence in aa_sequences:
            aa_output.append((header, sequence))
            if header.endswith("."):
                continue

            aa_count += 1

            start, end = find_index_pair(sequence, "-")
            nodes.append(
                NODE(
                    header=header,
                    base_header = header.split("|")[3],
                    score=this_gene_scores.get(header, 0),
                    sequence=sequence,
                    start=start,
                    end=end,
                )
            )

        ref_average_data_length = []
        ref_consensus = defaultdict(list)
        reference_seqs = [seq for header, seq in aa_output if header.endswith(".")]

        # Create a consensus using dumb_consensus from phymmr_tools
        ref_consensus_seq = dumb_consensus(reference_seqs, 0.5)
        

        # Create a flex consensus using the reference sequences
        for seq in reference_seqs:
            start, end = find_index_pair(seq, "-")
            for i in range(start, end):
                ref_consensus[i].append(seq[i])

            ref_average_data_length.append(len(seq) - seq.count("-"))

        # Calculate the average amount of data characters in the reference sequences
        ref_average_data_length = sum(ref_average_data_length) / len(
            ref_average_data_length
        )

        nodes, internal_log, internal_header_kicks = internal_filter_gene(nodes, args.debug, gene, args.min_overlap_internal)
        kicked_headers.update(internal_header_kicks)
        internal_kicks.extend(internal_log)

        x_cand_consensus, cand_coverage = do_consensus(
            nodes, args.consensus
        )

        for node in nodes:
            node.sequence = list(node.sequence)

        trimmed_pos = 0
        x_positions = defaultdict(set)

        has_x_before = 0
        has_x_after = 0

        for let in x_cand_consensus:
            if let == 'X': has_x_before += 1

        for node in nodes:
            trim_to = None
            for i, let in enumerate(x_cand_consensus):
                if let == 'X':
                    # Within 3 bp of start or end
                    cond_1 = i <= node.start + 2 and i >= node.start
                    cond_2 = i >= node.end - 2 and i <= node.end

                    if cond_1 or cond_2 and ref_consensus_seq[i] != node.sequence[i]:
                        trim_to = i
                        edge = "start" if cond_1 else "end"
                    elif trim_to is not None:
                        if edge == "start":
                            if trim_to == node.start + 2:
                                for i in range(node.start + 2, node.end):
                                    if x_cand_consensus[i] == "X":
                                        trim_to = i
                                    else:
                                        break

                            for x in range(node.start, trim_to):
                                node.sequence[x] = "-"
                                trimmed_pos += 1
                                x_positions[node.header].add(x)
                            node.start = trim_to + 1
                        elif edge == "end":
                            if trim_to == node.end - 2:
                                for i in range(node.end, node.start, -1):
                                    if x_cand_consensus[i] == "X":
                                        trim_to = i
                                    else:
                                        break

                            for x in range(trim_to, node.end):
                                node.sequence[x] = "-"
                                trimmed_pos += 1
                                x_positions[node.header].add(x)
                            node.end = trim_to
                        
                        trim_to = None

        x_cand_consensus, cand_coverage = do_consensus(
            nodes, args.consensus
        )

        for node in nodes:
            i = None
            for poss_i in range(node.start, node.start + 3):
                if node.sequence[poss_i] != x_cand_consensus[poss_i]:
                    i = poss_i

            if not i is None:
                for x in range(node.start , i + 1):
                    node.sequence[x] = "-"
                    trimmed_pos += 1
                    x_positions[node.header].add(x)
                node.start = i+1

            i = None
            for poss_i in range(node.end-3, node.end):
                if node.sequence[poss_i] != x_cand_consensus[poss_i]:
                    i = poss_i

            if not i is None:
                for x in range(i, node.end):
                    node.sequence[x] = "-"
                    trimmed_pos += 1
                    x_positions[node.header].add(x)
                node.end = i

        for node in nodes:
            node.sequence = "".join(node.sequence)

        x_cand_consensus, cand_coverage = do_consensus(
            nodes, args.consensus
        )

        for i, let in enumerate(x_cand_consensus):
            if let == "X":
                has_x_after += 1
                coverage = cand_coverage[i]
                reported_regions.append(f"{gene},{i},{coverage}")

        aa_output = [(header, del_cols(seq, x_positions[header])) for header, seq in aa_output if header not in kicked_headers]

        if aa_output:
            writeFasta(aa_out, aa_output, batch_args.compress)

        # Align kicks to the NT
        nt_sequences = [
            (header, del_cols(sequence, x_positions[header], True))
            for header, sequence in parseFasta(nt_in)
            if header not in kicked_headers
        ]
        if nt_sequences:
            writeFasta(nt_out, nt_sequences, batch_args.compress)

        total_x_before += has_x_before
        total_x_after += has_x_after

    return (
        consensus_kicks,
        internal_kicks,
        reported_regions,
        total_x_before,
        total_x_after,
    )

def do_folder(args: HmmfilterArgs, input_path: str):
    """
    Seperates the genes into batches and processes them in parallel

    Args:
    ----
        args (HmmfilterArgs): The arguments for the Hmmfilter
        input_path (str): The path to the input folder
    Returns:
    -------
        bool: True if the folder was processed successfully, False otherwise
    """
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    hmmfilter_path = path.join(input_path, "outlier", "hmmfilter")
    nt_out_path = path.join(hmmfilter_path, "nt")
    aa_out_path = path.join(hmmfilter_path, "aa")

    mkdir(hmmfilter_path)
    mkdir(nt_out_path)
    mkdir(aa_out_path)

    nt_db_path = path.join(input_path, "rocksdb", "sequences", "nt")
    gene_scores = {}
    if path.exists(nt_db_path):
        nt_db = RocksDB(nt_db_path)
        gene_scores = json.decode(nt_db.get("getall:hmm_gene_scores"), type=dict[str, dict[str, float]])
        del nt_db

    nt_input_path = path.join(input_path, "outlier", args.from_folder, "nt")
    aa_input_path = path.join(input_path, "outlier", args.from_folder, "aa")

    # Process NT
    genes = [
        gene
        for gene in listdir(nt_input_path)
        if gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]

    per_thread = ceil(len(genes) / args.processes)

    compress = not args.uncompress_intermediates or args.compress

    batched_arguments = []
    for i in range(0, len(genes), per_thread):
        this_batch_genes = genes[i : i + per_thread]
        this_batch_scores = {}

        for raw_gene in this_batch_genes:
            gene = raw_gene.split(".")[0]
            this_batch_scores[gene] = gene_scores.get(gene, {})

        batched_arguments.append(BatchArgs(
                args,
                this_batch_genes,
                nt_input_path,
                nt_out_path,
                aa_input_path,
                aa_out_path,
                compress,
                this_batch_scores
            ))


    if args.processes <= 1:
        results = []
        for batch in batched_arguments:
            results.append(process_batch(batch))
    else:
        with Pool(args.processes) as pool:
            results = pool.map(process_batch, batched_arguments)

    consensus_log = []
    internal_log = []
    reported_regions = []
    total_x_before = 0
    total_x_after = 0
    for consensus_kicks, internal_kicks, reported_regions, chunk_x_before, chunk_x_after in results:
        consensus_log.extend(consensus_kicks)
        internal_log.extend(internal_kicks)
        reported_regions.extend(reported_regions)
        total_x_before += chunk_x_before
        total_x_after += chunk_x_after

    print("Ambig columns before trim:", total_x_before, "Ambig columns after trim:", total_x_after)
    

    if args.debug:        
        print(f"Internal Kicks: {len(internal_log)}")

        with open(path.join(hmmfilter_path, "kicked_consensus.txt"), "w") as fp:
            fp.write("Kicked Gene,Header,Score,Start,End\n")
            fp.write("\n".join(consensus_log))

        with open(path.join(hmmfilter_path, "region_of_X.txt"), "w") as fp:
            fp.write("Gene,Position,Coverage\n")
            fp.write("\n".join(reported_regions))
        
        with open(path.join(hmmfilter_path, "internal_kicks.txt"), "w") as fp:
            fp.write("Kicked Gene,Header,Frame,Score,Start,End,Reason,Master Gene,Header,Frame,Score,Start,End\n")
            fp.write("\n".join(internal_log))

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose, 1)

    return True


def main(args, from_folder):

    this_args = HmmfilterArgs(
        compress=args.compress,
        uncompress_intermediates=args.uncompress_intermediates,
        processes=args.processes,
        verbose=args.verbose,
        debug=args.debug,
        consensus=args.consensus,
        from_folder=from_folder,
        min_overlap_internal=args.min_overlap_internal,
    )
    return do_folder(this_args, args.INPUT)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
