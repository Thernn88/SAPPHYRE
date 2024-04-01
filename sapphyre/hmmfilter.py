from collections import defaultdict
from itertools import combinations
from math import ceil
from multiprocessing import Pool
from os import listdir, mkdir, path
import re

from msgspec import Struct
from sapphyre_tools import (
    constrained_distance,
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    # del_cols,
    OverlapTree,
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
    score_diff_internal: float
    matching_consensus_percent: float
    add_hmmfilter_dupes: bool


class BatchArgs(Struct):
    args: HmmfilterArgs
    genes: list
    nt_input_path: str
    nt_out_path: str
    aa_input_path: str
    aa_out_path: str
    compress: bool
    gene_scores: dict
    prepare_dupe_counts: dict
    reporter_dupe_counts: dict
    has_dupes: bool
    add_hmmfilter_dupes: bool


class NODE(Struct):
    header: str
    bh_id: str
    score: float
    sequence: str | list
    start: int
    end: int
    index: int
    kicks: list
    saves: list
    has_been_saved: bool

def do_consensus(nodes, threshold, prepare_dupe_counts, reporter_dupe_counts):
    if not nodes:
        return "", False
    
    if prepare_dupe_counts or reporter_dupe_counts:
        bundle = [(node.header, node.sequence) for node in nodes]
        sequences = bundle_seqs_and_dupes(bundle, prepare_dupe_counts, reporter_dupe_counts)

        consensus_seq = dumb_consensus_dupe(sequences, threshold, 0)
        converted = convert_consensus([node.sequence for node in nodes], consensus_seq)
    else:
        sequences = [node.sequence for node in nodes]
        consensus_seq = dumb_consensus(sequences, threshold, 0)

        converted = convert_consensus(sequences, consensus_seq)
        
    start = len(converted)
    converted = converted.lstrip("X")
    left = start - len(converted)
    converted = converted.rstrip("X")
    right = start - left - len(converted)

    has_consensus = "X" in converted

    return (("-" * left) + converted.replace("?", "-") + ("-" * right)), has_consensus


class Leaf:

    __slots__ = "children", "length"

    def __init__(self, length):
        self.children = []
        self.length = length


def process_kicks(nodes, debug, gene, filtered_sequences_log):
    kicks = set()
    for hit_a in nodes:
        if hit_a is None:
            continue
        for b_index in hit_a.saves:
            if hit_b := nodes[b_index]:  # walrus operator introduced in version 3.8
                hit_b.has_been_saved = True
        for b_index in hit_a.kicks:
            hit_b = nodes[b_index]
            if hit_b is None:
                continue
            if hit_b.has_been_saved:
                continue
            kicks.add(hit_b.header)
            start, end = get_overlap(hit_a.start, hit_a.end, hit_b.start, hit_b.end, 0)
            kmer_a = hit_a.sequence[start:end]
            kmer_b = hit_b.sequence[start:end]
            if debug:
                filtered_sequences_log.append(
                    f"{gene},{hit_b.header},{hit_b.score},{hit_b.start},{hit_b.end},Internal Overlapped with Highest Score,{gene},{hit_a.header},{hit_a.score},{hit_a.start},{hit_a.end}\nHit A Kmer: {kmer_a}\nHit B Kmer: {kmer_b}\n"
                )
                filtered_sequences_log.append(
                    f"{gene},{hit_b.header},{hit_b.score},{hit_b.start},{hit_b.end},Internal Overlapped with Highest Score,{gene},{hit_a.header},{hit_a.score},{hit_a.start},{hit_a.end}\nHit A Kmer: {kmer_a}\nHit B Kmer: {kmer_b}\n"
                )
            nodes[b_index] = None
    return [node for node in nodes if node], kicks, filtered_sequences_log





# def compare_hit_to_leaf(hit_a, targets, overlap, score_diff, kicks, safe, debug, gene, filtered_sequences_log) -> None:
def compare_hit_to_leaf(hit_a, targets, overlap, score_diff) -> None:

    for hit_b in targets:
        # if hit_b is None:
        #     continue
        # if hit_b.index in kicks:
        #     continue
        if hit_a.bh_id == hit_b.bh_id:
            continue
        # if hit_a.score > hit_b.score:
        # higher, lower = hit_a, hit_b
        # else:
        #     higher, lower = hit_b, hit_a
        # if lower.index in safe:
        #     continue
        if hit_b.score == 0:
            continue
        if hit_a.score / hit_b.score < score_diff:
            continue
        kmer_a = hit_a.sequence[overlap[0]:overlap[1]]
        kmer_b = hit_b.sequence[overlap[0]:overlap[1]]
        if not is_same_kmer(kmer_b, kmer_a):
            hit_a.kicks.append(hit_b.index)
            # if debug:
            #     filtered_sequences_log.append(
            #         f"{gene},{hit_b.header},{hit_b.score},{hit_b.start},{hit_b.end},Internal Overlapped with Highest Score,{gene},{hit_a.header},{hit_a.score},{hit_a.start},{hit_a.end}\nHit A Kmer: {kmer_a}\nHit B Kmer: {kmer_b}\n"
            #     )
        else:
            hit_a.saves.append(hit_b.index)


def internal_filter_gene2(nodes, debug, gene, min_overlap_internal, score_diff_internal):
    intervals = {(node.start, node.end) for node in nodes}
    intervals = {tup: Leaf(tup[1] - tup[0]) for tup in intervals}
    tree = OverlapTree()
    tree.insert_vector(list(intervals.keys()))
    # tree = IntervalTree.from_tuples(intervals)
    for i, node in enumerate(nodes):
        node.index = i
        intervals[(node.start, node.end)].children.append(node)
    # kicks = [False] * len()
    # safe = set()
    filtered_sequence_log = []
    memoize = {}
    for node in nodes:
        if node.score == 0:
            continue
        # if node.index in kicks:
        #     continue
        index_tuple = (node.start, node.end)
        if index_tuple not in memoize:
            overlap = tree.query_overlap(index_tuple)
            node_length = node.end - node.start
            working = []
            # print(f"made for ({node.start}, {node.end})")
            for interval in overlap:
                interval_start, interval_end = interval[0], interval[1]
                i_length = interval[1] - interval[0]
                if i_length < node_length:
                    length = i_length
                else:
                    length = node_length
                length = length + 1
                coords = get_overlap(node.start, node.end, interval_start, interval_end, 0)
                start, end = coords
                if (end - start)/length < min_overlap_internal:
                    continue
                working.append(((interval[0], interval[1]), coords))
                memoize[(node.start, node.end)] = working
        # else:
        #     print(f"hit for ({node.start}, {node.end})")
        score_target = node.score / score_diff_internal
        for interval, coords in memoize[(node.start, node.end)]:
                targets = (hit_b for hit_b in intervals[(interval[0], interval[1])].children if
                           hit_b.score <= score_target)

                           # hit_b.score <= score_target and
                           # hit_b.index not in kicks and
                           # hit_b.index not in safe)
                # compare_hit_to_leaf(node, targets, coords, score_diff_internal, kicks, safe, debug, gene, filtered_sequence_log)
                compare_hit_to_leaf(node, targets, coords, score_diff_internal)


    # return [node for node in nodes if node.index not in kicks], [node.header for node in nodes if node.index in kicks], filtered_sequence_log
    return process_kicks(nodes, debug, gene, filtered_sequence_log)


def internal_filter_gene(nodes, debug, gene, min_overlap_internal, score_diff_internal):
    nodes.sort(key=lambda hit: hit.score, reverse=True)
    filtered_sequences_log = []
    kicks = set()
    safe = set()
    kicked_indices = set()

    for i, hit_a in enumerate(nodes):
        if i in kicked_indices:
            continue
        for j in range(len(nodes) - 1, i, -1):
            if j in kicked_indices:
                continue
            hit_b = nodes[j]
            if j in safe:
                continue
            if hit_a.bh_id == hit_b.bh_id:
                continue
            if hit_b.score == 0:
                continue
            if hit_a.score / hit_b.score < score_diff_internal:
                break

            overlap_coords = get_overlap(
                hit_a.start,
                hit_a.end,
                hit_b.start,
                hit_b.end,
                1,
            )
            amount_of_overlap = 0 if overlap_coords is None else overlap_coords[1] - overlap_coords[0]

            length = min((hit_b.end - hit_b.start), (hit_a.end - hit_a.start)) + 1  # Inclusive
            percentage_of_overlap = amount_of_overlap / length
            if percentage_of_overlap >= min_overlap_internal:

                kmer_a = hit_a.sequence[overlap_coords[0]: overlap_coords[1]]
                kmer_b = hit_b.sequence[overlap_coords[0]: overlap_coords[1]]

                if not is_same_kmer(kmer_a, kmer_b):
                    kicks.add(hit_b.header)
                    kicked_indices.add(j)
                    if debug:
                        filtered_sequences_log.append(
                            f"{gene},{hit_b.header},{hit_b.score},{hit_b.start},{hit_b.end},Internal Overlapped with Highest Score,{gene},{hit_a.header},{hit_a.score},{hit_a.start},{hit_a.end}\nHit A Kmer: {kmer_a}\nHit B Kmer: {kmer_b}\n"
                        )
                else:
                    safe.add(j)


    return [node for i, node in enumerate(nodes) if not i in kicked_indices], filtered_sequences_log, kicks


def kick_read_consensus(
    ref_consensus,
    match_percent,
    nodes,
    kicked_headers,
    consensus_kicks,
    debug,
    gene,
):
    """
    Generates a consensus sequence from the reference sequences and kicks any reads that have a score
    generated by average_match() below the given threshold.

    Args:
    ----
        aa_output (list): The list of reference sequences
        match_percent (float): The minimum percent of matching characters required for a kick
        nodes (list): The list of nodes to be kicked
        kicked_headers (set): The set of headers that have been kicked
        consensus_kicks (list): The list of kicks for debug
        debug (bool): Whether or not debug mode is enabled
        gene (str): The gene being processed
    Returns:
    -------
        float: The average length of the reference sequences
        list: The list of nodes to be kicked
        str: The consensus sequence formed from the reference sequences
    """

    # Kick reads with a score below the threshold
    for i, read in enumerate(nodes):
        average_matching_cols = average_match(
            read.sequence,
            ref_consensus,
            read.start,
            read.end,
        )

        if average_matching_cols < match_percent:
            if debug:
                length = read.end - read.start
                consensus_kicks.append(
                    f"{gene},{read.header},{average_matching_cols},{length}\n"
                )
            kicked_headers.add(read.header)
            nodes[i] = None

    return [i for i in nodes if i is not None]


def bundle_seqs_and_dupes(sequences: list, prepare_dupe_counts, reporter_dupe_counts):
    """
    Pairs each record object with its dupe count from prepare and reporter databases.
    Given dupe count dictionaries and a list of Record objects, makes tuples of the records
    and their dupe counts. Returns the tuples in a list.
    """
    output = []
    for header, seq in sequences:
        node = header.split("|")[3]
        dupes = prepare_dupe_counts.get(node, 1) + sum(
            prepare_dupe_counts.get(node, 1)
            for node in reporter_dupe_counts.get(node, [])
        )
        output.append((seq, dupes))
    return output


def average_match(seq_a, consensus, start, end):
    """
    Returns a score based on the matching percent of characters in the given sequence and the consensus.
    The score is calculated by taking the total matching characters where the sequence does not contain a gap
    and less than 50% of references contain a gap, and dividing it by the total number of characters in the sequence.

    Args:
    ----
        seq_a (str): The sequence to be checked
        consensus (dict): The consensus formed by all the candidates
        start (int): The start index of the sequence
        end (int): The end index of the sequence
    Returns:
    -------
        float: The matching score
    """
    match = 0
    total = 0
    for i in range(start, end):
        if consensus[i].count("-") / len(consensus[i]) > 0.5:
            continue

        total += 1

        if seq_a[i] == "-":
            match -= 1
            continue

        if seq_a[i] in consensus[i]:
            match += 1

    if total == 0:
        return 0

    return match / total

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
    overlap_kicks = []
    internal_kicks = []

    for input_gene in batch_args.genes:
        gene = input_gene.split('.')[0]
        prepare_dupe_count = batch_args.prepare_dupe_counts.get(gene, {})
        reporter_dupe_count = batch_args.reporter_dupe_counts.get(gene, {})
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
        base_to_id = {}
        for header, sequence in aa_sequences:
            aa_output.append((header, sequence))
            if header.endswith("."):
                continue

            aa_count += 1

            start, end = find_index_pair(sequence, "-")
            base = header.split("|")[3]
            score_key = "_".join(base.split("_")[0:2])
            nodes.append(
                NODE(
                    header=header,
                    bh_id=base_to_id.setdefault(base, len(nodes)),
                    score=this_gene_scores.get(score_key, 0),
                    sequence=sequence,
                    start=start,
                    end=end,
                    index=0,
                    kicks=[],
                    saves=[],
                    has_been_saved=False
                )
            )

        ref_average_data_length = []
        ref_consensus = defaultdict(list)
        
        # Create a consensus using dumb_consensus from sapphyre_tools
        reference_seqs = [seq for header, seq in aa_output if header.endswith(".")]
        # Create a flex consensus using the reference sequences
        for seq in reference_seqs:
            start, end = find_index_pair(seq, "-")
            for i in range(start, end):
                ref_consensus[i].append(seq[i])

            ref_average_data_length.append(len(seq) - seq.count("-"))

        nodes = kick_read_consensus(
                ref_consensus,
                args.matching_consensus_percent,
                nodes,
                kicked_headers,
                consensus_kicks,
                args.debug,
                gene,
            )

        # Calculate the average amount of data characters in the reference sequences
        ref_average_data_length = sum(ref_average_data_length) / len(
            ref_average_data_length
        )

        nodes.sort(key=lambda hit: hit.score, reverse=True)
        header_to_hits = defaultdict(list)
        for node in nodes:
            header_to_hits[node.bh_id].append(node)
        
        for node, hits in header_to_hits.items():
            for (i, hit), (j, hit_b) in combinations(enumerate(hits), 2):

                overlap_coords = get_overlap(hit.start, hit.end, hit_b.start, hit_b.end, 1)
                amount_of_overlap = 0 if overlap_coords is None else overlap_coords[1] - overlap_coords[0]
                distance = (hit_b.end - hit_b.start) + 1
                percentage_of_overlap = amount_of_overlap / distance

                if percentage_of_overlap >= 0.8:
                        kicked_headers.add(hit_b.header)
                        if args.debug:
                            overlap_kicks.append(
                                f"{hit_b.header},{hit_b.score},{hit_b.start},{hit_b.end},Same Header Overlap Lowest Score,{hit.header},{hit.score},{start},{end}"
                            )
        
        nodes = [i for i in nodes if i.header not in kicked_headers]
        nodes, internal_header_kicks, internal_log = internal_filter_gene2(nodes.copy(), args.debug, gene, args.min_overlap_internal, args.score_diff_internal)
        # nodes, internal_log, internal_header_kicks = internal_filter_gene(nodes, args.debug, gene, args.min_overlap_internal, args.score_diff_internal)
        # assert nodes == nodes2, "nodes have changed"
        # assert internal_header_kicks == internal_header_kicks2, "kicks have changed"
        kicked_headers.update(internal_header_kicks)
        internal_kicks.extend(internal_log)

        x_cand_consensus, _ = do_consensus(
            nodes, args.consensus, prepare_dupe_count, reporter_dupe_count
        )

        for i, let in enumerate(x_cand_consensus):
            if let == "X":
                reported_regions.append(f"{gene},{i}")

        aa_output = [(header, seq) for header, seq in aa_output if header not in kicked_headers]
        # aa_output = [(header, _del_cols(seq, x_positions[header])) for header, seq in aa_output if header not in kicked_headers]
        # assert aa_output == aa_output2, "aa_output has changed"
        # Align kicks to the NT
        nt_sequences = [
            (header, sequence)
            for header, sequence in parseFasta(nt_in)
            if header not in kicked_headers
        ]

        aa_out_dupes = []
        nt_out_dupes = []
        if batch_args.add_hmmfilter_dupes and batch_args.has_dupes:
            #insert aa dupes
            for header, seq in aa_output:
                node = header.split("|")[3]
                dupes = prepare_dupe_count.get(node, 1) + sum(
                    prepare_dupe_count.get(node, 1)
                    for node in reporter_dupe_count.get(node, [])
                )
                aa_out_dupes.append((header, seq))
                for i in range(dupes - 1):
                    aa_out_dupes.append((f"{header}_dupe_{i}", seq))
            #insert nt dupes
            for header, seq in nt_sequences:
                node = header.split("|")[3]
                dupes = prepare_dupe_count.get(node, 1) + sum(
                    prepare_dupe_count.get(node, 1)
                    for node in reporter_dupe_count.get(node, [])
                )
            
                nt_out_dupes.append((header, seq))
                for i in range(dupes - 1):
                
                    nt_out_dupes.append((f"{header}_dupe{i}", seq))
        else:
            aa_out_dupes = aa_output
            nt_out_dupes = nt_sequences

        aa_has_candidate = False
        for header, sequence in aa_output:
            if not header.endswith("."):
                aa_has_candidate = True
                break

        if aa_has_candidate:
            writeFasta(aa_out, aa_out_dupes, batch_args.compress)
            writeFasta(nt_out, nt_out_dupes, batch_args.compress)


    return (
        consensus_kicks,
        internal_kicks,
        reported_regions,
        overlap_kicks,
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
    has_dupes = False
    prepare_dupe_counts, reporter_dupe_counts = {}, {}
    if path.exists(nt_db_path):
        has_dupes = True
        nt_db = RocksDB(nt_db_path)
        prepare_dupe_counts = json.decode(
            nt_db.get("getall:gene_dupes"), type=dict[str, dict[str, int]]
        )
        reporter_dupe_counts = json.decode(
            nt_db.get("getall:reporter_dupes"), type=dict[str, dict[str, list]]
        )
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
                this_batch_scores,
                prepare_dupe_counts, 
                reporter_dupe_counts,
                has_dupes,
                args.add_hmmfilter_dupes,
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
    overlap_log = []
    for consensus_kicks, internal_kicks, reported_regions, overlap_kicks in results:
        consensus_log.extend(consensus_kicks)
        internal_log.extend(internal_kicks)
        reported_regions.extend(reported_regions)
        overlap_log.extend(overlap_kicks)

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

        with open(path.join(hmmfilter_path, "overlap_kicks.txt"), "w") as fp:
            fp.write("Header,Score,Start,End,Reason,Header,Frame,Score,Start,End\n")
            fp.write("\n".join(overlap_log))

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose, 1)

    return True


def main(args, from_folder):
    if not (0 < args.matching_consensus_percent < 1.0):
        if 0 < args.matching_consensus_percent <= 100:
            args.matching_consensus_percent = args.matching_consensus_percent / 100
        else:
            raise ValueError(
                "Cannot convert matching consensus percent threshold to a percent. Use a decimal or a whole number between 0 and 100"
            )
    this_args = HmmfilterArgs(
        compress=args.compress,
        uncompress_intermediates=args.uncompress_intermediates,
        processes=args.processes,
        verbose=args.verbose,
        debug=args.debug,
        consensus=args.hmmfilter_consensus,
        from_folder=from_folder,
        min_overlap_internal = args.min_overlap_internal,
        score_diff_internal = args.score_diff_internal,
        matching_consensus_percent = args.matching_consensus_percent,
        add_hmmfilter_dupes = args.add_hmmfilter_dupes
    )
    return do_folder(this_args, args.INPUT)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
