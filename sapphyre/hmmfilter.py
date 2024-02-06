from collections import defaultdict
from itertools import combinations
from math import ceil
from multiprocessing import Pool
from os import listdir, mkdir, path
import re

from msgspec import Struct
from sapphyre_tools import (
    constrained_distance,
    dumb_consensus,
    dumb_consensus_dupe,
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
    score_diff_internal: float
    matching_consensus_percent: float


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


class NODE(Struct):
    header: str
    base_header: str
    score: float
    sequence: str | list
    start: int
    end: int


def do_consensus(nodes, threshold, prepare_dupe_counts, reporter_dupe_counts):
    if not nodes:
        return "", {}

    length = len(nodes[0].sequence)
    consensus_sequence = ""
    cand_coverage = {}

    node_with_dupes = []
    for node in nodes:
        dupes = prepare_dupe_counts.get(node.base_header, 1) + sum(
            prepare_dupe_counts.get(node.base_header, 1)
            for node in reporter_dupe_counts.get(node.base_header, [])
        )

        node_with_dupes.append((node, dupes))
        

    for i in range(length):
        counts = {}

        for node, count in node_with_dupes:
            if i >= node.start and i < node.end:
                counts.setdefault(node.sequence[i], 0)
                counts[node.sequence[i]] += count

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

def internal_filter_gene(nodes, debug, gene, min_overlap_internal, score_diff_internal):
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
                        continue

                    

                    overlap_coords = get_overlap(
                        hit_a.start,
                        hit_a.end,
                        hit_b.start,
                        hit_b.end,
                        1,
                    )
                    amount_of_overlap = 0 if overlap_coords is None else overlap_coords[1] - overlap_coords[0]

                    distance = min((hit_b.end - hit_b.start), (hit_a.end - hit_a.start)) + 1  # Inclusive
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
        
        # Create a consensus using dumb_consensus from sapphyre_tools
        reference_seqs = [seq for header, seq in aa_output if header.endswith(".")]
        if batch_args.has_dupes:
            bundle = [(header, seq) for header, seq in aa_output if header.endswith(".")]
            sequences = bundle_seqs_and_dupes(bundle, batch_args.prepare_dupe_counts, batch_args.reporter_dupe_counts)

            ref_consensus_seq = dumb_consensus_dupe(sequences, 0.5, 1)
        else:
            ref_consensus_seq = dumb_consensus(reference_seqs, 0.5, 1)
        

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
            header_to_hits[node.base_header.split("_")[1]].append(node)
        
        kicked_headers = set()
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

        nodes, internal_log, internal_header_kicks = internal_filter_gene(nodes, args.debug, gene, args.min_overlap_internal, args.score_diff_internal)
        kicked_headers.update(internal_header_kicks)
        internal_kicks.extend(internal_log)

        x_cand_consensus, cand_coverage = do_consensus(
            nodes, args.consensus, batch_args.prepare_dupe_counts, batch_args.reporter_dupe_counts
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
            nodes, args.consensus, batch_args.prepare_dupe_counts, batch_args.reporter_dupe_counts
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
            for poss_i in range(node.end -1, node.end - 4, -1):
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
            nodes, args.consensus, batch_args.prepare_dupe_counts, batch_args.reporter_dupe_counts
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
    total_x_before = 0
    total_x_after = 0
    for consensus_kicks, internal_kicks, reported_regions, chunk_x_before, chunk_x_after, overlap_kicks in results:
        consensus_log.extend(consensus_kicks)
        internal_log.extend(internal_kicks)
        reported_regions.extend(reported_regions)
        overlap_log.extend(overlap_kicks)
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
        consensus=args.consensus,
        from_folder=from_folder,
        min_overlap_internal = args.min_overlap_internal,
        score_diff_internal = args.score_diff_internal,
        matching_consensus_percent = args.matching_consensus_percent,
    )
    return do_folder(this_args, args.INPUT)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
