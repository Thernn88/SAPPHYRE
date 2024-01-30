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


class CollapserArgs(Struct):
    compress: bool
    uncompress_intermediates: bool
    processes: int

    merge_overlap: int
    matching_percent: float
    overlap_percent: float

    verbose: int
    debug: int
    consensus: float
    matching_consensus_percent: float
    gross_diference_percent: float

    rolling_matching_percent: float
    rolling_consensus_percent: float
    rolling_window_size: int

    true_cluster_threshold: int

    from_folder: str


class BatchArgs(Struct):
    args: CollapserArgs
    genes: list
    nt_input_path: str
    nt_out_path: str
    aa_input_path: str
    aa_out_path: str
    compress: bool
    is_assembly_or_genome: bool
    gene_scores: dict


class NODE(Struct):
    header: str
    base_header: str
    score: float
    sequence: str | list
    start: int
    end: int
    internal_gaps: int
    length: int
    children: list
    is_contig: bool
    kick: bool

    def extend(self, node_2, overlap_coord):
        """
        Merges two nodes together, extending the current node to include the other node.
        Merges only occur if the overlapping kmer is a perfect match.

        Args:
        ----
            node_2 (NODE): The node to be merged into the current node
            overlap_coord (int): The index of the first matching character in the overlapping kmer
        Returns:
        -------
            None
        """
        # If node_2 is contained inside self
        if node_2.start >= self.start and node_2.end <= self.end:
            self.sequence = (
                self.sequence[:overlap_coord]
                + node_2.sequence[overlap_coord : node_2.end]
                + self.sequence[node_2.end :]
            )
        # If node_2 contains self
        elif self.start >= node_2.start and self.end <= node_2.end:
            self.sequence = (
                node_2.sequence[:overlap_coord]
                + self.sequence[overlap_coord : self.end]
                + node_2.sequence[self.end :]
            )

            self.start = node_2.start
            self.end = node_2.end
        # If node_2 is to the right of self
        elif node_2.start >= self.start:
            self.sequence = (
                self.sequence[:overlap_coord] + node_2.sequence[overlap_coord:]
            )

            self.end = node_2.end
        # If node_2 is to the left of self
        else:
            self.sequence = (
                node_2.sequence[:overlap_coord] + self.sequence[overlap_coord:]
            )

            self.start = node_2.start

        self.length = self.end - self.start

        # Save node_2 and the children of node_2 to self
        self.children.append(node_2.header)
        self.children.extend(node_2.children)

        # Self is now a contig
        self.is_contig = True

    def is_kick(self, node_2, overlap_coords, kick_percent, overlap_amount):
        """
        Compares two nodes to determine if the current node should be kicked.
        Kick occurs if the overlapping kmers match below a certain percent.

        Args:
        ----
            node_2 (NODE): The node to be compared to the current node
            overlap_coords (tuple): The start and end index of the overlapping kmer
            kick_percent (float): The minimum percent of matching characters in the overlapping kmer required
            overlap_amount (int): The length of the overlapping kmer
        """

        kmer_current = self.sequence[overlap_coords[0] : overlap_coords[1]]
        kmer_next = node_2.sequence[overlap_coords[0] : overlap_coords[1]]
        non_matching_chars = constrained_distance(kmer_current, kmer_next)
        non_mathching_percent = non_matching_chars / overlap_amount
        matching_percent = 1 - non_mathching_percent

        return matching_percent < kick_percent, matching_percent

    def contig_header(self):
        """
        Generates a header containg the contig node and all children nodes for debug

        Returns:
        -------
            str: The contig header
        """
        contig_node = self.header.split("|")[3]
        children_nodes = "|".join([i.split("|")[3] for i in self.children])
        return f"CONTIG_{contig_node}|{children_nodes}"


def rolling_window_consensus(
    seq,
    candidate_consensus,
    reference_consensus,
    start,
    end,
    window_size,
    min_match,
    step,
):
    """
    Performs an identity check between the given sequence, the reference consensus and
    the consensus formed by all the candidates over a rolling window of the given size.
    The consensus is performed on the candidate consensus where the consensus as a percentage
    of ambigous characters, X, is greater than 0.2. Otherwise the consensus is performed on the
    reference consensus. If the percent of ambigous characters is greater than 0.7, the scan skips
    the current window.

    Args:
    ----
        seq (str): The sequence to be checked
        candidate_consensus (str): The consensus formed by all the candidates
        reference_consensus (str): The reference consensus
        start (int): The start index of the sequence
        end (int): The end index of the sequence
        window_size (int): The size of the rolling window
        min_match (float): The minimum percent of matching characters required for a kick
        step (int): The step size of the rolling window
    Returns:
    -------
        bool: True if the sequence should be kicked, False otherwise
        int: The start index of the window
        float: The percent of matching characters in the window
        float: The percent of ambigous characters in the window
    """
    for i in range(start, end - window_size, step):
        window = seq[i : i + window_size]
        cand_window = candidate_consensus[i : i + window_size]
        ambig_percent = cand_window.count("X") / len(cand_window)

        if ambig_percent > 0.7:
            continue

        if ambig_percent > 0.2:
            ref_window = reference_consensus[i : i + window_size]
            matching = sum(1 for x, y in zip(window, ref_window) if x == y or y == "X")

        else:
            matching = window_size - constrained_distance(window, cand_window)

        if matching / window_size < min_match:
            return (
                True,
                i,
                matching / window_size,
                cand_window.count("X") / len(cand_window),
            )

    return False, None, None, None


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


def do_folder(args: CollapserArgs, input_path: str):
    """
    Seperates the genes into batches and processes them in parallel

    Args:
    ----
        args (CollapserArgs): The arguments for the collapser
        input_path (str): The path to the input folder
    Returns:
    -------
        bool: True if the folder was processed successfully, False otherwise
    """
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    collapsed_path = path.join(input_path, "outlier", "collapsed")
    nt_out_path = path.join(collapsed_path, "nt")
    aa_out_path = path.join(collapsed_path, "aa")

    mkdir(collapsed_path)
    mkdir(nt_out_path)
    mkdir(aa_out_path)

    nt_db_path = path.join(input_path, "rocksdb", "sequences", "nt")
    is_assembly_or_genome = False
    gene_scores = {}
    if path.exists(nt_db_path):
        nt_db = RocksDB(nt_db_path)
        dbis_assembly = nt_db.get("get:isassembly") == "True"
        dbis_genome = nt_db.get("get:isgenome") == "True"
        gene_scores = json.decode(nt_db.get("getall:hmm_gene_scores"), type=dict[str, dict[str, float]])
        is_assembly_or_genome = dbis_genome or dbis_assembly
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
                is_assembly_or_genome,
                this_batch_scores
            ))


    if args.processes <= 1:
        results = []
        for batch in batched_arguments:
            results.append(process_batch(batch))
    else:
        with Pool(args.processes) as pool:
            results = pool.map(process_batch, batched_arguments)

    all_passed = all(i[0] for i in results)

    total_kicks = sum(i[2] for i in results if i[2] != 0)
    total_sequences = sum(i[5] for i in results)

    before_total = sum(i[9] for i in results)
    after_total = sum(i[10] for i in results)
    

    print("Ambig columns before trim:", before_total, "Ambig columns after trim:", after_total)
    

    if args.debug:
        kicked_genes = "\n".join(["\n".join(i[3]) for i in results])
        genes_kicked_count = len(kicked_genes.split("\n"))
        kicked_consensus = "".join(["".join(i[4]) for i in results])
        rescued_kicks = "\n".join(["\n".join(i[6]) for i in results])
        reported_nodes = "True Node,Node,Overlap Percent\n" + "".join(["\n".join(i[7]) for i in results])
        xRegions = "Gene,Index,Candidate Coverage\n" + "\n".join(["\n".join(i[8]) for i in results])
        region_kicks =  "Triggering X Position,Master Header,Master Score,Kicked Header,Kicked Score\n" + "\n".join(["\n".join(i[11]) for i in results])
        int_kicks = "\n".join(["\n".join(i[12]) for i in results])
        internal_kicks = "Kicked Gene,Header,Frame,Score,Start,End,Reason,Master Gene,Header,Frame,Score,Start,End" + int_kicks
        kick_count = int_kicks.count('\n')
        print(f"Internal Kicks: {kick_count}")

        with open(path.join(collapsed_path, "kicked_genes.txt"), "w") as fp:
            fp.write(kicked_genes)

        with open(path.join(collapsed_path, "kicked_consensus.txt"), "w") as fp:
            fp.write(kicked_consensus)

        with open(path.join(collapsed_path, "rescued_kicks.txt"), "w") as fp:
            fp.write(rescued_kicks)

        with open(path.join(collapsed_path, "reported_nodes.txt"), "w") as fp:
            fp.write(reported_nodes)

        with open(path.join(collapsed_path, "region_of_X.txt"), "w") as fp:
            fp.write(xRegions)

        with open(path.join(collapsed_path, "kicks.txt"), "w") as fp:
            fp.write(f"Total Kicks: {total_kicks}\n")
            for data in results:
                kick_list = data[1]
                fp.write("".join(kick_list))

        with open(path.join(collapsed_path, "region_kicks.txt"), "w") as fp:
            fp.write(region_kicks)
        
        with open(path.join(collapsed_path, "internal_kicks.txt"), "w") as fp:
            fp.write(internal_kicks)
    else:
        genes_kicked_count = sum(len(i[3]) for i in results)

    printv(f"Kicked {genes_kicked_count} gene(s)", args.verbose, 1)
    printv(
        f"Kicked {total_kicks} sequences and wrote {total_sequences} sequences",
        args.verbose,
        1,
    )

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose, 1)

    return all_passed


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
    for read in nodes:
        average_matching_cols = average_match(
            read.sequence,
            ref_consensus,
            read.start,
            read.end,
        )

        if average_matching_cols < match_percent:
            if debug:
                consensus_kicks.append(
                    f"{gene},{read.header},{average_matching_cols},{read.length}\n"
                )
            kicked_headers.add(read.header)
            read.kick = True

    return nodes


def merge_overlapping_reads(nodes, minimum_overlap):
    """
    Forms contig sequences by merging reads with a perfect overlapping kmer of a
    length over the minimum_overlap.

    Args:
    ----
        nodes (list): The list of nodes to be merged
        minimum_overlap (int): The minimum length of the overlapping kmer required
    Returns:
    -------
        list: The list of left over reads and formed contigs after merging
    """
    # Rescurive scan
    for i, node in enumerate(nodes):
        # Node is none if it has been merged elsewhere. Node is possibly kicked in a previous function
        if node is None or node.kick:
            continue

        # Recursively merge into node until no merge occurs
        splice_occured = True
        while splice_occured:
            possible_extensions = []
            splice_occured = False
            for j, node_2 in enumerate(nodes):
                if node_2 is None or node_2.kick:
                    continue
                if i == j:
                    continue

                overlap_coords = get_overlap(
                    node.start, node.end, node_2.start, node_2.end, minimum_overlap
                )

                if overlap_coords:
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    overlap_coord = overlap_coords[0]
                    possible_extensions.append((overlap_amount, overlap_coord, j))

            # Merge the smallest overlap first
            for _, overlap_coord, j in sorted(
                possible_extensions, reverse=False, key=lambda x: x[0]
            ):
                node_2 = nodes[j]
                if node_2 is None:
                    continue
                # Confirm still overlaps
                overlap_coords = get_overlap(
                    node.start, node.end, node_2.start, node_2.end, minimum_overlap
                )
                if overlap_coords:
                    # Ensure previous merges still allow for this merge to occur
                    fail = False
                    for x in range(overlap_coords[0], overlap_coords[1]):
                        if node.sequence[x] != node_2.sequence[x]:
                            fail = True
                            break

                    # Merge node into node_2
                    if not fail:
                        splice_occured = True
                        node.extend(node_2, overlap_coords[0])
                        nodes[j] = None

    return list(filter(lambda x: x is not None, nodes))


def get_coverage(nodes, ref_average_data_length):
    """
    Returns the coverage of the given nodes. Coverage equals the amount of candidate
    data columns divided by the average amount of data columns in the reference sequences.

    Args:
    ----
        nodes (list): The list of nodes to be checked
        ref_average_data_length (float): The average amount of data columns in the reference sequences
    Returns:
    -------
        float: The coverage percentage calculated
    """
    read_alignments = [node.sequence for node in nodes]
    data_cols = 0

    for i in range(len(read_alignments[0])):
        if any(seq[i] != "-" for seq in read_alignments):
            data_cols += 1

    return data_cols / ref_average_data_length


def kick_overlapping_nodes(
    nodes,
    min_overlap_percent,
    required_matching_percent,
    gross_difference_percent,
    debug,
):
    """
    Kicks overlapping reads and contigs where the overlapping kmer does not match.

    Args:
    ----
        nodes (list): The list of nodes to be checked
        min_overlap_percent (float): The minimum percent of the overlapping kmer required for a kick
        required_matching_percent (float): The threshold percent of matching characters
        gross_difference_percent (float): The threshold percent of matching characters for sequences with a length difference of 85% or more
        debug (bool): Whether or not debug mode is enabled
    Returns:
    -------
        set: The set of headers that have been kicked
        list: The list of kicks for debug
    """
    kicked_headers = set()
    kicks = []
    for i, node_kick in enumerate(node for node in nodes if not node.kick):
        for j, node_2 in enumerate(node for node in nodes if not node.kick):
            if i == j:
                continue
            # Ensure the node being kicked is smaller than the node it is being compared to
            if node_2.length < node_kick.length:
                continue

            overlap_coords = get_overlap(
                node_2.start, node_2.end, node_kick.start, node_kick.end, 1
            )
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]
                percent = overlap_amount / (node_kick.length - node_kick.internal_gaps)

                if percent >= min_overlap_percent:
                    is_kick, matching_percent = node_kick.is_kick(
                        node_2,
                        overlap_coords,
                        required_matching_percent,
                        overlap_amount,
                    )

                    # If the difference in length is greater than 85%, the matching percent is calculated differently
                    # using the gross_difference_percent variable.
                    length_percent = 0
                    if not is_kick and node_kick.is_contig and node_2.is_contig:
                        length_percent = min(node_kick.length, node_2.length) / max(
                            node_kick.length, node_2.length
                        )

                        if (
                            not is_kick
                            and length_percent <= 0.15
                            and matching_percent < gross_difference_percent
                        ):
                            is_kick = True

                    if is_kick:
                        node_kick.kick = True
                        if debug:
                            kicks.append(
                                f"{node_kick.contig_header()},Kicked By,{node_2.contig_header()},{percent},{matching_percent},{length_percent}\n"
                            )
                        kicked_headers.add(node_kick.header)
                        if node_kick.is_contig:
                            kicked_headers.update(node_kick.children)
                        break

    return kicked_headers, kicks


def kick_rolling_consensus(
    nodes,
    ref_consensus_seq,
    kicked_headers,
    consensus_kicks,
    debug,
    gene,
    window_match,
    consensus_percent,
    window_size,
    step=1,
):
    """
    Creates a candidate consensus and kicks any nodes who fail the rolling_window_consensus check

    Args:
    ----
        nodes (list): The list of nodes to be checked
        ref_consensus_seq (str): The reference consensus sequence
        kicked_headers (set): The set of headers that have been kicked to add to
        consensus_kicks (list): The list of kicks for debug
        debug (bool): Whether or not debug mode is enabled
        gene (str): The gene being processed
    Returns:
    -------
        list: The filtered list of nodes after kicking out the fails
    """
    cand_consensus = dumb_consensus(
        [node.sequence for node in nodes], consensus_percent
    )

    for node in nodes:
        (
            is_kick,
            window_start,
            matching_percent,
            ambig_percent,
        ) = rolling_window_consensus(
            node.sequence,
            cand_consensus,
            ref_consensus_seq,
            node.start,
            node.end,
            window_size,
            window_match,
            step,
        )
        if is_kick:
            if debug:
                consensus_kicks.append(
                    f"{gene},{node.header},{window_start}:{window_start+window_size}\nCand: {node.sequence[window_start: window_start+window_size]}\nCons: {cand_consensus[window_start: window_start+window_size]}\nRefc: {ref_consensus_seq[window_start: window_start+window_size]}\nMatch: {matching_percent},Ambig: {ambig_percent}\n"
                )
            kicked_headers.add(node.header)
            node.kick = True

    nodes = list(filter(lambda x: not x.kick, nodes))
    return nodes


def report_overlaps(nodes, true_cluster_headers):
    THRESHOLD = 0
    true, not_true = [], []
    reported = []

    for node in nodes:
        if node.header in true_cluster_headers:
            true.append(node)
        else:
            not_true.append(node)

    for node_2 in not_true:
        for node in true:
            overlap_coords = get_overlap(
                node_2.start, node_2.end, node.start, node.end, 1
            )
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]
                percent = overlap_amount / min(node.length, node_2.length)

                if percent > THRESHOLD:
                    reported.append(f"{node.header},{node_2.header},{percent}")
                    break

    return reported


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

def find_regions_with_x(sequence, min_x_count=3, window_min_size=30):
    pattern = re.compile(r'[^-]{3,}')

    matches = pattern.finditer(sequence)
    matching_regions = [(match.start(), match.end() - 1) for match in matches if match.group().count('X') >= min_x_count and len(match.group()) >= window_min_size]

    return matching_regions

def get_nodes_in_region(nodes, region, minimum_overlap=15):
    overlapping_nodes = []

    for node in nodes:
        if region >= node.start and region <= node.end:
        # overlap = get_overlap(region[0], region[-1], node.start, node.end, minimum_overlap)
        # if overlap is not None:
            overlapping_nodes.append(node)

    overlapping_nodes.sort(key = lambda x: x.score, reverse= True)

    return [i.header for i in overlapping_nodes]


def internal_filter_gene(nodes, debug, gene, min_overlap_internal=0.9, score_diff_internal=1.5):
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
    kicked_genes = []
    consensus_kicks = []
    total = 0
    passed_total = 0

    has_x_genes = 0
    has_x_genes_after = 0 

    kicks = []
    reported = []
    rescues = []
    reported_regions = []
    internal_kicks = []
    region_kicks = []
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

        true_cluster_raw = []
        this_gene_scores = batch_args.gene_scores.get(gene, {})

        # make nodes out of the non reference sequences for processing
        aa_count = 0
        for header, sequence in aa_sequences:
            aa_output.append((header, sequence))
            if header.endswith("."):
                continue

            aa_count += 1

            start, end = find_index_pair(sequence, "-")

            internal_gaps = sequence[start:end].count("-")

            node_is_contig = batch_args.is_assembly_or_genome or "&&" in header

            id = header.split("|")[3].split("_")[1]
            true_cluster_raw.append((int(id), header))

            nodes.append(
                NODE(
                    header=header,
                    base_header = header.split("|")[3],
                    score=this_gene_scores.get(header, 0),
                    sequence=sequence,#list(sequence),
                    start=start,
                    end=end,
                    length=(end - start),
                    internal_gaps=internal_gaps,
                    children=[],
                    is_contig=node_is_contig,
                    kick=False,
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

        nodes, internal_log, kicks = internal_filter_gene(nodes, args.debug, gene)
        kicked_headers.update(kicks)
        internal_kicks.extend(internal_log)

        x_cand_consensus, cand_coverage = do_consensus(
            nodes, args.consensus
        )
        trimmed_pos = 0
        x_positions = defaultdict(set)

        has_x_before = 0
        has_x_after = 0

        for let in x_cand_consensus:
            if let == 'X': has_x_before += 1

        if False: #Toggle: trim edge regions
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

        # for node in nodes:
        #     i = None
        #     for poss_i in range(node.start, node.start + 3):
        #         if node.sequence[poss_i] != x_cand_consensus[poss_i]:
        #             i = poss_i

        #     if not i is None:
        #         for x in range(node.start , i + 1):
        #             node.sequence[x] = "-"
        #             trimmed_pos += 1
        #             x_positions[node.header].add(x)
        #         node.start = i+1

        #     i = None
        #     for poss_i in range(node.end-3, node.end):
        #         if node.sequence[poss_i] != x_cand_consensus[poss_i]:
        #             i = poss_i

        #     if not i is None:
        #         for x in range(i, node.end):
        #             node.sequence[x] = "-"
        #             trimmed_pos += 1
        #             x_positions[node.header].add(x)
        #         node.end = i

        # for node in nodes:
        #     node.sequence = "".join(node.sequence)

        x_cand_consensus, cand_coverage = do_consensus(
            nodes, args.consensus
        )

        for i, let in enumerate(x_cand_consensus):
            if let == "X":
                has_x_after += 1
                coverage = cand_coverage[i]
                reported_regions.append(f"{gene},{i},{coverage}")

        

        regions = find_regions_with_x(x_cand_consensus)
        x_nkicks = 0
        if False: #Toggle: Kick x regions of 60 bp
            threshold = 0.5
            for region in regions:
                this_x_positions = []
                for i, let in enumerate(x_cand_consensus[region[0]: region[1]], region[0]):
                    if let == "X":
                        this_x_positions.append(i)

                for x_pos in this_x_positions:
                    headers = get_nodes_in_region(nodes, x_pos)

                    # Grab master
                    master = None
                    for node in nodes:
                        if node.kick:
                            continue
                        if node.header in headers:
                            fail = False
                            for i in range(node.start, node.end):
                                if node.sequence[i] != x_cand_consensus[i] and x_cand_consensus[i] != "X":
                                    fail = True
                                    break
                            
                            if node.score < 30:
                                continue

                            if not fail:
                                master = node
                    if not master:
                        continue

                    for node in nodes:
                        if node.kick:
                            continue
                        if node == master:
                            continue
                        if not node.header in headers:
                            continue
                        if node.score >= master.score:
                            continue
                        overlap_coords = get_overlap(node.start, node.end, master.start, master.end, 1)

                        if overlap_coords is None:
                            continue

                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        percent = overlap_amount / min(node.length, master.length)

                        if percent < threshold:
                            continue

                        if node.sequence[i] != master.sequence[i]:
                            x_nkicks += 1
                            node.kick = True
                            kicked_headers.add(node.header)
                            region_kicks.append(f"{gene},{i},{master.header},{master.score},{node.header},{node.score}")
                            break

        nodes = [node for node in nodes if not node.kick]

        true_cluster_raw.sort(key = lambda x: x[0])
        before_true_clusters = []

        current_cluster = []
        for child_index, child_full in true_cluster_raw:
            if not current_cluster:
                current_cluster.append(child_full)
                current_index = child_index
            else:
                if child_index - current_index <= args.true_cluster_threshold:
                    current_cluster.append(child_full)
                    current_index = child_index
                else:
                    before_true_clusters.append(current_cluster)
                    current_cluster = [child_full]
                    current_index = child_index
        
        if current_cluster:
            before_true_clusters.append(current_cluster)

        if True: # Toggle: read consensus
            nodes = kick_read_consensus(
                ref_consensus,
                args.matching_consensus_percent,
                nodes,
                kicked_headers,
                consensus_kicks,
                args.debug,
                gene,
            )

        if not nodes:
            total += aa_count
            kicked_genes.append(
                f"No valid sequences after consensus: {gene}"
            )
            continue

        if False: # Toggle: rolling consensus
            if not batch_args.is_assembly_or_genome:
                nodes = kick_rolling_consensus(
                    nodes,
                    ref_consensus_seq,
                    kicked_headers,
                    consensus_kicks,
                    args.debug,
                    gene,
                    args.rolling_matching_percent,
                    args.rolling_consensus_percent,
                    args.rolling_window_size,
                )

        if not nodes:
            total += aa_count
            kicked_genes.append(
                f"No valid sequences after rolling candidate consensus: {gene}"
            )
            continue

        

        printv("Merging Overlapping Reads", args.verbose, 3)
        nodes = merge_overlapping_reads(nodes, args.merge_overlap)

        printv("Calculating Coverage", args.verbose, 3)
        coverage = get_coverage(nodes, ref_average_data_length)

        req_coverage = 0.3 if batch_args.is_assembly_or_genome else 0.01
        if coverage < req_coverage:
            total += aa_count
            kicked_genes.append(f"{gene} -> failed due to Coverage: {coverage}")
            continue

        if args.debug:
            kicks.append(
                f"Kicks for {gene}\nHeader B,,Header A,Overlap Percent,Matching Percent,Length Ratio\n"
            )
        nodes.sort(key=lambda x: x.length, reverse=True)

        if False: # Toggle: Kick overlapping assembly contigs
            printv("Kicking Overlapping Reads", args.verbose, 3)
            this_kicks, this_debug = kick_overlapping_nodes(
                nodes,
                args.overlap_percent,
                args.matching_percent,
                args.gross_diference_percent,
                args.debug,
            )

            if args.debug:
                kicks.extend(this_debug)

            kicked_headers.update(this_kicks)
        after_data = []
        for node in nodes:
            if node.header in kicked_headers:
                continue
            after_data.append((node.header.split("|")[3].split("_")[1], node.header))

        after_data.sort(key = lambda x: x[0])
        after_true_clusters = []

        current_cluster = []
        for child_index, child_full in true_cluster_raw:
            if not current_cluster:
                current_cluster.append(child_full)
                current_index = child_index
            else:
                if child_index - current_index <= args.true_cluster_threshold:
                    current_cluster.append(child_full)
                    current_index = child_index
                else:
                    after_true_clusters.append(current_cluster)
                    current_cluster = [child_full]
                    current_index = child_index
        
        if current_cluster:
            after_true_clusters.append(current_cluster)

        after_true_cluster = set(max(after_true_clusters, key=lambda x: len(x)))
        matches = []
        for i, cluster in enumerate(before_true_clusters):
            matches.append((len(after_true_cluster.intersection(set(cluster))) / len( cluster), i))
        matches.sort(key=lambda x: x[0], reverse= True)
        
        best_match = max(matches, key=lambda x: x[0])[1]
        this_rescues = [f"Rescued in gene: {gene}"]
        if False: # Toggle: Rescuing
            for header in before_true_clusters[best_match]:
                if header in kicked_headers:
                    kicked_headers.remove(header)
                    this_rescues.append(header)

        if len(this_rescues) > 1:
            rescues.append("\n".join(this_rescues))

        reported.extend(report_overlaps(nodes, before_true_clusters[best_match]))

        aa_output_after_kick = sum(
            1
            for i in aa_output
            if i[0] not in kicked_headers and not i[0].endswith(".")
        )

        if not aa_output_after_kick:
            kicked_genes.append(f"No valid sequences after kick: {gene}")
            continue

        # Output detailed debug information
        if args.debug == 2:
            nodes.sort(key=lambda x: x.start)
            with open(aa_out.strip(".gz"), "w") as f:
                for header, sequence in aa_output:
                    if header.endswith("."):
                        f.write(f">{header}\n{sequence}\n")

                f.write(f">Cand_Consensus\n{x_cand_consensus}\n")

                for node in nodes:
                    if node is None:
                        continue
                    is_kick = (
                        "_KICKED" if node.kick or node.header in kicked_headers else ""
                    )
                    f.write(f">{node.contig_header()}{is_kick}\n{node.sequence}\n")
        else:
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

        count = len(kicked_headers)
        if args.debug:
            kicks.append(f"Total Kicks: {count}\n")
        passed_total += aa_output_after_kick
        total += count

        if has_x_before:
            has_x_genes += has_x_before

        if has_x_after:
            has_x_genes_after += has_x_after

    if args.debug:
        return True, kicks, total, kicked_genes, consensus_kicks, passed_total, rescues, reported, reported_regions, has_x_genes, has_x_genes_after, region_kicks, internal_kicks

    return True, [], total, kicked_genes, consensus_kicks, passed_total, rescues, reported, reported_regions, has_x_genes, has_x_genes_after, region_kicks, internal_kicks


def main(args, from_folder):
    if not (0 < args.matching_percent < 1.0):
        if 0 < args.matching_percent <= 100:
            args.matching_percent = args.matching_percent / 100
        else:
            raise ValueError(
                "Cannot convert matching percent threshold to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.matching_consensus_percent < 1.0):
        if 0 < args.matching_consensus_percent <= 100:
            args.matching_consensus_percent = args.matching_consensus_percent / 100
        else:
            raise ValueError(
                "Cannot convert matching consensus percent threshold to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.gross_diference_percent < 1.0):
        if 0 < args.gross_diference_percent <= 100:
            args.gross_diference_percent = args.gross_diference_percent / 100
        else:
            raise ValueError(
                "Cannot convert gross difference percent threshold to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.rolling_matching_percent < 1.0):
        if 0 < args.rolling_matching_percent <= 100:
            args.rolling_matching_percent = args.rolling_matching_percent / 100
        else:
            raise ValueError(
                "Cannot convert rolling matching percent threshold to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.rolling_consensus_percent < 1.0):
        if 0 < args.rolling_consensus_percent <= 100:
            args.rolling_consensus_percent = args.rolling_consensus_percent / 100
        else:
            raise ValueError(
                "Cannot convert rolling consensus percent threshold to a percent. Use a decimal or a whole number between 0 and 100"
            )
            

    this_args = CollapserArgs(
        compress=args.compress,
        uncompress_intermediates=args.uncompress_intermediates,
        processes=args.processes,
        merge_overlap=args.merge_overlap,
        overlap_percent=args.kick_overlap,
        matching_percent=args.matching_percent,
        verbose=args.verbose,
        debug=args.debug,
        consensus=args.consensus,
        matching_consensus_percent=args.matching_consensus_percent,
        gross_diference_percent=args.gross_diference_percent,
        from_folder=from_folder,
        rolling_matching_percent=args.rolling_matching_percent,
        rolling_consensus_percent=args.rolling_consensus_percent,
        rolling_window_size=args.rolling_window_size,
        true_cluster_threshold=args.true_cluster_threshold,
    )
    return do_folder(this_args, args.INPUT)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
