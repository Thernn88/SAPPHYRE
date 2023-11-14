from collections import defaultdict
from itertools import combinations
from math import ceil
from multiprocessing import Pool
from os import path, mkdir, listdir

from msgspec import Struct
from phymmr_tools import (
    constrained_distance,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    dumb_consensus,
)
from wrap_rocks import RocksDB
from .timekeeper import KeeperMode, TimeKeeper
from .utils import writeFasta, parseFasta, printv


class CollapserArgs(Struct):
    compress: bool
    uncompress_intermediates: bool
    processes: int

    merge_overlap: int
    matching_percent: float
    overlap_percent: float

    verbose: int
    debug: int
    matching_consensus_percent: float
    gross_diference_percent: float

    from_folder: str


class BatchArgs(Struct):
    args: CollapserArgs
    genes: list
    nt_input_path: str
    nt_out_path: str
    aa_input_path: str
    aa_out_path: str
    compress: bool
    is_assembly: bool


class NODE(Struct):
    header: str
    sequence: str
    start: int
    end: int
    internal_gaps: int
    length: int
    children: list
    is_contig: bool
    kick: bool

    def get_extension(self, node_2, overlap_coord):
        if node_2.start >= self.start and node_2.end <= self.end:
            return 0
        elif self.start >= node_2.start and self.end <= node_2.end:
            return (overlap_coord - node_2.start) + (node_2.end - self.end)
        elif node_2.start >= self.start:
            return node_2.start - overlap_coord
        else:
            return overlap_coord - node_2.start

    def extend(self, node_2, overlap_coord):
        if node_2.start >= self.start and node_2.end <= self.end:
            self.sequence = (
                self.sequence[:overlap_coord]
                + node_2.sequence[overlap_coord : node_2.end]
                + self.sequence[node_2.end :]
            )

        elif self.start >= node_2.start and self.end <= node_2.end:
            self.sequence = (
                node_2.sequence[:overlap_coord]
                + self.sequence[overlap_coord : self.end]
                + node_2.sequence[self.end :]
            )

            self.start = node_2.start
            self.end = node_2.end

        elif node_2.start >= self.start:
            self.sequence = (
                self.sequence[:overlap_coord] + node_2.sequence[overlap_coord:]
            )

            self.end = node_2.end

        else:
            self.sequence = (
                node_2.sequence[:overlap_coord] + self.sequence[overlap_coord:]
            )

            self.start = node_2.start
        self.length = self.end - self.start

        self.children.append(node_2.header)
        self.children.extend(node_2.children)
        self.is_contig = True

    def is_kick(self, node_2, overlap_coords, kick_percent, overlap_amount):
        kmer_current = self.sequence[overlap_coords[0] : overlap_coords[1]]
        kmer_next = node_2.sequence[overlap_coords[0] : overlap_coords[1]]
        non_matching_chars = constrained_distance(kmer_current, kmer_next)
        non_mathching_percent = non_matching_chars / overlap_amount
        matching_percent = 1 - non_mathching_percent

        return matching_percent < kick_percent, matching_percent

    def contig_header(self):
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


def do_folder(args, input_path):
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    collapsed_path = path.join(input_path, "outlier", "collapsed")
    nt_out_path = path.join(collapsed_path, "nt")
    aa_out_path = path.join(collapsed_path, "aa")

    mkdir(collapsed_path)
    mkdir(nt_out_path)
    mkdir(aa_out_path)

    nt_db_path = path.join(input_path, "rocksdb", "sequences", "nt")
    is_assembly = False
    if path.exists(nt_db_path):
        nt_db = RocksDB(nt_db_path)
        dbis_assembly = nt_db.get("get:isassembly")

        if dbis_assembly and dbis_assembly == "True":
            is_assembly = True
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

    batched_arguments = [
        BatchArgs(
            args,
            genes[i : i + per_thread],
            nt_input_path,
            nt_out_path,
            aa_input_path,
            aa_out_path,
            compress,
            is_assembly,
        )
        for i in range(0, len(genes), per_thread)
    ]

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

    if args.debug:
        kicked_genes = "\n".join(["\n".join(i[3]) for i in results])
        genes_kicked_count = len(kicked_genes.split("\n"))
        kicked_consensus = "".join(["".join(i[4]) for i in results])

        with open(path.join(collapsed_path, "kicked_genes.txt"), "w") as fp:
            fp.write(kicked_genes)

        with open(path.join(collapsed_path, "kicked_consensus.txt"), "w") as fp:
            fp.write(kicked_consensus)

        with open(path.join(collapsed_path, "kicks.txt"), "w") as fp:
            fp.write(f"Total Kicks: {total_kicks}\n")
            for data in results:
                kick_list = data[1]
                fp.write("".join(kick_list))
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
    aa_output, match_percent, nodes, kicked_headers, consensus_kicks, debug, gene
):
    ref_average_data_length = []
    ref_consensus = defaultdict(list)
    reference_seqs = [seq for header, seq in aa_output if header.endswith(".")]
    ref_consensus_seq = dumb_consensus(reference_seqs, 0.5)
    for seq in reference_seqs:
        start, end = find_index_pair(seq, "-")
        for i in range(start, end):
            ref_consensus[i].append(seq[i])

        ref_average_data_length.append(len(seq) - seq.count("-"))

    ref_average_data_length = sum(ref_average_data_length) / len(
        ref_average_data_length
    )

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

    nodes = list(filter(lambda x: not x.kick, nodes))

    return ref_average_data_length, nodes, ref_consensus_seq


def merge_overlapping_reads(nodes, minimum_overlap):
    # Rescurive scan
    for i, node in enumerate(nodes):
        if node is None or node.kick:
            continue
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
                    # Get distance

                    fail = False
                    for x in range(overlap_coords[0], overlap_coords[1]):
                        if node.sequence[x] != node_2.sequence[x]:
                            fail = True
                            break

                    if not fail:
                        splice_occured = True
                        node.extend(node_2, overlap_coords[0])
                        nodes[j] = None

    return list(filter(lambda x: x is not None, nodes))


def get_coverage(nodes, ref_average_data_length):
    read_alignments = [node.sequence for node in nodes]
    data_cols = 0

    for i in range(len(read_alignments[0])):
        if any(seq[i] != "-" for seq in read_alignments):
            data_cols += 1

    return data_cols / ref_average_data_length


def kick_overlapping_reads(
    nodes,
    min_overlap_percent,
    required_matching_percent,
    gross_difference_percent,
    debug,
):
    kicked_headers = set()
    kicks = []
    for i, node_kick in enumerate(node for node in nodes if not node.kick):
        for j, node_2 in enumerate(node for node in nodes if not node.kick):
            if i == j:
                continue
            if node_2.length < node_kick.length:
                continue

            overlap_coords = get_overlap(
                node_2.start, node_2.end, node_kick.start, node_kick.end, 1
            )
            if overlap_coords:
                # this block can probably just be an overlap percent call
                overlap_amount = overlap_coords[1] - overlap_coords[0]
                percent = overlap_amount / (node_kick.length - node_kick.internal_gaps)
                # this block can probably just be an overlap percent call

                if percent >= min_overlap_percent:
                    is_kick, matching_percent = node_kick.is_kick(
                        node_2,
                        overlap_coords,
                        required_matching_percent,
                        overlap_amount,
                    )

                    length_percent = ""
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
    nodes, ref_consensus_seq, kicked_headers, consensus_kicks, debug, gene
):
    CONSENSUS_PERCENT = 0.5
    WINDOW_SIZE = 14
    WINDOW_MATCHING_PERCENT = 0.7
    STEP = 1
    cand_consensus = dumb_consensus(
        [node.sequence for node in nodes], CONSENSUS_PERCENT
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
            WINDOW_SIZE,
            WINDOW_MATCHING_PERCENT,
            STEP,
        )
        if is_kick:
            if debug:
                consensus_kicks.append(
                    f"{gene},{node.header},{window_start}:{window_start+WINDOW_SIZE}\nCand: {node.sequence[window_start: window_start+WINDOW_SIZE]}\nCons: {cand_consensus[window_start: window_start+WINDOW_SIZE]}\nRefc: {ref_consensus_seq[window_start: window_start+WINDOW_SIZE]}\nMatch: {matching_percent},Ambig: {ambig_percent}\n"
                )
            kicked_headers.add(node.header)
            node.kick = True

    nodes = list(filter(lambda x: not x.kick, nodes))
    return nodes


def process_batch(
    batch_args: BatchArgs,
):
    args = batch_args.args
    kicked_genes = []
    consensus_kicks = []
    total = 0
    passed_total = 0

    kicks = []
    for gene in batch_args.genes:
        kicked_headers = set()
        printv(f"Doing: {gene}", args.verbose, 2)

        nt_out = path.join(batch_args.nt_out_path, gene)
        aa_out = path.join(batch_args.aa_out_path, gene.replace(".nt.", ".aa."))

        nt_in = path.join(batch_args.nt_input_path, gene)
        aa_in = path.join(batch_args.aa_input_path, gene.replace(".nt.", ".aa."))

        aa_sequences = parseFasta(aa_in)

        aa_output = []

        nodes = []

        # make nodes out of nt_input for processing
        aa_count = 0
        for header, sequence in aa_sequences:
            aa_output.append((header, sequence))
            if header.endswith("."):
                continue

            aa_count += 1

            start, end = find_index_pair(sequence, "-")

            internal_gaps = sequence[start:end].count("-")

            node_is_contig = batch_args.is_assembly or "&&" in header

            nodes.append(
                NODE(
                    header=header,
                    sequence=sequence,
                    start=start,
                    end=end,
                    length=(end - start),
                    internal_gaps=internal_gaps,
                    children=[],
                    is_contig=node_is_contig,
                    kick=False,
                )
            )

        ref_average_data_length, nodes, ref_consensus_seq = kick_read_consensus(
            aa_output,
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
                f"No valid sequences after consensus: {gene.split('.')[0]}"
            )
            continue

        nodes = kick_rolling_consensus(
            nodes, ref_consensus_seq, kicked_headers, consensus_kicks, args.debug, gene
        )

        if not nodes:
            total += aa_count
            kicked_genes.append(
                f"No valid sequences after rolling candidate consensus: {gene.split('.')[0]}"
            )
            continue

        printv("Merging Overlapping Reads", args.verbose, 3)
        nodes = merge_overlapping_reads(nodes, args.merge_overlap)

        printv("Calculating Coverage", args.verbose, 3)
        coverage = get_coverage(nodes, ref_average_data_length)

        req_coverage = 0.3 if batch_args.is_assembly else 0.01
        if coverage < req_coverage:
            total += aa_count
            kicked_genes.append(f"{gene} -> failed due to Coverage: {coverage}")
            continue

        if args.debug:
            kicks.append(
                f"Kicks for {gene}\nHeader B,,Header A,Overlap Percent,Matching Percent,Length Ratio\n"
            )
        nodes.sort(key=lambda x: x.length, reverse=True)

        printv("Kicking Overlapping Reads", args.verbose, 3)
        this_kicks, this_debug = kick_overlapping_reads(
            nodes,
            args.overlap_percent,
            args.matching_percent,
            args.gross_diference_percent,
            args.debug,
        )

        if args.debug:
            kicks.extend(this_debug)

        kicked_headers.update(this_kicks)

        aa_output_after_kick = sum(
            1
            for i in aa_output
            if i[0] not in kicked_headers and not i[0].endswith(".")
        )

        if not aa_output_after_kick:
            kicked_genes.append(f"No valid sequences after kick: {gene.split('.')[0]}")
            continue

        if args.debug == 2:
            nodes.sort(key=lambda x: x.start)
            with open(aa_out.strip(".gz"), "w") as f:
                for header, sequence in aa_output:
                    if header.endswith("."):
                        f.write(f">{header}\n{sequence}\n")

                for node in nodes:
                    if node is None:
                        continue
                    is_kick = (
                        "_KICKED" if node.kick or node.header in kicked_headers else ""
                    )
                    f.write(f">{node.contig_header()}{is_kick}\n{node.sequence}\n")
        else:
            aa_output = [pair for pair in aa_output if pair[0] not in kicked_headers]
            writeFasta(aa_out, aa_output, batch_args.compress)

        nt_sequences = [
            (header, sequence)
            for header, sequence in parseFasta(nt_in)
            if header not in kicked_headers
        ]
        writeFasta(nt_out, nt_sequences, batch_args.compress)

        count = len(kicked_headers)
        if args.debug:
            kicks.append(f"Total Kicks: {count}\n")
        passed_total += aa_output_after_kick
        total += count

    if args.debug:
        return True, kicks, total, kicked_genes, consensus_kicks, passed_total

    return True, [], total, kicked_genes, consensus_kicks, passed_total


def main(args, from_folder):
    this_args = CollapserArgs(
        compress=args.compress,
        uncompress_intermediates=args.uncompress_intermediates,
        processes=args.processes,
        merge_overlap=args.merge_overlap,
        overlap_percent=args.kick_overlap,
        matching_percent=args.matching_percent,
        verbose=args.verbose,
        debug=args.debug,
        matching_consensus_percent=args.matching_consensus_percent,
        gross_diference_percent=args.gross_diference_percent,
        from_folder=from_folder,
    )
    return do_folder(this_args, args.INPUT)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
