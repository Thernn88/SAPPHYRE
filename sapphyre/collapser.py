from collections import Counter, defaultdict
from itertools import combinations
from math import ceil
from multiprocessing import Pool
from shutil import rmtree
import os

import blosum as bl

from msgspec import Struct
from phymmr_tools import constrained_distance, find_index_pair, get_overlap, is_same_kmer
from .timekeeper import KeeperMode, TimeKeeper
import wrap_rocks
from .utils import writeFasta, parseFasta, printv

class CollapserArgs(Struct):
    compress: bool
    uncompress_intermediates: bool
    processes: int

    merge_overlap: int
    read_percent: float
    contig_percent: float

    required_read_percent: float
    required_contig_percent: float
    keep_read_percent: float

    sub_percent: float

    verbose: int
    debug: int
    matching_consensus_percent: float


class BatchArgs(Struct):
    args: CollapserArgs
    genes: list
    nt_input_path: str
    nt_out_path: str
    aa_input_path: str
    aa_out_path: str
    compress: bool
    is_assembly: bool
    matching_consensus_percent: float

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


def average_match(seq_a, consensus, start, end):
    match = 0
    total = 0
    for i in range(start, end):
        if consensus[i].count("-") / len(consensus[i]) > 0.5:
            continue

        total += 1

        if seq_a[i] in consensus[i]:
            match += 1

    if total == 0:
        return 0

    return match / total

def do_folder(args, input_path):
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    nt_input_path = os.path.join(input_path, "outlier", "blosum", "nt")
    aa_input_path = os.path.join(input_path, "outlier", "blosum", "aa")

    collapsed_path = os.path.join(input_path, "outlier", "collapsed")
    nt_out_path = os.path.join(collapsed_path, "nt")
    aa_out_path = os.path.join(collapsed_path, "aa")

    # Reset folders
    if os.path.exists(collapsed_path):
        rmtree(collapsed_path)

    os.mkdir(collapsed_path)
    os.mkdir(nt_out_path)
    os.mkdir(aa_out_path)

    nt_db_path = os.path.join(input_path, "rocksdb", "sequences", "nt")
    if os.path.exists(nt_db_path):
        nt_db = wrap_rocks.RocksDB(nt_db_path)
        dbis_assembly = nt_db.get("get:isassembly")
        is_assembly = False
        if dbis_assembly and dbis_assembly == "True":
            is_assembly = True
        del nt_db

    # Process NT
    genes = [
        gene
        for gene in os.listdir(nt_input_path)
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
            args.matching_consensus_percent,
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
    

    if args.debug:
        total_kicks = sum(i[2] for i in results if i[2] != 0)
        kicked_genes = "\n".join(["\n".join(i[3]) for i in results])
        kicked_consensus = "".join(["".join(i[4]) for i in results])

        with open(os.path.join(collapsed_path, "kicked_genes.txt"), "w") as fp:
            fp.write(kicked_genes)

        with open(os.path.join(collapsed_path, "kicked_consensus.txt"), "w") as fp:
            fp.write(kicked_consensus)

        with open(os.path.join(collapsed_path, "kicks.txt"), "w") as fp:
            fp.write(f"Total Kicks: {total_kicks}\n")
            for data in results:
                kick_list = data[1]
                fp.write("".join(kick_list))

                

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose, 1)

    return all_passed


def blosum_sub_merge(mat, overlap_coords, overlap_amount, aa_sequence, aa_sequence_2):
    allow_sub = True
    BLOSUM_LIMIT = 2
    blosum_positions = 0
    for k in range(0, overlap_amount, 3):
        nt_pos = overlap_coords[0] + k
        aa_pos = (nt_pos) // 3
        aa_node_bp = aa_sequence[aa_pos]
        aa_other_bp = aa_sequence_2[aa_pos]

        if aa_node_bp == aa_other_bp:
            continue

        subs = mat[aa_node_bp][aa_other_bp]
        if subs >= 0:
            allow_sub = True
            blosum_positions += 1
            if blosum_positions >= BLOSUM_LIMIT:
                allow_sub = False
                break
        if subs < 0:
            allow_sub = False
            break
    
    return allow_sub


def get_score(seq, start, end, consensus_dict):
    total = 0
    score = 0
    for i in range(start, end):
        let = seq[i]
        if consensus_dict[i] == set():
            continue

        if let == "-":
            continue

        if seq[i] in consensus_dict[i]:
            score += 1
        total += 1

    percent = score / total
    return percent


def process_batch(
    batch_args: BatchArgs,
):
    args = batch_args.args
    kicked_genes = []
    consensus_kicks = []
    total = 0

    kicks = []
    for gene in batch_args.genes:
        kicked_headers = set()
        printv(f"Doing: {gene}", args.verbose, 2)

        nt_out = os.path.join(batch_args.nt_out_path, gene)
        aa_out = os.path.join(batch_args.aa_out_path, gene.replace(".nt.", ".aa."))

        nt_in = os.path.join(batch_args.nt_input_path, gene)
        aa_in = os.path.join(batch_args.aa_input_path, gene.replace(".nt.", ".aa."))

        aa_sequences = parseFasta(aa_in)
        
        aa_output = []

        nodes = []


        # make nodes out of nt_input for processing
        for header, sequence in aa_sequences:
            aa_output.append((header, sequence))
            if header.endswith("."):
                continue

            start,end = find_index_pair(sequence, "-")

            internal_gaps=  sequence[start:end].count("-")

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

        ref_alignments = [seq for header, seq in aa_output if header.endswith(".")]
        ref_consensus = {i: [seq[i] for seq in ref_alignments] for i in range(len(ref_alignments[0]))}
        del ref_alignments

        match_percent = args.matching_consensus_percent if batch_args.is_assembly else 0.6

        for read in nodes:
            average_matching_cols = average_match(
                read.sequence,
                ref_consensus,
                read.start,
                read.end,
            )

            if average_matching_cols < match_percent:
                if args.debug:
                    consensus_kicks.append(f"{gene},{read.header},{average_matching_cols},{read.length}\n")
                kicked_headers.add(read.header)
                read.kick = True

        # Rescurive scan
        splice_occured = True
        failed = defaultdict(dict)
        merges_occured = False
        while splice_occured:
            splice_occured = False
            for i, node in enumerate(nodes):
                if node is None or node.kick:
                    continue
                possible_extensions = []
                for j, node_2 in enumerate(nodes):
                    if node_2 is None or node_2.kick:
                        continue
                    if i == j:
                        continue
                    if failed[i].get(j, False):
                        continue

                    overlap_coords = get_overlap(node.start, node.end, node_2.start, node_2.end, args.merge_overlap)

                    if overlap_coords:
                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        overlap_coord = overlap_coords[0]
                        possible_extensions.append((overlap_amount, overlap_coord, j))
                for _, overlap_coord, j in sorted(possible_extensions, reverse=True, key = lambda x: x[0]):
                    if failed[i].get(j, False):
                        continue

                    node_2 = nodes[j]
                    if node_2 is None:
                        continue
                    #Confirm still overlaps
                    overlap_coords = get_overlap(node.start, node.end, node_2.start, node_2.end, args.merge_overlap)
                    if overlap_coords:
                        #Get distance
                        node_kmer = node.sequence[overlap_coords[0] : overlap_coords[1]]
                        other_kmer = node_2.sequence[overlap_coords[0] : overlap_coords[1]]
                        if is_same_kmer(node_kmer, other_kmer):
                            splice_occured = True
                            node.extend(node_2, overlap_coords[0])
                            merges_occured = True
                            nodes[j] = None
                            continue

                        failed[i][j] = True
                            
        valid_contigs = [node for node in nodes if node is not None and node.is_contig]

        if not merges_occured and not batch_args.is_assembly:
            read_alignments = [seq for header, seq in aa_output if not header.endswith(".")]
            read_consensus = {i: {seq[i] for seq in read_alignments if seq[i] != "-"} for i in range(len(read_alignments[0]))}

            total_cols = 0
            data_cols = 0

            for i, letters in read_consensus.items():
                total_cols += 1
                if letters:
                    data_cols += 1

            coverage = data_cols / total_cols

            del read_alignments
            del read_consensus

            if coverage < 0.1:
                kicked_genes.append(f"Coverage: {coverage} no contigs: {gene.split('.')[0]}")
                continue


 
        #contigs is a shallow copy of valid_contigs
        contigs = sorted(valid_contigs, key=lambda x: x.length, reverse=True)
        if args.debug:
            kicks.append(f"Kicks for {gene}\nHeader B,,Header A,Overlap Percent,Score A, Score B\n")
        if batch_args.is_assembly:
            reads = []
            for (i, contig_a), (j, contig_b) in combinations(enumerate(contigs), 2):
                if contig_a.kick or contig_b.kick:
                    continue

                overlap_coords = get_overlap(contig_a.start, contig_a.end, contig_b.start, contig_b.end, 1)
                if overlap_coords:
                    # this block can probably just be an overlap percent call
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    smaller = min(contig_a, contig_b, key=lambda x: x.length)
                    percent = overlap_amount / (smaller.length - smaller.internal_gaps)
                    # this block can probably just be an overlap percent call

                    if percent >= args.read_percent:

                        contig_a_score = get_score(contig_a.sequence, contig_a.start, contig_a.end, ref_consensus)
                        contig_b_score = get_score(contig_b.sequence, contig_b.start, contig_b.end, ref_consensus)

                        if contig_a_score > contig_b_score:
                            contig_b.kick = True
                            kicked_headers.add(contig_b.header)
                            kicked_headers.update(contig_b.children)
                            if args.debug:
                                kicks.append(
                                    f"{contig_b.contig_header()},Contig Kicked By,{contig_a.contig_header()},{percent},{contig_a_score},{contig_b_score}\n"
                                )
                        elif contig_b_score > contig_a_score:
                            contig_a.kick = True
                            kicked_headers.add(contig_a.header)
                            kicked_headers.update(contig_a.children)
                            if args.debug:
                                kicks.append(
                                    f"{contig_a.contig_header()},Contig Kicked By,{contig_b.contig_header()},{percent},{contig_a_score},{contig_b_score}\n"
                                )
        else:
            reads = [node for node in nodes if node is not None and not node.is_contig]

            for (i, contig_a), (j, contig_b) in combinations(enumerate(contigs), 2):
                if contig_a.kick or contig_b.kick:
                    continue
                # this block can probably just be an overlap percent call
                overlap_coords = get_overlap(contig_a.start, contig_a.end, contig_b.start, contig_b.end, 1)
                if overlap_coords:
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    percent = overlap_amount / (contig_b.length - contig_b.internal_gaps)
                # this block can probably just be an overlap percent call

                    if percent >= args.contig_percent:
                        is_kick, matching_percent = contig_a.is_kick(
                            contig_b,
                            overlap_coords,
                            args.required_contig_percent,
                            overlap_amount,
                        )
                        if is_kick:
                            contig_b.kick = True
                            if args.debug:
                                kicks.append(
                                    f"{contig_b.contig_header()},Contig Kicked By,{contig_a.contig_header()},{percent},{matching_percent}\n"
                                )
                            kicked_headers.add(contig_b.header)
                            kicked_headers.update(contig_b.children)
            contigs = [contig for contig in contigs if not contig.kick]

        
            if contigs:
                nodes = contigs + reads
            else:
                nodes = reads

            nodes.sort(key = lambda x: x.length, reverse=True)

            for i, node_kick in enumerate(nodes):
                kick = False
                keep = False
                for j, node_2 in enumerate(nodes):
                    if i == j:
                        continue
                    if node_kick.length > node_2.length:
                        continue
                    if node_2.kick:
                        continue

                    overlap_coords = get_overlap(node_2.start, node_2.end, node_kick.start, node_kick.end, 1)
                    if overlap_coords:
                        # this block can probably just be an overlap percent call
                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        percent = overlap_amount / (node_kick.length - node_kick.internal_gaps)
                        # this block can probably just be an overlap percent call

                        if percent >= args.read_percent:
                            is_kick, matching_percent = node_kick.is_kick(
                                node_2,
                                overlap_coords,
                                args.required_read_percent,
                                overlap_amount,
                            )

                            if is_kick:
                                kick_percent = percent
                                kick_matching_percent = matching_percent
                                kick_parent = node_2.contig_header()
                                kick = True

                            if matching_percent >= args.keep_read_percent:
                                keep = True
                                break
                if kick and not keep:
                    node_kick.kick = True
                    if args.debug:
                        kicks.append(
                            f"{node_kick.contig_header()},Kicked By,{kick_parent},{kick_percent},{kick_matching_percent}\n"
                        )
                    kicked_headers.add(node_kick.header)
                    if node_kick.is_contig:
                        kicked_headers.update(node_kick.children)

        aa_output_after_kick = sum(1 for i in aa_output if i[0] not in kicked_headers and not i[0].endswith(".")) > 0

        if not aa_output_after_kick:
            kicked_genes.append(f"No valid sequences after kick: {gene.split('.')[0]}")
            continue

        if args.debug == 2:
            output = valid_contigs + reads
            output.sort(key=lambda x: x.start)
            with open(aa_out, "w") as f:
                for header, sequence in aa_output:
                    if header.endswith("."):
                        f.write(f">{header}\n{sequence}\n")

                for node in output:
                    if node is None:
                        continue
                    is_kick = (
                        "_KICKED" if node.kick or node.header in kicked_headers else ""
                    )
                    f.write(f">{node.contig_header()}{is_kick}\n{node.sequence}\n")
        else:
            aa_output = [pair for pair in aa_output if pair[0] not in kicked_headers]
            writeFasta(
                aa_out, aa_output, batch_args.compress
            )

        nt_sequences = [
            (header, sequence)
            for header, sequence in parseFasta(nt_in)
            if header not in kicked_headers
        ]
        writeFasta(nt_out, nt_sequences, batch_args.compress)

        count = len(kicked_headers)
        if args.debug:
            kicks.append(f"Total Kicks: {count}\n")
        total += count

    if args.debug:
        return True, kicks, total, kicked_genes, consensus_kicks

    return True, [], 0, kicked_genes, consensus_kicks


def main(args):
    this_args = CollapserArgs(
        compress=args.compress,
        uncompress_intermediates=args.uncompress_intermediates,
        processes=args.processes,
        merge_overlap=args.merge_overlap,
        read_percent=args.read_overlap,
        contig_percent=args.contig_overlap,
        required_read_percent=args.read_matching_percent,
        required_contig_percent=args.contig_matching_percent,
        keep_read_percent=args.keep_read_percent,
        sub_percent=args.sub_percent,
        verbose=args.verbose,
        debug=args.debug,
        matching_consensus_percent = args.matching_consensus_percent,
    )
    return do_folder(this_args, args.INPUT)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
