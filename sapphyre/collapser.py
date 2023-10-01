from collections import defaultdict
from math import ceil
from multiprocessing import Pool
from shutil import rmtree
from os import path, mkdir, listdir

from msgspec import Struct
from phymmr_tools import constrained_distance, find_index_pair, get_overlap
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
    nt_input_path = path.join(input_path, "outlier", "internal", "nt")
    aa_input_path = path.join(input_path, "outlier", "internal", "aa")

    collapsed_path = path.join(input_path, "outlier", "collapsed")
    nt_out_path = path.join(collapsed_path, "nt")
    aa_out_path = path.join(collapsed_path, "aa")

    # Reset folders
    if path.exists(collapsed_path):
        rmtree(collapsed_path)

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
    printv(f"Kicked {total_kicks} sequences and wrote {total_sequences} sequences", args.verbose, 1)

                

    printv(f"Done! Took {time_keeper.differential():.2f} seconds", args.verbose, 1)

    return all_passed

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

        ref_average_data_length = []
        ref_consensus = defaultdict(list)
        for header, seq in aa_output:
            if header.endswith("."):
                start, end = find_index_pair(seq, "-")
                for i in range(start, end):
                    ref_consensus[i].append(seq[i])
                ref_average_data_length.append(len(seq) - seq.count("-"))

        ref_average_data_length = sum(ref_average_data_length) / len(ref_average_data_length)

        match_percent = args.matching_consensus_percent if not batch_args.is_assembly else 0.6

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

        read_alignments = [node.sequence for node in nodes if not node.kick]

        if not read_alignments:
            total += aa_count
            kicked_genes.append(f"No valid sequences after consensus: {gene.split('.')[0]}")
            continue

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

                    overlap_coords = get_overlap(node.start, node.end, node_2.start, node_2.end, args.merge_overlap)

                    if overlap_coords:
                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        overlap_coord = overlap_coords[0]
                        possible_extensions.append((overlap_amount, overlap_coord, j))
                for _, overlap_coord, j in sorted(possible_extensions, reverse=False, key = lambda x: x[0]):

                    node_2 = nodes[j]
                    if node_2 is None:
                        continue
                    #Confirm still overlaps
                    overlap_coords = get_overlap(node.start, node.end, node_2.start, node_2.end, args.merge_overlap)
                    if overlap_coords:
                        #Get distance

                        
                        fail = False
                        for x in range(overlap_coords[0], overlap_coords[1]):
                            if node.sequence[x] != node_2.sequence[x]:
                                fail = True
                                break

                        if not fail:
                            splice_occured = True
                            node.extend(node_2, overlap_coords[0])
                            nodes[j] = None
        nodes = [node for node in nodes if node is not None]

        data_cols = 0

        for i in range(len(read_alignments[0])):
            if any(seq[i] != "-" for seq in read_alignments):
                data_cols += 1
            

        coverage = data_cols / ref_average_data_length

        
        req_coverage = 0.3 if batch_args.is_assembly else 0.1
        if coverage < req_coverage:
            total += len(read_alignments)
            kicked_genes.append(f"{gene} -> failed due to Coverage: {coverage}, Ref average columns: {ref_average_data_length}, Data columns: {data_cols}")
            continue
        del read_alignments

        if args.debug:
            kicks.append(f"Kicks for {gene}\nHeader B,,Header A,Overlap Percent,Matching Percent,Length Ratio\n")
        nodes.sort(key = lambda x: x.length, reverse=True)

        for i, node_kick in enumerate(node for node in nodes if not node.kick):
            for j, node_2 in enumerate(node for node in nodes if not node.kick):
                if i == j:
                    continue
                if node_2.length < node_kick.length:
                    continue
                
                overlap_coords = get_overlap(node_2.start, node_2.end, node_kick.start, node_kick.end, 1)
                if overlap_coords:
                    # this block can probably just be an overlap percent call
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    percent = overlap_amount / (node_kick.length - node_kick.internal_gaps)
                    # this block can probably just be an overlap percent call

                    if percent >= args.overlap_percent:
                        is_kick, matching_percent = node_kick.is_kick(
                            node_2,
                            overlap_coords,
                            args.matching_percent,
                            overlap_amount,
                        )

                        length_percent = ""
                        if not is_kick and node_kick.is_contig and node_2.is_contig:
                            length_percent = min(node_kick.length, node_2.length) / max(node_kick.length, node_2.length)
                    
                            if not is_kick and length_percent <= 0.15 and matching_percent < args.gross_diference_percent:
                                is_kick = True

                        if is_kick:
                            node_kick.kick = True
                            if args.debug:
                                kicks.append(
                                    f"{node_kick.contig_header()},Kicked By,{node_2.contig_header()},{percent},{matching_percent},{length_percent}\n"
                                )
                            kicked_headers.add(node_kick.header)
                            if node_kick.is_contig:
                                kicked_headers.update(node_kick.children)
                            break

        aa_output_after_kick = sum(1 for i in aa_output if i[0] not in kicked_headers and not i[0].endswith("."))

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
        passed_total += aa_output_after_kick
        total += count

    if args.debug:
        return True, kicks, total, kicked_genes, consensus_kicks, passed_total

    return True, [], total, kicked_genes, consensus_kicks, passed_total


def main(args):
    this_args = CollapserArgs(
        compress=args.compress,
        uncompress_intermediates=args.uncompress_intermediates,
        processes=args.processes,
        merge_overlap=args.merge_overlap,
        overlap_percent=args.kick_overlap,
        matching_percent=args.matching_percent,
        verbose=args.verbose,
        debug=args.debug,
        matching_consensus_percent = args.matching_consensus_percent,
        gross_diference_percent = args.gross_diference_percent,
    )
    return do_folder(this_args, args.INPUT)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )
