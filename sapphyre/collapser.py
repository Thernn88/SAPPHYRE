from itertools import combinations
from math import ceil
from multiprocessing import Pool
from shutil import rmtree

import blosum as bl

from msgspec import Struct
from phymmr_tools import constrained_distance
from .timekeeper import KeeperMode, TimeKeeper
from .utils import writeFasta, parseFasta, printv
import os

class CollapserArgs(Struct):
    compress: bool
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

class NODE(Struct):
    header: str
    sequence: str
    start: int
    end: int
    length: int
    children: list
    is_contig: bool
    kick: bool
    splices: dict

    def get_overlap(self, node_2, min_overlap=0):
        overlap_start = max(self.start, node_2.start)
        overlap_end = min(self.end, node_2.end)
        if overlap_end - overlap_start >= min_overlap:
            return overlap_start, overlap_end
        
        return None
    
    def get_sequence_at_coord(self, position):
        seq_at_position = None
        for header, coord in self.splices.items():
            if position >= coord:
                seq_at_position = header
            else:
                break
        return seq_at_position

    def extend(self, node_2, overlap_coord):
        self.sequence = self.sequence[:overlap_coord] + node_2.sequence[overlap_coord:]
        self.splices[node_2.header] = overlap_coord
        self.end = node_2.end
        self.length = self.end - self.start
        self.children.append(node_2.header)
        self.children.extend(node_2.children)
        self.is_contig = True

    def is_kick(self, node_2, overlap_coords, kick_percent, overlap_amount):
        kmer_current = self.sequence[overlap_coords[0]:overlap_coords[1]]
        kmer_next = node_2.sequence[overlap_coords[0]:overlap_coords[1]]
        non_matching_chars = hamming_distance(kmer_current, kmer_next)
        non_mathching_percent = (non_matching_chars / overlap_amount)
        matching_percent = 1 - non_mathching_percent

        return matching_percent < kick_percent, matching_percent
    
    def contig_header(self):
        contig_node = self.header.split("|")[3]
        children_nodes = "|".join([i.split("|")[3] for i in self.children])
        return f"CONTIG_{contig_node}|{children_nodes}"

def hamming_distance(seq1, seq2):
    return sum([1 for i in range(len(seq1)) if seq1[i] != seq2[i]])

def get_start_end(seq):
    start = 0
    end = len(seq)
    for i in range(len(seq)):
        if seq[i] != "-":
            start = i
            break
    for i in range(len(seq)-1, -1, -1):
        if seq[i] != "-":
            end = i
            break
    return start, end

def do_folder(args, input_path):
    printv(f"Processing: {os.path.basename(input_path)}", args.verbose, 1)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    nt_input_path = os.path.join(input_path,"outlier","nt")
    aa_input_path = os.path.join(input_path,"outlier","aa")

    collapsed_path = os.path.join(input_path,"collapsed")
    nt_out_path = os.path.join(collapsed_path,"nt")
    aa_out_path = os.path.join(collapsed_path,"aa")

    #Reset folders
    if os.path.exists(collapsed_path):
        rmtree(collapsed_path)

    os.mkdir(collapsed_path)
    os.mkdir(nt_out_path)
    os.mkdir(aa_out_path)

    # Process NT
    genes = [
        gene
        for gene in os.listdir(nt_input_path)
        if gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]

    per_thread = ceil(len(genes) / args.processes)

    batched_arguments = [(args, genes[i : i + per_thread], nt_input_path, nt_out_path, aa_input_path, aa_out_path) for i in range(0, len(genes), per_thread)]

    if args.processes <= 1:
        results = []
        for batch in batched_arguments:
            results.append(process_batch(*batch))
    else:
        with Pool(args.processes) as pool:
            results = pool.starmap(process_batch, batched_arguments)

    all_passed = all([i[0] for i in results])

    total_kicks = sum([i[2] for i in results if i[2] != 0])
    
    if args.debug:
        with open(os.path.join(collapsed_path, "kicks.txt"), "w") as fp:
            fp.write(f"Total Kicks: {total_kicks}\n")
            for data in results:
                kick_list = data[1]
                fp.write("".join(kick_list))

    printv(f"Took {time_keeper.differential():.2f}s to process {len(genes)} genes.", args.verbose, 1)

    return all_passed

def process_batch(args, genes, nt_input_path, nt_out_path, aa_input_path, aa_out_path):
    kicks = []
    for gene in genes:
        printv(f"Doing: {gene}", args.verbose, 2)

        nt_out = os.path.join(nt_out_path, gene)
        aa_out = os.path.join(aa_out_path, gene.replace('.nt.','.aa.'))

        nt_in = os.path.join(nt_input_path, gene)
        aa_in = os.path.join(aa_input_path, gene.replace('.nt.','.aa.'))

        nt_sequences = parseFasta(nt_in)
        aa_sequences = {header:sequence for header,sequence in parseFasta(aa_in)}
        nt_output = []

        nodes = []

        for header, sequence in nt_sequences:
            nt_output.append((header, sequence))
            if header.endswith("."):
                continue

            start,end = get_start_end(sequence)

            nodes.append(NODE(header=header, sequence=sequence, start=start, end=end, length=(end-start), children=[], is_contig=False, kick=False, splices={header: start}))
        mat = bl.BLOSUM(62)
        #Rescurive scan
        splice_occured = True
        while splice_occured:
            splice_occured = False
            for i, node in enumerate(nodes):
                if node.kick:
                    continue
                for j, node_2 in enumerate(nodes):
                    if node_2.kick:
                        continue
                    if i == j:
                        continue
                    overlap_coords = node.get_overlap(node_2, args.merge_overlap)
                        
                    if overlap_coords:
                        node_kmer = node.sequence[overlap_coords[0]:overlap_coords[1]]
                        other_kmer = node_2.sequence[overlap_coords[0]:overlap_coords[1]]                        

                        distance = constrained_distance(node_kmer, other_kmer)

                        if distance == 0:
                            splice_occured = True
                            node.extend(node_2, overlap_coords[0])
                            nodes[j].kick = True
                            continue
                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        diff_percent = distance / overlap_amount

                        if diff_percent <= args.sub_percent and overlap_amount != node_2.length:
                            allow_sub = True
                            for k in range(0, overlap_amount, 3):
                                nt_pos = overlap_coords[0] + k
                                aa_pos = (nt_pos) // 3

                                aa_node_bp = aa_sequences[node.get_sequence_at_coord(nt_pos)][aa_pos]
                                aa_other_bp = aa_sequences[node_2.get_sequence_at_coord(nt_pos)][aa_pos]

                                if aa_node_bp == aa_other_bp:
                                    continue

                                subs = mat[aa_node_bp][aa_other_bp]
                                if subs >= 0:
                                    allow_sub = True
                                if subs < 0:
                                    allow_sub = False

                            if allow_sub:
                                splice_occured = True
                                node.extend(node_2, overlap_coords[0])
                                node_2.kick = True

        contigs = [node for node in nodes if not node.kick and node.is_contig]
        contigs.sort(key=lambda x: x.length, reverse=True)

        reads = [node for node in nodes if not node.kick and not node.is_contig]
        kicks.append(f"Kicks for {gene}\n")
        kicks.append("Header B,,Header A,Overlap Percent,Matching Percent\n")
        kicked_headers = set()
                            
        for (i, contig_a), (j, contig_b) in combinations(enumerate(contigs), 2):
            if contig_a.kick or contig_b.kick:
                continue
            overlap_coords = contig_a.get_overlap(contig_b)
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]
                percent = overlap_amount / contig_b.length
                if percent >= args.contig_percent:
                    is_kick, matching_percent = contig_a.is_kick(contig_b, overlap_coords, args.required_contig_percent, overlap_amount)
                    if is_kick:
                        contig_b.kick = True
                        kicks.append(f"{contig_b.contig_header()},Contig Kicked By,{contig_a.contig_header()},{percent},{matching_percent}\n")
                        kicked_headers.add(contig_b.header)
                        kicked_headers.update(contig_b.children)
        
        contigs = [contig for contig in contigs if not contig.kick]

        for read in reads:
            keep = False
            kick = False
            for contig in contigs:
                overlap_coords = read.get_overlap(contig)

                if overlap_coords:
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    percent = overlap_amount / read.length
                    if percent >= args.read_percent:
                        is_kick, matching_percent = read.is_kick(contig, overlap_coords, args.required_read_percent, overlap_amount)
                        
                        if is_kick:
                            kick = True
                        
                        if matching_percent >= args.keep_read_percent:
                            keep = True
                            break
            
            if kick and not keep:
                kicks.append(f"{read.header},Kicked By,{contig.contig_header()},{percent},{matching_percent}\n")
                kicked_headers.add(read.header)
            if keep:
                kicks.append(f"{read.header},Saved By,{contig.contig_header()},{percent},{matching_percent}\n")

        if args.debug == 2:
            output = contigs+reads
            output.sort(key=lambda x: x.start)
            with open(nt_out, "w") as f:
                for header, sequence in nt_output:
                    if header.endswith('.'):
                        f.write(f">{header}\n{sequence}\n")

                for node in output:
                    if node is None:
                        continue
                    is_kick = "_KICKED"  if node.kick or node.header in kicked_headers else ""
                    if node.is_contig:
                        f.write(f">{node.contig_header()}{is_kick}\n{node.sequence}\n")
                        continue
                    f.write(f">{node.header}{is_kick}\n{node.sequence}\n")
        else:     
            writeFasta(nt_out, [i for i in nt_output if i[0] not in kicked_headers], args.compress)

        aa_sequences = [(header,sequence) for header,sequence in aa_sequences.items() if header not in kicked_headers]
        writeFasta(aa_out, aa_sequences, args.compress)

    count = len(kicked_headers)
    kicks.append(f"Total Kicks: {count}\n")

    if args.debug:
        return True, kicks, count

    return True, [], 0
def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    results = []
    #args.processes
    this_args = CollapserArgs(
        compress=args.compress,
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
    )

    for input_path in args.INPUT:
        results.append(do_folder(this_args, input_path))
        
    if len(args.INPUT) >= 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return all(results)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )