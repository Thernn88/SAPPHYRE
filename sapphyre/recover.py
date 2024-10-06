from collections import defaultdict
from itertools import combinations
from math import ceil, floor
from multiprocessing import Pool
from os import listdir, makedirs, path
from pathlib import Path
import copy
import warnings
from msgspec import Struct
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    constrained_distance,
    join_with_exclusions,
    join_triplets_with_exclusions,
)
from Bio import BiopythonWarning

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta


class NODE(Struct):
    header: str
    count: int
    nt_sequence: str
    sequence: str
    start: int
    end: int
    children: list
    codename: str
    is_contig: bool

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
        before = len(self.nt_sequence) - self.nt_sequence.count("-")
        # If node_2 is contained inside self
        if node_2.start >= self.start and node_2.end <= self.end:
            self.nt_sequence = (
                self.nt_sequence[:overlap_coord]
                + node_2.nt_sequence[overlap_coord : node_2.end]
                + self.nt_sequence[node_2.end :]
            )

        # If node_2 contains self
        elif self.start >= node_2.start and self.end <= node_2.end:
            self.nt_sequence = (
                node_2.nt_sequence[:overlap_coord]
                + self.nt_sequence[overlap_coord : self.end]
                + node_2.nt_sequence[self.end :]
            )

            self.start = node_2.start
            self.end = node_2.end
        # If node_2 is to the right of self
        elif node_2.start >= self.start:
            self.nt_sequence = (
                self.nt_sequence[:overlap_coord] + node_2.nt_sequence[overlap_coord:]
            )

            self.end = node_2.end
        # If node_2 is to the left of self
        else:
            self.nt_sequence = (
                node_2.nt_sequence[:overlap_coord] + self.nt_sequence[overlap_coord:]
            )

            self.start = node_2.start

        # Save node_2 and the children of node_2 to self
        self.children.append(node_2.header)
        self.children.extend(node_2.children)
        if (len(self.nt_sequence) - self.nt_sequence.count("-")) > before:
            self.is_contig = True

    def get_children(self):
        """
        Returns the children of the current node

        Returns:
        -------
            list: The children of the current node
        """
        return [self.header.split("|")[3]] + [i.split("|")[3] for i in self.children]

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

def scan_extend(node, nodes, i, merged, min_overlap_percent=0.15, min_overlap_chars=10):
    merge_occured = True
    while merge_occured:
        merge_occured = False
        possible_merges = []
        for j, node_b in enumerate(nodes):  # When merge occurs start again at the beginning with highest count
            if j in merged:
                continue
            if i == j:
                continue
            overlap_coords = get_overlap(node.start, node.end, node_b.start, node_b.end, 1)

            
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]

                # Calculate percent overlap and compare to minimum overlap
                overlap_percent = overlap_amount / ((node.end - node.start))
                required_overlap = max(min_overlap_chars, min((node.end - node.start), (node_b.end - node_b.start)) * min_overlap_percent)
                
                # Use whichever is greater: percentage overlap or 10 characters
                if overlap_amount < required_overlap:
                    continue

                kmer_a = node.nt_sequence[overlap_coords[0]:overlap_coords[1]]
                kmer_b = node_b.nt_sequence[overlap_coords[0]:overlap_coords[1]]

                if not is_same_kmer(kmer_a, kmer_b):
                    continue

                merged.add(j)

                overlap_coord = overlap_coords[0]
                merge_occured = True
                node.extend(node_b, overlap_coord)
                break
    return node

def simple_assembly(nodes):
    nodes.sort(key=lambda x: x.count, reverse=True)
    merged = set()
    for i, node in enumerate(nodes):
        if i in merged:
            continue
        nodes[i] = scan_extend(node, nodes, i, merged)

    nodes = [node for node in nodes if node is not None]

    return nodes


def merge_regions(regions, buffer=10):
    # Sort regions by start position
    regions = sorted(regions, key=lambda x: x[0])
    
    merged = []
    current_start, current_end, current_seq = regions[0]
    current_seqs = [current_seq]

    for start, end, contig in regions[1:]:
        # Check if the current region overlaps or is within 10bp of the next region
        if start <= current_end + buffer:
            # Extend the current region
            current_end = max(current_end, end)
            current_seqs.append(contig)
        else:
            # Append the current region and start a new one
            merged.append((current_start, current_end, current_seqs))
            current_start, current_end, current_seq = start, end, contig
            current_seqs = [current_seq]
    
    # Add the last region
    merged.append((current_start, current_end, current_seqs))
    
    return merged


def do_gene(gene, blosum_folder, trimmed_folder, aa_input, nt_input, aa_output, nt_output, debug):

    internal_seqs = list(parseFasta(Path(nt_input, gene)))
    blosum_seqs = list(parseFasta(Path(blosum_folder, gene)))
    trimmed_seqs = list(parseFasta(Path(trimmed_folder, "nt", gene)))
    trimmed_aa_seqs = list(parseFasta(Path(trimmed_folder, "aa", gene.replace(".nt.", ".aa."))))
    
    blosum_headers = {header for header, _ in blosum_seqs}
    trimmed_headers = {header for header, _ in trimmed_seqs}
    kicked = trimmed_headers - blosum_headers
    
    internal_headers = {header for header, _ in internal_seqs}
    kicked = kicked - internal_headers
    kicked_sequences = [(header, sequence) for header, sequence in trimmed_seqs if header in kicked]
    
    # Form contigs with internal seqs
    internal_nodes = []
    for header, sequence in internal_seqs:
        start, end = find_index_pair(sequence, "-")
        internal_nodes.append(
            NODE(
                header=header,
                count=int(header.split("|")[5]),
                sequence=None,
                nt_sequence=sequence,
                start=start,
                end=end,
                children=[],
                codename=None,
                is_contig=False,
            )
        )
        
    internal_contigs = [contig for contig in simple_assembly(internal_nodes) if contig.is_contig]
    
    out_nt = internal_seqs.copy()
    out_aa = list(parseFasta(Path(aa_input, gene.replace(".nt.",".aa."))))
    aa_seqs = dict(trimmed_aa_seqs)
    recovered = 0
    if internal_contigs:
        kicked_nodes = []
        for header, sequence in kicked_sequences:
            start, end = find_index_pair(sequence, "-")
            kicked_nodes.append(
                NODE(
                    header=header,
                    count=int(header.split("|")[5]),
                    sequence=aa_seqs[header],
                    nt_sequence=sequence,
                    start=start,
                    end=end,
                    children=[],
                    codename=None,
                    is_contig=False,
                )
            )
        
        add_occured = True
        saved_headers = set()
        while add_occured:
            add_occured = False
            contig_regions = [(contig.start, contig.end, contig) for contig in internal_contigs]
            flattened_regions = merge_regions(contig_regions)
            nodes_outside_of_region = []
            for rstart, rend, contigs in flattened_regions:
                for node in kicked_nodes:
                    if node.header in saved_headers:
                        continue
                    overlap_coords = get_overlap(rstart, rend, node.start, node.end, 1)
                    overlap_amount = 0 if not overlap_coords else overlap_coords[1] - overlap_coords[0]
                    # If overlaps but is not contained in region
                    if overlap_amount == min((node.end - node.start), (rend - rstart)):
                        continue
                    
                    nodes_outside_of_region.append(node)
                
                for contig in contigs:
                    merged_indices = set()
                    this_possible_merges = sorted(copy.deepcopy(nodes_outside_of_region), key=lambda x: x.start)
                    scan_extend(contig, this_possible_merges, None, merged_indices)
                    for i in merged_indices:
                        saved_headers.add(this_possible_merges[i].header)
                    
            already_added = set()
            for node in kicked_nodes:
                if not node.header in saved_headers:
                    continue
                if node.header in already_added:
                    continue
                already_added.add(node.header)
                recovered += 1
                print("Recovered",node.header)
                out_nt.append((node.header, node.nt_sequence))
                out_aa.append((node.header, node.sequence))
                
        out_aa.sort(key=lambda x: find_index_pair(x[1], "-")[0])
        aa_order = [header for header, _ in out_aa]
        nt_sequences = dict(out_nt)
        out_nt = [(header, nt_sequences[header]) for header in aa_order if not header.endswith('.')]
                
        writeFasta(Path(nt_output, gene), out_nt)
        writeFasta(Path(aa_output, gene.replace(".nt.",".aa.")), out_aa)

    return recovered

def main(args, sub_dir):
    timer = TimeKeeper(KeeperMode.DIRECT)
    
    folder = args.INPUT
    outlier_folder = Path(folder, "outlier")
    blosum_folder = Path(outlier_folder, "blosum", "nt")
    trimmed_folder = Path(folder, "trimmed")
    input_folder = Path(outlier_folder, sub_dir)
    if not input_folder.exists():
        input_folder = Path(folder, sub_dir)

    printv(f"Processing: {folder}", args.verbose)

    output_folder = Path(folder, "outlier", "recovered")

    aa_output = output_folder.joinpath("aa")
    nt_output = output_folder.joinpath("nt")

    if not path.exists(aa_output):
        makedirs(str(aa_output), exist_ok=True)
        makedirs(str(nt_output), exist_ok=True)

    aa_input = input_folder.joinpath("aa")
    nt_input = input_folder.joinpath("nt")

    compress = not args.uncompress_intermediates or args.compress

    genes = [fasta for fasta in listdir(nt_input) if ".fa" in fasta]

    arguments = [(gene, blosum_folder, trimmed_folder, aa_input, nt_input, aa_output, nt_output, args.debug) for gene in genes]
    if args.processes > 1:
        with Pool(args.processes) as pool:
            results = pool.starmap(do_gene, arguments)
    else:
        results = []
        for arg in arguments:
            results.append(
                do_gene(
                    *arg
                )
            )
 
    print("Added", sum(results), "sequences")
 
    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)

    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
