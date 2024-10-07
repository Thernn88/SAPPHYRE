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


def make_duped_consensus(
    raw_sequences: list, threshold: float
) -> str:
    bundled_seqs = [(seq, int(header.split("|")[5])) for header, seq in raw_sequences if header[-1] != "."]
    return dumb_consensus_dupe(bundled_seqs, threshold, 0)

def find_non_depth_regions(sequence, buffer=10):
    regions = []
    start = None  # To track the start of a region

    # First pass: Identify all regions of non-'X' and non-'?'
    for i, char in enumerate(sequence):
        if char != 'X' and char != '?':
            if start is None:
                start = i  # Start of a new region
        else:
            if start is not None:
                regions.append((start, i))  # End of a region
                start = None

    # If the last character is not 'X' or '?', we need to close the final region
    if start is not None:
        regions.append((start, len(sequence)))

    # Second pass: Merge regions within the buffer distance
    merged_regions = []
    current_start, current_end = regions[0]

    for start, end in regions[1:]:
        if start <= current_end + buffer:
            current_end = max(current_end, end)  # Extend the current region
        else:
            merged_regions.append((current_start, current_end))  # Finalize the current region
            current_start, current_end = start, end  # Start a new region

    # Add the last merged region
    merged_regions.append((current_start, current_end))

    return merged_regions


def do_gene(gene, blosum_folder, trimmed_folder, aa_input, nt_input, aa_output, nt_output, debug):
    # print(gene)
    no_dupes = False

    internal_seqs = list(parseFasta(Path(nt_input, gene)))
    out_aa = list(parseFasta(Path(aa_input, gene.replace(".nt.", ".aa."))))
    blosum_seqs = list(parseFasta(Path(blosum_folder, gene)))
    trimmed_seqs = list(parseFasta(Path(trimmed_folder, "nt", gene)))
    trimmed_aa_seqs = list(parseFasta(Path(trimmed_folder, "aa", gene.replace(".nt.", ".aa."))))
    
    blosum_headers = {header for header, _ in blosum_seqs}
    trimmed_headers = {header for header, _ in trimmed_seqs}
    kicked = trimmed_headers - blosum_headers
    
    internal_headers = {header for header, _ in internal_seqs}
    kicked = kicked - internal_headers
    kicked_sequences = [(header, sequence) for header, sequence in trimmed_seqs if header in kicked]
    recovered = 0
    saved_headers = set()
    
    # added_sequences = set()
    new_seq_added = True
    while new_seq_added:
        new_seq_added = False
        internal_sequences = [x[1] for x in internal_seqs]
        if no_dupes:
            consensus_seq = dumb_consensus(internal_sequences, 0.5, 0)
        else:
            consensus_seq = make_duped_consensus(
                internal_seqs, 0.5
            )
            
        kicked_sequences = [(header, sequence) for header, sequence in kicked_sequences if header not in saved_headers]
        if not kicked_sequences:
            break
        
        consensus_seq = convert_consensus(internal_sequences, consensus_seq)
        regions = find_non_depth_regions(consensus_seq)
        to_check = []
        for header, sequence in kicked_sequences:
            start, end = find_index_pair(sequence, "-")
            
            doesnt_overlap = True
            best_small_overlap = 12
            bso_coords = None
            for rstart, rend in regions:
                overlap_coords = get_overlap(start, end, rstart, rend, 1)
                if overlap_coords:
                    amount = overlap_coords[1] - overlap_coords[0]
                    
                    if amount >= min((end - start), (rend - rstart)):
                        doesnt_overlap = False
                        break
                    
                    if amount >= best_small_overlap: # Greater than 12
                        best_small_overlap = amount
                        bso_coords = overlap_coords
                        region = (rstart, rend)
            
            if doesnt_overlap and bso_coords:
                to_check.append((header, sequence, bso_coords[0], bso_coords[1], region[0], region[1], start, end))
        
        if to_check:   
            for header, sequence, start, end, rstart, rend, seq_start, seq_end in to_check:
                
                kick_kmer = sequence[start:end]
                consensus_kmer = consensus_seq[start:end]
                if is_same_kmer(kick_kmer, consensus_kmer):
                    recovered += 1
                    print("Recovered", header)
                    saved_headers.add(header)
                    internal_seqs.append((header, sequence))
                    new_seq_added = True

    
    out_nt = internal_seqs.copy()
    if saved_headers:
        for header, sequence in trimmed_aa_seqs:
            if header in saved_headers:
                out_aa.append((header, sequence))
        
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
