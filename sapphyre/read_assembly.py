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
    sequence: str
    count: int
    nt_sequence: str
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
        before = len(self.sequence) - self.sequence.count("-")
        # If node_2 is contained inside self
        if node_2.start >= self.start and node_2.end <= self.end:
            self.sequence = (
                self.sequence[:overlap_coord // 3]
                + node_2.sequence[overlap_coord // 3 : node_2.end]
                + self.sequence[node_2.end :]
            )

            self.nt_sequence = (
                self.nt_sequence[:overlap_coord]
                + node_2.nt_sequence[overlap_coord : node_2.end * 3]
                + self.nt_sequence[node_2.end * 3 :]
            )

        # If node_2 contains self
        elif self.start >= node_2.start and self.end <= node_2.end:
            self.sequence = (
                node_2.sequence[:overlap_coord // 3]
                + self.sequence[overlap_coord // 3 : self.end]
                + node_2.sequence[self.end :]
            )

            self.nt_sequence = (
                node_2.nt_sequence[:overlap_coord]
                + self.nt_sequence[overlap_coord : self.end * 3]
                + node_2.nt_sequence[self.end * 3 :]
            )

            self.start = node_2.start
            self.end = node_2.end
        # If node_2 is to the right of self
        elif node_2.start >= self.start:
            self.sequence = (
                self.sequence[:overlap_coord // 3] + node_2.sequence[overlap_coord // 3:]
            )

            self.nt_sequence = (
                self.nt_sequence[:overlap_coord] + node_2.nt_sequence[overlap_coord:]
            )

            self.end = node_2.end
        # If node_2 is to the left of self
        else:
            self.sequence = (
                node_2.sequence[:overlap_coord // 3] + self.sequence[overlap_coord // 3:]
            )

            self.nt_sequence = (
                node_2.nt_sequence[:overlap_coord] + self.nt_sequence[overlap_coord:]
            )

            self.start = node_2.start

        # Save node_2 and the children of node_2 to self
        self.children.append(node_2.header)
        self.children.extend(node_2.children)
        if (len(self.sequence) - self.sequence.count("-")) > before:
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


def check_covered_bad_regions(nodes, consensus, min_ambiguous, max_distance, ambig_char):
    x_indices = []
    current_group = []
    regions = []

    start, stop = find_index_pair(consensus, "X")

    for i, base in enumerate(consensus[start:stop], start):
        if base == ambig_char:
            x_indices.append(i)

    for num in x_indices:
        if not current_group:
            current_group.append(num)
        elif num - current_group[-1] <= max_distance:
            current_group.append(num)
        else:
            if len(current_group) >= min_ambiguous:
                regions.append((current_group[0], current_group[-1] + 1))
            current_group = [num]

    if current_group:
        if len(current_group) >= min_ambiguous:
            regions.append((current_group[0], current_group[-1] + 1))

    # MERGE
    merge_occured = True
    while merge_occured:
        merge_occured = False
        for (i, region_a), (j, region_b) in combinations(enumerate(regions), 2):
            #gap_between = region_a[1], region_b[0]
            aa_gap_region = floor(region_a[1] / 3), ceil(region_b[0] / 3)
            seqs_in_gap = [(node.start, node.end) for node in nodes if get_overlap(node.start, node.end, aa_gap_region[0], aa_gap_region[1], 1)]
            if seqs_in_gap:
                all_overlap = True
                cr_start, cr_end = seqs_in_gap[0]
                for start, end in seqs_in_gap[1:]:
                    if get_overlap(cr_start, cr_end, start, end, 1) is None:       
                        all_overlap = False
                        break
                    else:
                        cr_start, cr_end = min(cr_start, start), max(cr_end, end)
                    
                if all_overlap:
                    regions[i] = (region_a[0], region_b[1])
                    regions[j] = None
                    merge_occured = True

        regions = [i for i in regions if i]
        
    if regions:
        return regions, consensus

    return [(None, None)], consensus

def scan_extend(node, nodes, i, merged, min_overlap_percent, min_overlap_chars):
    merge_occured = True
    while merge_occured:
        merge_occured = False
        
        for j, node_b in enumerate(nodes):  # When merge occurs start again at the beginning with highest count
            if j in merged:
                continue
            if i == j:
                continue
            overlap_coords = get_overlap(node.start * 3, node.end * 3, node_b.start * 3, node_b.end * 3, 1)

            
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]

                # Calculate percent overlap and compare to minimum overlap
                # overlap_percent = overlap_amount / ((node.end - node.start) * 3)
                required_overlap = max(min_overlap_chars, min((node.end - node.start), (node_b.end - node_b.start)) * 3 * min_overlap_percent)
                
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
                
                nodes[j] = None
                break
    return node

def simple_assembly(nodes, min_overlap_percent, min_overlap_chars):
    nodes.sort(key=lambda x: x.count, reverse=True)
    merged = set()
    for i, node in enumerate(nodes):
        if i in merged:
            continue
        nodes[i] = scan_extend(node, nodes, i, merged, min_overlap_percent, min_overlap_chars)

    nodes = [node for node in nodes if node is not None]

    return nodes

def contigs_that_resolve(possible_contigs, nodes_out_of_region, min_overlap = 0.25):
    contigs = []
    for contig in possible_contigs:
         # Check if contig is longer than 250bp and immediately pass it if so
        if (contig.end - contig.start) * 3 > 250:
            contigs.append(contig)
            continue
        for node in nodes_out_of_region:
            overlap_coords = get_overlap(contig.start * 3, contig.end * 3, node.start * 3, node.end * 3, 1)
            if overlap_coords is None:
                continue
            
            percent = (overlap_coords[1] - overlap_coords[0]) / min(contig.end - contig.start, node.end - node.start)
            if percent < min_overlap:
                continue

            # Overlap with kmer match
            # kmer_node = node.nt_sequence[overlap_coords[0]:overlap_coords[1]]
            # kmer_contig = contig.nt_sequence[overlap_coords[0]:overlap_coords[1]]

            # if is_same_kmer(kmer_node, kmer_contig):
            contigs.append(contig)
            break
    
    return contigs


def get_regions(seqeunces, nodes, threshold, no_dupes, minimum_ambig, return_regions, minimum_distance=30, ambig_char='X'):
    sequences = [x[1] for x in seqeunces]
    if no_dupes:
        consensus_seq = dumb_consensus(sequences, threshold, 0)
    else:
        consensus_seq = make_duped_consensus(
            seqeunces, threshold
        )

    consensus_seq = convert_consensus(sequences, consensus_seq)
    if return_regions:
        return check_covered_bad_regions(nodes, consensus_seq, minimum_ambig, minimum_distance, ambig_char)
    else:
        rstart, rend = find_index_pair(consensus_seq, "X")
        return consensus_seq[rstart:rend].find("X") != -1

def apply_positions(aa_nodes, x_positions, kicked_headers, log_output, debug, position_subset):
    for node in aa_nodes:
        positions = position_subset[node.header]
        if positions:
            node.sequence = del_cols(node.sequence, positions)
            node.start, node.end = find_index_pair(node.sequence, "-")
            if len(node.sequence) - node.sequence.count("-") < 15:
                kicked_headers.add(node.header)
                if debug:
                    log_output.append(f"Kicked >{node.header} due to low length after trimming (<15 AA)\n")
                    if debug > 1:
                        log_output.append(f"{node.sequence}\n")
                continue
            x_positions[node.header].update(positions)


def del_cols(sequence, columns, nt=False):
    if nt:
        return join_triplets_with_exclusions(sequence, set(), columns)

    return join_with_exclusions(sequence, columns)


def do_trim(aa_nodes, x_positions, ref_consensus, kicked_headers, no_dupes, excise_trim_consensus, log_output, debug):
    aa_sequences = [x.sequence for x in aa_nodes]
    if no_dupes:
        consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)
    else:
        current_raw_aa = [(node.header, node.sequence) for node in aa_nodes]
        consensus_seq = make_duped_consensus(
            current_raw_aa, excise_trim_consensus
        )

    cstart, cend = find_index_pair(consensus_seq, "X")
    first_positions = defaultdict(set)
    second_positions = defaultdict(set)

    for i, maj_bp in enumerate(consensus_seq[cstart:cend], cstart):
        if maj_bp != "X":
            continue

        in_region = []
        out_of_region = False
        for x, node in enumerate(aa_nodes):
            # within 3 bp
            if node.start <= i <= node.start + 3:
                in_region.append((x, node.sequence[i], False))
            elif node.end - 3 <= i <= node.end:
                in_region.append((x, node.sequence[i], True))
            elif i >= node.start and i <= node.end:
                out_of_region = True
                
        if not out_of_region and not in_region:
            continue

        if not out_of_region and in_region:
            for node_index, bp, on_end in in_region:
                node = aa_nodes[node_index]
                if bp in ref_consensus[i]:
                    continue
                
                if on_end:
                    for x in range(i, node.end):
                        first_positions[node.header].add(x * 3)
                else:
                    for x in range(node.start, i + 1):
                        first_positions[node.header].add(x * 3)

        if out_of_region and in_region:
            for node_index, bp, on_end in in_region:
                node = aa_nodes[node_index]
                if on_end:
                    for x in range(i, node.end):
                        first_positions[node.header].add(x * 3)
                else:
                    for x in range(node.start, i + 1):
                        first_positions[node.header].add(x * 3)


    #refresh aa
    apply_positions(aa_nodes, x_positions, kicked_headers, log_output, debug, first_positions)

    nodes = [node for node in aa_nodes if node.header not in kicked_headers]
    if not nodes:
        if debug:
            log_output.append("No nodes left after trimming\n")
        return

    if no_dupes:
        aa_sequences = [x.sequence for x in nodes]
        consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)
    else:
        current_raw_aa = [(node.header, node.sequence) for node in nodes]
        consensus_seq = make_duped_consensus(
            current_raw_aa, excise_trim_consensus
        )


    for node in aa_nodes:
        i = None
        for poss_i in range(node.start, node.start + 3):
            if node.sequence[poss_i] != consensus_seq[poss_i]:
                i = poss_i

        if not i is None:
            for x in range(node.start , i + 1):
                second_positions[node.header].add(x * 3)

        i = None
        for poss_i in range(node.end -1, node.end - 4, -1):
            if node.sequence[poss_i] != consensus_seq[poss_i]:
                i = poss_i

        if not i is None:
            for x in range(i, node.end):
                second_positions[node.header].add(x * 3)
    
    #refresh aa
    apply_positions(aa_nodes, x_positions, kicked_headers, log_output, debug, second_positions)


def get_score_difference(score_a: float, score_b: float) -> float:
    """Get the decimal difference between two scores.

    Args:
    ----
        score_a (float): The first score.
        score_b (float): The second score.

    Returns:
    -------
        float: The decimal difference between the two scores.
    """
    # If either score is zero return zero.
    if score_a == 0.0 or score_b == 0.0:
        return 0.0

    # Return the decimal difference between the largest score and the smallest score.
    return max(score_a, score_b) / min(score_a, score_b)


def get_score(contigs, flex_consensus, nodes, max_score, contig_depth_reward, log_output, debug):
    with_identity = []
    for node in contigs:
        matches = 0
        for i in range(node.start, node.end):
            matches += min(flex_consensus[i].count(node.sequence[i]), max_score)
        
        score = matches * (1 + (len(node.children) * contig_depth_reward))
        
        length = (node.end - node.start)
        children = node.get_children()
        if debug:
            log_output.append(f"{node.codename} with ({len(children)}: {', '.join(children)})\nhas a score of {score} over {length} AA\n{node.nt_sequence}\n")

        consists_of = []
        if debug > 1:
            for header in node.children:
                for node in nodes:
                    if node.header == header:
                        consists_of.append(node)
                        break
        
        consists_of.sort(key=lambda x: x.start)
        for node in consists_of:
            log_output.append(f">{node.header}\n{node.nt_sequence}")

        with_identity.append((node, score, length))

    return with_identity

def do_gene(gene, aa_input, nt_input, aa_output, nt_output, no_dupes, compress, excise_consensus, allowed_mismatches, region_min_ambig, debug):
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    kicks = 0
    
    max_score = 8 # Maximum count for matching columns in motif scoring
    min_children = 4 # Minimum children a contig must have to be considered
    contig_depth_reward = 0.001 # Reward for each child a contig has
    tie_breaker_percent = 1.15 # If top contig scores are within 15% trigger tie breaker which merges non ambig reads into the contigs to see if they can continue merging
    min_contig_bp = 100 # Minimum contig bp
    contig_read_overlap = 0.1 # Min overlap between reads and contigs for kick
    min_overlap_percent = 0.15 # Min overlap percent for reads in simple assembly
    min_overlap_chars = 10 # Min overlap chars for reads in simple assembly
    min_ambig_overlap_chars = 10 # Min overlap chars for reads into regions
    region_overlap = 0.15 # Min overlap percent for reads into regions

    kicked_nodes = set()
    unresolved = []
    raw_nodes = []
    log_output = []
    cull_positions = defaultdict(set)
    aa_gene = path.join(aa_input, gene)
    nt_gene = gene.replace(".aa.", ".nt.")
    raw_nodes = list(parseFasta(path.join(nt_input, nt_gene)))

    # Check bad region
    
    gene_has_mismatch = get_regions(raw_nodes, None, excise_consensus, no_dupes, region_min_ambig, False)
    if not gene_has_mismatch:
        writeFasta(path.join(aa_output, gene), parseFasta(aa_gene), compress)
        writeFasta(path.join(nt_output, nt_gene), raw_nodes, compress)
        return log_output, False, kicks, unresolved
    # Assembly
    if debug:
        log_output.append(f"Log output for {gene}\n")

    nodes = {header:
        NODE(header, "", int(header.split("|")[5]), sequence, None, None, [], None, False) for header, sequence in raw_nodes
    }
    
    
    flex_consensus = defaultdict(list)
    raw_aa = list(parseFasta(aa_gene))
    for header, sequence in raw_aa:
        if header.endswith("."):
            start, end = find_index_pair(sequence, "-")
            for i, bp in enumerate(sequence[start:end]):
                flex_consensus[i].append(bp)
            continue
        
        if header in nodes:
            parent = nodes[header]
            parent.sequence = sequence
            parent.start, parent.end = find_index_pair(sequence, "-")

    nodes = list(nodes.values())
    
    # do_trim(nodes, cull_positions, flex_consensus, kicked_nodes, no_dupes, excise_consensus, log_output, debug)
    
    nodes = [node for node in nodes if node.header not in kicked_nodes]
    if not nodes:
        return log_output, False, len(raw_nodes), unresolved
    for node in nodes:
        node.nt_sequence = del_cols(node.nt_sequence, cull_positions[node.header], True)

    recursion_limit = 5
    regions, consensus_seq = get_regions([(i.header, i.nt_sequence) for i in nodes], nodes, excise_consensus, no_dupes, region_min_ambig, True)
    changes_made = True
    had_region = False
    while regions[0][0] is not None and changes_made and recursion_limit >= 0:
        had_region = True
        recursion_limit -= 1
        changes_made = False
        for start, end in regions:
            nodes = [node for node in nodes if node.header not in kicked_nodes]
            if debug:
                log_output.append(f"Found region between {start} - {end}")
                log_output.append(f"{consensus_seq}\n")
            nodes_in_region = []
            nodes_out_of_region = []
            for node in nodes:
                coords = get_overlap(node.start * 3, node.end * 3, start, end, 1)
                overlap_amount = 0 if coords is None else coords[1] - coords[0]
                percent = 0 if coords is None else overlap_amount / (node.end - node.start)
                if percent < region_overlap and overlap_amount < min_ambig_overlap_chars:
                    nodes_out_of_region.append(node)
                    continue

                nodes_in_region.append(node)

            merged_nodes = simple_assembly(copy.deepcopy(nodes_in_region), min_overlap_percent, min_overlap_chars)
            
            possible_contigs = [node for node in merged_nodes if (len(node.children) >= min_children and 
                                                                  node.is_contig and # BP increase due to merge
                                                                  len(node.nt_sequence) - node.nt_sequence.count("-") >= min_contig_bp)]
            contigs = contigs_that_resolve(possible_contigs, nodes_out_of_region)
            for i, node in enumerate(contigs):
                node.codename = f"Contig{i}"

            if not contigs:
                changes_made = True
                for read in nodes_in_region:
                    kicked_nodes.add(read.header)
                    if debug:
                        log_output.append(f"Kicked >{read.header} due to no contig resolution in gene\n")
                        if debug > 1:
                            log_output.append(f"{read.nt_sequence}\n")
            else:
                with_identity = get_score(contigs, flex_consensus, nodes, max_score, contig_depth_reward, log_output, debug)
                    
                with_identity.sort(key=lambda x: x[1], reverse=True)
                best_contig, best_score, _ = with_identity[0]
                
                tie_breaker = []
                for contig, score, _ in with_identity[1:]:
                    if get_score_difference(best_score, score) < tie_breaker_percent:
                        tie_breaker.append(contig)
                        
                if tie_breaker and True: # Set True to False to disable tiebreaker
                    tie_breaker.insert(0, best_contig)
                    if debug:
                        log_output.append(f"\nBest contig: {best_contig.codename}\n\nComparing against reads")
                        log_output.append(f"Multiple contigs with similar scores, choosing the best contig by merging into unambig\n")
                    for i, contig in enumerate(tie_breaker):
                        tie_breaker[i] = scan_extend(contig, sorted(copy.deepcopy(nodes_out_of_region), key=lambda x: x.start), None, set(), min_overlap_percent, min_overlap_chars)
                    with_identity = get_score(tie_breaker, flex_consensus, nodes, max_score, contig_depth_reward, log_output, debug)
                    best_contig = max(with_identity, key=lambda x: x[1])[0]
                    
                if debug:
                    log_output.append(f"\nBest contig: {best_contig.codename}\n\nComparing against reads")

                for i, read in enumerate(nodes_in_region):
                    overlap_coords = get_overlap(best_contig.start * 3, best_contig.end * 3, read.start * 3, read.end * 3, 1)
                    if overlap_coords is None:
                        continue
                    
                    overlap_percent = (overlap_coords[1] - overlap_coords[0]) / min(best_contig.end - best_contig.start, read.end - read.start)
                    if overlap_percent < contig_read_overlap:
                        continue
                    
                    kmer_node = read.nt_sequence[overlap_coords[0]:overlap_coords[1]]
                    kmer_best = best_contig.nt_sequence[overlap_coords[0]:overlap_coords[1]]
                    distance = constrained_distance(kmer_node, kmer_best)

                    if distance <= allowed_mismatches:
                        if debug > 1:
                            log_output.append(f"Kept >{read.header} due to matching the best contig ({distance} mismatches)\n{read.nt_sequence}\n")
                        continue
                    
                    changes_made = True
                    kicked_nodes.add(read.header)
                    if debug:
                        log_output.append(f"Kicked >{read.header} due to distance from best contig ({distance} mismatches)\n")
                        if debug > 1:
                            log_output.append(f"{read.nt_sequence}\n")

        nodes = [node for node in nodes if node.header not in kicked_nodes]
        if not nodes:
            regions = [(None, None)]
            break
        regions, consensus_seq = get_regions([(i.header, i.nt_sequence) for i in nodes], nodes, excise_consensus, no_dupes, region_min_ambig, True)

    if regions[0][0] is not None:
        unresolved.append(f"Unresolved region in {gene} - {regions}")

    nt_out = []
    for header, seq in raw_nodes:
        if header in kicked_nodes:
            continue
        else:
            nt_out.append((header, seq))

    if nt_out:
        aa_out = []
        for header, seq in raw_aa:
            if header in kicked_nodes:
                kicks += 1
                continue
            else:
                aa_out.append((header, seq))
    
        writeFasta(path.join(aa_output, gene), aa_out, compress)
        writeFasta(path.join(nt_output, nt_gene), nt_out, compress)
            
    return log_output, had_region, kicks, unresolved


def main(args, sub_dir):
    timer = TimeKeeper(KeeperMode.DIRECT)
    
    folder = args.INPUT
    input_folder = Path(folder, "outlier", sub_dir)
    if not input_folder.exists():
        input_folder = Path(folder, sub_dir)

    printv(f"Processing: {folder}", args.verbose)

    output_folder = Path(folder, "outlier", "excise")

    aa_output = output_folder.joinpath("aa")
    nt_output = output_folder.joinpath("nt")

    if not path.exists(aa_output):
        makedirs(str(aa_output), exist_ok=True)
        makedirs(str(nt_output), exist_ok=True)

    aa_input = input_folder.joinpath("aa")
    nt_input = input_folder.joinpath("nt")

    genes = [fasta for fasta in listdir(aa_input) if ".fa" in fasta]

    arguments = [(gene, aa_input, nt_input, aa_output, nt_output, args.no_dupes, args.compress, args.excise_consensus, args.excise_allowed_distance, args.excise_minimum_ambig, args.debug) for gene in genes]
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
 
    log_final = []
    ambig_count = 0
    total_kicks = 0
    total_unresolved = []
    for log, has_ambig, kicks, unresolved in results:
        total_unresolved.extend(unresolved)
        ambig_count += 1 if has_ambig else 0
        total_kicks += kicks
        if len(log) > 1:
            log_final.extend(log)
    
    printv(f"{folder}: {ambig_count} ambiguous loci found. Kicked {total_kicks} sequences total. {len(total_unresolved)} unresolved ambiguous", args.verbose)

    if args.debug:
        with open(output_folder.joinpath("excise.log"), "w") as log_file:
            for line in log_final:
                log_file.write(line+"\n")
        with open(output_folder.joinpath("unresolved.log"), "w") as log_file:
            for line in total_unresolved:
                log_file.write(line+"\n")

    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)

    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
