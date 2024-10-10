from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from os import listdir, makedirs, path
from pathlib import Path
from msgspec import Struct
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
)
from .directional_cluster import cluster_ids, within_distance, node_to_ids, quick_rec
from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

def check_covered_bad_regions(consensus, min_ambiguous, ambig_char='X', max_distance=18):
    x_indices = []
    current_group = []

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
                return (current_group[0], current_group[-1] + 1)
            current_group = [num]

    if current_group:
        if len(current_group) >= min_ambiguous:
            return (current_group[0], current_group[-1] + 1)


    return None, None


def make_duped_consensus(
    raw_sequences: list, threshold: float
) -> str:
    bundled_seqs = [(seq, int(header.split("|")[5])) for header, seq in raw_sequences if header[-1] != "."]
    return dumb_consensus_dupe(bundled_seqs, threshold, 0)


class NODE(Struct):
    header: str
    frame: int
    sequence: str
    nt_sequence: str
    start: int
    end: int
    children: list

def del_cols(sequence, columns, nt=False):
    if nt:
        seq = [sequence[i: i+3] for i in range(0, len(sequence), 3)]
        for i in columns:
            seq[i] = "---"
        return "".join(seq)
    seq = list(sequence)
    for i in columns:
        seq[i] = "-"
    return "".join(seq)


def get_coverage(aa_seqs, ref_avg_len):
    data_cols = 0
    for i in range(len(aa_seqs[0])):
        for seq in aa_seqs:
            if seq[i] != "-":
                data_cols += 1
                break

    return data_cols / ref_avg_len


def cluster(cluster_raw, cluster_threshold):
    clusters = []

    current_cluster = []
    for child_index, child_full in cluster_raw:
        if not current_cluster:
            current_cluster.append(child_full)
            current_index = child_index
        else:
            if child_index - current_index <= cluster_threshold:
                current_cluster.append(child_full)
                current_index = child_index
            else:
                clusters.append(current_cluster)
                current_cluster = [child_full]
                current_index = child_index
    
    if current_cluster:
        clusters.append(current_cluster)

    return clusters


def calculate_split(node_a: str, node_b: str, overlapping_coords: tuple, ref_consensus: dict) -> int:
    """Iterates over each position in the overlap range of sequence A and sequence B and
    creates a frankenstein sequence of sequence A + Sequence B joined at each
    position in the overlap.

    Final split position = pos in overlap with the highest score.

    Score is determined by the amount of characters that are the same between each
    position in the frankenstein sequence and the comparison sequence.
    """
    overlap_start, overlap_end = overlapping_coords

    highest_score = -1
    highest_scoring_pos = overlap_start

    for i in range(overlap_start, overlap_end + 1):
        prev_kmer = node_a.sequence[overlap_start:i]
        next_kmer = node_b.sequence[i:overlap_end]

        window = prev_kmer + next_kmer

        this_window_score = 0
        for x, let in enumerate(window, overlap_start):
            if let in ref_consensus[x]:
                this_window_score += ref_consensus[x].count(let)
            else:
                this_window_score -= 1
            if let == "-":
                this_window_score -= 1
                
        if this_window_score >= highest_score:
            highest_score = this_window_score
            highest_scoring_pos = i

    return highest_scoring_pos


def do_trim(aa_nodes, cluster_sets, x_positions, ref_consensus, kicked_headers, excise_trim_consensus, no_dupes):
    for cluster_set in cluster_sets:
        sub_aa_nodes = [node for node in aa_nodes if node.header not in kicked_headers and within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0)]

        aa_sequences = [node.sequence for node in sub_aa_nodes]
        if aa_sequences:
            if not no_dupes:
                current_raw_aa = [(node.header, node.sequence) for node in sub_aa_nodes]
                consensus_seq = make_duped_consensus(
                    current_raw_aa, excise_trim_consensus
                )
            else:
                consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)

            cstart, cend = find_index_pair(consensus_seq, "X")
            cstart, cend = find_index_pair(consensus_seq, "X")

            for i, maj_bp in enumerate(consensus_seq[cstart:cend], cstart):
                if maj_bp != "X":
                    continue
            for i, maj_bp in enumerate(consensus_seq[cstart:cend], cstart):
                if maj_bp != "X":
                    continue

                in_region = []
                out_of_region = []
                for x, node in enumerate(sub_aa_nodes):
                    # within 3 bp
                    within_left = node.start <= i <= node.start + 3
                    within_right = node.end - 3 <= i <= node.end

                    if within_left or within_right:
                        in_region.append((x, node.sequence[i], within_right))
                    elif i >= node.start and i <= node.end:
                        out_of_region.append((x, node.sequence[i], within_right))
                    if within_left or within_right:
                        in_region.append((x, node.sequence[i], within_right))
                    elif i >= node.start and i <= node.end:
                        out_of_region.append((x, node.sequence[i], within_right))

                if not out_of_region and not in_region:
                    continue
                if not out_of_region and not in_region:
                    continue

                if not out_of_region and in_region:
                    for node_index, bp, on_end in in_region:
                        node = sub_aa_nodes[node_index]
                        if bp in ref_consensus[i]:
                            continue
                        
                        if on_end:
                            for x in range(i, node.end):
                                x_positions[node.header].add(x)
                        else:
                            for x in range(node.start, i + 1):
                                x_positions[node.header].add(x)

                if out_of_region and in_region:
                    for node_index, bp, on_end in in_region:
                        node = sub_aa_nodes[node_index]
                        if on_end:
                            for x in range(i, node.end):
                                x_positions[node.header].add(x)
                        else:
                            for x in range(node.start, i + 1):
                                x_positions[node.header].add(x)


            #refresh aa
            if x_positions:
                for node in sub_aa_nodes:
                    node.sequence = del_cols(node.sequence, x_positions[node.header])
                    node.start, node.end = find_index_pair(node.sequence, "-")

            if not no_dupes:
                current_raw_aa = [(node.header, node.sequence) for node in sub_aa_nodes if node.header not in kicked_headers]
                consensus_seq = make_duped_consensus(
                    current_raw_aa, excise_trim_consensus
                )
            else:
                aa_sequences = [x.sequence for x in sub_aa_nodes if x.header not in kicked_headers]
                consensus_seq = dumb_consensus(aa_sequences, excise_trim_consensus, 0)

            for node in sub_aa_nodes:
                i = None
                for poss_i in range(node.start, node.start + 3):
                    if node.sequence[poss_i] != consensus_seq[poss_i]:
                        i = poss_i

                if not i is None:
                    for x in range(node.start , i + 1):
                        x_positions[node.header].add(x)

                i = None
                for poss_i in range(node.end -1, node.end - 4, -1):
                    if node.sequence[poss_i] != consensus_seq[poss_i]:
                        i = poss_i

                if not i is None:
                    for x in range(i, node.end):
                        x_positions[node.header].add(x)


def finalize_cluster(current_cluster, ids, ref_coords, clusters, kicks, req_seq_coverage):
    if len(current_cluster) >= 2:
        cluster_data_cols = set()
        for _, _, index in current_cluster:
            cluster_data_cols.update(ids[index][1])
            
        cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)
        clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage))
    elif len(current_cluster) == 1:
        if current_cluster[0][1] > req_seq_coverage:
            cluster_coverage = current_cluster[0][1]
            clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage))
        else:
            kicks.add(current_cluster[0][0])

def determine_direction(start, end, current_start, current_end, current_direction):
    this_direction = None
    if start == current_start and end == current_end:
        this_direction = "bi"
    else:
        if start == current_start:
            if end >= current_end:
                this_direction = "forward"
            else:
                this_direction = "reverse"
        else:
            if start >= current_start:
                this_direction = "forward"
            else:
                this_direction = "reverse"
    
    if current_direction == "bi" or this_direction == "bi" or this_direction == current_direction:
        return this_direction
    return None


def do_cluster(ids, ref_coords, max_gap_size, id_chomp_distance=100):
    clusters = []
    kicks = set()
    ids.sort(key=lambda x: x[0])
    grouped_ids = defaultdict(list)
    for i, (child_index, seq_coords, start, end) in enumerate(ids):
        id = int(child_index.split("_")[0])
        grouped_ids[id].append((i, child_index, seq_coords, start, end))
        
    ids_ascending = sorted(grouped_ids.keys())
    req_seq_coverage = 0.5
    current_cluster = []

    for id in ids_ascending:
        seq_list = grouped_ids[id]
        if not current_cluster:
            current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list]
            current_index = id
            current_seqs = seq_list
            current_direction = "bi"
        else:
            passed = False
            passed_direction = None

            if id - current_index <= id_chomp_distance:
                for i, child_index, seq_coords, start, end in seq_list:
                    for _, _, _, current_start, current_end in current_seqs:
                        overlap = get_overlap(start, end, current_start, current_end, -max_gap_size)
                        if not overlap:
                            continue
                        
                        passed_direction = determine_direction(start, end, current_start, current_end, current_direction)
                        if passed_direction:
                            passed = True
                            break
                    if passed:
                        break
            
            if passed:
                current_cluster.extend([(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list])
                current_index = id
                current_seqs = seq_list
                if passed_direction != "bi":
                    current_direction = passed_direction
            else:
                finalize_cluster(current_cluster, ids, ref_coords, clusters, kicks, req_seq_coverage)
                current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list]
                current_index = id
                current_seqs = seq_list
                current_direction = "bi"
    
    if current_cluster:
        finalize_cluster(current_cluster, ids, ref_coords, clusters, kicks, req_seq_coverage)
                
    return clusters, kicks



def log_excised_consensus(
    gene: str,
    input_path: Path,
    output_path: Path,
    compress_intermediates: bool,
    excise_consensus,
    excise_minimum_ambig,
    excise_rescue_match,
    excise_trim_consensus,
    no_dupes,
    true_cluster_threshold = 24,
    min_overlap_percent = 0.15,
    min_overlap_amount = 10,
):
    """
    By default, this does non-dupe consensus. If you pass "dupes=True", it will call
    the dupe-weighted version of consensus. If dupes is enable, the program needs dupe counts
    in the standard rocksDB location.

    Consensus_threshold is a value between 0.0 and 1.0 which represents a ratio.
    At each location, a bp is chosen as the consensus bp if the ratio of occurrences/total_characters_at_location
    exceeds the consensus_threshold. Raising this value makes the selection more strict. If you want this to match
    the default in outlier.py, set this to 0.65

    After the consensus sequence is made, the tail is checked for excessive X placeholder characters.
    Excise_threshold is a value between 0.0 and 1.0, which represents the maximum allowable ratio of X/total_characters.
    The tail is checked in blocks of 16 characters. The block checks are cumulative, so the percentage at each block is
    affected directly by all the preceding blocks. For example, if there are 8 X in the first block, and 10 X in the
    second block, the percentage at the second block is calculated as (10+8)/(16+16) = 18/32 = 0.5625 .
    Lowering this value makes the check more strict. Small changes in the excise_threshold can result in
    disproportionately large increases in truncation. My best results happen around 0.40 so far.

    The first field in the log file is the gene name.
    The second field in the log file is the cut-index in regular slice notation. It starts at the index of the
    first removed character, inclusive. The end is the length of the string, which is exclusive.
    """
    log_output = []
    this_rescues = []

    aa_in = input_path.joinpath("aa", gene)
    aa_out = output_path.joinpath("aa", gene)

    nt_in = input_path.joinpath("nt", gene.replace(".aa.", ".nt."))
    nt_out = output_path.joinpath("nt", gene.replace(".aa.", ".nt."))

    x_positions = defaultdict(set)

    bp_count = lambda x: len(x) - x.count("-")

    raw_aa = list(parseFasta(str(aa_in)))

    ref_lens = []
    aa_nodes = []
    reference_cluster_data = set()
    ref_consensus = defaultdict(list)
    for header, seq in raw_aa:
        if header.endswith('.'):
            start, end = find_index_pair(seq, "-")
            for i, bp in enumerate(seq[start:end], start):
                if bp != "-":
                    reference_cluster_data.add(i)
                ref_consensus[i].append(bp)
                
            ref_lens.append(bp_count(seq))
            continue

        aa_nodes.append(NODE(header, int(header.split("|")[4]), seq, None, *find_index_pair(seq, "-"), []))

    ref_avg_len = sum(ref_lens) / len(ref_lens)
    kicked_headers = set()
    
    ids = []
    for node in aa_nodes:
        start, end = find_index_pair(node.sequence, "-")
        ids.append(quick_rec(node.header.split("|")[3], node.frame, node.sequence, start, end))

    max_gap_size = round(len(aa_nodes[0].sequence) * 0.3) # Half MSA length

    clusters, _ = cluster_ids(ids, 100, max_gap_size, reference_cluster_data) #TODO: Make distance an arg
    cluster_sets = [set(range(a, b+1)) for a, b, _ in clusters]

    head_to_node = {}
    for node in aa_nodes:
        head_to_node[node.header] = node

    raw_sequences = {header: del_cols(seq, x_positions[header], True) for header, seq in parseFasta(str(nt_in))}
    for node in aa_nodes:
        nt_seq = raw_sequences[node.header]
        node.nt_sequence = nt_seq
    
    true_cluster_raw = []

    for header in raw_sequences:
        true_cluster_raw.append((int(header.split("|")[3].split("&&")[0].split("_")[1]), header))

    true_cluster_raw.sort(key = lambda x: x[0])
    before_true_clusters = cluster(true_cluster_raw, true_cluster_threshold)


    # Search for slices of the consensus seq with a high ratio of 'X' to total characters
    has_region = True
    had_region = False
    recursion_max = 5 * len(cluster_sets)
    last_region = {}
    consensus_seq = ""
    while has_region:
        has_region = False

        for cluster_i, cluster_set in enumerate(cluster_sets):

            aa_subset = [node for node in aa_nodes if node.header not in kicked_headers and within_distance(node_to_ids(node.header.split("|")[3]), cluster_set, 0)]

            sequences = [node.nt_sequence for node in aa_subset]
            if not sequences:
                break

            if not no_dupes:
                nt_sequences = [(node.header, node.nt_sequence) for node in aa_subset]
                consensus_seq = make_duped_consensus(
                    nt_sequences, excise_consensus
                )
            else:
                consensus_seq = dumb_consensus(sequences, excise_consensus, 0)

            consensus_seq = convert_consensus(sequences, consensus_seq)
            region_start, region_end = check_covered_bad_regions(consensus_seq, excise_minimum_ambig)
            if region_start == last_region.get(cluster_i, -1):
                recursion_max -= 1
                if recursion_max <= 0:
                    region_start = None
                    continue
            last_region[cluster_i] = region_start

            if region_start:
                has_region = True
                had_region = True
                sequences_in_region = []
                sequences_out_of_region = []      

                for i, node in enumerate(aa_subset):

                    if node.header in kicked_headers:
                        continue

                    overlap_coords = get_overlap(region_start, region_end, node.start * 3, node.end * 3, 1)
                    if overlap_coords:
                        overlap_amount = overlap_coords[1] - overlap_coords[0]
                        overlap_percent = overlap_amount / (node.end - node.start)
                        
                        
                        
                        if overlap_percent >= min_overlap_percent or overlap_amount >= min_overlap_amount:
                            sequences_in_region.append(aa_subset[i])
                        else:
                            sequences_out_of_region.append(aa_subset[i])
                
                tagged_in_region = [(int(node.header.split("|")[3].split("&&")[0].split("_")[1]), node) for node in sequences_in_region]
                tagged_in_region.sort(key=lambda x: x[0])
                clusters = cluster(tagged_in_region, true_cluster_threshold)

                log_output.append(f">{gene}_ambig_{region_start}:{region_end}\n{consensus_seq}")

                for clust in clusters:
                    if len(clust) <= 1:
                        continue

                    clust.sort(key = lambda node: node.start)
                    for (_, prev_node), (i, node) in combinations(enumerate(clust), 2):
                        overlapping_coords = get_overlap(node.start, node.end, prev_node.start, prev_node.end, 1)
                        if overlapping_coords:
                            kmer = node.sequence[overlapping_coords[0]:overlapping_coords[1]]
                            prev_kmer = prev_node.sequence[overlapping_coords[0]:overlapping_coords[1]]

                            if is_same_kmer(kmer, prev_kmer):
                                continue

                            splice_index = calculate_split(prev_node, node, overlapping_coords, ref_consensus)
                            prev_positions = set()
                            node_positions = set()
                            
                            consecutive_match = True
                            for x in range(splice_index - 1, node.start - 1, -1):
                                if consecutive_match and node.sequence[x] == prev_node.sequence[x]:
                                    continue
                                
                                consecutive_match = False
                                node_positions.add(x)

                            for x in range(splice_index, prev_node.end):
                                prev_positions.add(x)

                            log_output.append(f">{prev_node.header} vs {node.header}")
                            log_output.append(f"Split at {splice_index}")

                            prev_node.sequence = del_cols(prev_node.sequence, prev_positions)
                            node.sequence = del_cols(node.sequence, node_positions)

                            if len(prev_node.sequence) - prev_node.sequence.count("-") < 15:
                                log_output.append(f"Kicking {prev_node.header} due to < 15 bp after splice")
                                kicked_headers.add(prev_node.header)
                                

                            if len(node.sequence) - node.sequence.count("-") < 15:
                                log_output.append(f"Kicking {node.header} due to < 15 bp after splice")
                                kicked_headers.add(node.header)

                            prev_node.nt_sequence = del_cols(prev_node.nt_sequence, prev_positions, True)
                            node.nt_sequence = del_cols(node.nt_sequence, node_positions, True)

                            node.start, node.end = find_index_pair(node.sequence, "-")
                            prev_node.start, prev_node.end = find_index_pair(prev_node.sequence, "-")

                            x_positions[node.header].update(node_positions)
                            x_positions[prev_node.header].update(prev_positions)
        if recursion_max <= 0:
            break

    do_trim(aa_nodes, cluster_sets, x_positions, ref_consensus, kicked_headers, excise_trim_consensus, no_dupes)

    if had_region:
        after_data = []
        for node in aa_nodes:
            if node.header in kicked_headers:
                continue
            after_data.append((node.header.split("|")[3].split("&&")[0].split("_")[1], node.header))

        after_data.sort(key = lambda x: x[0])
        after_true_clusters = []

        current_cluster = []
        for child_index, child_full in true_cluster_raw:
            if not current_cluster:
                current_cluster.append(child_full)
                current_index = child_index
            else:
                if child_index - current_index <= true_cluster_threshold:
                    current_cluster.append(child_full)
                    current_index = child_index
                else:
                    after_true_clusters.append(current_cluster)
                    current_cluster = [child_full]
                    current_index = child_index

        if current_cluster:
            after_true_clusters.append(current_cluster)

        after_true_cluster = set(max(after_true_clusters, key=len))
        matches = []
        for i, before_cluster in enumerate(before_true_clusters):
            matches.append((len(after_true_cluster.intersection(set(before_cluster))) / len(before_cluster), i))
        matches.sort(key=lambda x: x[0], reverse= True)
        
        best_match = max(matches, key=lambda x: x[0])[1]
        this_rescues = [f"Rescued in gene: {gene}"]
        for header in before_true_clusters[best_match]:
            if header in kicked_headers:

                node = head_to_node[header]
                
                if not node.nt_sequence:
                    continue
                
                if len(node.nt_sequence) - node.nt_sequence.count("-") < 15:
                    continue
                
                this_sequence = node.sequence
                start, end = find_index_pair(this_sequence, "-")
                matching_char = 0
                for i, let in enumerate(this_sequence[start: end], start):
                    if let in ref_consensus[i]:
                        matching_char += 1

                if matching_char / (end - start) < excise_rescue_match:
                    continue

                kicked_headers.remove(header)
                this_rescues.append(header)

    aa_output = [(header, del_cols(seq, x_positions[header])) for header, seq in raw_aa if header not in kicked_headers]

    aa_has_candidate = False
    for header, _ in aa_output:
        if not header.endswith("."):
            aa_has_candidate = True
            break

    if aa_has_candidate:
        gene_coverage = get_coverage([seq for header, seq in aa_output if header[-1] != "."], ref_avg_len)
        if gene_coverage < 0.4:
            log_output.append(f">{gene}_kicked_coverage_{gene_coverage}_of_0.4\n{consensus_seq}")
            return log_output, False, False, gene, False, len(aa_nodes), this_rescues
        
        writeFasta(aa_out, aa_output, compress_intermediates)
        nt_output = [(header, del_cols(seq, x_positions[header], True)) for header, seq in raw_sequences.items() if header not in kicked_headers]
        writeFasta(nt_out, nt_output, compress_intermediates)

        return log_output, had_region, False, False, gene, len(kicked_headers), this_rescues
    return log_output, had_region, gene, False, None, len(kicked_headers), this_rescues


def main(args, sub_dir):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if args.excise_consensus > 1.0:
        if 0 < args.excise_consensus <= 100:
            args.excise_consensus = args.excise_consensus / 100
        else:
            raise ValueError(
                "Cannot convert excise_consensus to a percent. Use a decimal or a whole number between 0 and 100"
            )
    folder = args.INPUT
    input_folder = Path(folder, "outlier", sub_dir)
    if not input_folder.exists():
        input_folder = Path(folder, sub_dir)

    output_folder = Path(folder, "outlier", "excise")

    output_aa_folder = output_folder.joinpath("aa")
    output_nt_folder = output_folder.joinpath("nt")

    if not path.exists(output_aa_folder):
        makedirs(str(output_aa_folder), exist_ok=True)
        makedirs(str(output_nt_folder), exist_ok=True)

    aa_input = input_folder.joinpath("aa")

    no_dupes = False

    genes = [fasta for fasta in listdir(aa_input) if ".fa" in fasta]
    log_path = Path(output_folder, "excise_regions.txt")
    gene_log_path = Path(output_folder, "excise_genes.txt")
    if args.processes > 1:
        arguments = [
            (
                gene,
                input_folder,
                output_folder,
                args.compress,
                args.excise_consensus,
                args.excise_minimum_ambig,
                args.excise_rescue_match,
                args.excise_trim_consensus,
                no_dupes,
            )
            for gene in genes
        ]
        with Pool(args.processes) as pool:
            results = pool.starmap(log_excised_consensus, arguments)
    else:
        results = []
        for gene in genes:
            results.append(
                log_excised_consensus(
                    gene,
                    input_folder,
                    output_folder,
                    args.compress,
                    args.excise_consensus,
                    args.excise_minimum_ambig,
                    args.excise_rescue_match,
                    args.excise_trim_consensus,
                    no_dupes,
                )
            )

    kicked_sequences = 0
    log_output = []
    loci_containing_bad_regions = 0
    kicked_no_resolve = []
    kicked_coverage = []
    has_resolution = []
    no_ambig = []
    rescues = []

    for glog, ghas_ambig, ghas_no_resolution, gcoverage_kick, g_has_resolution, gkicked_seq, grescues in results:
        log_output.extend(glog)
        rescues.extend(grescues)
        if ghas_ambig:
            loci_containing_bad_regions += 1
        kicked_sequences += gkicked_seq
        if ghas_no_resolution:
            kicked_no_resolve.append(ghas_no_resolution)
        if gcoverage_kick:
            kicked_coverage.append(gcoverage_kick)
        if g_has_resolution:
            if ghas_ambig:
                has_resolution.append(g_has_resolution)
            else:
                no_ambig.append(g_has_resolution)



    printv(
        f"{input_folder}: {loci_containing_bad_regions} ambiguous loci found. {len(kicked_coverage)} genes kicked due to low coverage. Kicked {kicked_sequences} sequences total.",
        args.verbose,
    )

    if args.debug:
        with open(log_path, "w") as f:
            f.write("Gene,Cut-Indices,Consensus Sequence\n")
            f.write("\n".join(log_output))
        with open(gene_log_path, "w") as f:
            f.write(f"{len(kicked_no_resolve)} gene(s) kicked due to no seqs with resolution:\n")
            f.write("\n".join(kicked_no_resolve))
            f.write("\n\n")
            f.write(f"{len(kicked_coverage)} gene(s) kicked due to low coverage:\n")
            f.write("\n".join(kicked_coverage))
            f.write("\n\n")
            f.write(f"{len(has_resolution)} gene(s) with resolution and ambig:\n")
            f.write("\n".join(has_resolution))
            f.write("\n\n")
            f.write(f"{len(no_ambig)} gene(s) resolve with no ambig:\n")
            f.write("\n".join(no_ambig))
            f.write("\n\n")
        with open(Path(output_folder, "rescues.txt"), "w") as f:
            f.write("\n".join(rescues))


    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)
    if len(genes) == 0:
        printv("WARNING: No genes in output.", args.verbose, 0)
        return True
    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )