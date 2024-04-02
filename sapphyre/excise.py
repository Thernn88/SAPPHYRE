from collections import defaultdict
from multiprocessing import Pool
from os import listdir, mkdir, path
from pathlib import Path
from shutil import move
import copy
from msgspec import Struct, json
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    constrained_distance,
)
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta


def min_aa_check(sequence: list, minimum: int) -> bool:
    """
    Checks if the given sequence has at least the minimum amount of data characters.
    '-' and 'X' are not counted as data. All other characters are counted as data.
    Given sequence must be ascii characters.
    Returns the answer as a boolean.
    If enough data is found, returns True. Otherwise, returns False.
    """
    valid_count = len(sequence) - sequence.count("-") - sequence.count("X")
    return valid_count >= minimum


def first_last_X(string: str) -> tuple:
    """
    Finds the first and last occurrence of 'X' in the given string.
    If found, returns a tuple containing the slice indices [first, last).
    First in inclusive, last is exclusive.
    If no 'X' is found in the string, returns (None, None).
    """
    first = None
    last = None
    for i, bp in enumerate(string):
        if bp == "X":
            first = i
            break
    for i in range(len(string) - 1, -1, -1):
        if string[i] == "X":
            last = i
            break
    if first is None or last is None:
        return None, None
    return first, last + 1

def check_bad_regions(
    consensus: str, limit: float, initial_window=16, offset=0, _rig=False
) -> list:
    """
    Checks a consensus sequence for regions containing too many 'X' characters.

    Iterates over the consensus using a sliding window.
    If the required number of 'X' are found, expands the end of the window until
    the ratio of X falls below the given threshold, then reduces the window to the
    indices of the first and last 'X' to make the final result.

    If the required number of 'X' are not found, the current window is dumped and the
    exclusive ending index of the current window becomes the inclusive starting index of the
    next window.

    Consensus is the region of the consensus sequence to be searched.

    Limit is a float between 0.0 and 1.0 that represents the minimum ratio of
    X to all characters to flag a region.

    Initial_window is an int that sets the starting length of the window.

    Offset is the index of the first character to be searched. By passing this to
    the function and adding it to the results, we can search an arbitrary slice of the
    consensus sequence and return indices that are usable without modification.

    Rig is a debugging flag which can be ignored by a regular user.
    """
    if not _rig:
        first, last = 0, len(consensus)
    else:
        first, last = find_index_pair(consensus, "X")
    if first is None or last is None:
        return []
    left, right = first, first + initial_window
    num_x = consensus[left:right].count("X")
    begin = None
    output = {}

    while right <= last:
        ratio = num_x / (right - left)
        if ratio > limit:  # if ratio is above limit, this is a bad region
            if begin is None:  # if begin is None, this is the start of a region
                begin = left
            if right == len(consensus):
                # This break stops an index error. The region will still be
                # caught by the last check in this function.
                break
            # increment the X count, then expand the window
            # without this we have to recount all the characters in the window
            num_x += consensus[right] == "X"
            right += 1
            continue
        # If the following block executes, it means the ratio of X is below the threshold.
        if begin is not None:  # end of region, so make output
            a, b = first_last_X(consensus[begin : right + 1])
            a, b = a + begin + offset, b + begin + offset
            if b not in output:
                output[b] = (a, b)
        # since the current window is below the threshold, discard it and make a new one
        left = right
        right = left + initial_window
        begin = None
        num_x = consensus[left:right].count("X")
    # this check catches bad regions at the end of the string
    if begin is not None:
        a, b = first_last_X(consensus[begin:last])
        a, b = a + begin + offset, b + begin + offset
        if a is not None and b is not None:
            if b not in output:
                output[b] = (a, b)
    return [*output.values()]

def check_covered_bad_regions(consensus, min_ambiguous, ambig_char='X', max_distance=18):
    x_regions = []
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
                x_regions.append((current_group[0], current_group[-1] + 1))
            current_group = [num]

    if current_group:
        if len(current_group) >= min_ambiguous:
            x_regions.append((current_group[0], current_group[-1] + 1))


    return x_regions

def bundle_seqs_and_dupes(sequences: list, prepare_dupe_counts, reporter_dupe_counts):
    output = []
    for header, seq in sequences:
        node = header.split("|")[3]
        dupes = prepare_dupe_counts.get(node, 1) + sum(
            prepare_dupe_counts.get(node, 1)
            for node in reporter_dupe_counts.get(node, [])
        )
        output.append((seq, dupes))
    return output


def make_duped_consensus(
    raw_sequences: list, prepare_dupes: dict, reporter_dupes: dict, threshold: float
) -> str:
    seqs = [(header, seq) for header, seq in raw_sequences if header[-1] != "."]
    bundled_seqs = bundle_seqs_and_dupes(seqs, prepare_dupes, reporter_dupes)
    return dumb_consensus_dupe(bundled_seqs, threshold, 0)


class NODE(Struct):
    header: str
    sequence: str
    start: int
    end: int
    children: list

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

        # Save node_2 and the children of node_2 to self
        self.children.append(node_2.header)
        self.children.extend(node_2.children)

    def clone(self):
        return NODE(self.header, self.sequence, self.start, self.end, self.children)
    
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


def simple_assembly(nodes, min_merge_overlap_percent):
    merged = set()
    for i, node in enumerate(nodes):
        if i in merged:
            continue
        merge_occured = True
        while merge_occured:
            merge_occured = False
            for j, node_b in enumerate(nodes):
                if j in merged:
                    continue
                if i == j:
                    continue

                overlap_coords=  get_overlap(node.start, node.end, node_b.start, node_b.end, 1)
                if overlap_coords:
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    overlap_percent = overlap_amount / (node.end - node.start)
                    if overlap_percent < min_merge_overlap_percent:
                        continue


                    kmer_a = node.sequence[overlap_coords[0]:overlap_coords[1]]
                    kmer_b = node_b.sequence[overlap_coords[0]:overlap_coords[1]]

                    if not is_same_kmer(kmer_a, kmer_b):
                        continue

                    merged.add(j)

                    overlap_coord = overlap_coords[0]
                    merge_occured = True
                    node.extend(node_b, overlap_coord)

                    nodes[j] = None
    
    nodes = [node for node in nodes if node is not None]

    return nodes
        

def identity(nodes_in_region, best_index, allowed_mismatches):
    indices = [i for i, node in enumerate(nodes_in_region) if i != best_index]
    best_node = nodes_in_region[best_index]
    best_length = len(best_node.sequence) - best_node.sequence.count("-")
    extend_region = set()
    for i in indices:
        node = nodes_in_region[i]
        distance = constrained_distance(best_node.sequence, node.sequence)
        node_length = len(node.sequence) - node.sequence.count("-")

        # if node_length < best_length * (1 - length_percentage):
        #     continue

        if distance < allowed_mismatches:
            extend_region.add(i)

    return extend_region


def indices_that_resolve(nodes, sequences_out_of_region, merge_percent):
    """
    This function is used to find the longest sequence in the list of sequences_out_of_region
    that can be merged into the list of nodes. If a sequence is found that can be merged,
    it is removed from the list of sequences_out_of_region and added to the list of nodes.
    """
    indices = []
    for i, node in enumerate(nodes):
        this_seq_out = sequences_out_of_region.copy()
        has_unambig_overlap = False
        for j, node_b in enumerate(this_seq_out):
            overlap_coords = get_overlap(node.start, node.end, node_b.start, node_b.end, 1)
            
            if overlap_coords:
                overlap_amount = overlap_coords[1] - overlap_coords[0]
                overlap_percent = overlap_amount / (node.end - node.start)
                
                if overlap_percent > merge_percent:
                    has_unambig_overlap = True

        if has_unambig_overlap:
            indices.append(i)

    return indices
        

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


def get_coverage(consensus_seq, aa_nodes, ref_avg_len):
    cons_start, cons_end = find_index_pair(consensus_seq, "-")
    data_cols = 0
    for i in range(cons_start, cons_end):
        for node in aa_nodes:
            if i >= node.start and i < node.end:
                data_cols += 1
                break

    return data_cols / ref_avg_len


def log_excised_consensus(
    gene: str,
    is_assembly_or_genome: bool,
    input_path: Path,
    output_path: Path,
    compress_intermediates: bool,
    excise_overlap_merge,
    excise_overlap_ambig,
    excise_region_overlap,
    excise_consensus,
    excise_maximum_depth,
    excise_minimum_ambig,
    allowed_distance,
    prepare_dupes: dict,
    reporter_dupes: dict,
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

    aa_in = input_path.joinpath("aa", gene)
    aa_out = output_path.joinpath("aa", gene)

    nt_in = input_path.joinpath("nt", gene.replace(".aa.", ".nt."))
    nt_out = output_path.joinpath("nt", gene.replace(".aa.", ".nt."))

    x_positions = defaultdict(set)

    bp_count = lambda x: len(x) - x.count("-")

    raw_aa = list(parseFasta(str(aa_in)))

    ref_lens = []
    for header, seq in raw_aa:
        if header.endswith('.'):
            ref_lens.append(bp_count(seq))

    ref_avg_len = sum(ref_lens) / len(ref_lens)

    aa_nodes = [NODE(header, seq, *find_index_pair(seq,"-"), []) for header, seq in raw_aa if header[-1] != "."]
    
    if prepare_dupes and reporter_dupes:
        consensus_seq = make_duped_consensus(
            raw_aa, prepare_dupes, reporter_dupes, excise_consensus
        )
    else:
        aa_sequences = [x.sequence for x in aa_nodes]
        consensus_seq = dumb_consensus(aa_sequences, excise_consensus, 0)
        
    
    gene_coverage = get_coverage(consensus_seq, aa_nodes, ref_avg_len)
    req_coverage = 0.4 if is_assembly_or_genome else 0.01
    if gene_coverage < req_coverage:
        log_output.append(f">{gene}_kicked_coverage_{gene_coverage}_of_{req_coverage}\n{consensus_seq}")
        return log_output, False, True, len(aa_nodes)

    TRIM_MAX = 6

    for node in aa_nodes:
        new_start = node.start
        new_end = node.end

        this_positions = set()

        for i, bp in enumerate(node.sequence):
            let = consensus_seq[i]
            if let == "X":
                if (i <= new_start + 2 and i >= new_start):
                    new_start = i
                    if new_start - node.start >= TRIM_MAX:
                        break
                if (i >= new_end - 2 and i <= new_end):
                    new_end = i
                    if node.end - new_end >= TRIM_MAX:
                        break

        if new_start != node.start or new_end != node.end:
            for i in range(node.start, new_start):
                this_positions.add(i)
                x_positions[node.header].add(i)
            for i in range(new_end, node.end):
                this_positions.add(i)
                x_positions[node.header].add(i)
            node.start = new_start
            node.end = new_end

        node.sequence = "".join(let if i not in this_positions else "-" for i, let in enumerate(node.sequence))

    for node in aa_nodes:
        this_positions = set()
        i = None
        for poss_i in range(node.start, node.start + 3):
            if node.sequence[poss_i] != consensus_seq[poss_i]:
                i = poss_i

        if not i is None:
            for x in range(node.start , i + 1):
                x_positions[node.header].add(x)
            node.start = i+1

        i = None
        for poss_i in range(node.end -1, node.end - 4, -1):
            if node.sequence[poss_i] != consensus_seq[poss_i]:
                i = poss_i

        if not i is None:
            for x in range(i, node.end):
                x_positions[node.header].add(x)
            node.end = i

    raw_sequences = [(header, del_cols(seq, x_positions[header], True)) for header, seq in parseFasta(str(nt_in))]
    sequences = [x[1] for x in raw_sequences if x[0][-1] != "."]

    nodes = []
    for header, seq in raw_sequences:
        if header[-1] == ".":
            continue
        start, end = find_index_pair(seq, "-")
        nodes.append(
            NODE(header, seq, start, end,  [])
        )

    if prepare_dupes and reporter_dupes:
        consensus_seq = make_duped_consensus(
            raw_sequences, prepare_dupes, reporter_dupes, excise_consensus
        )
    else:
        consensus_seq = dumb_consensus(sequences, excise_consensus, 0)

    consensus_seq = convert_consensus(sequences, consensus_seq)

    # Search for slices of the consensus seq with a high ratio of 'X' to total characters
    bad_regions = check_covered_bad_regions(consensus_seq, excise_minimum_ambig)

    kicked_headers = set()
    if bad_regions:
        for region in bad_regions:
            sequences_in_region = []
            sequences_out_of_region = []
            a, b = region         

            for i, node in enumerate(nodes):

                if node.header in kicked_headers:
                    continue

                overlap_coords = get_overlap(a, b, node.start, node.end, 1)
                if overlap_coords:
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    overlap_percent = overlap_amount / (node.end - node.start)
                    if overlap_percent >= excise_region_overlap: # Adjustable percent
                        sequences_in_region.append(nodes[i].clone())
                    else:
                        sequences_out_of_region.append(nodes[i].clone())

            if len(sequences_in_region) > excise_maximum_depth:
                continue

            nodes_in_region = simple_assembly(sequences_in_region, excise_overlap_ambig)

            copy_of_nodes_in_region = copy.deepcopy(nodes_in_region)

            node_indices = indices_that_resolve(copy_of_nodes_in_region, sequences_out_of_region, excise_overlap_merge)

            if not node_indices:
                continue
            
            best_index = max(node_indices, key=lambda x: len(nodes_in_region[x].sequence) - nodes_in_region[x].sequence.count("-"))

            keep_indices = set()
            if best_index is not None:
                keep_indices.update([best_index])

                similar_indices = identity(nodes_in_region, best_index, allowed_distance)
                keep_indices.update(similar_indices)

            for i, node in enumerate(nodes_in_region):
                if i in keep_indices:
                    continue

                kicked_headers.add(node.header)
                kicked_headers.update(node.children)

            if nodes_in_region:
                log_output.append(f">{gene}_ambig_{a}:{b}\n{consensus_seq}")
                log_output.extend([f">{node.contig_header()}_{'kept' if i in keep_indices else 'kicked'}\n{node.sequence}" for i, node in enumerate(nodes_in_region)])
                log_output.append("\n")

    aa_output = [(header, del_cols(seq, x_positions[header])) for header, seq in raw_aa if header not in kicked_headers]

    aa_has_candidate = False
    for header, _ in aa_output:
        if not header.endswith("."):
            aa_has_candidate = True
            break

    if aa_has_candidate:
        writeFasta(aa_out, aa_output, compress_intermediates)
        nt_output = [(header, del_cols(seq, x_positions[header], True)) for header, seq in raw_sequences if header not in kicked_headers]
        writeFasta(nt_out, nt_output, compress_intermediates)

    return log_output, bad_regions != [], False, len(kicked_headers)


def load_dupes(rocks_db_path: Path):
    rocksdb_db = RocksDB(str(rocks_db_path))
    prepare_dupe_counts = json.decode(
        rocksdb_db.get("getall:gene_dupes"), type=dict[str, dict[str, int]]
    )
    reporter_dupe_counts = json.decode(
        rocksdb_db.get("getall:reporter_dupes"), type=dict[str, dict[str, list]]
    )
    return prepare_dupe_counts, reporter_dupe_counts


### USED BY __main__
def do_move(from_, to_):
    move(from_, to_)


def move_flagged(to_move, processes):
    if processes <= 1:
        for from_, to_ in to_move:
            move(from_, to_)
    else:
        with Pool(processes) as pool:
            pool.starmap(do_move, to_move)


###


def main(args, override_cut, sub_dir, is_assembly_or_genome):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if not (0 < args.excise_overlap_merge < 1.0):
        if 0 < args.excise_overlap_merge <= 100:
            args.excise_overlap_merge = args.excise_overlap_merge / 100
        else:
            raise ValueError(
                "Cannot convert excise_overlap_merge to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.excise_overlap_ambig < 1.0):
        if 0 < args.excise_overlap_ambig <= 100:
            args.excise_overlap_ambig = args.excise_overlap_ambig / 100
        else:
            raise ValueError(
                "Cannot convert excise_overlap_ambig to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.excise_consensus < 1.0):
        if 0 < args.excise_consensus <= 100:
            args.excise_consensus = args.excise_consensus / 100
        else:
            raise ValueError(
                "Cannot convert excise_consensus to a percent. Use a decimal or a whole number between 0 and 100"
            )
    if not (0 < args.excise_region_overlap < 1.0):
        if 0 < args.excise_region_overlap <= 100:
            args.excise_region_overlap = args.excise_region_overlap / 100
        else:
            raise ValueError(
                "Cannot convert excise_region_overlap to a percent. Use a decimal or a whole number between 0 and 100"
            )

    folder = args.INPUT
    input_folder = Path(folder, "outlier", sub_dir)
    output_folder = Path(folder, "outlier", "excise")

    output_aa_folder = output_folder.joinpath("aa")
    output_nt_folder = output_folder.joinpath("nt")

    if not path.exists(output_folder):
        mkdir(str(output_folder))
        mkdir(str(output_aa_folder))
        mkdir(str(output_nt_folder))

    aa_input = input_folder.joinpath("aa")

    rocksdb_path = Path(folder, "rocksdb", "sequences", "nt")
    if args.no_dupes:
        prepare_dupes, reporter_dupes = {}, {}
    else:
        if not rocksdb_path.exists():
            err = f"cannot find rocksdb for {folder}"
            raise FileNotFoundError(err)
        prepare_dupes, reporter_dupes = load_dupes(rocksdb_path)

    compress = not args.uncompress_intermediates or args.compress

    genes = [fasta for fasta in listdir(aa_input) if ".fa" in fasta]
    log_path = Path(output_folder, "excise_regions.txt")
    if args.processes > 1:
        arguments = [
            (
                gene,
                is_assembly_or_genome,
                input_folder,
                output_folder,
                compress,
                args.excise_overlap_merge,
                args.excise_overlap_ambig,
                args.excise_region_overlap,
                args.excise_consensus,
                args.excise_maximum_depth,
                args.excise_minimum_ambig,
                args.excise_allowed_distance,
                prepare_dupes.get(gene.split(".")[0], {}),
                reporter_dupes.get(gene.split(".")[0], {}),
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
                    is_assembly_or_genome,
                    input_folder,
                    output_folder,
                    compress,
                    args.excise_overlap_merge,
                    args.excise_overlap_ambig,
                    args.excise_region_overlap,
                    args.excise_consensus,
                    args.excise_maximum_depth,
                    args.excise_minimum_ambig,
                    args.excise_allowed_distance,
                    prepare_dupes.get(gene.split(".")[0], {}),
                    reporter_dupes.get(gene.split(".")[0], {}),
                )
            )

    log_output = ["\n".join(x[0]) for x in results]
    loci_containing_bad_regions = len([1 for x in results if x[1]])
    kicked_coverage = len([1 for x in results if x[2]])
    kicked_sequences = sum(x[3] for x in results)
    printv(
        f"{input_folder}: {loci_containing_bad_regions} ambiguous loci found. {kicked_coverage} genes kicked due to low coverage. Kicked {kicked_sequences} sequences total.",
        args.verbose,
    )

    if args.debug:
        with open(log_path, "w") as f:
            f.write("Gene,Cut-Indices,Consensus Sequence\n")
            f.write("".join(log_output))

    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)
    if len(genes) == 0:
        printv("WARNING: No genes in output.", args.verbose, 0)
        return True, True
    return True, loci_containing_bad_regions / len(genes) >= args.majority_excise


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
