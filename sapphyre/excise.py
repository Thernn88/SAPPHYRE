from multiprocessing import Pool
from os import listdir, mkdir, path
from pathlib import Path
from shutil import move

from msgspec import json
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
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


def find_coverage_regions(consensus: str, start: int, stop: int, ambiguous="?") -> list:
    """
    Searches a given string slice for non-ambiguous regions.
    The slice is consensus[start:stop]. Start is inclusive, stop is exclusive.
    Returns a list of tuples. Each tuple is the [start,stop) indices for a region.
    """
    output = []
    last_seen_data = None
    for i, bp in enumerate(consensus[start:stop], start):
        # if last_seen_data is unset, we are in an ambiguous region
        # if it is set, we are in a data region
        if bp is not ambiguous:
            # this is the first character in a non-ambiguous region
            if last_seen_data is None:  # start of data region
                last_seen_data = i
        else:  # bp is ambiguous character
            if last_seen_data is not None:
                # this means the function just transitioned from a data region to
                # an ambiguous region, so output tuple and reset variables
                output.append((last_seen_data, i))
            last_seen_data = None
    # check if the string end in a non-ambiguous region
    if last_seen_data is not None:
        output.append((last_seen_data, stop))
    return output


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


def check_covered_bad_regions(
    consensus: str, limit: float, initial_window=16, ambiguous="?"
) -> list:
    """
    Creates a list of covered indices, then checks those slices for bad regions.

    Ignores leading and trailing 'X'. Makes a list of
    non-ambiguous slices, then checks each slice for bad regions.
    Returns a list of all found bad regions.

    Consensus is the entire consensus sequence.

    Limit is a float between 0.0 and 1.0 that represents the minimum ratio of
    X to all characters to flag a region.

    Initial_window is an int that sets the starting length of the window.

    Ambiguous is the character used for ambiguous locations. Defaults to '?'.
    """
    start, stop = find_index_pair(consensus, "X")
    bad_regions = []
    covered_indices = find_coverage_regions(consensus, start, stop, ambiguous=ambiguous)
    for begin, end in covered_indices:
        subregions = check_bad_regions(
            consensus[begin:end],
            limit,
            offset=begin,
            initial_window=initial_window,
        )
        if subregions:
            bad_regions.extend(subregions)
    return bad_regions


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


def log_excised_consensus(
    gene: str,
    input_path: Path,
    output_path: Path,
    compress_intermediates: bool,
    consensus_threshold,
    excise_threshold,
    prepare_dupes: dict,
    reporter_dupes: dict,
    debug,
    cut,
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

    raw_sequences = list(parseFasta(str(aa_in)))
    sequences = [x[1] for x in raw_sequences if x[0][-1] != "."]

    nodes = [(header, seq, *find_index_pair(seq, "-")) for header, seq in raw_sequences]

    # Make consensus sequence. Use dupe counts if available.
    if prepare_dupes and reporter_dupes:
        consensus_seq = make_duped_consensus(
            raw_sequences, prepare_dupes, reporter_dupes, consensus_threshold
        )
    else:
        consensus_seq = dumb_consensus(sequences, consensus_threshold, 0)
    consensus_seq = convert_consensus(sequences, consensus_seq)

    # Search for slices of the consensus seq with a high ratio of 'X' to total characters
    bad_regions = check_covered_bad_regions(consensus_seq, excise_threshold)
    kicked_headers = set()
    if bad_regions:
        for region in bad_regions:
            sequences_in_region = []
            a, b = region

            for header,seq,start,end in nodes:
                overlap_coords = get_overlap(a, b, start, end, 1)
                if overlap_coords:
                    overlap_amount = overlap_coords[1] - overlap_coords[0]
                    overlap_percent = overlap_amount / (end - start)

                    if overlap_percent > 0.5: # Adjustable percent
                        sequences_in_region.append((header, seq))

            if sequences_in_region:
                log_output.extend([f">{header}\n{seq}" for header, seq in sequences_in_region])
                log_output.append("\n")
    # if cut:
    #     if bad_regions:
    #         # log bad regions
    #         if len(bad_regions) == 1:
    #             a, b = bad_regions[0]
    #             if b - a != len(consensus_seq):
    #                 if debug:
    #                     log_output.append(f"{gene},{a}:{b},{consensus_seq}\n")
    #                 else:
    #                     log_output.append(f"{gene},{a}:{b}\n")
    #         else:
    #             for region in bad_regions:
    #                 a, b = region
    #                 if debug:
    #                     log_output.append(f"{gene},{a}:{b},{consensus_seq}\n")
    #                 else:
    #                     log_output.append(f"{gene},{a}:{b}\n")

    #         # make a list of locations that will be removed
    #         positions_to_cull = [i for a, b in bad_regions for i in range(a, b)]

    #         has_cand = False
    #         aa_output = []
    #         for header, sequence in raw_sequences:
    #             if header[-1] == ".":  # we don't want to alter reference sequences
    #                 aa_output.append((header, sequence))
    #                 continue
    #             # remove locations from the sequence, then remake string.
    #             sequence = list(sequence)
    #             for i in positions_to_cull:
    #                 sequence[i] = "-"

    #             sequence = "".join(sequence)

    #             # if the sequence no longer has the minimum amount of data characters,
    #             # kick the entire sequence instead of just excising parts of it
    #             if not min_aa_check(sequence, 20):
    #                 kicked_headers.add(header)
    #                 continue

    #             # sequence has enough data after the excision, so output it
    #             aa_output.append((header, sequence))
    #             has_cand = True
    #         if has_cand:
    #             writeFasta(str(aa_out), aa_output, compress_intermediates)

    #         # mirror the excisions in nt sequences
    #         nt_output = []
    #         for header, sequence in parseFasta(str(nt_in)):
    #             if header[-1] == ".":
    #                 nt_output.append((header, sequence))
    #                 continue

    #             if header in kicked_headers:
    #                 continue

    #             sequence = list(sequence)
    #             # sets bad locations to '-' instead of removing them
    #             for i in positions_to_cull:
    #                 sequence[i * 3 : i * 3 + 3] = ["-", "-", "-"]
    #             sequence = "".join(sequence)
    #             nt_output.append((header, sequence))
    #         if has_cand:
    #             writeFasta(str(nt_out), nt_output, compress_intermediates)
    #     else:
    #         if raw_sequences:
    #             writeFasta(str(aa_out), raw_sequences, compress_intermediates)
    #             writeFasta(str(nt_out), parseFasta(str(nt_in)), compress_intermediates)
    #         if debug:
    #             log_output.append(f"{gene},N/A,{consensus_seq}\n")

    return log_output, bad_regions != [], len(kicked_headers)


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


def main(args, override_cut, sub_dir):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if not (0 < args.consensus < 1.0):
        if 0 < args.consensus <= 100:
            args.consensus = args.consensus / 100
        else:
            raise ValueError(
                "Cannot convert consensus to a percent. Use a decimal or a whole number between 0 and 100"
            )

    if not (0 < args.excise < 1.0):
        if 0 < args.excise <= 100:
            args.excise = args.excise / 100
        else:
            raise ValueError(
                "Cannot convert excise to a percent. Use a decimal or a whole number between 0 and 100"
            )
    folder = args.INPUT
    if override_cut is None:
        cut = args.cut
    else:
        cut = override_cut

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
                input_folder,
                output_folder,
                compress,
                args.consensus,
                args.excise,
                prepare_dupes.get(gene.split(".")[0], {}),
                reporter_dupes.get(gene.split(".")[0], {}),
                args.debug,
                cut,
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
                    compress,
                    args.consensus,
                    args.excise,
                    prepare_dupes.get(gene.split(".")[0], {}),
                    reporter_dupes.get(gene.split(".")[0], {}),
                    args.debug,
                    cut,
                )
            )

    log_output = ["\n".join(x[0]) for x in results]
    loci_containing_bad_regions = len([x[1] for x in results if x[1]])
    kicked_sequences = sum(x[2] for x in results if x[1])
    printv(
        f"{input_folder}: {loci_containing_bad_regions} bad loci found. Kicked {kicked_sequences} sequences",
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
