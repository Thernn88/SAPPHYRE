from multiprocessing import Pool
import os
from shutil import move, rmtree
from pathlib import Path
import phymmr_tools as bd
import wrap_rocks

from msgspec import json
from .utils import parseFasta, writeFasta, printv
from .timekeeper import KeeperMode, TimeKeeper


def min_aa_check(sequence: list, minimum: int) -> bool:
    valid_count = len(sequence) - sequence.count("-") - sequence.count("X")
    return valid_count >= minimum


def find_quesion_marks(sequences: list, start: int, stop: int) -> set:
    def check_index(index: int):
        for seq in sequences:
            if seq[i] != "-":
                return False
        return True

    output = set()
    for i in range(start, stop):
        if check_index(i):
            output.add(i)
    return output


def replace_inner_gaps(seq: list, start, stop, question_marks) -> str:
    # start = None
    for i, bp in enumerate(seq[start:stop], start):
        if i not in question_marks:
            continue
        if bp == "?":
            if start is None:
                start = i
                continue
            if start is not None:
                if i - start >= 7:
                    for index in range(start, i + 1):
                        seq[index] = "-"
                start = None
    if start is not None:
        if len(seq) - start >= 7:
            for index in range(start, i + 1):
                seq[index] = "-"
    return seq


def first_last_X(string: str) -> tuple:
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


def sensitive_check_bad_regions(consensus: str, limit: float, window=16) -> list:
    # window = initial_window
    # total = window
    first, last = bd.find_index_pair(consensus, "X")
    if first is None or last is None:
        return []
    begin = None
    output = []
    for i in range(first, last, window):
        ratio = consensus[i : i + window].count("X") / window
        if ratio >= limit:
            if begin is None:
                begin = i
        else:
            if begin is not None:
                # the new first is the old window end
                a, b = first_last_X(consensus[begin:i])
                output.append((a + begin, b + begin))
                begin = None
    if begin is not None:
        a, b = first_last_X(consensus[begin:last])
        if a is not None and b is not None:
            output.append((a + begin, b + begin))

    return output


def find_coverage_regions(consensus: str, start: int, stop: int, ambiguous="?") -> list:
    output = []
    last_seen_data = None
    for i, bp in enumerate(consensus[start:stop], start):
        if bp is not ambiguous:
            if last_seen_data is None:  # start of data region
                last_seen_data = i
        else:  # bp is ambiguous character
            if last_seen_data is not None:
                output.append((last_seen_data, i))
            last_seen_data = None
    if last_seen_data is not None:
        output.append((last_seen_data, stop))
    return output


def check_bad_regions(
    consensus: str, limit: float, initial_window=16, offset=0, test=False, rig=False
) -> list:
    # if start is None and stop is None:
    # else:
    if not rig:
        first, last = 0, len(consensus)
    else:
        first, last = bd.find_index_pair(consensus, "X")
    if first is None or last is None:
        return []
    l, r = first, first + initial_window
    num_x = consensus[l:r].count("X")
    begin = None
    output = {}

    while r <= last:
        ratio = num_x / (r - l)
        if ratio > limit:
            if begin is None:
                begin = l
            if r == len(consensus):
                break
            num_x += consensus[r] == "X"
            r += 1
            continue
        if begin is not None:
            a, b = first_last_X(consensus[begin:last])
            a, b = a + begin + offset, b + begin + offset
            if b not in output:
                output[b] = (a, b)
        l = r
        r = l + initial_window
        begin = None
        num_x = consensus[l:r].count("X")
    if begin is not None:
        a, b = first_last_X(consensus[begin:last])
        a, b = a + begin + offset, b + begin + offset
        if a is not None and b is not None:
            if b not in output:
                output[b] = (a, b)
    return [*output.values()]


def check_covered_bad_regions(
    consensus: str, limit: float, initial_window=16, ambiguous="?", test=False
) -> list:
    """
    Creates a list of covered indices, then checks those slices for bad regions
    """
    start, stop = bd.find_index_pair(consensus, "X")
    bad_regions = []
    covered_indices = find_coverage_regions(consensus, start, stop, ambiguous=ambiguous)
    for begin, end in covered_indices:
        subregions = check_bad_regions(
            consensus[begin:end],
            limit,
            offset=begin,
            initial_window=initial_window,
            test=test,
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
    return bd.dumb_consensus_dupe(bundled_seqs, threshold)


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
    verbose,
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
    log_output = ""

    aa_in = input_path.joinpath("aa", gene)
    aa_out = output_path.joinpath("aa", gene)

    nt_in = input_path.joinpath("nt", gene.replace(".aa.", ".nt."))
    nt_out = output_path.joinpath("nt", gene.replace(".aa.", ".nt."))

    raw_sequences = list(parseFasta(str(aa_in)))
    sequences = [x[1] for x in raw_sequences if x[0][-1] != "."]
    if prepare_dupes and reporter_dupes:
        consensus_seq = make_duped_consensus(
            raw_sequences, prepare_dupes, reporter_dupes, consensus_threshold
        )
    else:
        consensus_seq = bd.dumb_consensus(sequences, consensus_threshold)
    consensus_seq = bd.convert_consensus(sequences, consensus_seq)
    # if verbose:
    #     print(f"{gene}\n{consensus_seq}\n")
    old_bad_regions = check_bad_regions(consensus_seq, excise_threshold, rig=True)
    bad_regions = check_covered_bad_regions(consensus_seq, excise_threshold, test=False)
    if old_bad_regions and not bad_regions:
        bad_regions = check_covered_bad_regions(
            consensus_seq, excise_threshold, test=True
        )
    if bad_regions:
        if len(bad_regions) == 1:
            a, b = bad_regions[0]
            if b - a != len(consensus_seq):
                log_output += (
                    f"{gene}\t{a}:{b}\t{consensus_seq}\n"
                    if debug
                    else f"{gene}\t{a}:{b}\n"
                )
        else:
            for region in bad_regions:
                a, b = region
                log_output += (
                    f"{gene}\t{a}:{b}\t{consensus_seq}\n"
                    if debug
                    else f"{gene}\t{a}:{b}\n"
                )

        if cut:
            positions_to_cull = [i for a, b in bad_regions for i in range(a, b)]

            aa_output = []
            aa_kicks = set()
            for header, sequence in raw_sequences:
                if header[-1] == ".":
                    aa_output.append((header, sequence))
                    continue
                sequence = list(sequence)
                for i in positions_to_cull:
                    sequence[i] = "-"

                sequence = "".join(sequence)
                if not min_aa_check(sequence, 20):
                    aa_kicks.add(header)
                    continue
                aa_output.append((header, sequence))
            writeFasta(str(aa_out), aa_output, compress_intermediates)

            nt_output = []
            for header, sequence in parseFasta(str(nt_in)):
                if header[-1] == ".":
                    nt_output.append((header, sequence))
                    continue

                if header in aa_kicks:
                    continue

                sequence = list(sequence)
                for i in positions_to_cull:
                    sequence[i * 3 : i * 3 + 3] = ["-", "-", "-"]
                sequence = "".join(sequence)
                nt_output.append((header, sequence))
            writeFasta(str(nt_out), nt_output, compress_intermediates)
    else:
        if debug:
            log_output += f"{gene}\tN/A\t{consensus_seq}\n"

    return log_output, bad_regions != []


def load_dupes(rocks_db_path: Path):
    rocksdb_db = wrap_rocks.RocksDB(str(rocks_db_path))
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


def main(args, override_cut=None):
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

    if cut:
        sub_dir = "collapsed"
    else:
        sub_dir = "blosum"

    input_folder = Path(folder, "outlier", sub_dir)
    output_folder = Path(folder, "outlier", "excise")

    output_aa_folder = output_folder.joinpath("aa")
    output_nt_folder = output_folder.joinpath("nt")

    if os.path.exists(output_folder):
        rmtree(output_folder)
    os.mkdir(str(output_folder))
    os.mkdir(str(output_aa_folder))
    os.mkdir(str(output_nt_folder))

    aa_input = input_folder.joinpath("aa")

    rocksdb_path = Path(folder, "rocksdb", "sequences", "nt")

    if args.dupes:
        if not rocksdb_path.exists():
            err = f"cannot find rocksdb for {folder}"
            raise FileNotFoundError(err)
        prepare_dupes, reporter_dupes = load_dupes(rocksdb_path)
    else:
        prepare_dupes, reporter_dupes = {}, {}

    compress = not args.uncompress_intermediates or args.compress

    genes = [fasta for fasta in os.listdir(aa_input) if ".fa" in fasta]
    log_path = Path(output_folder, "excise_indices.tsv")
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
                args.verbose,
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
                    args.verbose,
                    cut,
                )
            )

    log_output = [x[0] for x in results]
    loci_containing_bad_regions = len([x[1] for x in results if x[1]])
    printv(
        f"{input_folder}: {loci_containing_bad_regions} bad loci found", args.verbose
    )

    if args.debug:
        with open(log_path, "w") as f:
            f.write(f"Gene\tCut-Indices\n")
            f.write("".join(log_output))

    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)
    return True, loci_containing_bad_regions / len(genes) >= args.majority_excise


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
