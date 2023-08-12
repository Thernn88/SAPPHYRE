import os
import phymmr_tools as bd
import wrap_rocks

from msgspec import json
from pathlib import Path
from pyfastx import Fastx


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
        if i not in question_marks: continue
        if bp == "?":
            if start is None:
                start = i
                continue
            else:
                if start is not None:
                    if i - start >= 7:
                        for index in range(start, i+1):
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
    for i in range(len(string)-1, -1, -1):
        if string[i] == "X":
            last = i
            break
    if first is None or last is None:
        return None, None
    return first, last+1


def sensitive_check_bad_regions(consensus: str, limit: float, window=16) -> list:
    # window = initial_window
    # total = window
    first, last = bd.find_index_pair(consensus, "X")
    if first is None or last is None:
        return []
    num_x = consensus[first:first+window].count("X")
    in_region = False
    begin = None
    output = []
    for i in range(first, last, window):
        ratio = consensus[i:i+window].count("X") / window
        if ratio >= limit:
            if begin is None:
                begin = i
        else:
            if begin is not None:
                # the new first is the old window end
                a, b = first_last_X(consensus[begin:i])
                output.append((a+begin, b+begin))
                begin = None
    if begin is not None:
        a, b = first_last_X(consensus[begin:last])
        if a is not None and b is not None:
            output.append((a + begin, b + begin))

    return output


def check_bad_regions(consensus: str, limit: float, initial_window=16) -> list:
    first, last = bd.find_index_pair(consensus, "X")
    if first is None or last is None:
        return []
    l, r = first, first + initial_window
    num_x = consensus[l:r].count("X")
    begin = None
    output = {}

    while r <= last:
        ratio = num_x / (r-l)
        if ratio > limit:
            if begin is None:
                begin = l
            if r == len(consensus):
                break
            num_x += consensus[r] == "X"
            r += 1
            continue
        else:
            if begin is not None:
                a, b = first_last_X(consensus[begin:last])
                a, b = a + begin, b + begin
                if b not in output:
                    output[b] = (a,b)
                # if a is not None and b is not None:
                # output.append((a + begin, b + begin))
                # begin = None
            l = r
            r = l + initial_window
            begin = None
            num_x = consensus[l:r].count("X")
    if begin is not None:
        a, b = first_last_X(consensus[begin:last])
        a, b = a + begin, b + begin
        if a is not None and b is not None:
            if b not in output:
                output[b] = (a, b)
    return [*output.values()]

    # for i in range(first, last, window):
    #     ratio = consensus[i:i + window].count("X") / window
    #     if ratio >= limit:
    #         if begin is None:
    #             begin = i
    #     else:
    #         if begin is not None:
    #             # the new first is the old window end
    #             a, b = first_last_X(consensus[begin:i])
    #             output.append((a + begin, b + begin))
    #             begin = None
    # if begin is not None:
    #     a, b = first_last_X(consensus[begin:last])
    #     if a is not None and b is not None:
    #         output.append((a + begin, b + begin))

    # return output

    # while (first + window) <= last:
    #     ratio = num_x / total
    #     if ratio > limit:
    #         in_region = True
    #         num_x += consensus[window] == "X"
    #         total += 1
    #         window += 1
    #     else:
    #         if in_region:
    #             output.append(bd.find_index_pair(consensus[first:first+window], "X"))


def log_excised_consensus(path: Path, consensus_threshold=0.65, excise_threshold=0.40, dupes=False, log_path=None, debug=False, verbose=False, cut=False):
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

    sub_folder = path.parent
    folder = sub_folder.parent.parent
    if not log_path:
        log_path = Path(path.parent.parent.parent, "excise.tsv")
    gene = Path(path).name.split(".")[0]

    def load_dupes():
        rocks_db_path = Path(folder, "rocksdb", "sequences", "nt")
        if rocks_db_path.exists():

            rocksdb_db = wrap_rocks.RocksDB(str(rocks_db_path))
            prepare_dupe_counts = json.decode(
                rocksdb_db.get("getall:gene_dupes"), type=dict[str, dict[str, int]]
            )
            reporter_dupe_counts = json.decode(
                rocksdb_db.get("getall:reporter_dupes"), type=dict[str, dict[str, list]]
            )
        else:
            err = f"cannot find dupe databases for {folder}"
            raise FileNotFoundError(err)
        return prepare_dupe_counts, reporter_dupe_counts

    def bundle_seqs_and_dupes(sequences: list, prepare_dupe_counts, reporter_dupe_counts):
        output = []
        for header, seq in sequences:
            node = header.split("|")[3]
            dupes = prepare_dupe_counts.get(node, 1) + sum(
                prepare_dupe_counts.get(node, 1) for node in reporter_dupe_counts.get(node, []))
            output.append((seq, dupes))
        return output

    if not dupes:
        prepare_dupes, reporter_dupes = load_dupes()
        sequences = [(header, seq) for header,seq in Fastx(str(path)) if header[-1] != "."]
        sequences = bundle_seqs_and_dupes(sequences, prepare_dupes, reporter_dupes)
        consensus_func = bd.dumb_consensus_dupe
    else:
        sequences = [x[1] for x in Fastx(str(path)) if x[0][-1] != "."]
        consensus_func = bd.dumb_consensus
    consensus_seq = consensus_func(sequences, consensus_threshold)
    consensus_seq = bd.convert_consensus(sequences, consensus_seq)
    if verbose:
        print(f"{gene}\n{consensus_seq}\n")
    bad_regions = check_bad_regions(consensus_seq, excise_threshold)
    if bad_regions:
        if len(bad_regions) == 1:
            a, b = bad_regions[0]
            if b - a != len(consensus_seq):
                with open(log_path, "a") as f:
                    f.write(f"{gene}\t{a}:{b}")
                    if debug:
                        f.write("\t"+consensus_seq)
                    f.write("\n")
        else:
            with open(log_path, "a") as f:
                f.write(f"{gene}\t")
                for region in bad_regions:
                    a, b = region
                    f.write(f"{a}:{b}\t")
                    if debug:
                        f.write(consensus_seq)
                f.write("\n")

    # TODO:  if cut==True, cut and output here

    else:
        if debug:
            with open(log_path, "a") as f:
                f.write(f"{gene}\tN/A\t{consensus_seq}\n")


def main(args):
    if not (0 < args.consensus < 1.0):
        if 0 < args.consensus <= 100:
            args.consensus = args.consensus/100
        else:
            raise ValueError("Cannot convert consensus to a percent. Use a decimal or a whole number between 0 and 100")

    if not (0 < args.excise < 1.0):
        if 0 < args.excise <= 100:
            args.excise = args.excise/100
        else:
            raise ValueError("Cannot convert excise to a percent. Use a decimal or a whole number between 0 and 100")

    for folder in args.INPUT:
        aa_folder = Path(folder, "trimmed", "aa")
        fastas = [Path(aa_folder, fasta) for fasta in os.listdir(aa_folder) if ".fa" in fasta]
        log_path = Path(folder, "excise_indices.tsv")
        with open(log_path, "w") as f:
            f.write(f"Gene\tCut-Indices\n")
        for fasta in fastas:
            log_excised_consensus(fasta, args.consensus, args.excise, dupes=args.dupes, log_path=log_path, debug=args.debug,
                                  verbose=args.verbose, cut=args.cut)
    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
