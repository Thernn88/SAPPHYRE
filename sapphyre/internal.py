import os
from multiprocessing.pool import Pool
from shutil import rmtree

from pathlib import Path
from msgspec import Struct, json

import phymmr_tools as bd
import wrap_rocks
from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta, write2Line2Fasta

ALLOWED_EXTENSIONS = (".fa", ".fas", ".fasta", ".fa", ".gz", ".fq", ".fastq")

def folder_check(path: Path, debug: bool) -> None:
    """Create subfolders 'aa' and 'nt' to given path."""
    aa_folder = Path(path, "aa")
    nt_folder = Path(path, "nt")

    rmtree(aa_folder, ignore_errors=True)
    rmtree(nt_folder, ignore_errors=True)

    aa_folder.mkdir(parents=True, exist_ok=True)
    nt_folder.mkdir(parents=True, exist_ok=True)

    if debug:
        logs_folder = Path(path, "logs")
        logs_folder.mkdir(parents=True, exist_ok=True)


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


def load_dupes(folder):
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


def aa_internal(
    gene: str,
    consensus_threshold,
    distance_threshold,
    dupes,
    prepare_dupes,
    reporter_dupes,
):
    passing = {}
    failing = {}
    candidates, references = [], []
    for header, seq in parseFasta(gene):
        if header[-1] != ".":
            candidates.append((header, seq))
        else:
            references.append((header, seq))
    sequences = [tup[1] for tup in candidates]
    if dupes:
        consensus_func = bd.dumb_consensus_dupe
        sequences = bundle_seqs_and_dupes(sequences, prepare_dupes, reporter_dupes)
    else:
        consensus_func = bd.dumb_consensus
    consensus = consensus_func(sequences, consensus_threshold)
    for candidate in candidates:
        distance = bd.constrained_distance(consensus, candidate[1]) / len(candidate[1])
        if distance >= distance_threshold:
            failing[candidate[0]] = candidate[1]
        else:
            passing[candidate[0]] = candidate[1]
    return passing, failing, references


def mirror_nt(input_path, output_path, failing, gene):
    output_path = Path(output_path, gene)
    input_path = Path(input_path, gene)
    if not os.path.exists(input_path):
        return
    with open(output_path, "w") as f:
        for header, seq in parseFasta(input_path):
            if header in failing:
                continue
            f.write(f">{header}\n{seq}\n")


def run_internal(
    gene: str,
    nt_input: str,
    output_path,
    nt_output_path,
    consensus_threshold,
    distance_threshold,
    dupes,
    prepare_dupes,
    reporter_dupes,
):
    passing, failing, references = aa_internal(
        gene,
        consensus_threshold,
        distance_threshold,
        dupes,
        prepare_dupes,
        reporter_dupes,
    )
    aa_output = Path(output_path, "aa", gene.name)
    passing_lines = [(head, seq) for head, seq in passing.items()]
    writeFasta(aa_output, passing_lines)
    mirror_nt(nt_input, nt_output_path,passing, aa_output.name.replace(".aa.", ".nt."))




#             distance = constrained_distance(consensus, candidate.raw) / len(
#                 candidate.sequence
#             )
#             if distance >= internal_kick_threshold:
#                 candidate.grade = "Internal Fail"
#                 failing.append(candidate)
#                 passing[i] = None
#         passing = [x for x in passing if x is not None]
#     except KeyError:
#         print(f"key error for {gene}, skipping consensus")


def do_folder(folder, args):
    aa_input = Path(folder, "outlier", "aa")
    nt_input = Path(folder, "outlier", "nt")
    prepare_dupe_counts, reporter_dupe_counts = load_dupes(folder)
    file_inputs = [
        gene
        for gene in aa_input.iterdir()
        if ".aa" in gene.suffixes and gene.suffix in ALLOWED_EXTENSIONS
    ]
    output_path = Path(folder, "internal")
    nt_output_path = os.path.join(output_path, "nt")
    folder_check(output_path, False)
    file_inputs.sort(key=lambda x: x.stat().st_size, reverse=True)
    arguments = []

    if args.consensus_threshold > 100 or args.consensus_threshold <= 0:
        raise ValueError("cannot express given consensus threshold as a percent")
    if args.consensus_threshold > 1:
        args.consensus_threshold = args.consensus_thesold / 100
    if args.distance_threshold > 100 or args.distance_threshold <= 0:
        raise ValueError("cannot express given distance threshold as a percent")
    if args.distance_threshold > 1:
        args.distance_threshold = args.distance_thesold / 100

    for gene in file_inputs:
        gene_raw = gene.stem.split(".")[0]
        if args.dupes:
            prepare_dupes = (prepare_dupe_counts.get(gene_raw, {}),)
            reporter_dupes = (reporter_dupe_counts.get(gene_raw, {}),)
        else:
            prepare_dupes, reporter_dupes = None, None
        arguments.append(
            (
                gene,
                nt_input,
                output_path,
                nt_output_path,
                args.consensus_threshold,
                args.distance_threshold,
                args.dupes,
                prepare_dupes,
                reporter_dupes,
            ),
        )
        with Pool(args.processes) as pool:
            pool.starmap(run_internal, arguments, chunksize=1)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for folder in args.INPUT:
        do_folder(Path(folder), args)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre internal"
    raise RuntimeError(
        MSG,
    )
