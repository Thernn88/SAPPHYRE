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

class Record(Struct):
    id: str
    seq: str

    def __str__(self):
        return f">{self.id}\n{self.seq}\n"

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
    for rec in sequences:
        node = rec.id.split("|")[3]
        dupes = prepare_dupe_counts.get(node, 1) + sum(
            prepare_dupe_counts.get(node, 1)
            for node in reporter_dupe_counts.get(node, [])
        )
        output.append((rec.seq, dupes))
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
    failing = set()
    candidates, references = [], []
    for header, seq in parseFasta(gene):
        if header[-1] != ".":
            candidates.append(Record(header, seq))
        else:
            references.append(Record(header, seq))
    if not candidates:  # if no candidates are found, report to user
        # print(f"{gene}: No Candidate Sequences found in file. Returning.")
        return [], {}, []
    if dupes:
        consensus_func = bd.dumb_consensus_dupe
        sequences = bundle_seqs_and_dupes(candidates, prepare_dupes, reporter_dupes)
    else:
        consensus_func = bd.dumb_consensus
        sequences = [rec.seq for rec in candidates]

    consensus = consensus_func(sequences, consensus_threshold)
    for i, candidate in enumerate(candidates):
        distance = bd.constrained_distance(consensus, candidate.seq) / len(candidate.seq)
        if distance >= distance_threshold:
            # failing[candidate[0]] = candidate[1]
            failing.add(candidate.id)
            candidates[i] = None
    candidates = [cand for cand in candidates if cand != None]
    return candidates, failing, references


def mirror_nt(input_path, output_path, failing, gene):
    output_path = Path(output_path, gene)
    input_path = Path(input_path, gene)
    if not os.path.exists(input_path):
        return
    records = ((header, seq) for header, seq in parseFasta(input_path) if header not in failing)
    with open(output_path, "w") as f:
        f.writelines((f">{header}\n{seq}\n" for header, seq in records))


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
    if not passing:  # if no eligible candidates, don't create the output filegi
        return
    aa_output = Path(output_path, "aa", gene.name)
    with open(aa_output,"w") as f:
        f.writelines((str(rec) for rec in passing))
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


# def do_folder(folder, args):


def main(args):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if args.internal_consensus_threshold > 100 or args.internal_consensus_threshold <= 0:
        raise ValueError("cannot express given consensus threshold as a percent")
    if args.internal_consensus_threshold > 1:
        args.internal_consensus_threshold = args.consensus_thesold / 100
    if args.internal_distance_threshold > 100 or args.internal_distance_threshold <= 0:
        raise ValueError("cannot express given distance threshold as a percent")
    if args.internal_distance_threshold > 1:
        args.internal_distance_threshold = args.distance_thesold / 100

    with Pool(args.processes) as pool:
        folder = args.INPUT
        aa_input = Path(folder, "outlier", args.sub_directory, "aa")
        nt_input = Path(folder, "outlier", args.sub_directory, "nt")
        prepare_dupe_counts, reporter_dupe_counts = load_dupes(folder)
        file_inputs = [
            gene
            for gene in aa_input.iterdir()
            if ".aa" in gene.suffixes and gene.suffix in ALLOWED_EXTENSIONS
        ]
        output_path = Path(folder, "outlier", "internal")
        nt_output_path = os.path.join(output_path, "nt")
        folder_check(output_path, False)
        file_inputs.sort(key=lambda x: x.stat().st_size, reverse=True)
        arguments = []


        for gene in file_inputs:
            gene_raw = gene.stem.split(".")[0]
            if args.dupes:
                prepare_dupes = prepare_dupe_counts.get(gene_raw, {})
                reporter_dupes = reporter_dupe_counts.get(gene_raw, {})
            else:
                prepare_dupes, reporter_dupes = None, None
            arguments.append(
                (
                    gene,
                    nt_input,
                    output_path,
                    nt_output_path,
                    args.internal_consensus_threshold,
                    args.internal_distance_threshold,
                    args.dupes,
                    prepare_dupes,
                    reporter_dupes,
                ),
            )
        pool.starmap(run_internal, arguments, chunksize=1)
    printv(f"Done! Took {timer.differential():.2f}s", args.verbose)
    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre internal"
    raise RuntimeError(
        MSG,
    )
