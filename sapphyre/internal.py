from multiprocessing.pool import Pool
from os import path
from pathlib import Path

from msgspec import Struct, json
from phymmr_tools import (
    constrained_distance,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
)
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

ALLOWED_EXTENSIONS = (".fa", ".fas", ".fasta", ".fa", ".gz", ".fq", ".fastq")


class Record(Struct):
    id: str
    seq: str

    def __str__(self):
        return f">{self.id}\n{self.seq}\n"

    def get_pair(self):
        return (self.id, self.seq)


def folder_check(taxa_path: Path, debug: bool) -> None:
    """Create subfolders 'aa' and 'nt' to given path."""
    aa_folder = Path(taxa_path, "aa")
    nt_folder = Path(taxa_path, "nt")

    aa_folder.mkdir(parents=True, exist_ok=True)
    nt_folder.mkdir(parents=True, exist_ok=True)

    if debug:
        logs_folder = Path(taxa_path, "logs")
        logs_folder.mkdir(parents=True, exist_ok=True)


def bundle_seqs_and_dupes(sequences: list, prepare_dupe_counts, reporter_dupe_counts):
    """
    Pairs each record object with its dupe count from prepare and reporter databases.
    Given dupe count dictionaries and a list of Record objects, makes tuples of the records
    and their dupe counts. Returns the tuples in a list.
    """
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
    """
    Load dupe counts from prepare and reporter runs.
    """
    rocks_db_path = Path(folder, "rocksdb", "sequences", "nt")
    if rocks_db_path.exists():
        rocksdb_db = RocksDB(str(rocks_db_path))
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


def get_data_path(gene: Path) -> list:
    """
    Turn the given path of a collapsed file into an equivalent path in the excise folder.
    If the excise path exists, read it and replace any candidate sequences with the version
    from the excise file.
    """
    if "hmmfilter" in str(gene):
        excise_path = str(gene).replace("/hmmfilter/", "/hmmfilter/")
        if path.exists(excise_path):
            return [Record(header, seq) for header, seq in parseFasta(excise_path)]

    return [Record(header, seq) for header, seq in parseFasta(str(gene))]


def has_candidates(records: list) -> bool:
    """
    Checks if there are any candidates in a list of sequences. If a candidate header is
    found, returns True. Otherwise, returns False.
    """
    for rec in records:
        if rec.id[-1] != ".":
            return True
    return False


def aa_internal(
    gene: str,
    consensus_threshold,
    distance_threshold,
    no_dupes,
    prepare_dupes,
    reporter_dupes,
):
    """
    Compares the candidate sequences to their consensus sequence. If the hamming distance between
    a candidate and the consesnus seq is too high, kicks the candidate.
    Returns a tuple containing a list of passing candidates, a set of failing seqs, and a list of
    reference seqs.

    Gene is the path to the input file.

    Consensus threshold is a float between 0.0 and 1.0 that represents the minimum
    ratio a character must reach to become the consensus character at that location.

    Distance threshold is a float between 0.0 and 1.0 that represesnts the maximum allowable
    ratio between hamming distance and the length of candidate data region.

    No dupes is a flag that enables or disables use of dupe counts from prepare and reporter.

    Prepare Dupes and Reporter Dupes are dictionaries that map node names to their duplicate counts.
    """
    failing = set()

    # load sequences from file
    raws = get_data_path(gene)

    # split sequences into references and candidates
    candidates, references = [], []
    for record in raws:
        if record.id[-1] != ".":
            candidates.append(record)
        else:
            references.append(record)

    # if no candidates found, return from function
    if not candidates:
        return [], set(), references

    # make consensus sequences, use dupe counts if available
    if no_dupes:
        consensus_func = dumb_consensus
        sequences = [rec.seq for rec in candidates]
    else:
        consensus_func = dumb_consensus_dupe
        sequences = bundle_seqs_and_dupes(candidates, prepare_dupes, reporter_dupes)

    consensus = consensus_func(sequences, consensus_threshold)

    # compare the hamming distance between each candidate and the appropriate slice of the consensus
    # divides by length of the candidate's data region to adjust for length
    for i, candidate in enumerate(candidates):
        start, stop = find_index_pair(candidate.seq, "-")
        distance = constrained_distance(consensus, candidate.seq) / (stop - start)
        # if the ratio of hamming_distance to length is too high, kick the candidate
        if distance >= distance_threshold:
            failing.add(candidate.id)
            candidates[i] = None
    # failed records are represented by None values, so filter them from candidates
    candidates = [cand for cand in candidates if cand is not None]
    return candidates, failing, references


def mirror_nt(input_path, output_path, failing, gene, compression):
    """
    Mirrors aa kicks in the appropriate nt file.
    """
    output_path = Path(output_path, gene)
    input_path = Path(input_path, gene)
    if not path.exists(input_path):
        return

    records = get_data_path(input_path)
    records = [rec for rec in records if rec.id not in failing]
    # if no candidates are left after the internal check, don't make an output file
    if not has_candidates(records):
        return
    writeFasta(
        str(output_path), [rec.get_pair() for rec in records], compress=compression
    )


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
    decompress,
):
    """
    Given a gene, reads the aa file and compares the candidates to their consensus sequence.
    If a candidate has too many locations that disagree with the consensus, kicks the sequence.
    """
    compression = not decompress
    passing, failing, references = aa_internal(
        gene,
        consensus_threshold,
        distance_threshold,
        dupes,
        prepare_dupes,
        reporter_dupes,
    )
    if not passing:  # if no eligible candidates, don't create the output file
        return
    aa_output = Path(output_path, "aa", gene.name)
    writeFasta(
        str(aa_output),
        [rec.get_pair() for rec in references + passing],
        compress=compression,
    )
    mirror_nt(
        nt_input,
        nt_output_path,
        failing,
        aa_output.name.replace(".aa.", ".nt."),
        compression,
    )


def main(args, after_collapser, from_folder):
    timer = TimeKeeper(KeeperMode.DIRECT)
    if (
        args.internal_consensus_threshold > 100
        or args.internal_consensus_threshold <= 0
    ):
        raise ValueError("cannot express given consensus threshold as a percent")
    if args.internal_consensus_threshold > 1:
        args.internal_consensus_threshold = args.consensus_thesold / 100
    if args.internal_distance_threshold > 100 or args.internal_distance_threshold <= 0:
        raise ValueError("cannot express given distance threshold as a percent")
    if args.internal_distance_threshold > 1:
        args.internal_distance_threshold = args.distance_thesold / 100

    with Pool(args.processes) as pool:
        folder = args.INPUT
        if after_collapser:
            aa_input = Path(folder, "outlier", "hmmfilter", "aa")
            nt_input = Path(folder, "outlier", "hmmfilter", "nt")
        else:
            aa_input = Path(folder, "outlier", from_folder, "aa")
            nt_input = Path(folder, "outlier", from_folder, "nt")
        if not args.no_dupes:
            prepare_dupe_counts, reporter_dupe_counts = load_dupes(folder)
        file_inputs = [
            gene
            for gene in aa_input.iterdir()
            if ".aa" in gene.suffixes and gene.suffix in ALLOWED_EXTENSIONS
        ]

        output_path = Path(folder, "outlier", "internal")
        nt_output_path = path.join(output_path, "nt")
        folder_check(output_path, False)
        file_inputs.sort(key=lambda x: x.stat().st_size, reverse=True)
        arguments = []

        if not (0 < args.internal_consensus_threshold < 1.0):
            if 0 < args.internal_consensus_threshold <= 100:
                args.internal_consensus_threshold = args.internal_consensus_threshold / 100
            else:
                raise ValueError(
                    "Cannot convert internal consensus threshold to a percent. Use a decimal or a whole number between 0 and 100"
                )
        if not (0 < args.internal_distance_threshold < 1.0):
            if 0 < args.internal_distance_threshold <= 100:
                args.internal_distance_threshold = args.internal_distance_threshold / 100
            else:
                raise ValueError(
                    "Cannot convert internal distance threshold to a percent. Use a decimal or a whole number between 0 and 100"
                )
            
        for gene in file_inputs:
            gene_raw = gene.stem.split(".")[0]
            if not args.no_dupes:
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
                    args.internal_distance_threshold,
                    args.internal_distance_threshold,
                    args.no_dupes,
                    prepare_dupes,
                    reporter_dupes,
                    args.uncompress_intermediates,
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
