from multiprocessing.pool import Pool
from os import path
from pathlib import Path

from msgspec import Struct, json
from sapphyre_tools import (
    dumb_consensus,
    dumb_consensus_dupe,
    consensus_distance,
    find_index_pair,
)
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


def folder_check(taxa_path: Path) -> Path:
    """Create subfolders 'aa' and 'nt' to given path."""
    aa_folder = Path(taxa_path, "aa")
    nt_folder = Path(taxa_path, "nt")

    aa_folder.mkdir(parents=True, exist_ok=True)
    nt_folder.mkdir(parents=True, exist_ok=True)

    # if debug:
    logs_folder = Path(taxa_path, "logs")
    logs_folder.mkdir(parents=True, exist_ok=True)
    aa_log_path = Path(logs_folder, "aa_fail.log")
    with open(aa_log_path, "w") as log:
        log.write("id\tdistance\tthreshold\n")
    nt_log_path = Path(logs_folder, "nt_fail.log")
    with open(nt_log_path, "w") as log:
        log.write("id\tdistance\tthreshold\n")
    return aa_log_path, nt_log_path
    # else:
    #     return None

def bundle_seqs_and_dupes(sequences: list):
    """
    Pairs each record object with its dupe count from prepare and reporter databases.
    Given dupe count dictionaries and a list of Record objects, makes tuples of the records
    and their dupe counts. Returns the tuples in a list.
    """
    output = []
    for rec in sequences:
        dupe_count = int(rec.id.split("|")[5])
        output.append((rec.seq, dupe_count))
    return output

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


def do_internal(
    gene: str,
    consensus_threshold,
    distance_threshold,
    no_dupes,
    minimum_depth,
    minimum_length,
    minimum_overlap,
    existing_fails = set(),
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
    failing_dict = {}
    # load sequences from file
    raws = get_data_path(gene)

    # split sequences into references and candidates
    candidates, references = [], []
    for record in raws:
        if record.id in existing_fails:
            continue

        if record.id[-1] != ".":
            candidates.append(record)
        else:
            references.append(record)

    # if no candidates found, return from function
    if not candidates:
        return candidates, failing, references, failing_dict, ""

    # make consensus sequences, use dupe counts if available
    if no_dupes:
        consensus_func = dumb_consensus
        sequences = [rec.seq for rec in candidates]
    else:
        consensus_func = dumb_consensus_dupe
        sequences = bundle_seqs_and_dupes(candidates)

    consensus = consensus_func(sequences, consensus_threshold, minimum_depth)

    # compare the hamming distance between each candidate and the appropriate slice of the consensus
    # divides by length of the candidate's data region to adjust for length
    for i, candidate in enumerate(candidates):
        # start, stop = find_index_pair(candidate.seq, "-")
        total_distance, total_length = consensus_distance(consensus, candidate.seq, minimum_length, minimum_overlap)
        if not total_length:
            continue
        distance = total_distance / total_length
        # if the ratio of hamming_distance to length is too high, kick the candidate
        if distance >= distance_threshold:
            failing.add(candidate.id)
            start, stop = find_index_pair(candidate.seq, '-')
            failing_dict[candidate.id] = (candidate, start, stop, distance, distance_threshold)
            candidates[i] = None
    # failed records are represented by None values, so filter them from candidates
    candidates = [cand for cand in candidates if cand is not None]
    return candidates, failing, references, failing_dict, consensus


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
    consensus_threshold_nt,
    distance_threshold_nt,
    no_dupes,
    compression,
    aa_log_path,
    nt_log_path,
    minimum_depth,
    minimum_length,
    minimum_overlap,
):
    """
    Given a gene, reads the aa file and compares the candidates to their consensus sequence.
    If a candidate has too many locations that disagree with the consensus, kicks the sequence.
    """
    aa_passing, aa_failing, aa_references, aa_fail_dict, aa_consensus = do_internal(
        gene,
        consensus_threshold,
        distance_threshold,
        no_dupes,
        minimum_depth,
        minimum_length,
        minimum_overlap,
    )

    nt_passing, nt_failing, nt_references, nt_fail_dict, nt_consensus = do_internal(
        Path(nt_input, gene .name.replace(".aa.", ".nt.")),
        consensus_threshold_nt,
        distance_threshold_nt,
        no_dupes,
        minimum_depth,
        minimum_length,
        minimum_overlap,
        existing_fails = aa_failing
    )

    aa_passing = [rec for rec in aa_passing if rec.id not in nt_failing]
    # log_path = Path(log_folder_path, "fail.log")
    if aa_fail_dict:
        with open(aa_log_path, 'a+') as f:
            f.write(f"{gene.name}\n")
            for tup in aa_fail_dict.values():
                candidate, start, stop, distance, threshold = tup
                f.write(f'{candidate.id}\tstart:{start}\tstop:{stop}distance:{distance}\tthreshold:{threshold}\n')
                f.write(f'cand seq:{candidate.seq}\ncons seq:{aa_consensus}\n')
    if nt_fail_dict:
        with open(nt_log_path, 'a+') as f:
            f.write(f"{gene.name}\n")
            for tup in nt_fail_dict.values():
                candidate, start, stop, distance, threshold = tup
                f.write(f'{candidate.id}\tstart:{start}\tstop:{stop}\tdistance:{distance}\tthreshold:{threshold}\n')
                f.write(f'cand seq:{candidate.seq[start:stop]}\ncons seq:{nt_consensus[start:stop]}\n')
            
    if not aa_passing:  # if no eligible candidates, don't create the output file
        return
    
    
    aa_out = [rec.get_pair() for rec in aa_references + aa_passing]
    nt_out = [rec.get_pair() for rec in nt_references + nt_passing]

    aa_output = Path(output_path, "aa", gene.name)
    writeFasta(
        str(aa_output),
        aa_out,
        compress=compression,
    )
    nt_output = Path(nt_output_path, gene.name.replace(".aa.", ".nt."))
    writeFasta(
        str(nt_output),
        nt_out,
        compress=compression,
    )


def main(args, from_folder):
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
        args.internal_distance_threshold = args.internal_distance_threshold / 100
    if args.internal_consensus_threshold_nt > 100 or args.internal_consensus_threshold_nt <= 0:
        raise ValueError("cannot express given consensus threshold as a percent")
    if args.internal_consensus_threshold_nt > 1:
        args.internal_consensus_threshold_nt = args.consensus_thesold_nt / 100
    if args.internal_distance_threshold_nt > 100 or args.internal_distance_threshold_nt <= 0:
        raise ValueError("cannot express given distance threshold as a percent")
    if args.internal_distance_threshold_nt > 1:
        args.internal_distance_threshold_nt = args.internal_distance_threshold_nt / 100

    with Pool(args.processes) as pool:
        folder = args.INPUT

        aa_input = Path(folder, "outlier", from_folder, "aa")
        nt_input = Path(folder, "outlier", from_folder, "nt")

        file_inputs = [
            gene
            for gene in aa_input.iterdir()
            if ".aa" in gene.suffixes and gene.suffix in ALLOWED_EXTENSIONS
        ]

        output_path = Path(folder, "outlier", "internal")
        nt_output_path = path.join(output_path, "nt")
        aa_log_path, nt_log_path = folder_check(output_path)
        file_inputs.sort(key=lambda x: x.stat().st_size, reverse=True)
        arguments = []
            
        for gene in file_inputs:
            arguments.append(
                (
                    gene,
                    nt_input,
                    output_path,
                    nt_output_path,
                    args.internal_consensus_threshold,
                    args.internal_distance_threshold,
                    args.internal_consensus_threshold_nt,
                    args.internal_distance_threshold_nt,
                    args.no_dupes,
                    args.compress,
                    aa_log_path,
                    nt_log_path,
                    args.minimum_depth,
                    args.minimum_candidate_length,
                    args.minimum_candidate_overlap,
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
