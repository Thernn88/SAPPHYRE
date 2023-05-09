from __future__ import annotations
from math import ceil
import os
from collections import defaultdict, namedtuple
from multiprocessing.pool import Pool
from shutil import rmtree
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
import xxhash
from .utils import printv, gettempdir, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode

KMER_LEN = 15
KMER_PERCENT = 0.15
SUBCLUSTER_AT = 1000
CLUSTER_EVERY = 500  # Aim for x seqs per cluster
SAFEGUARD_BP = 15000
SINGLETON_THRESHOLD = 5
UPPER_SINGLETON_TRESHOLD = 20000
SINGLETONS_REQUIRED = (
    50  # Amount of singleton clusters required to continue checking threshold
)


def find_kmers(fasta: dict) -> dict[str, set]:
    """
    Returns the kmers of length KMER_LEN for each sequence in the fasta dict.

    Args:
        fasta (dict): A dictionary of header -> sequence
    Returns:
        dict[str, set]: A dictionary of header -> set of kmers
    """
    kmers = defaultdict(set)
    for header, sequence in fasta.items():
        for i in range(0, len(sequence) - KMER_LEN):
            kmer = sequence[i : i + KMER_LEN]
            if "*" not in kmer and "-" not in kmer:
                kmers[header].add(kmer)
    return kmers


ALIGN_FOLDER = "align"
AA_FOLDER = "aa"
ALN_FOLDER = "aln"


def get_start(sequence: str) -> int:
    """
    Returns the index of the first non-gap character in the sequence.

    Args:
        sequence (str): The sequence to check
    Returns:
        int: The index of the first non-gap character
    """
    for i, char in enumerate(sequence):
        if char != "-":
            return i
    return -1


def process_genefile(
    path: str,
) -> tuple[int, dict[str, str], dict[str, str], dict[str, list[str]], dict[str, str]]:
    """
    Returns the number of sequences, the sequences, the targets, the reinsertions and the
    trimmed header to full header.

    Reinsertions are dupes where the sequences are exactly the same. To save on computation we
    insert these duplicates into the alignment after it has been generated right next to their
    duplicate parents with the same alignment.

    Args:
        path (str): The path to the gene file
    Returns:
        tuple[int, dict[str, str], dict[str, str], dict[str, list[str]], dict[str, str]]:
        The number of sequences, the sequences, the targets, the reinsertions and the
        trimmed header to full header.
    """
    data = {}
    targets = {}
    reinsertions = defaultdict(list)

    seq_to_first_header = {}
    trimmed_header_to_full = {}
    for header, sequence in parseFasta(path):
        if not header.endswith("."):
            seq_hash = xxhash.xxh3_64(sequence).hexdigest()
            trimmed_header_to_full[header[:127]] = header
            if seq_hash not in seq_to_first_header:
                seq_to_first_header[seq_hash] = header
            else:
                reinsertions[seq_to_first_header[seq_hash]].append(header)

                continue
            data[header] = sequence.replace("-","") # Delete gaps from previous alignment
        else:
            targets[header.split("|")[2]] = header

    return len(data), data, targets, reinsertions, trimmed_header_to_full


def delete_empty_cols(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """
    Deletes columns that are empty in all sequences.

    Args:
        records (list[tuple[str, str]]): A list of tuples of header, sequence
    Returns:
        list[tuple[str, str]]: A list of tuples of header, sequence with removed columns
    """
    cols_to_keep = set()
    output = []
    for _, sequence in records:
        for x, let in enumerate(sequence):
            if x not in cols_to_keep:
                if let != "-":
                    cols_to_keep.add(x)
    for header, sequence in records:
        sequence = [let for x, let in enumerate(sequence) if x in cols_to_keep]
        output.append((header, "".join(sequence)))

    return output


def compare_cluster(
    cluster_a: list[str], cluster_b: list[str], kmers: dict[str, set]
) -> float:
    """
    Calculates the average kmer similarity between children of two clusters.

    Args:
        cluster_a (list[str]): A list of headers in cluster a
        cluster_b (list[str]): A list of headers in cluster b
        kmers (dict[str, set]): A dictionary of header -> set of kmers
    Returns:
        float: The average similarity between children of two clusters
    """
    average_identity = []
    for head_a in cluster_a:
        for head_b in cluster_b:
            if head_a == head_b:
                continue
            akmers = kmers[head_a]
            bkmers = kmers[head_b]

            similar = akmers.intersection(bkmers)
            average_identity.append(len(similar) / min(len(akmers), len(bkmers)))

    return sum(average_identity) / len(average_identity)


def generate_clusters(data: dict[str, str]) -> list[list[str]]:
    """
    Generates clusters from a dictionary of header -> sequence.

    Clusters are based on the average identity of kmers between sequences.

    Args:
        data (dict[str, str]): A dictionary of header -> sequence
    Returns:
        list[list[str]]: A list of each clusters' headers
    """
    cluster_children = {header: [header] for header in data}
    kmers = find_kmers(data)
    gene_headers = list(kmers.keys())

    # Add a dictionary to store child sets for each primary set
    child_sets = {header: set() for header in data}

    for iteration in range(2):
        processed_headers = set()
        merge_occured = True
        while merge_occured:
            merge_occured = False
            for i in range(len(gene_headers) - 1, -1, -1):  # reverse iteration
                master_header = gene_headers[i]
                master = kmers[master_header]
                if master:
                    for candidate_header in gene_headers[:i]:
                        if candidate_header in processed_headers:
                            continue

                        candidate = kmers[candidate_header]
                        if candidate:
                            # Check similarity against both parent set and child sets
                            matched = False
                            for header_to_check in [master_header] + list(
                                child_sets[master_header]
                            ):
                                set_to_check = kmers[header_to_check]
                                if set_to_check is None:
                                    continue

                                similar = set_to_check.intersection(candidate)

                                if len(similar) != 0:
                                    if (
                                        len(similar)
                                        / min(len(set_to_check), len(candidate))
                                        >= KMER_PERCENT
                                    ):
                                        if (
                                            iteration == 0
                                            or iteration == 1
                                            and len(cluster_children[master_header])
                                            != 1
                                        ):
                                            # Add the candidate set as a child set of the primary set
                                            child_sets[master_header].add(
                                                candidate_header
                                            )
                                            cluster_children[master_header].extend(
                                                cluster_children[candidate_header]
                                            )

                                            # Remove candidate
                                            cluster_children[candidate_header] = None
                                            kmers[candidate_header] = None
                                            processed_headers.add(candidate_header)

                                            matched = True
                                            break

                            if matched:
                                merge_occured = True

    return cluster_children.values()


def seperate_into_clusters(
    cluster_children: list[list[str]], parent_tmpdir: str, data: dict[str, str]
) -> list[list[str]]:
    """
    Seperates sequence records into clusters and subclusters
    if they are larger than SUBCLUSTER_AT.

    Args:
        cluster_children (list[list[str]]): A list of each clusters' headers
        parent_tmpdir (str): The parent temporary directory
        data (dict[str, str]): A dictionary of header -> sequence
    Returns:
        list[list[str]]: A list of each clusters' headers
    """
    clusters = []
    for this_cluster in cluster_children:
        if this_cluster is not None:
            if len(this_cluster) > SUBCLUSTER_AT:
                clusters_to_create = ceil(len(this_cluster) / CLUSTER_EVERY)
                with NamedTemporaryFile(
                    mode="w+", dir=parent_tmpdir, suffix=".fa"
                ) as this_tmp:
                    writeFasta(
                        this_tmp.name,
                        [(header, data[header]) for header in this_cluster],
                    )
                    with NamedTemporaryFile("r", dir=gettempdir()) as this_out:
                        subprocess.run(
                            f"SigClust/SigClust -k 8 -c {clusters_to_create} {this_tmp.name} > {this_out.name}",
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                            check=True,
                            shell=True,
                        )
                        sig_out = this_out.read()
                    sub_clusters = defaultdict(list)
                    for line in sig_out.split("\n"):
                        line = line.strip()
                        if line:
                            seq_index, clust_index = line.strip().split(",")
                            sub_clusters[clust_index].append(
                                this_cluster[int(seq_index)]
                            )
                    clusters.extend(list(sub_clusters.values()))
                    continue
            clusters.append(this_cluster)

    return clusters


def generate_tmp_aln(
    aln_file: str,
    targets: dict[str, str],
    tmp: NamedTemporaryFile,
    debug: float,
    this_intermediates: str,
) -> None:
    """
    Grabs target reference sequences and removes empty columns from the alignment.

    Args:
        aln_file (str): The path to the alignment file
        targets (dict[str, str]): A dictionary of target -> header
        tmp (NamedTemporaryFile): The temporary file to write to
        debug (float): Whether to write debug files
        this_intermediates (str): The path to the intermediates debug folder
    Returns:
        None
    """
    sequences = []
    for header, sequence in parseFasta(aln_file, True):
        if header in targets:
            sequences.append((targets[header], sequence))

    empty_columns = None
    for header, sequence in sequences:
        if empty_columns is None:
            empty_columns = [True] * len(sequence)

        for col, let in enumerate(sequence):
            if let != "-":
                empty_columns[col] = False

    to_write = []
    for header, sequence in sequences:
        to_write.append(
            (
                header,
                "".join(
                    [let for col, let in enumerate(sequence) if not empty_columns[col]]
                ),
            )
        )

    writeFasta(tmp.name, to_write)
    if debug:
        writeFasta(os.path.join(this_intermediates, "references.fa"), to_write)
    tmp.flush()


CmdArgs = namedtuple(
    "CmdArgs",
    [
        "string",
        "gene_file",
        "result_file",
        "gene",
        "verbose",
        "compress",
        "aln_path",
        "debug",
        "only_singletons",
        "add_fragments"
    ],
)


def run_command(args: CmdArgs) -> None:
    keeper = TimeKeeper(KeeperMode.DIRECT)
    debug = args.debug
    printv(f"Doing: {args.gene} ", args.verbose, 2)

    temp_dir = gettempdir()

    aligned_ingredients = []

    this_intermediates = os.path.join("intermediates", args.gene)
    if debug:
        if not os.path.exists(this_intermediates):
            os.mkdir(this_intermediates)

    aln_file = os.path.join(args.aln_path, args.gene + ".aln.fa")
    to_merge = []

    with TemporaryDirectory(dir=temp_dir) as parent_tmpdir, TemporaryDirectory(
        dir=parent_tmpdir
    ) as raw_files_tmp, TemporaryDirectory(dir=parent_tmpdir) as aligned_files_tmp:
        printv(
            f"Cleaning NT file. Elapsed time: {keeper.differential():.2f}",
            args.verbose,
            3,
        )  # Debug

        (
            seq_count,
            data,
            targets,
            reinsertions,
            trimmed_header_to_full,
        ) = process_genefile(args.gene_file)

        only_singletons = args.only_singletons
        if len(data) > UPPER_SINGLETON_TRESHOLD:
            only_singletons = True

        cluster_time = 0
        align_time = 0
        merge_time = 0

        if seq_count == 1:
            printv(
                f"Outputting singleton alignment. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            aligned_file = os.path.join(aligned_files_tmp, f"{args.gene}_cluster_1_0")
            sequences = list(data.items())
            writeFasta(aligned_file, sequences)
            if args.debug:
                writeFasta(os.path.join(this_intermediates, aligned_file), sequences)
            aligned_ingredients.append((aligned_file, len(sequences)))
        else:
            printv(
                f"Generating Cluster. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug

            if args.add_fragments:
                clusters = [[header] for header in data]
            else:
                cluster_children = generate_clusters(data)

                clusters = seperate_into_clusters(cluster_children, parent_tmpdir, data)

            cluster_time = keeper.differential()
            printv(
                f"Found {seq_count} sequences over {len(clusters)} clusters. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            printv(
                f"Aligning Clusters. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug

            items = os.listdir(raw_files_tmp)
            items.sort(key=lambda x: int(x.split("_cluster")[-1]))

            under_threshold = 0
            for cluster in clusters:
                if 1 < len(cluster) < SINGLETON_THRESHOLD:
                    under_threshold += 1

            singleton_allowed = under_threshold >= SINGLETONS_REQUIRED

            cluster_i = 0
            for cluster in clusters:
                cluster_i += 1
                cluster_seqs = [(header, data[header]) for header in cluster]
                cluster_length = len(cluster)

                if (
                    only_singletons
                    or cluster_length == 1
                    or cluster_length < SINGLETON_THRESHOLD
                    and singleton_allowed
                ):
                    printv(
                        f"Storing singleton cluster {cluster_i}. Elapsed time: {keeper.differential():.2f}",
                        args.verbose,
                        3,
                    )  # Debug
                    to_merge.extend(cluster_seqs)
                    continue

                printv(
                    f"Aligning cluster {cluster_i}. Elapsed time: {keeper.differential():.2f}",
                    args.verbose,
                    3,
                )  # Debug

                this_clus_align = f"{args.gene}_cluster_{cluster_i}_length_{len(cluster)}_index_aligned"
                aligned_cluster = os.path.join(aligned_files_tmp, this_clus_align)

                raw_cluster = os.path.join(
                    raw_files_tmp, f"{args.gene}_cluster{cluster_i}"
                )
                writeFasta(raw_cluster, cluster_seqs)
                if debug:
                    writeFasta(
                        os.path.join(this_intermediates, this_clus_align),
                        cluster_seqs,
                    )

                command = args.string.format(
                    in_file=raw_cluster, out_file=aligned_cluster
                )
                os.system(command)

                aligned_sequences = list(parseFasta(aligned_cluster, True))

                if debug:
                    printv(command, args.verbose, 3)
                    writeFasta(
                        os.path.join(
                            this_intermediates, os.path.split(aligned_cluster)[-1]
                        ),
                        aligned_sequences,
                    )

                aligned_ingredients.append((aligned_cluster, len(aligned_sequences)))

        align_time = keeper.differential() - cluster_time
        has_singleton_merge = len(to_merge) > 0
        if has_singleton_merge:
            merged_singleton_final = os.path.join(
                aligned_files_tmp, f"{args.gene}_cluster_merged_singletons"
            )
            writeFasta(merged_singleton_final, to_merge)
            if args.debug:
                writeFasta(
                    os.path.join(
                        this_intermediates, f"{args.gene}_cluster_merged_singletons"
                    ),
                    to_merge,
                )
        if aligned_ingredients or has_singleton_merge:
            aligned_ingredients = [
                i[0] for i in sorted(aligned_ingredients, key=lambda x: x[1])
            ]
            if has_singleton_merge:
                aligned_ingredients.insert(0, merged_singleton_final)
            printv(
                f"Merging Alignments. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            with NamedTemporaryFile(dir=parent_tmpdir, mode="w+") as tmp_aln:
                generate_tmp_aln(aln_file, targets, tmp_aln, debug, this_intermediates)

                prev_file = tmp_aln.name

                for i, file in enumerate(aligned_ingredients):
                    printv(
                        f"Merging alignment {i+1} of {len(aligned_ingredients)}.",
                        args.verbose,
                        3,
                    )
                    out_file = os.path.join(parent_tmpdir, f"part_{i}.fa")
                    if len(aligned_ingredients) - 1 == i:
                        out_file = args.result_file

                    if has_singleton_merge and i == 0:
                        os.system(
                            f"mafft --anysymbol --jtt 1 --quiet --addfragments {file} --thread 1 {prev_file} > {out_file}"
                        )
                        if args.debug:
                            printv(
                                f"mafft --anysymbol --jtt 1 --quiet --addfragments {file} --thread 1 {prev_file} > {out_file}",
                                args.verbose,
                                3,
                            )
                    else:
                        os.system(
                            f"clustalo --p1 {prev_file} --p2 {file} -o {out_file} --threads=1 --full --is-profile --force"
                        )
                        if debug:
                            printv(
                                f"clustalo --p1 {prev_file} --p2 {file} -o {out_file} --threads=1  --full --is-profile --force",
                                args.verbose,
                                3,
                            )
                    if debug:
                        writeFasta(
                            os.path.join(this_intermediates, f"part_{i}.fa"),
                            parseFasta(out_file, True),
                        )

                    prev_file = out_file

    merge_time = keeper.differential() - align_time - cluster_time

    # Reinsert and sort
    to_write = []
    references = []
    inserted = 0
    for header, sequence in parseFasta(args.result_file, True):
        if header.endswith("."):
            references.append((header, sequence))
        else:
            header = trimmed_header_to_full[header[:127]]
            if header in reinsertions:
                for insertion_header in reinsertions[header]:
                    inserted += 1
                    to_write.append((insertion_header, sequence))
            to_write.append((header, sequence))

    to_write.sort(key=lambda x: get_start(x[1]))

    writeFasta(args.result_file, references + to_write)

    printv(f"Done. Took {keeper.differential():.2f}", args.verbose, 3)  # Debug

    return args.gene, cluster_time, align_time, merge_time, keeper.differential()


def do_folder(folder, args):
    printv(f"Processing: {os.path.basename(folder)}", args.verbose, 0)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    align_path = os.path.join(folder, ALIGN_FOLDER)
    aa_path = os.path.join(folder, AA_FOLDER)
    if not os.path.exists(aa_path):
        printv(f"ERROR: Can't find aa ({aa_path}) folder. Abort", args.verbose, 0)
        return False
    rmtree(align_path, ignore_errors=True)
    os.mkdir(align_path)

    genes = [
        (gene, os.stat(os.path.join(aa_path, gene)).st_size)
        for gene in os.listdir(aa_path)
        if gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]
    genes.sort(key=lambda x: x[1], reverse=True)
    orthoset_path = os.path.join(args.orthoset_input, args.orthoset)
    aln_path = os.path.join(orthoset_path, ALN_FOLDER)
    only_singletons = set()
    for gene, _ in genes:
        for _, seq in parseFasta(
            os.path.join(aln_path, gene.split(".")[0] + ".aln.fa"), True
        ):
            if len(seq) >= SAFEGUARD_BP:
                printv(f"{gene} will be using Singletons only", args.verbose, 3)
                only_singletons.add(gene)
                # break
    if not os.path.exists(orthoset_path):
        printv("ERROR: Orthoset path not found.", args.verbose, 0)
        return False
    if not os.path.exists(aln_path):
        printv("ERROR: Aln folder not found.", args.verbose, 0)
        return False

    command = "clustalo -i {in_file} -o {out_file} --threads=1 --full"
    if args.second_run:
        command += " --full-iter --iter=1"
    # command = f"mafft --maxiterate 2 --anysymbol --quiet --thread 1 {{in_file}} > {{out_file}}"

    intermediates = "intermediates"
    if not os.path.exists(intermediates):
        os.mkdir(intermediates)

    func_args = []
    for file, _ in genes:
        gene = file.split(".")[0]
        gene_file = os.path.join(aa_path, file)
        result_file = os.path.join(align_path, file.rstrip(".gz"))
        func_args.append(
            (
                CmdArgs(
                    command,
                    gene_file,
                    result_file,
                    gene,
                    args.verbose,
                    args.compress,
                    aln_path,
                    args.debug,
                    file in only_singletons,
                    args.add_fragments,
                ),
            )
        )

    if args.processes > 1:
        with Pool(args.processes) as pool:
            pool.starmap(run_command, func_args, chunksize=1)
    else:
        [run_command(arg[0]) for arg in func_args]

    # if args.debug:
    # with open("mafft_times.csv", "w") as fp:
    # fp.write("\n".join(["Gene,Clustering,Aligning,Merging,Total"]+[",".join(map(str, i)) for i in times]))

    printv(f"Done! Took {time_keeper.differential():.2f}s overall", args.verbose, 0)
    return True


def main(args):
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    for folder in args.INPUT:
        success = do_folder(folder, args)
        if not success:
            return False

    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {time_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception("Cannot be called directly, please use the module:\nsapphyre align")
