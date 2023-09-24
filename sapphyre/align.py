from __future__ import annotations

from os import path, mkdir, system, listdir, stat
from subprocess import DEVNULL, run
from collections import defaultdict, namedtuple
from math import ceil
from multiprocessing.pool import Pool
from shutil import rmtree
from tempfile import NamedTemporaryFile, TemporaryDirectory

from xxhash import xxh3_64
from phymmr_tools import find_index_pair
from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta

def find_kmers(fasta: dict) -> dict[str, set]:
    """Returns the kmers of length KMER_LEN for each sequence in the fasta dict.

    Args:
    ----
        fasta (dict): A dictionary of header -> sequence
    Returns:
        dict[str, set]: A dictionary of header -> set of kmers
    """
    KMER_LEN = 15
    kmers = defaultdict(set)
    for header, sequence in fasta.items():
        for i in range(0, len(sequence) - KMER_LEN):
            kmer = sequence[i : i + KMER_LEN]
            if "*" not in kmer and "-" not in kmer:
                kmers[header].add(kmer)
    return kmers


def get_aln_path(orthoset_dir: str) -> str:
    """
    Checks if /trimmed subdirectory exists. If so, returns a path to /trimmed.
    Otherwise, returns a path to /aln dir.
    """
    ALN_FOLDER = "aln"
    TRIMMED_FOLDER = "trimmed"

    trimmed = path.join(orthoset_dir, TRIMMED_FOLDER)
    if path.exists(trimmed):
        return trimmed
    else:
        return path.join(orthoset_dir, ALN_FOLDER)


def get_start(sequence: str) -> int:
    """Returns the index of the first non-gap character in the sequence.

    Args:
    ----
        sequence (str): The sequence to check
    Returns:
        int: The index of the first non-gap character
    """
    for i, char in enumerate(sequence):
        if char != "-":
            return i
    return -1


def process_genefile(
    gene_path: str,
) -> tuple[int, dict[str, str], dict[str, str], dict[str, list[str]], dict[str, str]]:
    """Returns the number of sequences, the sequences, the targets, the reinsertions and the
    trimmed header to full header.

    Reinsertions are dupes where the sequences are exactly the same. To save on computation we
    insert these duplicates into the alignment after it has been generated right next to their
    duplicate parents with the same alignment.

    Args:
    ----
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
    for header, sequence in parseFasta(gene_path):
        if not header.endswith("."):
            seq_hash = xxh3_64(sequence).hexdigest()
            trimmed_header_to_full[header[:127]] = header
            if seq_hash not in seq_to_first_header:
                seq_to_first_header[seq_hash] = header
            else:
                reinsertions[seq_to_first_header[seq_hash]].append(header)

                continue
            data[header] = sequence.replace(
                "-", ""
            )  # Delete gaps from previous alignment
        else:
            targets[header.split("|")[2]] = header

    return len(data), data, targets, reinsertions, trimmed_header_to_full


def compare_cluster(
    cluster_a: list[str],
    cluster_b: list[str],
    kmers: dict[str, set],
) -> float:
    """Calculates the average kmer similarity between children of two clusters.

    Args:
    ----
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
    """Generates clusters from a dictionary of header -> sequence.

    Clusters are based on the average identity of kmers between sequences.

    Args:
    ----
        data (dict[str, str]): A dictionary of header -> sequence
    Returns:
        list[list[str]]: A list of each clusters' headers
    """
    KMER_PERCENT = 0.15
    
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
                            for header_to_check in [
                                master_header,
                                *list(child_sets[master_header]),
                            ]:
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
                                                candidate_header,
                                            )
                                            cluster_children[master_header].extend(
                                                cluster_children[candidate_header],
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
    cluster_children: list[list[str]],
    parent_tmpdir: str,
    data: dict[str, str],
) -> list[list[str]]:
    """Seperates sequence records into clusters and subclusters
    if they are larger than SUBCLUSTER_AT.

    Args:
    ----
        cluster_children (list[list[str]]): A list of each clusters' headers
        parent_tmpdir (str): The parent temporary directory
        data (dict[str, str]): A dictionary of header -> sequence
    Returns:
        list[list[str]]: A list of each clusters' headers
    """
    SUBCLUSTER_AT = 1000
    CLUSTER_EVERY = 500  # Aim for x seqs per cluster
    clusters = []
    for this_cluster in cluster_children:
        if this_cluster is not None:
            if len(this_cluster) > SUBCLUSTER_AT:
                clusters_to_create = ceil(len(this_cluster) / CLUSTER_EVERY)
                with NamedTemporaryFile(
                    mode="w+",
                    dir=parent_tmpdir,
                    suffix=".fa",
                ) as this_tmp:
                    writeFasta(
                        this_tmp.name,
                        [(header, data[header]) for header in this_cluster],
                    )
                    with NamedTemporaryFile("r", dir=gettempdir()) as this_out:
                        run(
                            f"sigclust/SigClust -k 8 -c {clusters_to_create} {this_tmp.name} > {this_out.name}",
                            stdout=DEVNULL,
                            stderr=DEVNULL,
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
                                this_cluster[int(seq_index)],
                            )
                    clusters.extend(list(sub_clusters.values()))
                    continue
            clusters.append(this_cluster)

    return clusters


def generate_tmp_aln(
    aln_file: str,
    targets: dict[str, str],
    dest: NamedTemporaryFile,
    parent_tmpdir: str,
    debug: float,
    this_intermediates: str,
    align_method: str,
) -> None:
    """Grabs target reference sequences and removes empty columns from the alignment.

    Args:
    ----
        aln_file (str): The path to the alignment file
        targets (dict[str, str]): A dictionary of target -> header
        tmp (NamedTemporaryFile): The temporary file to write to
        debug (float): Whether to write debug files
        this_intermediates (str): The path to the intermediates debug folder
    Returns:
        None
    """
    with NamedTemporaryFile(dir=parent_tmpdir, mode="w+", prefix="References_") as tmp_prealign:
        sequences = []
        for header, sequence in parseFasta(aln_file, True):
            if header in targets:
                sequences.append((targets[header], sequence))

        if align_method in {"base", "frags"}:
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
                            [let for col, let in enumerate(sequence) if not empty_columns[col]],
                        ),
                    ),
                )
        else:
            to_write = []
            for header, sequence in sequences:
                to_write.append(
                    (
                        header,
                        sequence.replace("-",""),
                    ),
                )

        if len(to_write) <= 1 or align_method in {"base", "frags"}:
            writeFasta(dest.name, to_write)
        else:
            writeFasta(tmp_prealign.name, to_write)
            tmp_prealign.flush()

            if align_method == "mafft":
                system(f"mafft-linsi --quiet --thread 1 --anysymbol '{tmp_prealign.name}' > '{dest.name}'")
            elif align_method == "clustal":
                system(f"clustalo -i '{tmp_prealign.name}' -o '{dest.name}' --thread=1 --full --force")

    if debug:
        writeFasta(path.join(this_intermediates, "references.fa"), parseFasta(dest.name, True))

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
        "align_method",
    ],
)


def run_command(args: CmdArgs) -> None:
    keeper = TimeKeeper(KeeperMode.DIRECT)
    debug = args.debug
    printv(f"Doing: {args.gene} ", args.verbose, 2)

    temp_dir = gettempdir()

    aligned_ingredients = []

    this_intermediates = path.join("intermediates", args.gene)
    if debug and not path.exists(this_intermediates):
        mkdir(this_intermediates)
    aln_file = path.join(args.aln_path, args.gene + ".aln.fa")

    with TemporaryDirectory(dir=temp_dir) as parent_tmpdir, TemporaryDirectory(
        dir=parent_tmpdir,
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

        cluster_time = 0
        align_time = 0
        merge_time = 0

        if seq_count == 1:
            printv(
                f"Outputting singleton alignment. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            aligned_file = path.join(aligned_files_tmp, f"{args.gene}_cluster_1_0")
            sequences = list(data.items())
            writeFasta(aligned_file, sequences)
            if args.debug:
                writeFasta(path.join(this_intermediates, aligned_file), sequences)
                unaligned_path = path.join(this_intermediates, f"unaligned_cluster_singleton")
                writeFasta(unaligned_path, sequences)
            aligned_ingredients.append((aligned_file, len(sequences), 0))
        elif args.align_method == "frags":
            aligned_file = path.join(aligned_files_tmp, f"{args.gene}_sequences")

            sequences = list(data.items())

            if debug:
                writeFasta(
                    path.join(this_intermediates, f"{args.gene}_sequences"),
                    sequences,
                )

            
            writeFasta(aligned_file, sequences)
            
            aligned_ingredients = [(aligned_file, len(sequences), 0)]
        else:
            printv(
                f"Generating Cluster. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug

            cluster_children = generate_clusters(data)
            clusters = seperate_into_clusters(cluster_children, parent_tmpdir, data)
            cluster_time = keeper.differential()
            if debug:
                for i, cluster in enumerate(clusters):
                    unaligned_path = path.join(this_intermediates, f"unaligned_cluster_{i}")
                    test_output = [(header, data[header]) for header in cluster]
                    writeFasta(unaligned_path, test_output)
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

            items = listdir(raw_files_tmp)
            items.sort(key=lambda x: int(x.split("_cluster")[-1]))

            for cluster_i, cluster in enumerate(clusters):
                cluster_seqs = [(header, data[header]) for header in cluster]
                cluster_length = len(cluster_seqs)

                printv(
                    f"Aligning cluster {cluster_i}. Elapsed time: {keeper.differential():.2f}",
                    args.verbose,
                    3,
                )  # Debug

                # this_clus_align = f"{args.gene}_cluster_{cluster_i}_length_{len(cluster)}_index_aligned"
                this_clus_align = f"aligned_cluster_{cluster_i}"
                aligned_cluster = path.join(aligned_files_tmp, this_clus_align)

                raw_cluster = path.join(
                    raw_files_tmp,
                    f"{args.gene}_cluster{cluster_i}",
                )
                if (cluster_length == 1):
                    writeFasta(aligned_cluster, cluster_seqs)
                    aligned_ingredients.append((aligned_cluster, len(cluster_seqs), cluster_i))
                    if debug:
                        writeFasta(
                            path.join(
                                this_intermediates,
                                this_clus_align,
                            ),
                            cluster_seqs,
                        )
                    continue
                    
                writeFasta(raw_cluster, cluster_seqs)
                if debug:
                    writeFasta(
                        path.join(this_intermediates, this_clus_align),
                        cluster_seqs,
                    )

                command = args.string.format(
                    in_file=raw_cluster,
                    out_file=aligned_cluster,
                )
                system(command)

                aligned_sequences = list(parseFasta(aligned_cluster, True))

                if debug:
                    printv(command, args.verbose, 3)
                    writeFasta(
                        path.join(
                            this_intermediates,
                            this_clus_align,
                        ),
                        aligned_sequences,
                    )

                aligned_ingredients.append((aligned_cluster, len(aligned_sequences), cluster_i))
        align_time = keeper.differential() - cluster_time

        if aligned_ingredients:
            aligned_ingredients = [
                i for i in sorted(aligned_ingredients, key=lambda x: x[1])
            ]
            printv(
                f"Merging Alignments. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            with NamedTemporaryFile(dir=parent_tmpdir, mode="w+", prefix="References_") as tmp_aln:
                generate_tmp_aln(aln_file, targets, tmp_aln, parent_tmpdir, debug, this_intermediates, args.align_method)

                for i, (file, seq_count, cluster_i) in enumerate(aligned_ingredients):
                    printv(
                        f"Creating reference subalignment {i+1} of {len(aligned_ingredients)}.",
                        args.verbose,
                        3,
                    )
                    if args.align_method == "frags":
                        out_file = path.join(parent_tmpdir, f"aligned.fa")
                        command = f"mafft --anysymbol --quiet --jtt 1 --addfragments {file} --thread 1 {tmp_aln.name} > {out_file}"
                        
                        system(command)

                        final_sequences = list(parseFasta(out_file, True))
                        continue

                    out_file = path.join(parent_tmpdir, f"part_{cluster_i}.fa")

                    if seq_count >= 1:
                        is_profile = "--is-profile"
                    else:
                        is_profile = ""

                    system(
                        f"clustalo --p1 {tmp_aln.name} --p2 {file} -o {out_file} --threads=1 --full {is_profile} --force",
                    )

                    if debug:
                        printv(
                            f"clustalo --p1 {tmp_aln.name} --p2 {file} -o {out_file} --threads=1 --full {is_profile} --force",
                            args.verbose,
                            3,
                        )
                    if debug:
                        writeFasta(
                            path.join(this_intermediates, f"reference_subalignment_{cluster_i}.fa"),
                            parseFasta(out_file, True),
                        )

                if args.align_method != "frags":
                    data = []

                    alignment_insertion_coords = []
                    subalignments = {}
                    insertion_coords_global = set()

                    for item in listdir(parent_tmpdir):
                        if item.startswith("References_"):
                            refs = list(parseFasta(path.join(parent_tmpdir, item), True))
                            continue

                        if not item.startswith("part_"):
                            continue

                        references = []
                        this_seqs = []
                        for header, seq in parseFasta(path.join(parent_tmpdir, item), True):
                            if header.endswith("."):
                                references.append(seq)
                                continue

                            this_seqs.append((header, list(seq)))
                        
                        data_coords = set()
                        for seq in references:
                            for i, char in enumerate(seq):
                                if i in data_coords:
                                    continue
                                if char != "-":
                                    data_coords.add(i)

                        insertion_coords = set(range(len(seq))) - data_coords

                        insertion_after_overlap_filter = sorted([i for i in insertion_coords if i not in insertion_coords_global])
                        
                        if insertion_after_overlap_filter:
                            alignment_insertion_coords.append(insertion_after_overlap_filter)
                            insertion_coords_global.update(insertion_coords)

                        subalignments[item] = this_seqs, insertion_coords

                    alignment_insertion_coords.sort(key=lambda x: x[0])

                    final_refs = []
                    for header,seq in refs:
                        seq = list(seq)
                        adj = 0
                        for x, coords in enumerate(alignment_insertion_coords):
                            if x != 0:
                                adj += len(alignment_insertion_coords[x-1])
                            
                            for gap in coords:
                                seq.insert(gap+adj, "-")
                        
                        final_refs.append((header, "".join(seq)))

                    for item, (seqs, this_msa_insertions) in subalignments.items():
                        adj = 0
                        for x, coords in enumerate(alignment_insertion_coords):
                            if x != 0:
                                adj += len(alignment_insertion_coords[x-1])

                            for coord in coords:
                                if coord in this_msa_insertions:
                                    continue

                                for _, seq in seqs:
                                    seq.insert(coord+adj, "-")

                        subalignments[item] = [(i[0], "".join(i[1])) for i in seqs]

                    sequences = []
                    for group in subalignments.values():
                        sequences.extend(group)

                    sequences.sort(key=lambda x: find_index_pair(x[1], "-")[0])

                    final_sequences = final_refs+sequences

    merge_time = keeper.differential() - align_time - cluster_time

    # Reinsert and sort
    to_write = []
    references = []
    inserted = 0
    for header, sequence in final_sequences:
        if header.endswith("."):
            references.append((header, sequence))
        else:
            header = trimmed_header_to_full[header[:127]]
            # if not "NODE_15784382" in header:
            #     continue
            if header in reinsertions:
                for insertion_header in reinsertions[header]:
                    inserted += 1
                    to_write.append((insertion_header, sequence))
            to_write.append((header, sequence))

    to_write.sort(key=lambda x: get_start(x[1]))

    writeFasta(args.result_file, references + to_write, compress=args.compress)
    printv(f"Done. Took {keeper.differential():.2f}", args.verbose, 3)  # Debug

    return args.gene, cluster_time, align_time, merge_time, keeper.differential()


def do_folder(folder, args):
    ALIGN_FOLDER = "align"
    AA_FOLDER = "aa"

    printv(f"Processing: {path.basename(folder)}", args.verbose, 0)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    align_path = path.join(folder, ALIGN_FOLDER)
    aa_path = path.join(folder, AA_FOLDER)
    if not path.exists(aa_path):
        printv(f"ERROR: Can't find aa ({aa_path}) folder. Abort", args.verbose, 0)
        return False
    rmtree(align_path, ignore_errors=True)
    mkdir(align_path)

    genes = [
        (gene, stat(path.join(aa_path, gene)).st_size)
        for gene in listdir(aa_path)
        if gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]
    genes.sort(key=lambda x: x[1], reverse=True)
    orthoset_path = path.join(args.orthoset_input, args.orthoset)
    aln_path = get_aln_path(orthoset_path)
    if not path.exists(orthoset_path):
        printv("ERROR: Orthoset path not found.", args.verbose, 0)
        return False
    if not path.exists(aln_path):
        printv("ERROR: Aln folder not found.", args.verbose, 0)
        return False

    command = "clustalo -i {in_file} -o {out_file} --threads=1 --full"
    if args.second_run:
        command += " --full-iter --iter=1"

    intermediates = "intermediates"
    if not path.exists(intermediates):
        mkdir(intermediates)

    func_args = []
    for file, _ in genes:
        gene = file.split(".")[0]
        gene_file = path.join(aa_path, file)
        result_file = path.join(align_path, file.rstrip(".gz"))
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
                    args.align_method.lower(),
                ),
            ),
        )

    if args.processes > 1:
        with Pool(args.processes) as pool:
            pool.starmap(run_command, func_args, chunksize=1)
    else:
        [run_command(arg[0]) for arg in func_args]

    # if args.debug:
    # with open("mafft_times.csv", "w") as fp:

    printv(f"Done! Took {time_keeper.differential():.2f}s overall", args.verbose, 0)
    return True


def main(args):
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    if path.exists("intermediates"):
        rmtree("intermediates")
    mkdir("intermediates")

    for folder in args.INPUT:
        success = do_folder(folder, args)
        if not success:
            return False

    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {time_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre align"
    raise Exception(msg)
