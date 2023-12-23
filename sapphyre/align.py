from __future__ import annotations

from collections import Counter, defaultdict, namedtuple
from math import ceil
from multiprocessing.pool import Pool
from os import listdir, mkdir, path, stat, system
from shutil import rmtree
from tempfile import NamedTemporaryFile, TemporaryDirectory

from phymmr_tools import find_index_pair, sigclust
from xxhash import xxh3_64

from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta


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
    """Returns the number of sequences, the sequences, the target references, the reinsertions and the
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
        if header.endswith("."):
            targets[header.split("|")[2]] = header
            continue

        seq_hash = xxh3_64(sequence).hexdigest()
        trimmed_header_to_full[header[:127]] = header
        if seq_hash not in seq_to_first_header:
            seq_to_first_header[seq_hash] = header
        else:
            reinsertions[seq_to_first_header[seq_hash]].append(header)
            continue
        # Delete gaps from previous alignment
        data[header] = sequence.replace("-", "")

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


def generate_clusters(data: dict[str, str], second_run) -> list[list[str]]:
    """Generates clusters from a dictionary of header -> sequence.

    Clusters are based on the average identity of kmers between sequences.

    Args:
    ----
        data (dict[str, str]): A dictionary of header -> sequence
    Returns:
        list[list[str]]: A list of each clusters' headers
    """

    with NamedTemporaryFile(
        dir=gettempdir(), mode="w+", suffix=".fa"
    ) as tmp_in, NamedTemporaryFile(
        dir=gettempdir(), mode="w+", suffix=".txt"
    ) as tmp_result:
        writeFasta(tmp_in.name, data.items())
        tmp_in.flush()

        if second_run:
            system(
                f"diamond cluster -d {tmp_in.name} -o {tmp_result.name} --approx-id 85 --member-cover 65 --threads 1 --ignore-warnings --quiet"
            )
        else:
            system(
                f"diamond cluster -d {tmp_in.name} -o {tmp_result.name} --approx-id 85 --member-cover 70 --threads 1 --ignore-warnings --quiet"
            )

        cluster_children = defaultdict(list)

        for line in tmp_result.read().split("\n"):
            if line.strip():
                master, child = line.split("\t")
                cluster_children[master].append(child)

    return cluster_children.values()


def seperate_into_clusters(
    cluster_children: list[list[str]],
    data: dict[str, str],
) -> list[list[str]]:
    """Seperates sequence records into clusters and subclusters
    if they are larger than SUBCLUSTER_AT.

    Args:
    ----
        cluster_children (list[list[str]]): A list of each clusters' headers
        data (dict[str, str]): A dictionary of header -> sequence
    Returns:
        list[list[str]]: A list of each clusters' headers
    """
    SUBCLUSTER_AT = 1000  # Subcluster threshold
    CLUSTER_EVERY = 500  # Sequences per subcluster
    clusters = []
    for this_cluster in cluster_children:
        if this_cluster is not None:
            if len(this_cluster) > SUBCLUSTER_AT:
                clusters_to_create = ceil(len(this_cluster) / CLUSTER_EVERY)
                records = [(header, data[header]) for header in this_cluster]
                sub_clusters = sigclust(records, 8, clusters_to_create)
                clusters.extend(sub_clusters)
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
    """Grabs target reference sequences. If the current align method is equal to base, then delete empty columns.
    For all other methods, return the sequences without gaps and realign them.

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
    with NamedTemporaryFile(
        dir=parent_tmpdir, mode="w+", prefix="References_"
    ) as tmp_prealign:
        sequences = []
        for header, sequence in parseFasta(aln_file, True):
            if header in targets:
                if len(sequence) != sequence.count("-"):
                    sequences.append((targets[header], sequence))

        if align_method == "base":
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
                            [
                                let
                                for col, let in enumerate(sequence)
                                if not empty_columns[col]
                            ],
                        ),
                    ),
                )
        else:
            to_write = []
            for header, sequence in sequences:
                to_write.append(
                    (
                        header,
                        sequence.replace("-", ""),
                    ),
                )
        # Skip realignment if there is only one sequence
        if len(to_write) <= 1 or align_method == "base":
            writeFasta(dest.name, to_write)
        else:
            writeFasta(tmp_prealign.name, to_write)
            tmp_prealign.flush()

            if align_method == "clustal":
                system(
                    f"clustalo -i '{tmp_prealign.name}' -o '{dest.name}' --thread=1 --full --force"
                )
            else:
                system(
                    f"mafft --localpair --quiet --thread 1 --anysymbol '{tmp_prealign.name}' > '{dest.name}'"
                )

    if debug:
        writeFasta(
            path.join(this_intermediates, "references.fa"), parseFasta(dest.name, True)
        )


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
        "second_run",
    ],
)


def insert_gaps(input_string, positions, existing_gaps={}):
    """
    Inserts gaps into a sequence at the given positions skipping any existing gaps in the profile.

    Args:
    ----
        input_string (str): The sequence to insert gaps into
        positions (list[int]): The positions to insert gaps at
        existing_gaps (dict[int, int]): The existing gaps in the profile
    Returns:
        str: The sequence with gaps inserted
    """
    input_string = list(input_string)
    gap_offset = 0

    for coord in positions:
        if existing_gaps.get(coord, 0) > 0:
            gap_offset += 1
            existing_gaps[coord] -= 1
            continue

        input_string.insert(coord + gap_offset, "-")
        gap_offset += 1

    return "".join(input_string)


def align_cluster(
    cluster_seqs,
    cluster_i,
    gene,
    raw_files_tmp,
    aligned_files_tmp,
    debug,
    this_intermediates,
    command_string,
):
    """
    Performs a pairwise alignment on a cluster of sequences.

    Args:
    ----
        cluster_seqs (list[tuple[str, str]]): A list of sequence records
        cluster_i (int): And index for the current cluster
        gene (str): The name of the current gene
        raw_files_tmp (str): The path to the raw files temporary directory
        aligned_files_tmp (str): The path to the aligned files temporary directory
        debug (bool): Whether to write debug files
        this_intermediates (str): The path to the intermediates debug folder
        command_string (str): The command to run for alignment
    Returns:
    -------
        tuple[str, int, int]: The path to the aligned file, the number of sequences in the cluster and the cluster index
    """
    cluster_length = len(cluster_seqs)
    this_clus_align = f"aligned_cluster_{cluster_i}"
    aligned_cluster = path.join(aligned_files_tmp, this_clus_align)

    raw_cluster = path.join(
        raw_files_tmp,
        f"{gene}_cluster{cluster_i}",
    )
    if cluster_length == 1:
        writeFasta(aligned_cluster, cluster_seqs)

        if debug:
            writeFasta(
                path.join(
                    this_intermediates,
                    this_clus_align,
                ),
                cluster_seqs,
            )
        return (aligned_cluster, len(cluster_seqs), cluster_i)

    writeFasta(raw_cluster, cluster_seqs)
    if debug:
        writeFasta(
            path.join(this_intermediates, this_clus_align),
            cluster_seqs,
        )

    command = command_string.format(
        in_file=raw_cluster,
        out_file=aligned_cluster,
    )
    system(command)

    aligned_sequences = list(parseFasta(aligned_cluster, True))

    if debug:
        print(command)
        writeFasta(
            path.join(
                this_intermediates,
                this_clus_align,
            ),
            aligned_sequences,
        )

    return (aligned_cluster, len(aligned_sequences), cluster_i)


def create_subalignment(
    align_method: str,
    parent_tmpdir: str,
    file: str,
    tmp_aln: NamedTemporaryFile,
    cluster_i: int,
    seq_count: int,
    debug: bool,
    this_intermediates: str,
) -> list[tuple[str, str]]:
    """
    Creates a subalignment by aligning an aligned cluster with the reference alignment.

    Args:
    ----
        align_method (str): The alignment method to use
        parent_tmpdir (str): The path to the parent temporary directory
        file (str): The path to the aligned cluster
        tmp_aln (NamedTemporaryFile): The temporary file containing the reference alignment
        cluster_i (int): The index of the current cluster
        seq_count (int): The number of sequences in the cluster
        debug (bool): Whether to write debug files
        this_intermediates (str): The path to the intermediates debug folder
    Returns:
    -------
        list[tuple[str, str]]: The sequences in the subalignment
    """
    if align_method == "frags":
        out_file = path.join(parent_tmpdir, f"aligned.fa")
        command = f"mafft --anysymbol --quiet --jtt 1 --addfragments {file} --thread 1 {tmp_aln.name} > {out_file}"

        system(command)

        return list(parseFasta(out_file, True))

    out_file = path.join(parent_tmpdir, f"part_{cluster_i}.fa")

    if seq_count >= 1:
        is_profile = "--is-profile"
    else:
        is_profile = ""

    system(
        f"clustalo --p1 {tmp_aln.name} --p2 {file} -o {out_file} --threads=1 --full {is_profile} --force",
    )

    if debug:
        print(
            f"clustalo --p1 {tmp_aln.name} --p2 {file} -o {out_file} --threads=1 --full {is_profile} --force"
        )
        writeFasta(
            path.join(
                this_intermediates,
                f"reference_subalignment_{cluster_i}.fa",
            ),
            parseFasta(out_file, True),
        )

    return []


def get_insertions(parent_tmpdir: str) -> tuple[list, list, list]:
    """
    Gets the insertion coords for each subalignment

    Args:
    ----
        parent_tmpdir (str): The path to the parent temporary directory containing all the subalignmnets
    Returns:
    -------
        tuple[list, list, list]: The insertion coords, the subalignments and the reference sequences
    """
    subalignments = {}

    global_insertions = Counter()

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

            this_seqs.append((header, seq))

        insertion_coords = []
        for i in range(len(references[0])):
            if all(seq[i] == "-" for seq in references):
                insertion_coords.append(i - len(insertion_coords))

        this_counter = Counter(insertion_coords).items()

        for key, value in this_counter:
            if global_insertions.get(key, 0) < value:
                global_insertions[key] = value

        subalignments[item] = this_seqs, this_counter

    alignment_insertion_coords = []
    for key, value in global_insertions.items():
        alignment_insertion_coords.extend([key] * value)

    return alignment_insertion_coords, subalignments.values(), refs


def insert_refs(refs, alignment_insertion_coords):
    """
    Inserts insertion gaps into the reference sequences at the given positions.

    Args:
    ----
        refs (list[tuple[str, str]]): The reference sequences
        alignment_insertion_coords (list[int]): The insertion coords
    Returns:
    -------
        list[tuple[str, str]]: The reference sequences with gaps inserted
    """
    final_refs = []

    alignment_insertion_coords.sort()
    for header, seq in refs:
        if alignment_insertion_coords:
            seq = insert_gaps(seq, alignment_insertion_coords)

        final_refs.append((header, seq))

    return final_refs


def insert_sequences(
    subalignments: list, alignment_insertion_coords: list[int]
) -> list[tuple[str, str]]:
    """
    Inserts insertion gaps into the candidate sequences while ignoring and masking existing gaps.

    Args:
    ----
        subalignments (list): The subalignments
        alignment_insertion_coords (list[int]): The insertion coords
    Returns:
    -------
        list[tuple[str, str]]: The candidate sequences with gaps inserted
    """
    out_seqs = []
    for seqs, this_msa_insertions in subalignments:
        for header, seq in seqs:
            if alignment_insertion_coords:
                existing_gaps = {k: v for k, v in this_msa_insertions if v > 0}
                seq = insert_gaps(seq, alignment_insertion_coords, existing_gaps)
            out_seqs.append((header, seq))

    out_seqs.sort(key=lambda x: find_index_pair(x[1], "-")[0])

    return out_seqs


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
            f"Cleaning AA file. Elapsed time: {keeper.differential():.2f}",
            args.verbose,
            3,
        )  # Debug

        # Grab sequences, reinsertions, trimmed header to full header dict and targets
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

        # If only one sequence is present we can just output the singleton
        if seq_count == 1:
            printv(
                f"Outputting singleton alignment. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            aligned_file = path.join(aligned_files_tmp, f"{args.gene}_cluster_1_0")
            sequences = list(data.items())
            writeFasta(aligned_file, sequences)
            if debug:
                writeFasta(path.join(this_intermediates, aligned_file), sequences)
                unaligned_path = path.join(
                    this_intermediates, f"unaligned_cluster_singleton"
                )
                writeFasta(unaligned_path, sequences)
            aligned_ingredients.append((aligned_file, len(sequences), 0))
        # Otherwise if align method is frags
        elif args.align_method == "frags":
            aligned_file = path.join(aligned_files_tmp, f"{args.gene}_sequences")

            sequences = list(data.items())

            if debug:
                writeFasta(
                    path.join(this_intermediates, f"{args.gene}_sequences"),
                    sequences,
                )

            # Write the aligned ingredients to a temp file and set the aligned ingredients
            # to equal just this file
            writeFasta(aligned_file, sequences)

            aligned_ingredients = [(aligned_file, len(sequences), 0)]
        # Else any other method is selected
        else:
            printv(
                f"Generating Cluster. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug

            # Seperate into clusters if more than 5 sequences are present
            if len(data) < 5:
                cluster_children = [[i] for i in data]
            else:
                cluster_children = generate_clusters(data, args.second_run)
            clusters = seperate_into_clusters(cluster_children, data)
            cluster_time = keeper.differential()

            if debug:
                for i, cluster in enumerate(clusters):
                    unaligned_path = path.join(
                        this_intermediates, f"unaligned_cluster_{i}"
                    )
                    test_output = [(header, data[header]) for header in cluster]
                    writeFasta(unaligned_path, test_output)

            printv(
                f"Found {seq_count} sequences over {len(clusters)} clusters. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )
            printv(
                f"Aligning Clusters. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )

            # Align each cluster
            for cluster_i, cluster in enumerate(clusters):
                cluster_seqs = [(header, data[header]) for header in cluster]

                printv(
                    f"Aligning cluster {cluster_i}. Elapsed time: {keeper.differential():.2f}",
                    args.verbose,
                    3,
                )  # Debug

                aligned_ingredients.append(
                    align_cluster(
                        cluster_seqs,
                        cluster_i,
                        args.gene,
                        raw_files_tmp,
                        aligned_files_tmp,
                        debug,
                        this_intermediates,
                        args.string,
                    )
                )

        align_time = keeper.differential() - cluster_time

        if aligned_ingredients:
            # Merge subalignments in order by seq count
            aligned_ingredients = [
                i for i in sorted(aligned_ingredients, key=lambda x: x[1])
            ]
            printv(
                f"Merging Alignments. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            with NamedTemporaryFile(
                dir=parent_tmpdir, mode="w+", prefix="References_"
            ) as tmp_aln:
                generate_tmp_aln(
                    aln_file,
                    targets,
                    tmp_aln,
                    parent_tmpdir,
                    debug,
                    this_intermediates,
                    args.align_method,
                )

                for i, (file, seq_count, cluster_i) in enumerate(aligned_ingredients):
                    printv(
                        f"Creating reference subalignment {i+1} of {len(aligned_ingredients)}.",
                        args.verbose,
                        3,
                    )
                    # If frags we want to save the final_sequences from the aligned result
                    final_sequences = create_subalignment(
                        args.align_method,
                        parent_tmpdir,
                        file,
                        tmp_aln,
                        cluster_i,
                        seq_count,
                        debug,
                        this_intermediates,
                    )

                # Otherwise we want to do insertion logic
                if args.align_method != "frags":
                    # Grab insertions in each subalignment
                    alignment_insertion_coords, subalignments, refs = get_insertions(
                        parent_tmpdir
                    )

                    # Insert into refs
                    final_refs = insert_refs(refs, alignment_insertion_coords)

                    # Insert into each subalignment
                    sequences = insert_sequences(
                        subalignments, alignment_insertion_coords
                    )

                    # Consolidate output
                    final_sequences = final_refs + sequences

    merge_time = keeper.differential() - align_time - cluster_time

    # Reinsert and sort by start position
    to_write = []
    references = []
    inserted = 0
    for header, sequence in final_sequences:
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
        printv(f"Please make sure Reporter finished succesfully", args.verbose, 0)
        return False
    if args.overwrite:
        rmtree(align_path, ignore_errors=True)
    if not path.exists(align_path):
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

    intermediates = "intermediates"
    if not path.exists(intermediates):
        mkdir(intermediates)

    func_args = []
    printv(
        f"Aligning AA Files. Elapsed time: {time_keeper.differential():.2f}s",
        args.verbose,
    )
    for file, _ in genes:
        gene = file.split(".")[0]
        gene_file = path.join(aa_path, file)
        compress_result_file = path.join(align_path, file)
        result_file = compress_result_file.rstrip(".gz")

        exst_file = result_file if path.exists(result_file) else None
        if not exst_file:
            exst_file = (
                compress_result_file if path.exists(compress_result_file) else None
            )

        if not exst_file or stat(exst_file).st_size == 0:
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
                        args.second_run,
                    ),
                ),
            )

    if args.processes > 1:
        with Pool(args.processes) as pool:
            pool.starmap(run_command, func_args, chunksize=1)
    else:
        [run_command(arg[0]) for arg in func_args]

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
