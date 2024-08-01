from __future__ import annotations
from collections import Counter, defaultdict, namedtuple
import itertools
import subprocess
from math import ceil
from multiprocessing.pool import Pool
from os import listdir, mkdir, path, stat, system
from shutil import rmtree
from statistics import median
from numpy import percentile
from subprocess import Popen, run
from tempfile import NamedTemporaryFile, TemporaryDirectory
from time import time
from psutil import Process
from sapphyre_tools import blosum62_distance, find_index_pair, sigclust, delete_empty_columns
from xxhash import xxh3_64
from wrap_rocks import RocksDB
from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta

def get_aln_path(orthoset_dir: str) -> str:
    """
    Checks if /trimmed subdirectory exists. If so, returns a path to /trimmed.
    Otherwise, returns a path to /aln dir.
    """
    ALN_FOLDER = "aln"
    TRIMMED_FOLDER = "trimmed"
    CLEAN_FOLDER = "cleaned"
    FINAL_FOLDER = "final"
    RAW_FOLDER = "raw"

    final = path.join(orthoset_dir, FINAL_FOLDER)
    if path.exists(final):
        return final
    
    cleaned = path.join(orthoset_dir, CLEAN_FOLDER)
    if path.exists(cleaned):
        return cleaned

    trimmed = path.join(orthoset_dir, TRIMMED_FOLDER)
    if path.exists(trimmed):
        return trimmed
    
    aln = path.join(orthoset_dir, ALN_FOLDER)
    if path.exists(aln):    
        return aln
    
    return path.join(orthoset_dir, RAW_FOLDER)


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

        # if second_run:
        #     system(
        #         f"diamond cluster -d {tmp_in.name} -o {tmp_result.name} --approx-id 85 --member-cover 65 --threads 1 --ignore-warnings --quiet"
        #     )
        # else:
        #     system(
        #         f"diamond cluster -d {tmp_in.name} -o {tmp_result.name} --approx-id 85 --member-cover 70 --threads 1 --ignore-warnings --quiet"
        #     )
        terminal_args = [
            "diamond", "cluster",
            "-d", str(tmp_in.name),
            "-o", str(tmp_result.name),
            "--approx-id", "85",
            "--member-cover", "70",
            "--threads", "1",
            "--ignore-warnings", "--quiet"]
        if second_run:
            terminal_args[9] = "65"

        with Popen(terminal_args) as diamond_run:
            diamond_process = Process(diamond_run.pid)
            first_zero_time = None
            DIAMOND_WAIT = 10
            while True:
                try:
                    cpu_percent = diamond_process.cpu_percent()
                    if diamond_run.poll() is not None:  # if process has a finished, break
                        break
                    if cpu_percent == 0.0:  # if process is not using cpu, break
                        if first_zero_time is None:
                            first_zero_time = time()
                        else:
                            if time() - first_zero_time > DIAMOND_WAIT:
                                diamond_process.terminate()
                                print(f"Zombie diamond process {diamond_run.pid} found. Breaking and rerunning")
                                diamond_process.terminate()
                                break
                    else:
                        first_zero_time = None
                    #sleep(1)  # sleep to avoid spamming the cpu with polling ops
                except KeyboardInterrupt:
                    break
                
            if diamond_run.returncode != 0:
                print(f"Non-zero exit code for diamond, rerunning process {diamond_run.pid} on {tmp_in.name}")
                run(terminal_args, check=False)
            cluster_children = defaultdict(list)

        for line in tmp_result.read().splitlines():
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
            if header.split(" ")[0] in targets:
                if len(sequence) != sequence.count("-"):
                    sequences.append((targets[header.split(" ")[0]], sequence))

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
                
                cmd = ["clustalo", "-i", tmp_prealign.name, "-o", dest.name, "--thread=1", "--full", "--force"]
                subprocess.run(cmd, stdout=subprocess.DEVNULL)
                #system(
                    #f"clustalo -i '{tmp_prealign.name}' -o '{dest.name}' --thread=1 --full --force"
                #)
            else:
                system(
                    f"mafft --localpair --quiet --thread 1 --anysymbol '{tmp_prealign.name}' > '{dest.name}'"
                )

            recs = list(parseFasta(dest.name, True))
        
            del_columns = set()
            for i in range(len(recs[0][1])):
                if all(j[1][i] == "-" or j[1][i] == "X" for j in recs):
                    del_columns.add(i)

            out = [(header, "".join([let for i, let in enumerate(seq) if i not in del_columns])) for header, seq in recs]

            writeFasta(dest.name, out)
    dest.flush()
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
        "top_folder",
        "chomp_max_distance",
        "is_genome",
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
    tmp_aln: str,
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
        out_file = path.join(parent_tmpdir, "aligned.fa")
        command = f"mafft --anysymbol --quiet --jtt 1 --addfragments {file} --thread 1 {tmp_aln} > {out_file}"

        system(command)

        return list(parseFasta(out_file, True))

    out_file = path.join(parent_tmpdir, f"part_{cluster_i}.fa")

    if seq_count >= 1:
        is_profile = "--is-profile"
    else:
        is_profile = ""

    cmd = ["clustalo", "--p1", tmp_aln, "--p2", file, "-o", out_file, "--threads=1", "--full", is_profile, "--force"]
    subprocess.run(cmd)
    #system(
        #f"clustalo --p1 {tmp_aln} --p2 {file} -o {out_file} --threads=1 --full {is_profile} --force",
    #)

    if debug:
        print(
            f"clustalo --p1 {tmp_aln} --p2 {file} -o {out_file} --threads=1 --full {is_profile} --force"
        )
        writeFasta(
            path.join(
                this_intermediates,
                f"reference_subalignment_{cluster_i}.fa",
            ),
            parseFasta(out_file, True),
        )

    return []


def get_insertions(parent_tmpdir: str, target_dict, ref_path) -> tuple[list, list, list]:
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

    refs = []
    for header, seq in parseFasta(ref_path, True):
        if header in target_dict:
            refs.append((target_dict[header], seq))
        elif header.endswith('.'):
            refs.append((header, seq))

    for item in listdir(parent_tmpdir):
        if not item.startswith("part_"):
            continue

        references = []
        this_seqs = []
        for header, seq in parseFasta(path.join(parent_tmpdir, item), True):
            if "|" not in header or header.endswith("."):
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


def do_cluster(ids, ref_coords, id_chomp_distance=100):
    clusters = []
    ids.sort(key = lambda x: x[0])
    grouped_ids = defaultdict(list)
    for i, (child_index, seq_coords, start, end) in enumerate(ids):
        id = int(child_index.split("_")[0])
        grouped_ids[id].append((i, child_index, seq_coords, start, end))
        
    ids_ascending = sorted(grouped_ids.keys())
    

    req_seq_coverage = 0.5

    current_cluster = []
    
    for id in ids_ascending:
        seq_list = grouped_ids[id]
        if not current_cluster:
            current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list]
            current_index = id
            current_seqs = seq_list
            current_direction = "bi"
        else:
            passed = False
            passing_direction = None
            if abs(id - current_index) <= id_chomp_distance:
                for i, child_index, seq_coords, start, end in seq_list:
                    for _, _, _, current_start, current_end in current_seqs:
                        this_direction = None
                        if start == current_start and end == current_end:
                            this_direction = "bi"
                        else:
                            if start == current_start:
                                if end >= current_end:
                                    this_direction = "forward"
                                else:
                                    this_direction = "reverse"
                            else:
                                if start >= current_start:
                                    this_direction = "forward"
                                else:
                                    this_direction = "reverse"

                     
                        if current_direction == "bi" or this_direction == "bi" or this_direction == current_direction:
                            passed = True
                            passing_direction = this_direction
                            break
                    if passed:
                        break
                
            if passed:
                current_cluster.extend([(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list])
                current_index = id
                current_seqs = seq_list
                if passing_direction != "bi":
                    current_direction = passing_direction
                
            else:
                if len(current_cluster) >= 2:
                    cluster_data_cols = set()
                    for _, _, index in current_cluster:
                        cluster_data_cols.update(ids[index][1])
                        
                    cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

                    clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage))
                elif len(current_cluster) == 1:
                    if current_cluster[0][1] > req_seq_coverage:
                        cluster_coverage = current_cluster[0][1]
                        clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage))
                        
                current_cluster = [(id, len(seq_coords.intersection(ref_coords)) / len(ref_coords), i) for i, _, seq_coords, _, _ in seq_list]
                current_index = id
                current_seqs = seq_list
                current_direction = "bi"
    
    if current_cluster:
        if len(current_cluster) >= 2:
            cluster_data_cols = set()
            for _, _, index in current_cluster:
                cluster_data_cols.update(ids[index][1])
                
            cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

            clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage))
        elif len(current_cluster) == 1:
            if current_cluster[0][1] > req_seq_coverage:
                cluster_coverage = current_cluster[0][1]
                clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage))
                
    return clusters


def cull_reference_outliers(reference_records: list, debug: int) -> list:
    """
    Removes reference sequences which have an unusually large mean
    blosum distance. Finds the constrained blosum distance between
    each reference pull any reference with a mean 1.5x higher than
    the group mean. Returns the remaining references in a list.
    """
    needs_to_check = True
    filtered = []
    loop = 0
    while needs_to_check:
        loop += 1
        needs_to_check = False
        distances_by_index = defaultdict(list)
        all_distances = []
        indices = {i:find_index_pair(reference_records[i][1], '-') for i in range(len(reference_records))}
        # generate reference distances
        for i, ref1 in enumerate(reference_records[:-1]):
            start1, stop1 = indices[i]
            for j, ref2 in enumerate(reference_records[i+1:],i+1):
                start2, stop2 = indices[j]
                start = max(start1, start2)
                stop = min(stop1, stop2)
                # avoid the occasional rust-nan result
                if start >= stop:
                    distances_by_index[i].append(1)
                    distances_by_index[j].append(1)
                    continue
                dist = blosum62_distance(ref1[1][start:stop], ref2[1][start:stop]) ** 2
                distances_by_index[i].append(dist)
                distances_by_index[j].append(dist)
                all_distances.append(dist)

        if not all_distances:
            return list(itertools.chain(*reference_records)), filtered, 0, 0, 0

        total_median = median(all_distances)
        q3, q1 = percentile(all_distances, [75, 25])
        iqr = q3 - q1
        iqr_coeff = 1

        allowable = max(total_median + iqr_coeff * iqr, 0.02)

        # if a record's mean is too high, cull it
        for index, distances in distances_by_index.items():
            this_median = median(distances)
            if this_median > allowable:# or mean > 1:
                distances_by_index[index] = None
                filtered.append( (reference_records[index], this_median, "kicked") )
                needs_to_check = True

            if debug == 2:
                filtered.append( (reference_records[index], this_median, "") )
        reference_records = [reference_records[i] for i in range(len(reference_records)) if distances_by_index[i] is not None]
            # else:
            #     distances_by_index[index] = mean
# get all remaining records
    output = list(itertools.chain(*(reference_records[i] for i in range(len(reference_records)) if distances_by_index[i] is not None)))

    return output, filtered, total_median, allowable, iqr

def run_command(args: CmdArgs) -> None:
    keeper = TimeKeeper(KeeperMode.DIRECT)
    debug = args.debug
    printv(f"Doing: {args.gene} ", args.verbose, 2)

    temp_dir = gettempdir()

    aligned_ingredients = []
    if args.is_genome:
        get_id = lambda header: header.split("|")[3].split("&&")[0].replace("NODE_","")
    else:
        get_id = lambda header: header

    if not path.exists(args.result_file) or stat(args.result_file).st_size == 0:
        this_intermediates = path.join("intermediates", args.gene)
        if debug and not path.exists(this_intermediates):
            mkdir(this_intermediates)
        aln_file = path.join(args.aln_path, args.gene + ".fa" if "raw" in args.aln_path else args.gene + ".aln.fa")

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
                        this_intermediates, "unaligned_cluster_singleton"
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


            if aligned_ingredients:
                # Merge subalignments in order by seq count
                aligned_ingredients.sort(key=lambda x: x[1])
                printv(
                    f"Merging Alignments. Elapsed time: {keeper.differential():.2f}",
                    args.verbose,
                    3,
                )  # Debug

                # Grab target reference sequences

                top_aln_path = path.join(args.top_folder, args.gene + ".aln.fa")
                realign_rec = args.second_run # Realign if required
                if not path.exists(top_aln_path):
                    realign_rec = True # Force realignment if the top alignment doesn't exist
                    top_aln_path = path.join(args.aln_path, args.gene + ".aln.fa")
                # cull reference outliers
                references = list(parseFasta(top_aln_path, has_interleave=True))
                references, filtered_refs, ref_total_median, ref_allowable, ref_iqr = cull_reference_outliers(references, args.debug)
                references, _ = delete_empty_columns(references)
                references = [(references[i], references[i + 1]) for i in range(0, len(references), 2)]
                
                culled_references = []
                if filtered_refs:
                    culled_references.append(f'{args.gene} total median: {ref_total_median}\n')
                    culled_references.append(f'{args.gene} threshold: {ref_allowable}\n')
                    culled_references.append(f'{args.gene} standard deviation: {ref_iqr}\n')
                    for ref_kick, ref_median, kick in filtered_refs:
                        culled_references.append(f'{ref_kick[0]},{ref_median},{kick}\n')
                tmp_aln = NamedTemporaryFile(dir=parent_tmpdir, mode="w+", prefix="References_")

                if realign_rec:
                    print(aln_file)
                    generate_tmp_aln(
                        aln_file,
                        # references,
                        targets,
                        tmp_aln,
                        parent_tmpdir,
                        debug,
                        this_intermediates,
                        args.align_method,
                    )
                else:
                    writeFasta(tmp_aln.name, references)
                    tmp_aln.flush()
                    
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
                        tmp_aln.name, #if realign_rec else top_aln_path,
                        cluster_i,
                        seq_count,
                        debug,
                        this_intermediates,
                    )

                # Otherwise we want to do insertion logic
                if args.align_method != "frags":
                    # Grab insertions in each subalignment
                    alignment_insertion_coords, subalignments, refs = get_insertions(
                        parent_tmpdir, targets, tmp_aln.name
                    )

                    # Insert into refs
                    final_refs = insert_refs(refs, alignment_insertion_coords)

                    # Insert into each subalignment
                    sequences = insert_sequences(
                        subalignments, alignment_insertion_coords
                    )

                    # Consolidate output
                    final_sequences = final_refs + sequences
                
                # if realign_rec:
                del tmp_aln

        # Reinsert and sort by start position
        to_write = []
        references = []
        inserted = 0
        for header, sequence in final_sequences:
            if "|" not in header:
                if header.split(" ")[0] not in targets:
                    continue
                
                header = targets[header.split(" ")[0]]
            elif not header.endswith("."):
                header = trimmed_header_to_full[header[:127]]

            if header.endswith('.'):
                references.append((header, sequence))
                continue

            if header in reinsertions:
                for insertion_header in reinsertions[header]:
                    inserted += 1
                    to_write.append((insertion_header, sequence))
            to_write.append((header, sequence))

        if args.is_genome:
            reordered = []
            group = defaultdict(list)
            for header, seq in to_write:
                group[int(header.split("|")[3].split("&&")[0].split("_")[1])].append((header, seq))
                
            for key in sorted(group.keys()):
                seqs = group[key]
                
                if len(seqs) == 1:
                    reordered.append(seqs[0])
                else:
                    strand = int(seqs[0][0].split("|")[4])
                    if strand < 0:
                        seqs.reverse()
                    
                    reordered.extend(seqs)
                
            to_write = reordered
        else:
            to_write.sort(key=lambda x: find_index_pair(x[1], "-")[0])
        writeFasta(args.result_file, references + to_write, compress=args.compress)
        
        ids = []
        for header, seq in to_write:
            start, end = find_index_pair(seq, "-")
            data_cols = {i for i, let in enumerate(seq[start:end], start) if let != "-"}
            ids.append((get_id(header), data_cols, start, end))
        ref_coords = set()
        for _, seq in references:
            start, end = find_index_pair(seq, "-")
            for i, let in enumerate(seq[start:end], start):
                if let != "-":
                    ref_coords.add(i)
    else:
        ids = []
        ref_coords = set()
        for header, seq in parseFasta(args.result_file):
            if header.endswith("."):
                start, end = find_index_pair(seq, "-")
                for i, let in enumerate(seq[start:end], start):
                    if let != "-":
                        ref_coords.add(i)
                continue

            start, end = find_index_pair(seq, "-")
            data_cols = {i for i, let in enumerate(seq[start:end], start) if let != "-"}
            ids.append((get_id(header), data_cols, start, end))

    clusters = []
    if args.is_genome:
        clusters = do_cluster(ids, ref_coords, args.chomp_max_distance)

        cluster_string = ", ".join([f"{cluster[0]}-{cluster[1]} {(cluster[2]*100):.2f}%" for cluster in clusters])         
            
    printv(f"Done. Took {keeper.differential():.2f}", args.verbose, 3)  # Debug

    if not args.is_genome:
        return
    return (args.gene, f"{args.gene},{len(ids)},{len(clusters)},{cluster_string}")


def do_folder(folder, args):
    ALIGN_FOLDER = "align"
    AA_FOLDER = "aa"

    printv(f"Processing: {path.basename(folder)}", args.verbose, 0)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    align_path = path.join(folder, ALIGN_FOLDER)
    aa_path = path.join(folder, AA_FOLDER)
    if args.use_miniprot:
        aa_path = path.join(folder, "miniprot", "aa")
        align_path = path.join(folder, "miniprot", "align")
    if not path.exists(aa_path):
        printv(f"ERROR: Can't find aa ({aa_path}) folder. Abort", args.verbose, 0)
        printv("Please make sure Reporter finished succesfully", args.verbose, 0)
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

    top_folder = path.join(folder, "top")

    intermediates = "intermediates"
    if not path.exists(intermediates) and args.debug:
        mkdir(intermediates)

    func_args = []
    printv(
        f"Aligning AA Files. Elapsed time: {time_keeper.differential():.2f}s",
        args.verbose,
    )

    rocks_db_path = path.join(folder, "rocksdb", "sequences", "nt")
    is_genome = False
    if args.second_run:
        pass
    elif not path.exists(rocks_db_path):
        printv(f"WARNING: Can't find rocksdb folder in {folder}. Unable to determine if datsets is genomic", args.verbose, 0)
    else:
        rocksdb_db = RocksDB(str(rocks_db_path))
        is_genome = rocksdb_db.get("get:isgenome")
        is_genome = is_genome == "True"
        del rocksdb_db
    
    for file, _ in genes:
        gene = file.split(".")[0]
        gene_file = path.join(aa_path, file)
        compress_result_file = path.join(align_path, file)
        result_file = compress_result_file.rstrip(".gz")
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
                    top_folder,
                    args.chomp_max_distance,
                    is_genome,
                ),
            ),
        )

    if args.processes > 1:
        with Pool(args.processes) as pool:
            cluster_logs = pool.starmap(run_command, func_args, chunksize=1)
    else:
        cluster_logs = [run_command(arg[0]) for arg in func_args]

    if any(cluster_logs):
        cluster_logs.sort(key=lambda x: x[0])
        cluster_logs = [x[1] for x in cluster_logs]
            
        with open(path.join(folder, "align_clusters.csv"), "w") as f:
            f.write("Gene,Seq count,Cluster count,Cluster ranges\n")
            f.write("\n".join(cluster_logs))

    printv(f"Done! Took {time_keeper.differential():.2f}s overall", args.verbose, 0)
    return True


def main(args):
    if not all(path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    if path.exists("intermediates"):
        rmtree("intermediates")
    if args.debug:
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