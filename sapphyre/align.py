from __future__ import annotations
from itertools import combinations
from math import ceil
import os
from collections import defaultdict, namedtuple
from multiprocessing.pool import Pool
from shutil import rmtree
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
import numpy as np
from .utils import printv, gettempdir, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode

KMER_LEN = 12
KMER_PERCENT = 0.70
SUBCLUSTER_AT = 1000
CLUSTER_EVERY = 500  # Aim for x seqs per cluster
SAFEGUARD_BP = 15000
SINGLETON_THRESHOLD = 5
IDENTITY_THRESHOLD = 80
SUBCLUSTER_AMOUNT = 2

def find_kmers(fasta):
    kmers = {}
    for header, sequence in fasta.items():
        kmers[header] = set()
        for i in range(0, len(sequence) - KMER_LEN):
            kmer = sequence[i : i + KMER_LEN]
            if "*" not in kmer and "-" not in kmer:
                kmers[header].add(kmer)
    return kmers


ALIGN_FOLDER = "align"
AA_FOLDER = "aa"
ALN_FOLDER = "aln"


def get_start(sequence):
    for i, char in enumerate(sequence):
        if char != "-":
            return i
    return -1


def process_genefile(fileread):
    data = {}
    targets = {}
    reinsertions = {}
    seq_to_first_header = {}
    trimmed_header_to_full = {}
    for header, sequence in parseFasta(fileread):
        if not header.endswith("."):
            seq_hash = hash(sequence)
            trimmed_header_to_full[header[:127]] = header
            if seq_hash not in seq_to_first_header:
                seq_to_first_header[seq_hash] = header
            else:
                reinsertions.setdefault(seq_to_first_header[seq_hash], []).append(
                    header
                )
                continue
            data[header] = sequence.replace("-", "")
        else:
            targets[header.split("|")[2]] = header

    return len(data), data, targets, reinsertions, trimmed_header_to_full


def get_identity(aligned_cluster):
    identity_out = subprocess.run(
        f"identity/identity {aligned_cluster}",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True
    )
    identity_out = identity_out.stdout.decode("utf-8")
    identity = float(identity_out.replace("Average pairwise identity: ","").replace("%",""))
    return identity

def subcluster(aligned_sequences):
    sequences = [np.array(list(seq)) for _, seq in aligned_sequences]
    distances = np.zeros((len(sequences), len(sequences)))
    for i, j in combinations(range(len(sequences)), 2):

        distances[i, j] = np.count_nonzero(sequences[i] != sequences[j])
        distances[j, i] = distances[i, j]

    # Perform hierarchical clustering using complete linkage
    Z = linkage(squareform(distances), method='complete')

    # Cut the dendrogram to create subclusters
    subclusters = fcluster(Z, SUBCLUSTER_AMOUNT, criterion='maxclust')

    

    # Save each subcluster to a separate FASTA file
    for i in range(1, SUBCLUSTER_AMOUNT+1):
        subcluster_indices = [j for j, c in enumerate(subclusters) if c == i]
        subcluster_records = [aligned_sequences[j] for j in subcluster_indices]
        yield subcluster_records


def delete_empty_cols(records):
     #Delete empty cols
    cols_to_keep = set()
    output = []
    for _, sequence in records:
        for x,let in enumerate(sequence):
            if x not in cols_to_keep:
                if let != "-":
                    cols_to_keep.add(x)
    for header, sequence in records:
        sequence = [let for x, let in enumerate(sequence) if x in cols_to_keep]
        output.append((header, "".join(sequence)))

    return output


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
    ],
)


def run_command(args: CmdArgs) -> None:
    keeper = TimeKeeper(KeeperMode.DIRECT)
    debug = args.debug
    printv(f"Doing: {args.gene} ", args.verbose, 2)
    # print(args.only_singletons)
    temp_dir = gettempdir()

    aligned_ingredients = []

    if debug:
        this_intermediates = os.path.join("intermediates", args.gene)
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
            clusters = []
            cluster_children = {header: [header] for header in data}
            kmers = find_kmers(data)
            gene_headers = list(kmers.keys())
            merge_occured = True
            while merge_occured:
                merge_occured = False
                for i in range(len(gene_headers) - 1, -1, -1):  # reverse iteration
                    master_header = gene_headers[i]
                    master = kmers[master_header]
                    if master:
                        for candidate_header in gene_headers[:i]:
                            candidate = kmers[candidate_header]
                            if candidate:
                                similar = master.intersection(candidate)
                                if len(similar) != 0:
                                    if (
                                        len(similar) / min(len(master), len(candidate))
                                        >= KMER_PERCENT
                                    ):
                                        # Merge kmeans and merge cluster headers
                                        master.update(candidate)
                                        cluster_children[master_header].extend(cluster_children[candidate_header])

                                        # Remove candidate
                                        cluster_children[candidate_header] = None
                                        kmers[candidate_header] = None

                                        

                                        merge_occured = True

            for this_cluster in cluster_children.values():
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
                            sig_out = subprocess.run(
                                f"SigClust/SigClust -c {clusters_to_create} {this_tmp.name}",
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=True
                            )
                            sig_out = sig_out.stdout.decode("utf-8")
                            sub_clusters = defaultdict(list)
                            for line in sig_out.split("\n"):
                                line = line.strip()
                                if line:
                                    seq_index, clust_index = line.strip().split(",")
                                    sub_clusters[clust_index].append(
                                        this_cluster[int(seq_index)]
                                    )
                            clusters = list(sub_clusters.values())
                            continue
                    clusters.append(this_cluster)
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

            cluster_i = 0
            for cluster in clusters:
                
                cluster_i += 1
                cluster_seqs = [(header, data[header]) for header in cluster]

                if args.only_singletons or len(cluster_seqs) < SINGLETON_THRESHOLD:
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
                aligned_cluster = os.path.join(
                    aligned_files_tmp, this_clus_align
                )

                raw_cluster = os.path.join(raw_files_tmp, f"{args.gene}_cluster{cluster_i}")
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
                        aligned_sequences
                    )
                identity = get_identity(aligned_cluster)
                if identity <= IDENTITY_THRESHOLD:
                    printv(f"{args.gene} cluster {cluster_i} has identity {identity}. Subclustering", args.verbose, 1)
                    # Calculate the pairwise distances between sequences using Hamming distance
                    subcluster_i = 0
                    to_subcluster = [aligned_sequences]
                    done = set()
                    while to_subcluster:
                        this_sequences = to_subcluster[0]
                        removed = False
                        for subcluster_records in subcluster(this_sequences):
                            subcluster_i += 1
                            if subcluster_records:
                                output = delete_empty_cols(subcluster_records)

                                this_id = hash("".join([i[0] for i in output]))
                                

                                cluster = f"{args.gene}_cluster_{cluster_i}_length_{len(cluster)}_subcluster{subcluster_i}_aligned"
                                aligned_cluster = os.path.join(
                                    aligned_files_tmp, cluster
                                )
                                if debug:
                                    writeFasta(
                                        os.path.join(
                                            this_intermediates, cluster
                                        ),
                                        output
                                    )
                                

                                writeFasta(aligned_cluster, output)
                                if len(output) <= 1:
                                    identity = 100.0
                                else:
                                    identity = get_identity(aligned_cluster)
                                print(aligned_cluster, "has identity",identity )
                                if identity <= IDENTITY_THRESHOLD and this_id not in done:
                                    to_subcluster.append(output)
                                else:
                                    if len(output) < SINGLETON_THRESHOLD:
                                        print(aligned_cluster, "is a singleton")
                                        to_merge.extend(output)
                                    else:
                                        aligned_ingredients.append((aligned_cluster, len(output)))

                                if not removed:
                                    to_subcluster.pop(0)
                                    removed = True

                                done.add(this_id)
                            else:
                                print(f"Subcluster {i} for cluster {cluster_i} is empty")
                else:
                    aligned_ingredients.append((aligned_cluster, len(aligned_sequences)))


        align_time = keeper.differential() - cluster_time
        if to_merge:
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
        if aligned_ingredients or to_merge:
            aligned_ingredients = [i[0] for i in sorted(aligned_ingredients, key = lambda x: x[1])]
            if to_merge:
                aligned_ingredients.insert(0, merged_singleton_final)
            printv(
                f"Merging Alignments. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            printv(
                f"Merging 1 of {len(aligned_ingredients)}. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            with NamedTemporaryFile(dir=parent_tmpdir, mode="w+") as tmp:
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
                                [
                                    let
                                    for col, let in enumerate(sequence)
                                    if not empty_columns[col]
                                ]
                            ),
                        )
                    )

                writeFasta(tmp.name, to_write)
                if debug:
                    writeFasta(
                        os.path.join(this_intermediates, "references.fa"), to_write
                    )
                tmp.flush()

                out_file = os.path.join(parent_tmpdir, "part_0.fa")
                if len(aligned_ingredients) == 1:
                    out_file = args.result_file

                if to_merge:
                    os.system(
                        f"mafft --anysymbol --jtt 1 --quiet --addfragments {aligned_ingredients[0]} --thread 1 {tmp.name} > {out_file}"
                    )
                    if args.debug:
                        printv(
                            f"mafft --anysymbol --jtt 1 --quiet --addfragments {aligned_ingredients[0]} --thread 1 {tmp.name} > {out_file}", args.verbose, 3
                        )
                else:
                    os.system(
                        f"clustalo --p1 {tmp.name} --p2 {aligned_ingredients[0]} -o {out_file} --threads=1 --iter=3 --full --full-iter --is-profile --force"
                    )
                    if debug:
                        printv(
                            f"clustalo --p1 {tmp.name} --p2 {aligned_ingredients[0]} -o {out_file} --threads=1 --iter=3 --full --full-iter --is-profile --force",
                            args.verbose,
                            3,
                        )
                        writeFasta(
                            os.path.join(
                                this_intermediates, os.path.basename(out_file)
                            ),
                            parseFasta(out_file, True),
                        )

            for i, file in enumerate(aligned_ingredients[1:]):
                printv(
                    f"Merging {i+2} of {len(aligned_ingredients)}. Elapsed time: {keeper.differential():.2f}",
                    args.verbose,
                    3,
                )  # Debug
                out_file = os.path.join(parent_tmpdir, f"part_{i+1}.fa")
                in_file = os.path.join(parent_tmpdir, f"part_{i}.fa")
                if file == aligned_ingredients[-1]:
                    out_file = args.result_file

                os.system(
                    f"clustalo --p1 {in_file} --p2 {file} -o {out_file} --threads=1 --iter=3 --full --full-iter --is-profile --force"
                )
                if debug:
                    printv(
                        f"clustalo --p1 {in_file} --p2 {file} -o {out_file} --threads=1 --iter=3 --full --full-iter --is-profile --force",
                        args.verbose,
                        3,
                    )
                    writeFasta(
                        os.path.join(this_intermediates, os.path.basename(out_file)),
                        parseFasta(out_file, True),
                    )
    merge_time = keeper.differential() - align_time - cluster_time

    # Reinsert and sort
    to_write = []
    references = []
    for header, sequence in parseFasta(args.result_file, True):
        if header.endswith("."):
            references.append((header, sequence))
        else:
            header = trimmed_header_to_full[header[:127]]
            if header in reinsertions:
                for insertion_header in reinsertions[header]:
                    to_write.append((insertion_header, sequence))
            to_write.append((header, sequence))

    to_write.sort(key=lambda x: get_start(x[1]))

    writeFasta(args.result_file, references + to_write)

    printv(f"Done. Took {keeper.differential():.2f}", args.verbose, 3)  # Debug

    return args.gene, cluster_time, align_time, merge_time, keeper.differential()


def do_folder(folder, args):
    printv(f"Processing: {os.path.basename(folder)}", args.verbose)
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
            os.path.join(aln_path, gene.split(".")[0] + ".aln.fa"),
            True
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

    command = f"clustalo -i {{in_file}} -o {{out_file}} --threads=1 --full --iter=1 --full-iter"
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
            (CmdArgs(
                    command,
                    gene_file,
                    result_file,
                    gene,
                    args.verbose,
                    args.compress,
                    aln_path,
                    args.debug,
                    file in only_singletons,
                ),)
            )
        
    if args.processes > 1:
        with Pool(args.processes) as pool:
            times = pool.starmap(run_command, func_args, chunksize=1)
    else:
        times = [run_command(arg[0]) for arg in func_args]

    # if args.debug:
    # with open("mafft_times.csv", "w") as fp:
    # fp.write("\n".join(["Gene,Clustering,Aligning,Merging,Total"]+[",".join(map(str, i)) for i in times]))

    printv(f"Done! Took {time_keeper.differential():.2f}s overall", args.verbose)
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
