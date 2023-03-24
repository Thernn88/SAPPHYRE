from __future__ import annotations
from math import ceil
import os
from collections import namedtuple
from multiprocessing.pool import ThreadPool
from shutil import rmtree
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
from threading import Lock
from .utils import printv, gettempdir, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode

KMER_LEN = 7
KMER_PERCENT = 0.55

SUBCLUSTER_AT = 500
CLUSTER_EVERY = 250 # Aim for x seqs per cluster

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


CmdArgs = namedtuple(
    "CmdArgs",
    [
        "string",
        "gene_file",
        "result_file",
        "gene",
        "lock",
        "verbose",
        "compress",
        "aln_path",
        "debug",
    ],
)


def run_command(args: CmdArgs) -> None:
    keeper = TimeKeeper(KeeperMode.DIRECT)
    debug = args.debug
    if args.lock is not None:
        with args.lock:
            printv(f"Doing: {args.gene}", args.verbose, 2)
    else:
        printv(f"Doing: {args.gene}", args.verbose, 2)
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
            writeFasta(aligned_file, data.items())
            if args.debug:
                writeFasta(os.path.join(this_intermediates, aligned_file), data.items())
            aligned_ingredients.append(aligned_file)
        else:
            printv(
                f"Generating Cluster. Elapsed time: {keeper.differential():.2f}",
                args.verbose,
                3,
            )  # Debug
            clusters = []
            cluster_children = {header: [] for header in data}
            kmers = find_kmers(data)
            master_headers = list(kmers.keys())
            for i in range(len(master_headers) - 1, -1, -1):  # reverse iteration
                master_header = master_headers[i]
                master = kmers[master_header]
                if master:
                    for key in master_headers[:i]:
                        candidate = kmers[key]
                        if candidate:
                            similar = master.intersection(candidate)
                            if len(similar) != 0:
                                if (
                                    len(similar) / min(len(master), len(candidate))
                                    >= KMER_PERCENT
                                ):
                                    candidate.update(master)
                                    children = [
                                        master_header.strip(">")
                                    ] + cluster_children[master_header.strip(">")]
                                    cluster_children[key.strip(">")].extend(
                                        children
                                    )
                                    cluster_children[
                                        master_header.strip(">")
                                    ] = None
                                    kmers[master_header] = None
                                    break

            for master, children in cluster_children.items():
                if children is not None:
                    this_cluster = [master] + children
                    if len(this_cluster) > SUBCLUSTER_AT:
                        clusters_to_create = ceil(len(this_cluster) / CLUSTER_EVERY)
                        with NamedTemporaryFile(mode="w+", dir=parent_tmpdir, suffix=".fa") as this_tmp:
                            writeFasta(this_tmp.name, [(header, data[header]) for header in this_cluster])  
                            sig_out = subprocess.run(f"SigClust/SigClust -c {clusters_to_create} {this_tmp.name}", shell=True, stdout=subprocess.PIPE)
                            sig_out = sig_out.stdout.decode("utf-8")
                            sub_clusters = {}
                            for line in sig_out.split("\n"):
                                line = line.strip()
                                if line:
                                    seq_index, clust_index = line.strip().split(',')
                                    sub_clusters.setdefault(clust_index, []).append(this_cluster[int(seq_index)])
                            for sub_cluster in sub_clusters.values():
                                clusters.append(sub_cluster)
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

            for i, cluster in enumerate(clusters):
                printv(
                    f"Aligning cluster {i}. Elapsed time: {keeper.differential():.2f}",
                    args.verbose,
                    3,
                )  # Debug
                aligned_cluster = os.path.join(
                    aligned_files_tmp, f"{args.gene}_cluster_{len(cluster)}_{i}"
                )

                cluster_seqs = [
                    (header, data[header]) for header in cluster
                ]

                if len(cluster_seqs) == 1:
                    to_merge.extend(cluster_seqs)

                else:
                    aligned_ingredients.append(aligned_cluster)
                    raw_cluster = os.path.join(raw_files_tmp, f"{args.gene}_cluster{i}")
                    writeFasta(raw_cluster, cluster_seqs)
                    if debug:
                        writeFasta(
                            os.path.join(this_intermediates, f"{args.gene}_cluster{i}"),
                            cluster_seqs,
                        )

                    command = args.string.format(
                        in_file=raw_cluster, out_file=aligned_cluster
                    )

                    os.system(command)
                    if debug:
                        printv(command, args.verbose, 3)
                        writeFasta(
                            os.path.join(
                                this_intermediates, f"{args.gene}_cluster{i}_aligned"
                            ),
                            parseFasta(aligned_cluster),
                        )
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
            aligned_ingredients.sort(key=lambda x: int(x.split("_")[-2]), reverse=True)
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
                for header, sequence in parseFasta(aln_file):
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
                        print(
                            f"mafft --anysymbol --jtt 1 --quiet --addfragments {aligned_ingredients[0]} --thread 1 {tmp.name} > {out_file}"
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
                            parseFasta(out_file),
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
                        parseFasta(out_file),
                    )
    merge_time = keeper.differential() - align_time - cluster_time

    # Reinsert and sort
    to_write = []
    references = []
    for header, sequence in parseFasta(args.result_file):
        if header.endswith("."):
            references.append((header, sequence))
        else:
            header = trimmed_header_to_full[header]
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
        gene
        for gene in os.listdir(aa_path)
        if gene.split(".")[-1] in ["fa", "gz", "fq", "fastq", "fasta"]
    ]
    orthoset_path = os.path.join(args.orthoset_input, args.orthoset)
    aln_path = os.path.join(orthoset_path, ALN_FOLDER)
    if not os.path.exists(orthoset_path):
        printv("ERROR: Orthoset path not found.", args.verbose, 0)
        return False
    if not os.path.exists(aln_path):
        printv("ERROR: Aln folder not found.", args.verbose, 0)
        return False

    command = f"clustalo -i {{in_file}} -o {{out_file}} --threads=1 --full --iter=1 --full-iter"
    #command = f"mafft --maxiterate 2 --anysymbol --quiet --thread 1 {{in_file}} > {{out_file}}"

    intermediates = "intermediates"
    if not os.path.exists(intermediates):
        os.mkdir(intermediates)

    if args.processes > 1:
        arguments = []
        func = arguments.append
        lock = Lock()
    else:
        func = run_command
        lock = None

    times = []
    for file in genes:
        gene = file.split(".")[0]
        gene_file = os.path.join(aa_path, file)
        result_file = os.path.join(align_path, file.rstrip(".gz"))
        if func == run_command:
            times.append(
                run_command(
                    CmdArgs(
                        command,
                        gene_file,
                        result_file,
                        gene,
                        lock,
                        args.verbose,
                        args.compress,
                        aln_path,
                        args.debug,
                    )
                )
            )
        else:
            func(
                (
                    CmdArgs(
                        command,
                        gene_file,
                        result_file,
                        gene,
                        lock,
                        args.verbose,
                        args.compress,
                        aln_path,
                        args.debug,
                    ),
                )
            )

    if args.processes > 1:
        with ThreadPool(args.processes) as pool:
            times = pool.starmap(run_command, arguments, chunksize=1)

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
