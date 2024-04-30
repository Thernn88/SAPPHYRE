from collections import defaultdict
from itertools import combinations
from math import ceil
from multiprocessing import Pool
from os import listdir, makedirs, mkdir, path
from shutil import rmtree

from msgspec import Struct
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    # del_cols,
    OverlapTree,
)
from wrap_rocks import RocksDB
from msgspec import json

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

def do_cluster(ids, ref_coords, max_distance=100):
    clusters = []
    ids.sort(key = lambda x: x[0])

    req_seq_coverage = 0.5

    current_cluster = []
    for i, (child_index, seq_coords) in enumerate(ids):

        coverage = len(seq_coords.intersection(ref_coords)) / len(ref_coords)


        if not current_cluster:
            current_cluster.append((child_index, coverage, i))
            current_index = child_index
        else:
            if child_index - current_index <= max_distance:
                current_cluster.append((child_index, coverage, i))
                current_index = child_index
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
                        
                current_cluster = [(child_index, coverage, i)]
                current_index = child_index

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

def within_highest_coverage(logs, clusters, within=0.1): # kick 10% less coverage
    highest_coverage_cluster = max([x[2] for x in clusters])
    req_coverage = highest_coverage_cluster - within

    within = []
    for cluster in clusters:
        if cluster[2] >= req_coverage:
            within.append(cluster)
        else:
            logs.append(f"Kicked due to not within 10% of max coverage:\n{cluster[0]}-{cluster[1]}\nCluster coverage: {cluster[2]:.2f}\nHighest coverage: {highest_coverage_cluster:.2f}")

    return within


def do_gene(gene: str, aa_gene_input_path: str, nt_gene_input_path: str, aa_gene_output_path: str, nt_gene_output_path: str):
    reference_data_cols = set()
    ids = []
    logs = []
    cluster_consensi = {}

    get_id = lambda header: int(header.split("|")[3].split("_")[1])
    ref_consensus = defaultdict(set)
    sequences = defaultdict(list)

    for header, seq in parseFasta(path.join(aa_gene_input_path, gene)):
        if header.endswith("."):
            start, end = find_index_pair(seq, "-")
            for i,let in enumerate(seq[start:end], start):
                if let != "-":
                    ref_consensus[i].add(let)
                    reference_data_cols.add(i)
                    
            continue

        this_id = get_id(header)
        start, end = find_index_pair(seq, "-")
        data_cols = {i for i, let in enumerate(seq[start:end], start) if let != "-"}

        ids.append((this_id, data_cols))
        sequences[this_id].append((header, seq))

    clusters = do_cluster(ids, reference_data_cols)

    if not clusters:
        if logs:
            logs.insert(0, f"Log for {gene}")
        return logs
    

    clusters = within_highest_coverage(logs, clusters)

    for i, cluster in enumerate(clusters):
        cluster_consensi[i] = None
   
    if len(clusters) >= 2: 
        cluster_percents = []

        for x, cluster in enumerate(clusters):
            this_cluster_sequences = []
            for i in range(cluster[0], cluster[1] + 1):
                if i in sequences:
                    for seq in sequences[i]:
                        this_cluster_sequences.append(seq[1])

            cluster_consensus = defaultdict(set)
            for seq in this_cluster_sequences:
                start, end = find_index_pair(seq, "-")
                
                for i, let in enumerate(seq[start:end], start):
                    if let != "-":
                        cluster_consensus[i].add(let)

            cluster_match = 0
            cluster_cols = 0
            for i, clets in cluster_consensus.items():
                cluster_cols += 1
                if clets.intersection(ref_consensus[i]):
                    cluster_match += 1
            
            cluster_consensi[x] = cluster_consensus
            cluster_percents.append((x, cluster_match / cluster_cols))

        max_cluster = max([i[1] for i in cluster_percents])

        for i, cluster_percent in cluster_percents:
            if abs(cluster_percent - max_cluster) > 0.1: # 10% difference to reference matching percent
                logs.append(f"Kicked due to not within 10% of top reference match:\n{clusters[i][0]}-{clusters[i][1]}\nCluster coverage: {clusters[i][2]:.2f}\nHighest matching cluster to references: {max_cluster:.2f}\nThis cluster matched by {cluster_percent:.2f}")
                cluster_consensi.pop(i)

        cluster_percents = [(index, percent) for index, percent in cluster_percents if index in cluster_consensi]
        cluster_percents.sort(key = lambda x: x[1], reverse = True)

   
        if len(cluster_percents) > 1:
            kicked = True
            while kicked:
                kicked = False
                cluster_percents = [(index, percent) for index, percent in cluster_percents if index in cluster_consensi]
                for cluster_a, cluster_b in combinations(cluster_percents, 2):
                    consensus_a = cluster_consensi[cluster_a[0]]
                    consensus_b = cluster_consensi[cluster_b[0]]

                    start_a = min(consensus_a.keys())
                    start_b = min(consensus_b.keys())

                    end_a = max(consensus_a.keys())
                    end_b = max(consensus_b.keys())

                    end_start = max(start_a, start_b)
                    first_end = min(end_a, end_b)

                    matching = 0
                    checked = 0
                    for i in range(end_start, first_end + 1):
                        checked += 1

                        if consensus_a[i].intersection(consensus_b[i]):
                            matching += 1

                    matching_percent = matching / checked

                    if matching_percent <= (1 - 0.1): #Matching percent
                        logs.append(f"Kicked due to cluster compare:\n{clusters[cluster_b[0]][0]}-{clusters[cluster_b[0]][1]}\nCluster coverage: {clusters[cluster_b[0]][2]:.2f}\nKicked by:\n{clusters[cluster_a[0]][0]}-{clusters[cluster_a[0]][1]}\nCluster coverage: {clusters[cluster_a[0]][2]:.2f}\nClusters matched by {matching_percent:.2f}")
                        cluster_consensi.pop(cluster_b[0])

                        kicked = True
                        break
    allowed = set()
    for i in cluster_consensi:
        cluster = clusters[i]
        id_start = cluster[0]
        id_end = cluster[1]

        for i in range(id_start, id_end + 1):
            allowed.add(i)

    aa_out = []
    nt_out = []

    for header, seq in parseFasta(path.join(aa_gene_input_path, gene)):
        if header.endswith("."):
            aa_out.append((header, seq))
            continue

        if get_id(header) in allowed:
            aa_out.append((header, seq))

    for header, seq in parseFasta(path.join(nt_gene_input_path, gene.replace('.aa.','.nt.'))):
        if header.endswith("."):
            nt_out.append((header, seq))
            continue

        if get_id(header) in allowed:
            nt_out.append((header, seq))

    writeFasta(path.join(aa_gene_output_path, gene), aa_out)
    writeFasta(path.join(nt_gene_output_path, gene.replace('.aa.','.nt.')), nt_out)

    if logs:
        logs.insert(0, f"Log for {gene}")
    return logs

def do_folder(args, folder: str, from_folder: str):
    arguments = []
    gene_input_path = path.join(folder, "outlier", from_folder)

    aa_gene_input_path = path.join(gene_input_path, "aa")
    nt_gene_input_path = path.join(gene_input_path, "nt")

    gene_output_path = path.join(folder, "outlier", "clusters")

    aa_gene_output_path = path.join(gene_output_path, "aa")
    nt_gene_output_path = path.join(gene_output_path, "nt")

    makedirs(aa_gene_output_path, exist_ok=True)
    makedirs(nt_gene_output_path, exist_ok=True)

    for gene in listdir(aa_gene_input_path):
        arguments.append((gene, aa_gene_input_path, nt_gene_input_path, aa_gene_output_path, nt_gene_output_path))
    

    if args.processes >= 1:
        with Pool(args.processes) as p:
            logs = p.starmap(do_gene, arguments)
    else:
        logs = []
        for arg in arguments:
            logs.append(do_gene(*arg))

    logs = [log for log in logs if log]
    log_output = path.join(gene_output_path, "logs.txt")
    with open(log_output, "w") as f:
        f.write("Cluster Consensus Logs:\n")
        for log in logs:
            f.write("\n\n".join(log) + "\n#############\n")


    return True


def main(args, from_folder):
    return do_folder(args, args.INPUT, from_folder)


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )