from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from os import listdir, makedirs, path
from msgspec import Struct
from sapphyre_tools import (
    find_index_pair,
)
from .directional_cluster import cluster_ids, quick_rec
from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

class cluster_obj(Struct):
    cluster_set: set
    start: int
    end: int
    coverage: float
    seq_data: list
    seq_data_coords: set
    

def within_highest_coverage(logs, clusters, overlap_perc, within=0.1): # kick 10% less coverage
    highest_cluster = max(clusters, key = lambda x: x.coverage)
    highest_coverage_cluster = highest_cluster.coverage

    highest_data_coords = highest_cluster.seq_data_coords

    req_coverage = highest_coverage_cluster - within

    within = []
    for cluster in clusters:

        cluster_data_coords = cluster.seq_data_coords
        overlap = highest_data_coords.intersection(cluster_data_coords)

        percent = 0
        if overlap:
            amount = len(overlap)
            start, end = min(cluster_data_coords), max(cluster_data_coords)
            highest_start, highest_end = min(highest_data_coords), max(highest_data_coords)
            percent = amount / (min((end-start), (highest_end-highest_start)) + 1)

        if percent >= overlap_perc:
            if cluster.coverage < req_coverage:
                logs.append(f"Kicked due to not within 10% of max coverage:\n{cluster.start}-{cluster.end}\nOverlap percent: {percent:.2f}\nCluster coverage: {cluster.coverage:.2f}\nHighest coverage: {highest_coverage_cluster:.2f}")
                continue
        
        within.append(cluster)

    return within


def do_gene(gene: str, aa_gene_input_path: str, nt_gene_input_path: str, aa_gene_output_path: str, nt_gene_output_path: str, verbose: int, cluster_overlap_requirement: float):
    printv(f"Processing {gene}", verbose, 2)
    reference_data_cols = set()
    ids = []
    logs = []
    cluster_consensi = {}

    ref_consensus = defaultdict(set)
    sequences = []

    for header, seq in parseFasta(path.join(aa_gene_input_path, gene)):
        if header.endswith("."):
            start, end = find_index_pair(seq, "-")
            for i,let in enumerate(seq[start:end], start):
                if let != "-":
                    ref_consensus[i].add(let)
                    reference_data_cols.add(i)
                    
            continue

        start, end = find_index_pair(seq, "-")

        frame = int(header.split("|")[4])
        ids.append(quick_rec(header.split("|")[3], frame, seq, start, end))
        sequences.append((header, seq))
        msa_length = len(seq)

    max_gap_size = round(msa_length * 0.3)
    clusters, _ = cluster_ids(ids, 100, max_gap_size, reference_data_cols, req_seq_coverage=0)

    final_clusters = []
    for cluster_set, cluster_range, cluster_coverage in clusters:
        this_cluster_seq = []

        cluster_min_start = None
        cluster_min_end = None
        for header, seq in sequences:
            if header.split("|")[3] in cluster_set:
                this_cluster_seq.append((header, seq, start, end))
                start, end = find_index_pair(seq, "-")
                if cluster_min_start is None:
                    cluster_min_start = start
                    cluster_min_end = end
                else:
                    cluster_min_start = min(cluster_min_start, start)
                    cluster_min_end = max(cluster_min_end, end)
                    
        this_cluster_coords = set(range(cluster_min_start, cluster_min_end))
        
        final_clusters.append(cluster_obj(cluster_set, cluster_range[0], cluster_range[1], cluster_coverage, this_cluster_seq, this_cluster_coords))

    if not clusters:
        if logs:
            logs.insert(0, f"Log for {gene}")
        return logs
    

    final_clusters = within_highest_coverage(logs, final_clusters, cluster_overlap_requirement)

    for i, cluster in enumerate(final_clusters):
        cluster_consensi[i] = None
   
    if len(final_clusters) >= 2: 
        cluster_percents = []

        for x, cluster in enumerate(final_clusters):
            cluster_consensus = defaultdict(set)
            for _, seq, start, end in cluster.seq_data:
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
            percent = 0 if cluster_cols == 0 else cluster_match / cluster_cols
            cluster_percents.append((x, cluster.seq_data_coords, percent))

        max_cluster = max(cluster_percents, key = lambda x: x[2])
        max_percent = max_cluster[2]
        max_data_coords = max_cluster[1]
        max_end = max(max_data_coords)
        max_start = min(max_data_coords)

        for i, data_coords, cluster_percent in cluster_percents:
            cluster_overlap = max_data_coords.intersection(data_coords)
            if not cluster_overlap:
                continue

            amount = len(cluster_overlap)
            percent = amount / min((end - start), (max_end - max_start))

            if percent < cluster_overlap_requirement:
                continue

            if abs(cluster_percent - max_percent) > 0.1: # 10% difference to reference matching percent
                logs.append(f"Kicked due to not within 10% of top reference match:\n{clusters[i][0]}-{clusters[i][1]}\nCluster coverage: {clusters[i][2]:.2f}\nHighest matching cluster to references: {max_percent:.2f}\nThis cluster matched by {cluster_percent:.2f}")
                cluster_consensi.pop(i)

        cluster_percents = [i for i in cluster_percents if i[0] in cluster_consensi]
        cluster_percents.sort(key = lambda x: x[1], reverse = True)

   
        if len(cluster_percents) > 1:
            kicked = True
            while kicked:
                kicked = False
                cluster_percents = [i for i in cluster_percents if i[0] in cluster_consensi]
                for cluster_a, cluster_b in combinations(cluster_percents, 2):
                    consensus_a = cluster_consensi[cluster_a[0]]
                    consensus_b = cluster_consensi[cluster_b[0]]

                    data_cols_a = cluster_a[1]

                    data_cols_b = cluster_b[1]
  
                    clusters_overlap = data_cols_a.intersection(data_cols_b)
                    if not clusters_overlap:
                        continue

                    start_a = min(data_cols_a)
                    end_a = max(data_cols_a)
                    
                    start_b = min(data_cols_b)
                    end_b = max(data_cols_b)

                    amount_overlap = len(clusters_overlap)
                    cluster_overlap_percent = amount_overlap / min((end_a - start_a), (end_b - start_b))

                    if cluster_overlap_percent < cluster_overlap_requirement:
                        continue

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
        cluster = final_clusters[i]
        for header, _, _, _ in cluster.seq_data:
            allowed.add(header)

    aa_out = []
    nt_out = []

    for header, seq in parseFasta(path.join(aa_gene_input_path, gene)):
        if header.endswith("."):
            aa_out.append((header, seq))
            continue

        if header in allowed:
            aa_out.append((header, seq))

    for header, seq in parseFasta(path.join(nt_gene_input_path, gene.replace('.aa.','.nt.'))):
        if header.endswith("."):
            nt_out.append((header, seq))
            continue

        if header in allowed:
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
        arguments.append((gene, aa_gene_input_path, nt_gene_input_path, aa_gene_output_path, nt_gene_output_path, args.verbose, args.cluster_overlap_requirement))
    

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
    global_time = TimeKeeper(KeeperMode.DIRECT)
    # print(f"Processing: {path.basename(args.INPUT)}")
    success = do_folder(args, args.INPUT, from_folder)
    printv(f"Done! Took {global_time.differential():.2f} seconds", args.verbose, 0)
    return success


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )