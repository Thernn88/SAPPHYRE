from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from os import listdir, makedirs, path
from sapphyre_tools import (
    find_index_pair,
    get_overlap,
)

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

def do_cluster(ids, ref_coords, max_distance=100):
    clusters = []
    ids.sort(key = lambda x: x[0])

    req_seq_coverage = 0

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

                    start = min(cluster_data_cols)
                    end = max(cluster_data_cols)
                        
                    cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

                    clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage, start, end))
                elif len(current_cluster) == 1:
                    if current_cluster[0][1] > req_seq_coverage:
                        cluster_coverage = current_cluster[0][1]

                        data_cols = ids[current_cluster[0][2]][1]
                        start = min(data_cols)
                        end = max(data_cols)

                        clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage, start, end))
                        
                current_cluster = [(child_index, coverage, i)]
                current_index = child_index

    if current_cluster:
        if len(current_cluster) >= 2:
            cluster_data_cols = set()
            for _, _, index in current_cluster:
                cluster_data_cols.update(ids[index][1])

            start = min(cluster_data_cols)
            end = max(cluster_data_cols)
                
            cluster_coverage = len(cluster_data_cols.intersection(ref_coords)) / len(ref_coords)

            clusters.append((current_cluster[0][0], current_cluster[-1][0], cluster_coverage, start, end))
        elif len(current_cluster) == 1:
            if current_cluster[0][1] > req_seq_coverage:
                cluster_coverage = current_cluster[0][1]

                data_cols = ids[current_cluster[0][2]][1]
                start = min(data_cols)
                end = max(data_cols)
                
                clusters.append((current_cluster[0][0], current_cluster[0][0], cluster_coverage, start, end))
                
    return clusters

def within_highest_coverage(logs, clusters, overlap_perc, within=0.1): # kick 10% less coverage
    highest_cluster = max(clusters, key = lambda x: x[2])
    highest_coverage_cluster = highest_cluster[2]

    highest_start = highest_cluster[3]
    highest_end = highest_cluster[4]

    req_coverage = highest_coverage_cluster - within

    within = []
    for cluster in clusters:

        start = cluster[3]
        end = cluster[4]

        overlap = get_overlap(start, end, highest_start, highest_end, 0)
        percent = 0
        if overlap:
            amount = (overlap[1] - overlap[0])
            percent = amount / min((end-start), (highest_end-highest_start))

        if percent >= overlap_perc:
            if cluster[2] < req_coverage:
                logs.append(f"Kicked due to not within 10% of max coverage:\n{cluster[0]}-{cluster[1]}\nCluster coverage: {cluster[2]:.2f}\nHighest coverage: {highest_coverage_cluster:.2f}")
                continue
        
        within.append(cluster)

    return within


def do_gene(gene: str, aa_gene_input_path: str, nt_gene_input_path: str, aa_gene_output_path: str, nt_gene_output_path: str, verbose: int, cluster_overlap_requirement: float):
    printv(f"Processing {gene}", verbose, 2)
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
    

    clusters = within_highest_coverage(logs, clusters, cluster_overlap_requirement)

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
            start = min(cluster_consensus.keys())
            end = max(cluster_consensus.keys())
            cluster_percents.append((x, start, end, cluster_match / cluster_cols))

        max_cluster = max(cluster_percents, key = lambda x: x[3])
        max_percent = max_cluster[3]
        max_start, max_end = max_cluster[1], max_cluster[2]

        for i, start, end, cluster_percent in cluster_percents:
            cluster_overlap = get_overlap(start, end, max_start, max_end, 0)
            if not cluster_overlap:
                continue

            amount = cluster_overlap[1] - cluster_overlap[0]
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

                    start_a = cluster_a[1]
                    end_a = cluster_a[2]

                    start_b = cluster_b[1]
                    end_b = cluster_b[2]

                    clusters_overlap = get_overlap(start_a, end_a, start_b, end_b, 0)
                    if not clusters_overlap:
                        continue

                    amount_overlap = clusters_overlap[1] - clusters_overlap[0]
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
