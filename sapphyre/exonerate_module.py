from itertools import combinations
from math import ceil
from pathlib import Path
from shutil import rmtree
import subprocess
from tempfile import NamedTemporaryFile
from wrap_rocks import RocksDB
from os import listdir, mkdir, path
from .directional_cluster import cluster_ids, within_distance, node_to_ids, quick_rec
from .utils import gettempdir, parseFasta, printv, writeFasta
from sapphyre_tools import (
    find_index_pair,
    get_overlap,
)
from Bio.Seq import Seq
from multiprocessing import Pool
from msgspec import json
from .pal2nal import worker
import xxhash
def get_head_to_seq(nt_db):
    """Get a dictionary of headers to sequences.

    Args:
    ----
        nt_db (RocksDB): The NT rocksdb database.
        recipe (list[str]): A list of each batch index.
    Returns:
    -------
        dict[str, str]: A dictionary of headers to sequences.
    """
    recipe = nt_db.get("getall:batches")
    if recipe:
        recipe = recipe.split(",")
        
    head_to_seq = {}
    for i in recipe:
        lines = nt_db.get_bytes(f"ntbatch:{i}").decode().splitlines()
        head_to_seq.update(
            {
                int(lines[i][1:]): lines[i + 1]
                for i in range(0, len(lines), 2)
                if lines[i] != ""
            },
        )

    return head_to_seq


class Node:
    def __init__(self, head, seq, start, end, score, frame) -> None:
        self.head = head
        self.seq = seq
        self.start = start
        self.end = end
        self.score = score
        self.frame = frame



class exonerate:
    def __init__(self, aa_input, nt_input, folder, chomp_max_distance, orthoset_raw_path, exonerate_path, max_extend) -> None:
        self.aa_input = aa_input
        self.nt_input = nt_input
        self.folder = folder
        self.chomp_max_distance = chomp_max_distance
        self.orthoset_raw_path = orthoset_raw_path
        self.exonerate_path = exonerate_path
        self.max_extend = max_extend
        
    def run(self, genes, original_coords, temp_source_file):
        head_to_seq = dict(parseFasta(temp_source_file))
        for gene in genes:
            gene_name = gene.split(".")[0]
            print(gene)
            aa_path = path.join(self.aa_input, gene)
            nt_path = path.join(self.nt_input, gene.replace(".aa.", ".nt."))
            
            ids = []
            ref_coords = set()
            raw_seq_hash = set()
            for (header, seq) in parseFasta(aa_path):
                if header.endswith('.'):
                    for i, let in enumerate(seq):
                        if let != "-":
                            ref_coords.add(i)
                    continue
                            
                frame = int(header.split("|")[4])
                start, end = find_index_pair(seq, "-")
                raw_seq_hash.add(xxhash.xxh64(seq[start:end].replace("-","")).hexdigest())
                ids.append(quick_rec(header.split("|")[3], frame, seq, start, end))
                
            msa_length = len(seq)
            gap_distance = round(msa_length * 0.3)
            
            clusters, _ = cluster_ids(ids, self.chomp_max_distance, gap_distance, ref_coords, 0)
            cluster_sets = [set(range(a, b+1)) for a, b, _ in clusters]
            raw_path = path.join(self.orthoset_raw_path, gene_name+".fa")
            # max_cluster = max(clusters, key=lambda x: x[1] - x[0])
            # cluster = max_cluster
            final_output = []
            raw_nt_final_output = []
            index = 0
            for cluster_i, cluster in enumerate(clusters):
                cluster_seq = ""
                
                outside_cluster = set()
                for x, cset in enumerate(cluster_sets):
                    if x != cluster_i:
                        outside_cluster.update(cset)
                
                cluster_start = cluster[0]
                for i in range(self.max_extend):
                    if cluster_start - i in outside_cluster:
                        break
                    cluster_start -= i
                    
                cluster_end = cluster[1]
                for i in range(self.max_extend):
                    if cluster_end + i in outside_cluster:
                        break
                    cluster_end += i
                
                for x, i in enumerate(range(cluster_start, cluster_end + 1)):
                    if x == 0:
                        cluster_seq += head_to_seq[str(i)]
                    else:
                        cluster_seq += head_to_seq[str(i)][250:]

                cluster_name = path.join(self.exonerate_path, f"{gene_name}_{cluster[0]}-{cluster[1]}.fa")
                with NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()) as f, NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".txt", dir=gettempdir()) as result:
                    writeFasta(f.name, [(f"{gene_name}_{cluster[0]}-{cluster[1]}", cluster_seq)])
                    
                    command = [
                        "exonerate",
                        "--geneticcode", "1",
                        "--ryo", '>%ti|%tcb-%tce|%s\n%tcs\n',
                        #"--model", "protein2genome",
                        #"--showcigar", "no",
                        "--showvulgar", "no",
                        "--showalignment", "no",
                        "--verbose", "0",
                        raw_path,
                        f.name
                    ]
                                        
                    subprocess.run(command, stdout=result)

                    if path.getsize(result.name) == 0:
                        continue
                    
                    this_nodes = []
                    for header, sequence in parseFasta(result.name, True):
                        header, coords, score = header.split("|")
                        start, end = map(int, coords.split("-"))
                        if start > end:
                            start, end = end, start
                            frame = -((start % 3) + 1)
                        else:
                            frame = (start % 3) + 1
                        index += 1
                        this_nodes.append(Node(f"{index}", sequence, start, end, int(score), frame))
                        
                    kicked_ids = set()
                    overlap_min = 0.1
                    for node_a, node_b in combinations(this_nodes, 2):
                        if node_a.head in kicked_ids or node_b.head in kicked_ids:
                            continue
                        overlap = get_overlap(node_a.start, node_a.end, node_b.start, node_b.end, 1)
                        if overlap:
                            amount = overlap[1] - overlap[0]
                            percent = amount / min((node_a.end - node_a.start), (node_b.end - node_b.start))

                            
                            if percent > overlap_min:
                                if node_a.score > node_b.score:
                                    kicked_ids.add(node_b.head)
                                else:
                                    kicked_ids.add(node_a.head)
                    
                    this_nodes = [node for node in this_nodes if node.head not in kicked_ids]

                    raw_nt_final_output.extend([(f"{gene_name}|placeholder|taxa|CONTIG_{node.head}|{node.frame}|1", node.seq) for node in this_nodes])
                    final_output.extend([(f"{gene_name}|placeholder|taxa|CONTIG_{node.head}|{node.frame}|1", str(Seq(node.seq).translate())) for node in this_nodes])

            kick_dupes = set()
            for node, seq in final_output:
                if xxhash.xxh64(seq).hexdigest() in raw_seq_hash:
                    kick_dupes.add(node)

            final_output = [(node, seq) for node, seq in final_output if node not in kick_dupes]
            raw_nt_final_output = [(node, seq) for node, seq in raw_nt_final_output if node not in kick_dupes]

            if final_output:
                result_path = path.join(self.exonerate_path, "aa", f"{gene_name}.aa.fa")
                result_nt_path = path.join(self.exonerate_path, "nt", f"{gene_name}.nt.fa")
                with NamedTemporaryFile(
                    prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()
                    ) as f, NamedTemporaryFile(
                        prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()
                    ) as result_file, NamedTemporaryFile(
                        prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()
                    ) as aln_file, NamedTemporaryFile(
                        prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()
                    ) as result_nt_file:
                        
                    aa_seqs = parseFasta(aa_path, False)
                    writeFasta(aln_file.name, aa_seqs)
                        
                    writeFasta(f.name, final_output)
                    command = ["mafft", "--anysymbol", "--quiet", "--jtt", "1", "--addfragments", f.name, "--thread", "1", aln_file.name]
                    subprocess.run(command, stdout=result_file)
                    result = list(parseFasta(result_file.name, True))
                    writeFasta(result_path, result)
                    
                    nt_sequences = [(header, sequence.replace("-","")) for header, sequence in parseFasta(nt_path, False)] + raw_nt_final_output
                    writeFasta(result_nt_file.name, nt_sequences)
                    worker(result_file.name, result_nt_file.name, Path(result_nt_path), 1, False, False)
            else:
                aa_seqs = parseFasta(aa_path, False)
                writeFasta(result_path, aa_seqs)
                
                nt_sequences = [(header, sequence.replace("-","")) for header, sequence in parseFasta(nt_path, False)]
                with NamedTemporaryFile(
                    prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()
                    ) as result_nt_file:
                    writeFasta(result_nt_file.name, nt_sequences)
                    worker(result_path, result_nt_file.name, Path(result_nt_path), 1, False, False)

            

def do_folder(folder, args):
    exonerate_path = path.join(folder, "exonerate")
    if path.exists(exonerate_path):
        rmtree(exonerate_path)
    mkdir(exonerate_path)
        
    exonerate_aa_path = path.join(exonerate_path, "aa")
    mkdir(exonerate_aa_path)
        
    exonerate_nt_path = path.join(exonerate_path, "nt")
    mkdir(exonerate_nt_path)
    
    aa_input = path.join(folder, "align")# "outlier", "excise", "aa")
    nt_input = path.join(folder, "nt")#"outlier", "excise", "nt")

    genes = []
    for gene in listdir(aa_input):
        genes.append(gene)
    
    nt_db_path = path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(nt_db_path)
    head_to_seq = get_head_to_seq(nt_db)

    temp_source_file = NamedTemporaryFile(dir=gettempdir(), prefix="seqs_", suffix=".fa")
    writeFasta(temp_source_file.name, head_to_seq.items())

    raw_data = nt_db.get("getall:original_coords")
    original_coords = json.decode(raw_data, type = dict[str, tuple[str, int, int, int, int]])
    
    orthoset_path = path.join(args.orthoset_input, args.orthoset)
    orthoset_raw_path = path.join(orthoset_path, "raw")

    
    per_batch = ceil(len(genes) / args.processes)
    batches = [(genes[i : i + per_batch], original_coords, temp_source_file.name) for i in range(0, len(genes), per_batch)]
    
    if args.processes <= 1:
        exonerate_obj = exonerate(aa_input, nt_input, folder, args.chomp_max_distance, orthoset_raw_path, exonerate_path, args.max_extend)
        for batch in batches:
            exonerate_obj.run(*batch)
    else:
        with Pool(args.processes) as pool:
            pool.starmap(
                exonerate(aa_input, nt_input, folder, args.chomp_max_distance, orthoset_raw_path, exonerate_path, args.max_extend).run,
                batches,
            )

    del temp_source_file

    return True    

def main(args):
    success = False
    if isinstance(args.INPUT, list):
        success = all(do_folder(folder, args) for folder in args.INPUT)
    elif isinstance(args.INPUT, str):
        success = do_folder(args.INPUT, args)
    return success

if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Outlier"
    raise Exception(
        msg,
    )