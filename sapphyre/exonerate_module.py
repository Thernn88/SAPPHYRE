from math import ceil
import subprocess
from tempfile import NamedTemporaryFile
from wrap_rocks import RocksDB
from os import listdir, mkdir, path
from .directional_cluster import cluster_ids, within_distance, node_to_ids, quick_rec
from .utils import gettempdir, parseFasta, printv, writeFasta
from sapphyre_tools import (
    find_index_pair,
)
from multiprocessing import Pool
from msgspec import json

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


class exonerate:
    def __init__(self, aa_input, nt_input, folder, chomp_max_distance, orthoset_raw_path, exonerate_path) -> None:
        self.aa_input = aa_input
        self.nt_input = nt_input
        self.folder = folder
        self.chomp_max_distance = chomp_max_distance
        self.orthoset_raw_path = orthoset_raw_path
        self.exonerate_path = exonerate_path
        
    def run(self, genes, head_to_seq, original_coords):
        for gene in genes:
            gene_name = gene.split(".")[0]
            print(gene)
            aa_path = path.join(self.aa_input, gene)
            nt_path = path.join(self.nt_input, gene.replace(".aa.", ".nt."))
            
            ids = []
            ref_coords = set()
            for (header, seq) in parseFasta(aa_path):
                if header.endswith('.'):
                    for i, let in enumerate(seq):
                        if let != "-":
                            ref_coords.add(i)
                    continue
                            
                frame = int(header.split("|")[4])
                start, end = find_index_pair(seq, "-")
                ids.append(quick_rec(header.split("|")[3], frame, seq, start, end))
                
            msa_length = len(seq)
            gap_distance = round(msa_length * 0.3)
            
            clusters, _ = cluster_ids(ids, self.chomp_max_distance, gap_distance, ref_coords, 0)
            raw_path = path.join(self.orthoset_raw_path, gene_name+".fa")
            # max_cluster = max(clusters, key=lambda x: x[1] - x[0])
            # cluster = max_cluster
            for cluster in clusters:
                cluster_seq = ""
                test = []
                for x, i in enumerate(range(cluster[0], cluster[1] + 1)):
                    if x == 0:
                        cluster_seq += head_to_seq[i]
                    else:
                        cluster_seq += head_to_seq[i][250:]
                    test.append(head_to_seq[i])
                cluster_name = path.join(self.exonerate_path, f"{gene_name}_{cluster[0]}-{cluster[1]}.txt")
                with NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()) as f, open(cluster_name, "w") as result:#NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".txt", dir=gettempdir()) as result:
                    writeFasta(f.name, [(f"{gene_name}_{cluster[0]}-{cluster[1]}", cluster_seq)])
                    
                    command = [
                        "exonerate",
                        "--geneticcode", "1",
                        "--model", "protein2genome",
                        "--showcigar", "no",
                        "--showvulgar", "true",
                        "--verbose", "0",
                        raw_path,
                        f.name
                    ]
                                        
                    subprocess.run(command, stdout=result)
                    # print(result.name)
                    # input("called")
                #writeFasta(f"{cluster_name}.test", [(i, test[i]) for i in range(len(test))])
                


def do_folder(folder, args):
    exonerate_path = path.join(folder, "Exonerate")
    if not path.exists(exonerate_path):
        mkdir(exonerate_path)
    
    aa_input = path.join(folder, "outlier", "excise", "aa")
    nt_input = path.join(folder, "outlier", "excise", "nt")
    
    genes = []
    for gene in listdir(aa_input):
        genes.append(gene)
    
    nt_db_path = path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(nt_db_path)
    head_to_seq = get_head_to_seq(nt_db)
    raw_data = nt_db.get("getall:original_coords")
    original_coords = json.decode(raw_data, type = dict[str, tuple[str, int, int, int, int]])
    
    orthoset_path = path.join(args.orthoset_input, args.orthoset)
    orthoset_raw_path = path.join(orthoset_path, "raw")

    
    per_batch = ceil(len(genes) / args.processes)
    batches = [(genes[i : i + per_batch], head_to_seq, original_coords) for i in range(0, len(genes), per_batch)]
    
    if args.processes <= 1:
        exonerate_obj = exonerate(aa_input, nt_input, folder, args.chomp_max_distance, orthoset_raw_path, exonerate_path)
        for batch in batches:
            exonerate_obj.run(*batch)
    else:
        with Pool(args.processes) as pool:
            pool.starmap(
                exonerate(aa_input, nt_input, folder, args.chomp_max_distance, orthoset_raw_path, exonerate_path).run,
                batches,
            )

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