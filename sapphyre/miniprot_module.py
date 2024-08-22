from collections import Counter, defaultdict
from itertools import combinations
from math import ceil
from pathlib import Path
from shutil import rmtree
import subprocess
from tempfile import NamedTemporaryFile
from time import time
from wrap_rocks import RocksDB
from os import listdir, mkdir, path
from .utils import gettempdir, parseFasta, printv, writeFasta
from sapphyre_tools import (
    find_index_pair,
    get_overlap,
    entropy_filter,
    bio_revcomp,
)
from .directional_cluster import node_to_ids
from Bio.Seq import Seq
from multiprocessing import Pool
from msgspec import json
from .pal2nal import worker
import xxhash
from .diamond import ReferenceHit, ReporterHit as Hit

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
    def __init__(self, head, seq, start, end, frame, ref_start, ref_end, ref_name) -> None:
        self.head = head
        self.seq = seq
        self.start = start
        self.end = end
        self.frame = frame
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ref_name = ref_name



def generate_sequence(ids, head_to_seq):
    ids_to_coords = {}
    prev_og = head_to_seq[ids[0]]
    start, end = 0, len(prev_og)
    ids_to_coords[ids[0]] = (start, end)
    for i, child in enumerate(ids):
        if i == 0:
            continue
        prev_og += head_to_seq[child][250:]
        start = end - 250
        end = len(prev_og)
        ids_to_coords[child] = (start, end)

    return ids_to_coords, prev_og


class miniprot:
    def __init__(self, folder, chomp_max_distance, top_path, miniprot_path, max_extend, target_to_taxon, entropy_percent) -> None:
        self.folder = folder
        self.chomp_max_distance = chomp_max_distance
        self.top_path = top_path
        self.miniprot_path = miniprot_path
        self.max_extend = max_extend
        self.target_to_taxon = target_to_taxon
        self.entropy_percent = entropy_percent
        
    def run(self, batches, temp_source_file):
        head_to_seq = dict(parseFasta(temp_source_file))
        batch_result = []
        additions = 0
        for gene_name, fasta_input in batches:
            print(gene_name)
            
            gene_ids = []
            done_ids = set()
            raw_references = []
            id_count = Counter()
            for header, sequence in parseFasta(fasta_input):
                if header.endswith("."):
                    raw_references.append((header, sequence.replace("-", "")))
                    continue
                for id in node_to_ids(header.split("|")[3]):
                    if id in done_ids:
                        continue
                    strand = "-" if int(header.split("|")[4]) < 0 else "+"
                    gene_ids.append((id, strand))

            diamond_ids = sorted(gene_ids, key=lambda x: x[0])
            clusters = []
            current_cluster = []

            for child_index, strand in diamond_ids:
                if not current_cluster:
                    current_cluster.append(child_index)
                    current_index = child_index
                    current_strand = strand
                else:
                    if strand == current_strand and (child_index - current_index <= self.chomp_max_distance):
                        current_cluster.append(child_index)
                        current_index = child_index
                        current_strand = strand
                    else:
                        if len(current_cluster) > 2:
                            clusters.append((current_cluster[0], current_cluster[-1], current_strand))
                        current_cluster = [child_index]
                        current_index = child_index
            
            if current_cluster:
                if len(current_cluster) > 2:
                    clusters.append((current_cluster[0], current_cluster[-1], current_strand))
            
            cluster_sets = [(set(range(a, b+1)), strand) for a, b, strand in clusters]
            raw_path = path.join(self.top_path, gene_name+".aln.fa")
            # max_cluster = max(clusters, key=lambda x: x[1] - x[0])
            # cluster = max_cluster
            
            aa_path = path.join(self.miniprot_path, "aa", gene_name+".aa.fa")
            nt_path = path.join(self.miniprot_path, "nt", gene_name+".nt.fa")
            final_aa = []
            final_nt = []

            for cluster_i, (cluster_start, cluster_end, strand) in enumerate(clusters):

                outside_cluster = set()
                for x, (cset, cstrand) in enumerate(cluster_sets):
                    if x != cluster_i:
                        outside_cluster.update(cset)

                for i in range(self.max_extend):
                    if cluster_start - i in outside_cluster:
                        break
                    cluster_start -= i
                    
                for i in range(self.max_extend):
                    if cluster_end + i in outside_cluster:
                        break
                    cluster_end += i
                    
                ids_to_coords, cluster_seq = generate_sequence([str(i) for i in range(cluster_start, cluster_end + 1)], head_to_seq)

                cluster_name = path.join(self.miniprot_path, f"{gene_name}_{cluster_start}-{cluster_end}.txt")
                cluster_input = path.join(self.miniprot_path, f"{gene_name}_{cluster_start}-{cluster_end}.fa")
                #NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()) as f
                with open(cluster_input, "w") as f, NamedTemporaryFile(prefix=f"{gene_name}_", suffix=".fa", dir=gettempdir()) as temp_raw, open(cluster_name, "w") as result:
                    writeFasta(f.name, [("cluster", cluster_seq)])
                    f.flush()
                    writeFasta(temp_raw.name, [(header, sequence.replace("-","")) for header, sequence in parseFasta(raw_path)])
                    command = [
                        "miniprot/miniprot",
                        f.name,
                        temp_raw.name,
                        "--gff",
                    ]
                                        
                    subprocess.run(command, stdout=result)
                    if path.getsize(result.name) == 0:
                        continue

                    this_nodes = defaultdict(list)
                    with open(result.name) as result_file:
                        for line in result_file:
                            if not line.strip():
                                continue
                            if line.startswith("#"):
                                if "PAF" in line:
                                    alignment_against = line.split("\t")[1]
                                    alignment_score = int(line.split("\t")[13].replace("AS:i:",""))
                                continue
                            
                            fields = line.strip().split("\t")
                            
                            if fields[2] != "CDS":
                                continue
                            
                            
                            
                            data = fields[8].split(" ")
                            data_fields = data[0].split(";")
                            ref_id, identity = None, None
                            for field in data_fields:
                                if "Target=" in field:
                                    ref_id = field.replace("Target=", "")
                                if "Identity=" in field:
                                    identity = float(field.replace("Identity=", ""))
                            
                            ref_start = int(data[1])
                            ref_end = int(data[2])
                            offset = int(fields[7])
                            
                            start = int(fields[3])
                            end = int(fields[4])
                            strand = fields[6]
                            
                            sequence = cluster_seq[start-1:end]
                            
                            
                            if strand == "-":
                                
                                sequence = bio_revcomp(cluster_seq[start-1:end])
                            else:
                                sequence = cluster_seq[start-1:end]
                                
                            sequence = sequence[offset:]
                            frame = (offset % 3 + 1)
                            if strand == "-":
                                frame = -frame  
                            if len(sequence) % 3 != 0:
                                sequence = sequence[:-(len(sequence) % 3)]

                            this_nodes[(alignment_score, alignment_against)].append(Node(None, sequence, start, end, frame, ref_start, ref_end, ref_id))
                    if not this_nodes:
                        continue
                    this_nodes = max(this_nodes.items(), key=lambda x: x[0])[1]

                    # sequence_template = ">{}|{}\n{}\n"
                    # as_fasta = [sequence_template.format(node.head, node.frame, node.seq) for node in this_nodes]
                    # passing_fasta = entropy_filter(as_fasta, self.entropy_percent)
                    # passing_headers = set()
                    # for header, _ in passing_fasta:
                    #     passing_headers.add(header.strip()[1:])
                            

                    for node in this_nodes:
                        ids_in_qstart = [int(id) for id, range in ids_to_coords.items() if node.start >= range[0] and node.start <= range[1] or node.end >= range[0] and node.end <= range[1]]
    
                        final_ids = []
                        for id in ids_in_qstart:
                            count = id_count[id]
                            if count == 0:
                                final_ids.append(str(id))
                            else:
                                final_ids.append(f"{id}_{count}")
                                
                            id_count[id] += 1
                            
                        node.head = "&&".join(final_ids)
                            
                        target = node.ref_name
                        key = f"{gene_name}|{target}"
                        
                        _, ref, _ = self.target_to_taxon[key]
                        additions += 1
                        
                        header = f"{gene_name}|{target}|{ref}|NODE_{node.head}|{node.frame}|1"
                        nt_sequence = node.seq
                        aa_sequence = Seq(node.seq).translate()
                        
                        final_aa.append((header, str(aa_sequence)))
                        final_nt.append((header, nt_sequence))
                        
            if final_aa:
                writeFasta(aa_path, raw_references+final_aa)
                writeFasta(nt_path, final_nt)
                        
                        

        return None, additions

def do_folder(folder, args):
    miniprot_path = path.join(folder, "miniprot")
    if path.exists(miniprot_path):
        rmtree(miniprot_path)
    mkdir(miniprot_path)
    
    miniprot_aa_path = path.join(miniprot_path, "aa")
    mkdir(miniprot_aa_path)
    miniprot_nt_path = path.join(miniprot_path, "nt")
    mkdir(miniprot_nt_path)
        
    nt_db_path = path.join(folder, "rocksdb", "sequences", "nt")
    nt_db = RocksDB(nt_db_path)
    head_to_seq = get_head_to_seq(nt_db)

    temp_source_file = NamedTemporaryFile(dir=gettempdir(), prefix="seqs_", suffix=".fa")
    writeFasta(temp_source_file.name, head_to_seq.items())
    top_path = path.join(folder, "top")

    orthoset_db_path = path.join(args.orthoset_input, args.orthoset, "rocksdb")
    orthoset_db = RocksDB(orthoset_db_path)
    
    for p_folder in ["excise", "clusters", "blosum"]:
        p_path = path.join(folder, "outlier", p_folder)
        if not path.exists(p_path):
            p_path = None
        else:
            break
    
    if p_path is None:
        printv("ERROR: Outlier folder not found.", args.verbose, 0)
        return False
    
    aa_input = path.join(p_path, "aa")

    genes = [(path.basename(f).split(".")[0], path.join(aa_input, f)) for f in listdir(aa_input) if ".fa" in f]
    per_batch = ceil(len(genes) / args.processes)
    batches = [(genes[i:i+per_batch], temp_source_file.name) for i in range(0, len(genes), per_batch)]
    
    target_to_taxon_raw = orthoset_db.get_bytes("getall:targetreference")
    target_to_taxon = json.decode(
        target_to_taxon_raw,
        type=dict[str, list[str | int]],
    )
    batch_result = []
    if args.processes <= 1:
        miniprot_obj = miniprot(folder, args.chomp_max_distance, top_path, miniprot_path, args.max_extend, target_to_taxon, args.entropy_percent)
        for batch in batches:
            batch_result.append(miniprot_obj.run(*batch))
    else:
        with Pool(args.processes) as pool:
            batch_result.extend(pool.starmap(
                miniprot(folder, args.chomp_max_distance, top_path, miniprot_path, args.max_extend, target_to_taxon, args.entropy_percent).run,
                batches,
            ))
            

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
