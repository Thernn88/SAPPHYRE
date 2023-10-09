from collections import defaultdict
import os

from msgspec import json
from phymmr_tools import find_index_pair, is_same_kmer
from Bio.Seq import Seq
from multiprocessing import Pool

from .utils import parseFasta, writeFasta
from . import rocky
from .timekeeper import KeeperMode, TimeKeeper
from .diamond import ReporterHit as Hit

taxa_paths = []

#input
#processes

def main(args):
    if not all(os.path.exists(i) for i in args.INPUT):
        print("ERROR: All folders passed as argument must exist.")#, args.verbose, 0)
        return False
    
    for input_path in args.INPUT:
        extend_super_dir(args, input_path)

    return True

def extend_super_dir(args, super_dir):
    tk = TimeKeeper(KeeperMode.DIRECT)

    taxa_paths, genes = get_taxa_paths(super_dir)

    for taxa in taxa_paths:
        rocky.create_pointer(
            taxa,
            os.path.join(super_dir, taxa, 'rocksdb', 'hits'),
        )

    arguments = [(gene, taxa_paths, super_dir) for gene in genes]

    with Pool(args.processes) as p:
        total = sum(p.starmap(find_extensions, arguments))

    print(f"Took {tk.differential():.2f}s to extend an additional {total} bp")

def get_taxa_paths(super_dir):
    for i, taxa in enumerate(os.listdir(super_dir)):
        taxa_dir = os.path.join(super_dir, taxa)
        # Is folder
        if os.path.isfile(taxa_dir):
            continue


        rocksdb_path = os.path.join(taxa_dir, 'rocksdb', 'hits')

        if not os.path.exists(rocksdb_path):
            continue

        # db = RocksDB(rocksdb_path)
        aa_path = os.path.join(taxa_dir, 'aa_merged')
        if i == 0:
            genes = list(os.listdir(aa_path))

        taxa_paths.append(taxa)

    return taxa_paths, genes

def find_extensions(gene, taxa_paths, super_dir):
    print("Doing:",gene)
    total_extension = 0
    gene_taxon_consensus = defaultdict(lambda: defaultdict(list))
    for taxa in taxa_paths:
        aa_path = os.path.join(super_dir, taxa, 'aa_merged', gene)

        if not os.path.exists(aa_path):
            continue

        this_gene_seqs = [(header, seq) for header, seq in parseFasta(aa_path) if not header.endswith(".")]

        for header, seq in this_gene_seqs:
            taxon = header.split("|")[1]
            for i, let in enumerate(seq):
                gene_taxon_consensus[taxon][i].append(let)

    this_formed_consensus = {}
    for taxa in taxa_paths:
        
        out_path = os.path.join(super_dir, taxa, 'aa_extended', gene)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        aa_path = os.path.join(super_dir, taxa, 'aa_merged', gene)

        if not os.path.exists(aa_path):
            continue

        db = rocky.get_rock(taxa)

        this_db_entry = json.decode(db.get(f"gethits:{gene.split('.')[0]}"), type = list[Hit])
        this_db_sequences = {}

        for entry in this_db_entry:
            this_db_sequences.setdefault(entry.node, {})[entry.frame] = entry.seq

        this_gene_last_component = json.decode(db.get(f"get_last:{gene.split('.')[0]}"), type=dict[str, tuple])

        out = []
        for header, seq in parseFasta(aa_path):
            if header.endswith("."):
                out.append((header, seq))
                continue

            new_seq = seq

            last_component = this_gene_last_component[header]

            lc_header, lc_frame, lc_start = last_component

            if lc_header.count("_") == 2:
                lc_header = "_".join(lc_header.split("_")[:-1])

            lc_start = int(lc_start)
            lc_frame = int(lc_frame)

            taxon = header.split("|")[1] 

            if taxon in this_formed_consensus:
                consensus_seq = this_formed_consensus[taxon]
            else:
                consensus_seq = []
                consensus = gene_taxon_consensus[taxon]
                for i, letters in consensus.items():
                    non_gap_letters = [let for let in letters if let != "-"]
                    if not non_gap_letters:
                        consensus_seq.append("-")
                        continue

                    consensus_seq.append(
                        max(set(non_gap_letters), key=non_gap_letters.count)
                    )
                
                consensus_seq = "".join(consensus_seq)
                this_formed_consensus[taxon] = consensus_seq

            # if this_db_sequences[lc_header]:
            #     out.append((header, seq))
            #     continue

            new_seq = list(seq)
            last_component_seq = this_db_sequences[lc_header][lc_frame]
            
            this_aa = Seq(last_component_seq).translate()

            end = lc_start + len(this_aa)
            this_aa = "-" * lc_start + str(this_aa)
            this_aa += "-" * (len(consensus_seq) - len(this_aa))

            target_region = (lc_start, end)

            query_region = find_index_pair(seq, "-")
        
            target_extension = this_aa[target_region[0]:target_region[1]]
            target_query = seq[target_region[0]:target_region[1]]

            if target_region[0] < query_region[0]:
                out.append((header, seq))
                continue

            extension = target_region[1] - query_region[1]
            

            if target_query.count("-") == len(target_query):
                out.append((header, seq))
                continue

            target_consensus = consensus_seq[target_region[0]:target_region[1]]

            if len(target_extension) != len(target_consensus):
                print("WARNING: Different length bug TBF\n")
                continue

            if (
                    target_region[1] - target_region[0] <= query_region[1] - query_region[0]
                ) or (
                    not is_same_kmer(target_extension, target_consensus)
                ) or ():
                out.append((header, seq))
                continue

            

            for i in range(query_region[1], target_region[1]):
                new_seq[i] = this_aa[i]
            
            print(header)
            total_extension += extension

            new_seq = "".join(new_seq)
            out.append((header, new_seq))
        writeFasta(out_path, out)

    return total_extension