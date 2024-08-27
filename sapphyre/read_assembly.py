from collections import defaultdict
from functools import cached_property
from itertools import combinations, product
from multiprocessing import Pool
from os import listdir, makedirs, path, stat
from pathlib import Path
from shutil import move
import copy
import subprocess
from tempfile import TemporaryDirectory
from msgspec import Struct, json
from sapphyre_tools import (
    convert_consensus,
    dumb_consensus,
    dumb_consensus_dupe,
    find_index_pair,
    get_overlap,
    is_same_kmer,
    constrained_distance,
    bio_revcomp,
    join_with_exclusions,
    join_triplets_with_exclusions,
)
from Bio.Seq import Seq
from .directional_cluster import cluster_ids, within_distance, node_to_ids, quick_rec
from wrap_rocks import RocksDB

from .timekeeper import KeeperMode, TimeKeeper
from .utils import gettempdir, parseFasta, printv, writeFasta

def make_duped_consensus(
    raw_sequences: list, threshold: float
) -> str:
    bundled_seqs = [(seq, int(header.split("|")[5])) for header, seq in raw_sequences if header[-1] != "."]
    return dumb_consensus_dupe(bundled_seqs, threshold, 0)


def check_covered_bad_regions(consensus, min_ambiguous, ambig_char='X', max_distance=18):
    x_indices = []
    current_group = []

    start, stop = find_index_pair(consensus, "X")

    for i, base in enumerate(consensus[start:stop], start):
        if base == ambig_char:
            x_indices.append(i)

    regions = []
    if not x_indices:
        return regions

    for num in x_indices:
        if not current_group:
            current_group.append(num)
        elif num - current_group[-1] <= max_distance:
            current_group.append(num)
        else:
            if len(current_group) >= min_ambiguous:
                regions.append((current_group[0], current_group[-1] + 1))
            current_group = [num]

    if current_group:
        if len(current_group) >= min_ambiguous:
            regions.append((current_group[0], current_group[-1] + 1))


    return regions


def has_region(nodes, threshold, no_dupes, minimum_ambig):
    sequences = [x[1] for x in nodes]
    if no_dupes:
        consensus_seq = dumb_consensus(sequences, threshold, 0)
    else:
        consensus_seq = make_duped_consensus(
            nodes, threshold
        )

    consensus_seq = convert_consensus(sequences, consensus_seq)
    regions = check_covered_bad_regions(consensus_seq, minimum_ambig)

    return regions

def do_gene(gene, aa_input, nt_input, aa_output, nt_output, no_dupes, compress):

    min_identity = 0.8
    region_threshold = 0.5
    region_min_ambig = 9
    min_ambig_bp_overlap = 6

    nodes = []
    log_output = []
    aa_gene = path.join(aa_input, gene)
    nt_gene = gene.replace(".aa.", ".nt.")
    for header, sequence in parseFasta(path.join(nt_input, nt_gene)):
        nodes.append((header, sequence))

    # Check bad region
    
    regions = has_region(nodes, region_threshold, no_dupes, region_min_ambig)

    if not regions:
        writeFasta(path.join(aa_output, gene), parseFasta(aa_gene), compress)
        writeFasta(path.join(nt_output, nt_gene), nodes, compress)
        return log_output

    # Assembly

    with TemporaryDirectory(dir=gettempdir(), prefix=f"{gene}_") as temp:
        clean_out = [
            (node[0].split("|")[3], node[1].replace("-","")) for node in nodes
        ]

        temp_gene = path.join(temp, "reads.fa")
        writeFasta(temp_gene, clean_out)

        log_file = path.join(temp, "cap3.log")
        with open(log_file, "w") as log:
            command = [
                "./cap3",
                temp_gene,
            ]
                                
            subprocess.run(command, stdout = log)
        contig_file = temp_gene + ".cap.contigs"
        contigs = defaultdict(list)
        with open(log_file) as fp:
            for line in fp:
                if  "*******************" in line:
                    contig = line.split("******************* ")[1].split(" *")[0].replace(" ", "")
                else:
                    if "NODE" in line:
                        contigs[contig].append(line.replace(" ","").split("+")[0])
                    elif ":" in line:
                        break

        if stat(contig_file).st_size == 0:
            log_output.append(f"{gene} - Ambig with no contigs found, Kicking ambig region")
            kicked = set()
            for region in regions:
                r_start, r_end = region
                for node in nodes:
                    n_start, n_end = find_index_pair(node[1], "-")
                    if get_overlap(r_start, r_end, n_start, n_end, min_ambig_bp_overlap):
                        log_output.append(f"{gene} - Kicked: {node[0]}")
                        kicked.add(node[0])
            nt_out = [i for i in nodes if i[0] not in kicked]
            if nt_out:
                writeFasta(path.join(aa_output, gene), [i for i in parseFasta(aa_gene) if i[0] not in kicked], compress)
                writeFasta(path.join(nt_output, nt_gene), nt_out, compress)
            return log_output
        
        # Translate and Align
        translate_file = path.join(temp, "translate.fa")
        writeFasta(translate_file, [(header, str(Seq(seq).translate())) for header, seq in parseFasta(contig_file, True)])

        temp_aa = path.join(temp, "aa.fa")
        writeFasta(temp_aa, parseFasta(aa_gene))

        mafft = ["mafft","--anysymbol","--quiet","--jtt","1","--addfragments",translate_file,"--thread","1",temp_aa]
        aligned_file = path.join(temp, "aligned.fa")
        with open(aligned_file, "w") as aln_tmp:
            subprocess.run(mafft, stdout=aln_tmp)

        out = list(parseFasta(aligned_file, True))
        writeFasta(aligned_file, out)

        flex_consensus = defaultdict(set)
        formed_contigs = []
        for header, seq in out:
            if "Contig" in header:
                formed_contigs.append((header, seq))
            elif header.endswith("."):
                start,end = find_index_pair(seq, "-")
                for i, base in enumerate(seq[start:end],start):
                    flex_consensus[i].add(base)
        
        kicked_nodes = set()
        for header, seq in formed_contigs:
            start, end = find_index_pair(seq, "-")
            matches = 0
            for i in range(start, end):
                if seq[i] in flex_consensus[i]:
                    matches += 1
            length = end - start
            identity = matches / length
            if identity < min_identity:
                kicked_nodes.update(contigs[header])
                log_output.append(
                    f"{gene} - {header} has {','.join(contigs[header])}\nwith {matches} matches over {length} length equals {identity:.2f}/{min_identity} identity Kicked -"
                )
            else:
                log_output.append(
                    f"{gene} - {header} has {','.join(contigs[header])}\nwith {matches} matches over {length} length equals {identity:.2f}/{min_identity} identity Kept +"
                )
    nt_out = []
    for header, seq in nodes:
        if header in kicked_nodes:
            continue
        else:
            nt_out.append((header, seq))

    aa_out = []
    for header, seq in parseFasta(aa_gene):
        if header in kicked_nodes:
            continue
        else:
            aa_out.append((header, seq))

    if nt_out:
        writeFasta(path.join(aa_output, gene), aa_out, compress)
        writeFasta(path.join(nt_output, nt_gene), nt_out, compress)
            
    return log_output


def main(args, sub_dir):
    timer = TimeKeeper(KeeperMode.DIRECT)
    
    folder = args.INPUT
    input_folder = Path(folder, "outlier", sub_dir)
    if not input_folder.exists():
        input_folder = Path(folder, sub_dir)

    printv(f"Processing: {folder}", args.verbose)

    output_folder = Path(folder, "outlier", "excise")

    aa_output = output_folder.joinpath("aa")
    nt_output = output_folder.joinpath("nt")

    if not path.exists(aa_output):
        makedirs(str(aa_output), exist_ok=True)
        makedirs(str(nt_output), exist_ok=True)

    aa_input = input_folder.joinpath("aa")
    nt_input = input_folder.joinpath("nt")

    compress = not args.uncompress_intermediates or args.compress

    genes = [fasta for fasta in listdir(aa_input) if ".fa" in fasta]

    arguments = [(gene, aa_input, nt_input, aa_output, nt_output, args.no_dupes, compress) for gene in genes]
    if args.processes > 1:
        with Pool(args.processes) as pool:
            results = pool.starmap(do_gene, arguments)
    else:
        results = []
        for arg in arguments:
            results.append(
                do_gene(
                    *arg
                )
            )
 
    log_final = []
    for log in results:
        log_final.extend(log)
    
    with open(output_folder.joinpath("excise.log"), "w") as log_file:
        log_file.write("\n".join(log_final))

    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)

    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
