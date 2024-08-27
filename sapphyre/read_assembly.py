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
import warnings
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
from Bio import BiopythonWarning
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


def get_regions(nodes, threshold, no_dupes, minimum_ambig):
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

def do_gene(gene, aa_input, nt_input, aa_output, nt_output, no_dupes, compress, excise_consensus):
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    # within_identity = 0.9
    max_score = 8
    min_difference = 0.05
    min_contig_overlap = 0.5
    region_min_ambig = 9
    min_ambig_bp_overlap = 6

    has_region = False
    kicks = 0

    nodes = []
    log_output = []
    aa_gene = path.join(aa_input, gene)
    nt_gene = gene.replace(".aa.", ".nt.")
    for header, sequence in parseFasta(path.join(nt_input, nt_gene)):
        nodes.append((header, sequence))

    # Check bad region
    
    regions = get_regions(nodes, excise_consensus, no_dupes, region_min_ambig)

    if not regions:
        writeFasta(path.join(aa_output, gene), parseFasta(aa_gene), compress)
        writeFasta(path.join(nt_output, nt_gene), nodes, compress)
        return log_output, has_region, kicks

    has_region = True
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
                        kicks += 1
            nt_out = [i for i in nodes if i[0] not in kicked]
            if nt_out:
                writeFasta(path.join(aa_output, gene), [i for i in parseFasta(aa_gene) if i[0] not in kicked], compress)
                writeFasta(path.join(nt_output, nt_gene), nt_out, compress)
            return log_output, has_region, kicks
        
        # Translate and Align
        translate_file = path.join(temp, "translate.fa")
        writeFasta(translate_file, [(header, str(Seq(seq).translate())) for header, seq in parseFasta(contig_file, True)])

        temp_aa = path.join(temp, "aa.fa")
        writeFasta(temp_aa, [i for i in parseFasta(aa_gene) if i[0].endswith(".")])

        mafft = ["mafft","--anysymbol","--quiet","--jtt","1","--addfragments",translate_file,"--thread","1",temp_aa]
        aligned_file = path.join(temp, "aligned.fa")
        with open(aligned_file, "w") as aln_tmp:
            subprocess.run(mafft, stdout=aln_tmp)

        out = list(parseFasta(aligned_file, True))
        writeFasta(aligned_file, out)

        flex_consensus = defaultdict(list)
        formed_contigs = []
        for header, seq in out:
            if "Contig" in header:
                formed_contigs.append((header, seq))
            elif header.endswith("."):
                start,end = find_index_pair(seq, "-")
                for i, base in enumerate(seq[start:end],start):
                    flex_consensus[i].append(base)
        
        with_identity = []
        for header, seq in formed_contigs:
            start, end = find_index_pair(seq, "-")
            matches = 0
            for i in range(start, end):
                matches += min(flex_consensus[i].count(seq[i]), max_score)
            length = end - start
            with_identity.append((header, start, end, seq, matches, length))


        kicked_nodes = set()
        kicked_contigs = set()
        with_identity.sort(key=lambda x: x[5], reverse=True)

        for contig_a, contig_b in combinations(with_identity, 2):
            if contig_a[0] in kicked_contigs or contig_b[0] in kicked_contigs:
                continue
            coords = get_overlap(contig_a[1], contig_a[2], contig_b[1], contig_b[2], 1)
            if coords:
                percent = ((coords[1] - coords[0]) / min((contig_a[2] - contig_a[1]), (contig_b[2] - contig_b[1])))
                if percent < min_contig_overlap:
                    continue

                contig_a_kmer = contig_a[3][coords[0]: coords[1]]
                contig_b_kmer = contig_b[3][coords[0]: coords[1]]

                distance = constrained_distance(contig_a_kmer, contig_b_kmer)
                difference = distance / (coords[1] - coords[0])

                if difference <= min_difference:
                    log_output.append(
                        f"{gene} - {contig_a[0]} has ({', '.join(contigs[contig_a[0]])})\nwith {contig_a[4]} score over {contig_a[5]} length + Kept\nvs ({percent * 100:.2f}% Overlap) ({difference * 100:.2f}% Difference)\n{gene} - {contig_b[0]} has ({', '.join(contigs[contig_b[0]])})\nwith {contig_b[4]} score over {contig_b[5]} length + Kept (Within 95% similar)\n"
                    )
                    continue

                if coords:
                    if contig_a[4] > contig_b[4]:
                        kicked_contigs.add(contig_b[0])
                        kicked_nodes.update(contigs[contig_b[0]])
                        log_output.append(
                            f"{gene} - {contig_a[0]} has ({', '.join(contigs[contig_a[0]])})\nwith {contig_a[4]} score over {contig_a[5]} length + Kept\nvs ({percent * 100:.2f}% Overlap) ({difference * 100:.2f}% Difference)\n{gene} - {contig_b[0]} has ({', '.join(contigs[contig_b[0]])})\nwith {contig_b[4]} matches over {contig_b[5]} length - Kicked\n"
                        )
                    else:
                        kicked_contigs.add(contig_a[0])
                        kicked_nodes.update(contigs[contig_a[0]])
                        log_output.append(
                            f"{gene} - {contig_b[0]} has ({', '.join(contigs[contig_b[0]])})\nwith {contig_b[4]} score over {contig_b[5]} length + Kept\nvs ({percent * 100:.2f}% Overlap) ({difference * 100:.2f}% Difference)\n{gene} - {contig_a[0]} has ({', '.join(contigs[contig_a[0]])})\nwith {contig_a[4]} matches over {contig_a[5]} length - Kicked\n"
                        )

    nt_out = []
    for header, seq in nodes:
        if header.split("|")[3] in kicked_nodes:
            continue
        else:
            nt_out.append((header, seq))

    aa_out = []
    for header, seq in parseFasta(aa_gene):
        if header.split("|")[3] in kicked_nodes:
            kicks += 1
            continue
        else:
            aa_out.append((header, seq))

    if nt_out:
        writeFasta(path.join(aa_output, gene), aa_out, compress)
        writeFasta(path.join(nt_output, nt_gene), nt_out, compress)
            
    return log_output, has_region, kicks


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

    arguments = [(gene, aa_input, nt_input, aa_output, nt_output, args.no_dupes, compress, args.excise_consensus) for gene in genes]
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
    ambig_count = 0
    total_kicks = 0
    for log, has_ambig, kicks in results:
        ambig_count += 1 if has_ambig else 0
        total_kicks += kicks
        log_final.extend(log)
    
    printv(f"{folder}: {ambig_count} ambiguous loci found. Kicked {total_kicks} sequences total.", args.verbose)

    with open(output_folder.joinpath("excise.log"), "w") as log_file:
        log_file.write("\n".join(log_final))

    printv(f"Done! Took {timer.differential():.2f} seconds", args.verbose)

    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre excise"
    raise RuntimeError(
        MSG,
    )
