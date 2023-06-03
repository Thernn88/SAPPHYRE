from collections import defaultdict
import os
from dataclasses import dataclass
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import rmtree

from .timekeeper import KeeperMode, TimeKeeper
from .utils import printv, parseFasta, writeFasta

AA_FOLDER = "aa_merged"
NT_FOLDER = "nt_merged"

STOPCODON = "*"
AA_REPLACE = "X"
NT_REPLACE = "N"


@dataclass
class GeneConfig:
    gene: str
    kick_columns: bool
    kick_percentage: float
    minimum_bp: int
    stopcodon: bool
    rename: bool
    sort: bool
    taxa_folder: Path
    aa_file: Path
    nt_file: Path
    to_kick: set
    taxa_to_taxon: dict
    target: set
    verbose: int
    count_taxa: bool


def kick_taxa(content: list[tuple, tuple], to_kick: set) -> list:
    out = []
    for header, sequence in content:
        taxon = header.split("|")[1]
        if not taxon in to_kick:
            out.append((header, sequence))
    return out


def align_kick_nt(
    fasta_content: list,
    to_kick: set,
    aa_kicks: set,
):
    result = []
    for header, sequence in fasta_content:
        if header not in aa_kicks:
            sequence = "".join(
                [
                    sequence[i : i + 3]
                    for i in range(0, len(sequence), 3)
                    if i not in to_kick
                ],
            )

            result.append((header, sequence))
    return result


def kick_empty_columns(
    fasta_content: list,
    kick_percent: float,
    minimum_bp: int,
) -> list:
    
    kicks = set()
    fasta_data = defaultdict(list)
    out = {}

    for (header, sequence) in fasta_content:
        for i, let in enumerate(sequence):
            fasta_data[i].append(let)

        out[header] = sequence

    to_kick = set()
    for i, characters in fasta_data.items():
        if characters.count("-") / len(characters) > kick_percent:
            to_kick.add(i * 3)

    result = []
    for header in out:
        sequence = "".join(
            [let for i, let in enumerate(out[header]) if i * 3 not in to_kick],
        )
        if len(sequence) - sequence.count("-") >= minimum_bp:
            result.append((header, sequence))
        else:
            kicks.add(header)

    return result, kicks, to_kick


def stopcodon(aa_content: list, nt_content: list) -> tuple:
    aa_out = []
    nt_out = []

    aa_refs = []
    aa_seqs = []
    for header, sequence in aa_content:
        if header.endswith("."):
            aa_refs.append((header, sequence))
        else:
            aa_seqs.append((header, sequence))

    for (aa_header, aa_line), (nt_header, nt_line) in zip(aa_seqs, nt_content):
        
        if aa_header != nt_header:
            print(
                "Warning STOPCODON: Nucleotide order doesn't match Amino Acid order",
            )
        elif STOPCODON in aa_line:
            aa_line = list(aa_line)
            nt_line = list(nt_line)
            for pos, char in enumerate(aa_line):
                if char == STOPCODON:
                    aa_line[pos] = AA_REPLACE

                    nt_line[pos * 3] = NT_REPLACE
                    nt_line[(pos * 3) + 1] = NT_REPLACE
                    nt_line[(pos * 3) + 2] = NT_REPLACE

            aa_line = "".join(aa_line)
            nt_line = "".join(nt_line)

        aa_out.append((aa_header, aa_line))
        nt_out.append((nt_header, nt_line))

    return aa_refs + aa_out, nt_out


def rename_taxon(aa_content: list, nt_content: list, taxa_to_taxon: dict) -> tuple:
    new_aa_content, new_nt_content = [], []
    for (aa_header, aa_line), (nt_header, nt_line) in zip(aa_content, nt_content):
        if not aa_header.endswith("."):
            if aa_header != nt_header:
                print("Warning RENAME: Nucleotide order doesn't match Amino Acid order")
            aa_components = aa_header.split("|")
            nt_components = nt_header.split("|")
            taxa = aa_components[2]
            if taxa not in taxa_to_taxon:
                print(
                    f"Error: Taxa ID, {taxa}, not found in names csv file",
                )

            taxon = taxa_to_taxon[taxa].rstrip("_SPM")
            # print(taxon)
            aa_components[1] = taxon
            nt_components[1] = taxon

            aa_header = "|".join(aa_components)
            nt_header = "|".join(nt_components)

        new_aa_content.append((aa_header, aa_line))
        new_nt_content.append((nt_header, nt_line))

    return new_aa_content, new_nt_content


def clean_gene(gene_config: GeneConfig):
    printv(f"Doing: {gene_config.gene}", gene_config.verbose, 2)
    aa_content = parseFasta(str(gene_config.aa_file))
    nt_content = parseFasta(str(gene_config.nt_file))

    if gene_config.stopcodon:
        aa_content, nt_content = stopcodon(aa_content, nt_content)

    if gene_config.rename:
        aa_content, nt_content = rename_taxon(
            aa_content,
            nt_content,
            gene_config.taxa_to_taxon,
        )

    if gene_config.to_kick:
        aa_content = kick_taxa(aa_content, gene_config.to_kick)
        nt_content = kick_taxa(nt_content, gene_config.to_kick)

    if gene_config.kick_columns:
        aa_content, aa_kicks, cols_to_kick = kick_empty_columns(
            aa_content,
            gene_config.kick_percentage,
            gene_config.minimum_bp,
        )
        nt_content = align_kick_nt(nt_content, cols_to_kick, aa_kicks)

    processed_folder = gene_config.taxa_folder.joinpath("Processed")

    on_target = Path(processed_folder).joinpath("Target")
    on_target_aa = Path(on_target).joinpath("aa")
    on_target_nt = Path(on_target).joinpath("nt")

    off_target = Path(processed_folder).joinpath("Off_Target")
    off_target_aa = Path(off_target).joinpath("aa")
    off_target_nt = Path(off_target).joinpath("nt")

    aa_target_content = []
    # nt_target_content = []
    taxa_count = {}

    if gene_config.gene in gene_config.target or not gene_config.sort:
        if gene_config.count_taxa:
            taxa_names = set(gene_config.taxa_to_taxon.keys())
            taxa_count = taxa_present(aa_content, taxa_names)






        aa_target_content.extend(aa_content)
        writeFasta(str(on_target_aa.joinpath(gene_config.aa_file.name)), aa_content)
        writeFasta(str(on_target_nt.joinpath(gene_config.nt_file.name)), nt_content)

    else:
        writeFasta(str(off_target_aa.joinpath(gene_config.aa_file.name)), aa_content)
        writeFasta(str(off_target_nt.joinpath(gene_config.nt_file.name)), nt_content)

    taxa_local = get_taxa_local(aa_target_content)

    return gene_config.gene, taxa_local, aa_target_content, taxa_count  # , nt_target_content


def get_taxa_local(aa_content: list) -> set:
    taxa_local = set()
    for header, _ in aa_content:
        taxa_local.add(header.split("|")[1])

    return taxa_local

def taxa_present(aa_content: list, names: list) -> dict:
    taxac_present = {i: 0 for i in names}
    for header, _ in aa_content:
        if not header.endswith('.'):
            taxa = header.split("|")[2]
            if taxa in taxac_present:
                taxac_present[taxa] += 1
            else:
                print(f"WARNING Taxa Present: Taxa, {taxa}, not found in names file")

    return taxac_present


def process_folder(args, input_path):
    tk = TimeKeeper(KeeperMode.DIRECT)
    taxa_folder = Path(input_path)
    basename = os.path.basename(taxa_folder)
    no_suffix = basename.split(".")[0]
    print(f"Processing: {basename}")
    aa_folder = taxa_folder.joinpath(AA_FOLDER)
    nt_folder = taxa_folder.joinpath(NT_FOLDER)

    def makent(x):
        return x + ".nt.fa"

    printv("Grabbing necessary files and directories", args.verbose)
    processed_folder = taxa_folder.joinpath("Processed")
    rmtree(processed_folder, ignore_errors=True)
    processed_folder.mkdir(exist_ok=True)

    on_target = Path(processed_folder).joinpath("Target")
    on_target_aa = Path(on_target).joinpath("aa")
    on_target_nt = Path(on_target).joinpath("nt")

    off_target = Path(processed_folder).joinpath("Off_Target")
    off_target_aa = Path(off_target).joinpath("aa")
    off_target_nt = Path(off_target).joinpath("nt")

    for dir in [
        on_target,
        on_target_aa,
        on_target_nt,
        off_target,
        off_target_aa,
        off_target_nt,
    ]:
        dir.mkdir(exist_ok=True)

    to_kick = set()
    if args.kick_taxa and args.kick_file:
        with open(taxa_folder.joinpath(args.kick_file), encoding="utf-8-sig") as fp:
            for line in fp:
                to_kick.add(line.strip())

    target = set()
    if args.sort:
        with open(taxa_folder.joinpath(args.target_file), encoding="utf-8-sig") as fp:
            for line in fp:
                target.add(line.strip())

    taxa_to_taxon = {}
    if args.rename or args.count:
        with open(taxa_folder.joinpath(args.names_csv), encoding="utf-8-sig") as fp:
            for line in fp:
                if line != "\n":
                    parse = line.strip().split(",")

                    id = parse[0]
                    name = parse[1]

                    taxa_to_taxon[id] = name
        

    arguments = []
    for aa_file in aa_folder.glob("*.fa"):
        gene = aa_file.name.split(".")[0]
        nt_file = nt_folder.joinpath(makent(gene))
        this_config = GeneConfig(
            gene,
            args.kick_columns,
            args.kick_percentage,
            args.minimum_bp,
            args.stopcodon,
            args.rename,
            args.sort,
            taxa_folder,
            aa_file,
            nt_file,
            to_kick,
            taxa_to_taxon,
            target,
            args.verbose,
            args.count,
        )
        arguments.append((this_config,))

    with Pool(args.processes) as pool:
        to_write = pool.starmap(clean_gene, arguments, chunksize=1)

    if args.count and not args.concat:
        total = {i: 0 for i in taxa_to_taxon.keys()}
        for _, _, _, taxa_count in to_write:
            for taxa, count in taxa_count.items():
                total[taxa] += count
 
        
    if args.concat:
        sequences = defaultdict(dict)
        gene_lengths = {}
        taxa_sequences_global = {}
        taxa_global = set()
        log = {}
        total = {i: 0 for i in taxa_to_taxon.keys()}
        for gene, taxa_local, aa_content, taxa_count in to_write:
            for taxa, count in taxa_count.items():

                total[taxa] += count
            taxa_global.update(taxa_local)
            this_gene_global_length = 0

            for i, (header, sequence) in enumerate(aa_content):
                if i == 0:
                    this_gene_global_length = len(sequence)
                    gene_lengths[gene] = this_gene_global_length

                taxon = header.split("|")[1]
                sequences[gene][taxon] = sequence

        for gene in sequences:
            this_sequences = sequences[gene]
            for taxa in taxa_global:
                if taxa not in this_sequences:
                    seq = "-" * gene_lengths[gene]
                else:
                    seq = this_sequences[taxa]
                if taxa not in taxa_sequences_global:
                    start = 1
                    taxa_sequences_global[taxa] = seq
                    end = len(taxa_sequences_global[taxa])
                else:
                    start = len(taxa_sequences_global[taxa]) + 1
                    taxa_sequences_global[taxa] += seq
                    end = len(taxa_sequences_global[taxa])

            log[gene] = (start, end)

        output_fas = processed_folder.joinpath(no_suffix + ".fas")
        output_nex = processed_folder.joinpath(no_suffix + ".nex")

        with open(output_fas, "w", encoding="utf-8-sig") as fp:
            for taxa in taxa_sequences_global:
                fp.write(">" + taxa + "\n")
                fp.write(taxa_sequences_global[taxa] + "\n")

        with open(output_nex, "w", encoding="utf-8-sig") as fp:
            fp.write("#nexus\nbegin sets;\n")
            for gene in log:
                start, end = log[gene]
                fp.write(f"CHARSET {gene} = {start}-{end} ;\n")

            fp.write("end;\n")
    if args.count:
        out = ["Taxa,Taxon,Total"]

        for taxa, total_taxa in total.items():

            out.append(f"{taxa},{taxa_to_taxon[taxa]},{total_taxa}")

        with open(str(processed_folder.joinpath("TaxaPresent.csv")), "w") as fp:

            fp.write("\n".join(out))
            
    printv(f"Done! Took {tk.lap():.2f}s", 1)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for input_path in args.INPUT:
        process_folder(args, input_path)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre finalize"
    raise Exception(
        msg,
    )
