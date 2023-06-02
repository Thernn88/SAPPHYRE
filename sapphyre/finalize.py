import os
from dataclasses import dataclass
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import rmtree

from .timekeeper import KeeperMode, TimeKeeper
from .utils import printv

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


def grab_gene(fp, to_kick: set) -> list:
    out = []
    for line in fp:
        if line.startswith(">"):
            taxon = line.split("|")[1]
            if taxon in to_kick:
                next(fp)  # Skip this row and sequence row
            else:
                out.append(line.strip())
        else:
            if line.strip() != "":
                out.append(line.strip())
    return out


def align_kick_nt(
    fasta_content: list,
    to_kick: set,
    aa_kicks: set,
):
    result = []
    for i in range(0, len(fasta_content), 2):
        header = fasta_content[i]
        sequence = fasta_content[i + 1]
        if header not in aa_kicks:
            sequence = "".join(
                [
                    sequence[i : i + 3]
                    for i in range(0, len(sequence), 3)
                    if i not in to_kick
                ],
            )

            result.append(header)
            result.append(sequence)

    return result


def kick_empty_columns(
    fasta_content: list,
    kick_percent: float,
    minimum_bp: int,
) -> list:
    kicks = set()
    fasta_data = {}
    out = {}

    for i in range(0, len(fasta_content), 2):
        header = fasta_content[i]
        sequence = fasta_content[i + 1]

        for j in range(0, len(sequence), 1):
            let = sequence[j : j + 1]
            fasta_data.setdefault(j, []).append(let)

        out[header] = sequence

    to_kick = set()
    if not to_kick:
        for i, characters in fasta_data.items():
            if characters.count("-") / len(characters) > kick_percent:
                to_kick.add(i * 3)

    result = []
    for header in out:
        sequence = "".join(
            [let for i, let in enumerate(out[header]) if i * 3 not in to_kick],
        )
        if len(sequence) - sequence.count("-") >= minimum_bp:
            result.append(header)
            result.append(sequence)
        else:
            kicks.add(header)

    return result, kicks, to_kick


def stopcodon(aa_content: list, nt_content: list) -> tuple:
    aa_out = []
    nt_out = []

    aa_refs = []
    aa_seqs = []
    for i in range(0, len(aa_content), 2):
        if aa_content[i].endswith("."):
            aa_refs.append(aa_content[i])
            aa_refs.append(aa_content[i + 1])
        else:
            aa_seqs.append(aa_content[i])
            aa_seqs.append(aa_content[i + 1])

    for aa_line, nt_line in zip(aa_seqs, nt_content):
        if aa_line.startswith(">"):
            if aa_line != nt_line:
                print(
                    "Warning STOPCODON: Nucleotide line doesn't match Amino Acid line",
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

        aa_out.append(aa_line)
        nt_out.append(nt_line)
    return aa_refs + aa_out, nt_out


def rename_taxon(aa_content: list, nt_content: list, taxa_to_taxon: dict) -> tuple:
    for i, (aa_line, nt_line) in enumerate(zip(aa_content, nt_content)):
        if aa_line.startswith(">") and not aa_line.endswith("."):
            if aa_line != nt_line:
                print("Warning RENAME: Nucleotide line doesn't match Amino Acid line")
            aa_components = aa_line.split("|")
            nt_components = nt_line.split("|")
            if aa_components[2] not in taxa_to_taxon:
                print(
                    f"Error: Taxa ID, {aa_components[2]}, not found in names csv file",
                )
            taxon = taxa_to_taxon[aa_components[2]].strip("_SPM")
            aa_components[1] = taxon
            nt_components[1] = taxon

            aa_line = "|".join(aa_components)
            nt_line = "|".join(nt_components)

            aa_content[i] = aa_line
            nt_content[i] = nt_line

    return aa_content, nt_content


def clean_gene(gene_config):
    printv(f"Doing: {gene_config.gene}", gene_config.verbose, 2)
    with gene_config.aa_file.open() as aa_fp, gene_config.nt_file.open() as nt_fp:
        aa_content = grab_gene(aa_fp, gene_config.to_kick)
        nt_content = grab_gene(nt_fp, gene_config.to_kick)

        if gene_config.kick_columns:
            aa_content, aa_kicks, cols_to_kick = kick_empty_columns(
                aa_content,
                gene_config.kick_percentage,
                gene_config.minimum_bp,
            )
            nt_content = align_kick_nt(nt_content, cols_to_kick, aa_kicks)

        if gene_config.stopcodon:
            aa_content, nt_content = stopcodon(aa_content, nt_content)

        if gene_config.rename:
            aa_content, nt_content = rename_taxon(
                aa_content,
                nt_content,
                gene_config.taxa_to_taxon,
            )

        processed_folder = gene_config.taxa_folder.joinpath("Processed")

        on_target = Path(processed_folder).joinpath("Target")
        on_target_aa = Path(on_target).joinpath("aa")
        on_target_nt = Path(on_target).joinpath("nt")

        off_target = Path(processed_folder).joinpath("Off_Target")
        off_target_aa = Path(off_target).joinpath("aa")
        off_target_nt = Path(off_target).joinpath("nt")

        aa_target_content = []
        nt_target_content = []

        if gene_config.gene in gene_config.target or not gene_config.sort:
            with on_target_aa.joinpath(gene_config.aa_file.name).open(mode="w") as fp:
                for line in aa_content:
                    fp.write(line + "\n")
                    aa_target_content.append(line)
            with on_target_nt.joinpath(gene_config.nt_file.name).open(mode="w") as fp:
                for line in nt_content:
                    fp.write(line + "\n")
                    nt_target_content.append(line)
        else:
            with off_target_aa.joinpath(gene_config.aa_file.name).open(mode="w") as fp:
                for line in aa_content:
                    fp.write(line + "\n")
            with off_target_nt.joinpath(gene_config.nt_file.name).open(mode="w") as fp:
                for line in nt_content:
                    fp.write(line + "\n")

        taxa_local = get_taxa_local(aa_target_content)

        return gene_config.gene, taxa_local, aa_target_content  # , nt_target_content


def get_taxa_local(aa_content: list) -> set:
    taxa_local = set()
    for line in aa_content:
        if line.startswith(">"):
            taxa_local.add(line.split("|")[1])

    return taxa_local


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
    if args.rename:
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
        )
        arguments.append((this_config,))

    with Pool(args.processes) as pool:
        to_write = pool.starmap(clean_gene, arguments, chunksize=1)

    if args.concat:
        sequences = {}
        gene_lengths = {}
        taxa_sequences_global = {}
        taxa_global = set()
        log = {}
        for gene, taxa_local, aa_content in to_write:
            taxa_global.update(taxa_local)
            this_gene_global_length = 0

            for i in range(0, len(aa_content), 2):
                header = aa_content[i]
                sequence = aa_content[i + 1]

                if i == 0:
                    this_gene_global_length = len(sequence)
                    gene_lengths[gene] = this_gene_global_length

                taxon = header.split("|")[1]
                sequences.setdefault(gene, {})
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
