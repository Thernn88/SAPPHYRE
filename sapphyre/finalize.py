import json
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from glob import glob
from math import ceil
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import rmtree

import requests
from bs4 import BeautifulSoup

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta

AA_FOLDER = "aa_merged"
NT_FOLDER = "nt_merged"

STOPCODON = "*"
AA_REPLACE = "X"
NT_REPLACE = "N"


@dataclass
class GeneConfig:
    gene: str
    kick_columns: float
    minimum_bp: int
    stopcodon: bool
    rename: bool
    taxa_folder: Path
    aa_file: Path
    nt_file: Path
    to_kick: set
    taxa_to_taxon: dict
    target: set
    verbose: int
    count_taxa: bool
    generating_names: bool
    no_references: bool
    compress: bool


def kick_taxa(content: list[tuple, tuple], to_kick: set) -> list:
    out = []
    for header, sequence in content:
        taxon = header.split("|")[1]
        if taxon not in to_kick:
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
    kick_percent = 1 - kick_percent
    for header, sequence in fasta_content:
        for i, let in enumerate(sequence):
            fasta_data[i].append(let)

        out[header] = sequence

    to_kick = set()
    for i, characters in fasta_data.items():
        if characters.count("-") / len(characters) > kick_percent:
            to_kick.add(i * 3)

    result = []
    for header, sequence in out.items():
        sequence = "".join(
            [let for i, let in enumerate(sequence) if i * 3 not in to_kick],
        )
        if len(sequence) - sequence.count("-") >= minimum_bp:
            result.append((header, sequence))
        else:
            kicks.add(header)

    return result, kicks, to_kick


def stopcodon(aa_content: list, nt_content: list) -> tuple:
    aa_out = []
    nt_out = []

    for (aa_header, aa_line), (nt_header, nt_line) in zip(aa_content, nt_content):
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

    return aa_out, nt_out


def rename_taxon(aa_content: list, nt_content: list, taxa_to_taxon: dict) -> tuple:
    new_aa_content = []
    new_nt_content = []
    aa_candidates = []
    for header, seq in aa_content:
        if header.endswith("."):
            new_aa_content.append((header, seq))
        else:
            aa_candidates.append((header, seq))
    nt_candidates = []
    for header, seq in nt_content:
        if header.endswith("."):
            new_nt_content.append((header, seq))
        else:
            nt_candidates.append((header, seq))
    for (aa_header, aa_line), (nt_header, nt_line) in zip(aa_candidates, nt_candidates):
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

            taxon = taxa_to_taxon[taxa]
            # print(taxon)
            aa_components[1] = taxon
            nt_components[1] = taxon

            aa_header = "|".join(aa_components)
            nt_header = "|".join(nt_components)

        new_aa_content.append((aa_header, aa_line))
        new_nt_content.append((nt_header, nt_line))

    return new_aa_content, new_nt_content


def get_taxon_total(aa_content):
    taxon_set = set()
    for header, _ in aa_content:
        if not header.endswith("."):
            taxon = header.split("|")[1]
            taxon_set.add(taxon)

    return len(taxon_set)


def taxon_only(content):
    return [(header.split("|")[1], sequence) for header, sequence in content]


def kick_gene(present_taxa, minimum_percentage, global_total_taxon):
    if len(present_taxa) == 0:
        return 0 <= minimum_percentage

    return (len(present_taxa) / len(global_total_taxon)) <= minimum_percentage


def clean_gene(gene_config: GeneConfig):
    printv(f"Doing: {gene_config.gene}", gene_config.verbose, 2)
    if not gene_config.no_references:
        aa_content = parseFasta(str(gene_config.aa_file))
        nt_content = parseFasta(str(gene_config.nt_file))
    else:
        aa_content = (
            pair for pair in parseFasta(str(gene_config.aa_file)) if pair[0][-1] != "."
        )
        nt_content = (
            pair for pair in parseFasta(str(gene_config.nt_file)) if pair[0][-1] != "."
        )

    if gene_config.stopcodon:
        aa_content, nt_content = stopcodon(aa_content, nt_content)

    if gene_config.rename and not gene_config.generating_names:
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
            gene_config.kick_columns,
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
    nt_target_content = []
    taxon_count = {}
    gene_taxon_to_taxa = {}
    gene_taxa_present = {}

    column_stats = {}
    candidate_count = 0

    if gene_config.gene in gene_config.target or not gene_config.sort:
        if gene_config.count_taxa or gene_config.generating_names:
            taxon_count, gene_taxon_to_taxa, gene_taxa_present = taxon_present(
                aa_content
            )

        col_dict = defaultdict(list)
        for header, sequence in aa_content:
            if not header.endswith("."):
                candidate_count += 1

            for i, char in enumerate(sequence):
                col_dict[i].append(char.replace("X", "-"))

        for col, chars in col_dict.items():
            this_most_common = Counter(chars).most_common()
            most_common_inclusive = this_most_common[0][1]
            if len(this_most_common) == 1:
                if this_most_common[0][0] == "-":
                    most_common_AA_count = 0
                else:
                    most_common_AA_count = this_most_common[0][1]
            else:
                most_common_AA_count = (
                    this_most_common[0][1]
                    if this_most_common[0][0] == "-"
                    else this_most_common[1][1]
                )
            column_stats[col] = (
                len(chars),
                len(chars) - chars.count("-"),
                most_common_inclusive,
                most_common_AA_count,
            )

        aa_target_content.extend(aa_content)
        nt_target_content.extend(nt_content)

        aa_path = str(on_target_aa.joinpath(gene_config.aa_file.name))
        nt_path = str(on_target_nt.joinpath(gene_config.nt_file.name))

    else:
        aa_path = str(off_target_aa.joinpath(gene_config.aa_file.name))
        nt_path = str(off_target_nt.joinpath(gene_config.nt_file.name))

    if gene_config.rename and not gene_config.generating_names:
        aa_content = taxon_only(aa_content)
        nt_content = taxon_only(nt_content)

    if nt_content:
        writeFasta(aa_path, aa_content, gene_config.compress)
        writeFasta(nt_path, nt_content, gene_config.compress)

    path_to = (aa_path, nt_path)

    taxa_local = get_taxa_local(aa_target_content)

    return (
        gene_config.gene,
        taxa_local,
        aa_target_content,
        nt_target_content,
        taxon_count,
        gene_taxon_to_taxa,
        path_to,
        column_stats,
        candidate_count,
        gene_taxa_present,
    )


def get_taxa_local(aa_content: list) -> set:
    taxa_local = set()
    for header, _ in aa_content:
        taxa_local.add(header.split("|")[1])

    return taxa_local


def taxon_present(aa_content: list) -> dict:
    taxonc_present = defaultdict({"p": 0, "bp": 0}.copy)
    gene_taxon_to_taxa = {}
    taxa_present = set()
    for header, sequence in aa_content:
        taxa = header.split("|")[2]
        taxon = header.split("|")[1]

        if not header.endswith("."):
            gene_taxon_to_taxa[taxon] = taxa
            taxa_present.add(taxa)

        taxonc_present[taxon]["p"] += 1
        taxonc_present[taxon]["bp"] += len(sequence) - sequence.count("-")

    return taxonc_present, gene_taxon_to_taxa, taxa_present


def scrape_taxa(taxas):
    result = {}
    for taxa in taxas:
        got_res = False
        try:
            req = requests.get(f"https://www.ncbi.nlm.nih.gov/sra/{taxa}[accn]")
            soup = BeautifulSoup(req.content, "html.parser")
            for anchor in soup.find("div", {"id": "maincontent"}).find_all("a"):
                if anchor.get("href", "").startswith("/Taxonomy/Browser/wwwtax.cgi?"):
                    result[taxa] = anchor.contents[0]
                    got_res = True
                    break
        except:
            pass

        if got_res:
            continue

        try:
            req = requests.get(
                f"https://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&view=all&search={taxa}"
            )
            soup = BeautifulSoup(req.content, "html.parser")
            for anchor in soup.find("table", {"class": "geo_zebra"}).find_all("a"):
                if anchor.get("href", "").lower() == taxa.lower():
                    for anchor_2 in anchor.parent.parent.find_all("a"):
                        if anchor_2.get("href", "").startswith(
                            "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?"
                        ):
                            result[taxa] = anchor_2.contents[0]
                    break
        except:
            pass

    return json.dumps(result)


def process_folder(args, input_path):
    tk = TimeKeeper(KeeperMode.DIRECT)
    taxa_folder = Path(input_path)
    basename = os.path.basename(taxa_folder)
    no_suffix = basename.split(".")[0]
    print(f"Processing: {basename}")
    aa_folder = taxa_folder.joinpath(AA_FOLDER)
    nt_folder = taxa_folder.joinpath(NT_FOLDER)

    def makent(gene, aa_path):
        name = gene + ".nt.fa"
        if ".gz" in str(aa_path):
            name += ".gz"
        return name

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
    if args.kick:
        with open(taxa_folder.joinpath(args.kick), encoding="utf-8-sig") as fp:
            for line in fp:
                to_kick.add(line.strip())

    target = set()
    if args.target_file:
        with open(taxa_folder.joinpath(args.target_file), encoding="utf-8-sig") as fp:
            for line in fp:
                target.add(line.strip())

    taxa_to_taxon = {}
    generate_names = False
    if args.rename or args.count:
        if args.names_csv and os.path.exists(taxa_folder.joinpath(args.names_csv)):
            with open(taxa_folder.joinpath(args.names_csv), encoding="utf-8-sig") as fp:
                for line in fp:
                    if line != "\n":
                        parse = line.strip().split(",")

                        id = parse[0]
                        name = parse[1]

                        taxa_to_taxon[id] = name
        else:
            printv(
                "No Names.csv file provided. Will continue and generate blank names.csv\nWARNING: Rename and Taxa Count will not run.",
                args.verbose,
            )
            generate_names = True

    arguments = []
    for aa_file in aa_folder.glob("*.fa*"):
        gene = aa_file.name.split(".")[0]
        nt_file = nt_folder.joinpath(makent(gene, aa_file))
        this_config = GeneConfig(
            gene,
            args.kick_columns if args.kick_columns <= 1 else args.kick_columns / 100,
            args.minimum_bp,
            args.stopcodon,
            args.rename,
            taxa_folder,
            aa_file,
            nt_file,
            to_kick,
            taxa_to_taxon,
            target,
            args.verbose,
            args.count,
            generate_names,
            args.no_references,
            args.compress,
        )
        arguments.append((this_config,))

    if args.processes > 1:
        with Pool(args.processes) as pool:
            to_write = pool.starmap(clean_gene, arguments, chunksize=1)
    else:
        to_write = [clean_gene(argument[0]) for argument in arguments]

    taxa_global = set()
    to_scrape = set()
    max_candidate_count = 0
    for _, taxa_local, _, _, _, _, _, _, candidate_count, taxa_present in to_write:
        taxa_global.update(taxa_local)
        if generate_names:
            to_scrape.update(taxa_present)

        max_candidate_count = max(max_candidate_count, candidate_count)

    if generate_names:
        printv("Attempting to scrape names.csv", args.verbose)
        to_scrape = list(to_scrape)
        per_thread = ceil(len(to_scrape) / args.processes)
        scrape_args = [
            to_scrape[i : i + per_thread] for i in range(0, len(to_scrape), per_thread)
        ]

        if args.processes > 1:
            with Pool(args.processes) as pool:
                scrape_components = pool.map(scrape_taxa, scrape_args)
        else:
            scrape_components = [scrape_taxa(i) for i in scrape_args]

        scrape = {}
        for component in scrape_components:
            scrape.update(json.loads(component))

        with open(taxa_folder.joinpath("names.csv"), "w") as fp:
            for taxa in to_scrape:
                organism_name = scrape.get(taxa, "").replace(" ", "_")
                fp.write(f"{taxa},{organism_name}\n")

    gene_kick = args.gene_kick if args.gene_kick <= 1 else (args.gene_kick / 100)

    sequences = {"aa": defaultdict(dict), "nt": defaultdict(dict)}
    gene_lengths = defaultdict(dict)
    total = defaultdict({"p": 0, "bp": 0}.copy)
    taxon_to_taxa = {}

    positions = None
    if args.position != {1, 2, 3}:
        positions = [i for i in range(3) if str(i + 1) not in set(str(args.position))]

    total_lax = 0
    total_strict = 0
    total_inform = 0
    total_inform_lax = 0
    column_stats_res = []
    MISMATCHES = 2

    for (
        gene,
        taxa_local,
        aa_content,
        nt_content,
        taxon_count,
        gene_taxon_to_taxa,
        path_to,
        column_stats,
        _,
        _,
    ) in to_write:
        this_lax = 0
        this_strict = 0
        this_inform = 0
        this_inform_lax = 0
        for (
            _,
            non_gap,
            most_common_inclusive,
            most_common_AA_count,
        ) in column_stats.values():
            if non_gap >= max_candidate_count:
                this_strict += 1
                if most_common_inclusive < max_candidate_count - MISMATCHES:
                    this_inform += 1

            if non_gap != 0:
                if non_gap >= 4:
                    if most_common_AA_count < non_gap - MISMATCHES:
                        this_inform_lax += 1

                this_lax += 1

        total_lax += this_lax
        total_strict += this_strict
        total_inform += this_inform
        total_inform_lax += this_inform_lax

        column_stats_res.append(
            (
                this_strict,
                f"{gene},{str(this_lax)},{str(this_strict)},{str(this_inform)},{str(this_inform_lax)}",
            )
        )

        if args.gene_kick:
            if kick_gene(taxa_local, gene_kick, taxa_global):
                aa_path, nt_path = path_to
                aa_glob = aa_path.replace(".gz", "") + "*"
                nt_glob = nt_path.replace(".gz", "") + "*"
                for fasta in glob(aa_glob):
                    os.remove(fasta)
                for fasta in glob(nt_glob):
                    os.remove(fasta)
                continue
        if args.count:
            taxon_to_taxa.update(gene_taxon_to_taxa)
            for taxon, count in taxon_count.items():
                total[taxon]["p"] += count["p"]
                total[taxon]["bp"] += count["bp"]

        if args.concat:
            this_gene_global_length = 0
            for i, (header, sequence) in enumerate(aa_content):
                if i == 0:
                    this_gene_global_length = len(sequence)
                    gene_lengths["aa"][gene] = this_gene_global_length

                taxon = header.split("|")[1]
                sequences["aa"][gene][taxon] = sequence

            for i, (header, sequence) in enumerate(nt_content):
                if positions:
                    sequence = list(sequence)
                    for ai in range(0, len(sequence), 3):
                        for ti in positions:
                            sequence[ai + ti] = ""
                    sequence = "".join(sequence)

                if i == 0:
                    this_gene_global_length = len(sequence)
                    gene_lengths["nt"][gene] = this_gene_global_length

                taxon = header.split("|")[1]
                sequences["nt"][gene][taxon] = sequence

    if args.concat:
        for type_ in ["aa", "nt"]:
            log = {}
            taxa_sequences_global = defaultdict(list)

            for gene_i, gene in enumerate(sequences[type_]):
                this_sequences = sequences[type_][gene]

                if gene_i == 0:
                    start = 1
                    end = gene_lengths[type_][gene]
                else:
                    start = end + 1
                    end = start + (gene_lengths[type_][gene] - 1)

                for taxa in taxa_global:
                    if taxa not in this_sequences:
                        seq = "-" * gene_lengths[type_][gene]
                    else:
                        seq = this_sequences[taxa]

                    taxa_sequences_global[taxa].append(seq)

                log[gene] = (start, end)

            output_fas = processed_folder.joinpath(no_suffix + f".{type_}.fas")
            output_nex = processed_folder.joinpath(no_suffix + f".{type_}.nex")

            taxa_present = sorted(taxa_sequences_global.keys())

            with open(output_fas, "w", encoding="UTF-8") as fp:
                for taxa in taxa_present:
                    taxa_contig_sequence = taxa_sequences_global[taxa]
                    fp.write(">" + taxa + "\n")
                    fp.write("".join(taxa_contig_sequence) + "\n")

            with open(output_nex, "w", encoding="UTF-8") as fp:
                fp.write("#nexus\nbegin sets;\n")
                for gene, (start, end) in log.items():
                    fp.write(f"CHARSET {gene} = {start}-{end} ;\n")

                fp.write("end;\n")

    column_stats_res.sort(key=lambda x: x[0], reverse=True)
    column_stats_res = [
        "Gene,Lax Count,Strict Count,Informative Count,Lax Informative Count"
    ] + [i[1] for i in column_stats_res]

    column_stats_res.append(
        f"Total,{str(total_lax)},{str(total_strict)},{str(total_inform)},{str(total_inform_lax)}"
    )

    with open(str(processed_folder.joinpath("ColumnStats.csv")), "w") as fwcsv:
        fwcsv.write("\n".join(column_stats_res))

    if args.count:
        out = ["Taxa,Taxon,Present,Total AA"]

        for taxon, total_taxa in total.items():
            out.append(
                f"{taxon_to_taxa.get(taxon, '')},{taxon},{total_taxa['p']},{total_taxa['bp']}"
            )

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
