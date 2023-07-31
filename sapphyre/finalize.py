from collections import Counter, defaultdict
import json
from math import ceil
import os
from dataclasses import dataclass
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import rmtree
from bs4 import BeautifulSoup

import requests

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
    kick_columns: float
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
    generating_names: bool


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

    for header, sequence in fasta_content:
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
    return [(header.split("|")[1],sequence) for header, sequence in content]


def kick_gene(present_taxa, minimum_percentage, global_total_taxon):
    if len(present_taxa) == 0:
        return 0 <= minimum_percentage

    return (len(present_taxa) / len(global_total_taxon)) <= minimum_percentage

def clean_gene(gene_config: GeneConfig):
    printv(f"Doing: {gene_config.gene}", gene_config.verbose, 2)
    aa_content = parseFasta(str(gene_config.aa_file))
    nt_content = parseFasta(str(gene_config.nt_file))

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
    consensus = set()

    if (gene_config.gene in gene_config.target or not gene_config.sort):
        if gene_config.count_taxa or gene_config.generating_names:
            taxon_count, gene_taxon_to_taxa, consensus = taxon_present(aa_content)

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

    writeFasta(aa_path, aa_content)
    writeFasta(nt_path, nt_content)

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
        consensus,
    )


def get_taxa_local(aa_content: list) -> set:
    taxa_local = set()
    for header, _ in aa_content:
        taxa_local.add(header.split("|")[1])

    return taxa_local


def taxon_present(aa_content: list) -> dict:
    taxonc_present = defaultdict({"p": 0, "bp": 0}.copy)
    gene_taxon_to_taxa = {}
    consensus = set()
    for header, sequence in aa_content:
        taxa = header.split("|")[2]
        taxon = header.split("|")[1]

        consensus.add(taxa)
        
        if not header.endswith("."):
            gene_taxon_to_taxa[taxon] = taxa

        taxonc_present[taxon]["p"] += 1
        taxonc_present[taxon]["bp"] += len(sequence) - sequence.count("-")

    return taxonc_present, gene_taxon_to_taxa, consensus


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
            req = requests.get(f"https://www.ncbi.nlm.nih.gov/nuccore/?term={taxa}", allow_redirects=True)
            if req.url.startswith("https://www.ncbi.nlm.nih.gov/nuccore/"):
                if "term=" in req.url:
                    soup = BeautifulSoup(req.content, "html.parser")
                    for anchor in soup.find_all("a", {"class":"dblinks"}):
                        if anchor.get("href", "").startswith("/taxonomy?LinkName=nuccore_taxonomy&from_uid="):
                            req = requests.get(f"https://www.ncbi.nlm.nih.gov{anchor.get('href')}")
                            break
                else:
                    id = req.url.split("/")[-1]
                    req = requests.get(f"https://www.ncbi.nlm.nih.gov/taxonomy?LinkName=nuccore_taxonomy&from_uid={id}")
            
            got_res = False
            soup = BeautifulSoup(req.content, "html.parser")
            for anchor in soup.find_all("a"):
                if anchor.get("href", "").startswith("/Taxonomy/Browser/wwwtax.cgi?"):
                    result[taxa] = anchor.contents[0]
                    got_res = True
                    break

            if not got_res:
                #Try again
                req = requests.get(req.url, allow_redirects=True)
                soup = BeautifulSoup(req.content, "html.parser")
                for anchor in soup.find_all("a"):
                    if anchor.get("href", "").startswith("/Taxonomy/Browser/wwwtax.cgi?"):
                        result[taxa] = anchor.contents[0]
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
            printv("No Names.csv file provided. Will continue and generate blank names.csv\nWARNING: Rename and Taxa Count will not run.", args.verbose)
            generate_names = True

    arguments = []
    for aa_file in aa_folder.glob("*.fa"):
        gene = aa_file.name.split(".")[0]
        nt_file = nt_folder.joinpath(makent(gene))
        this_config = GeneConfig(
            gene,
            args.kick_columns if args.kick_columns <= 1 else args.kick_columns / 100,
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
            generate_names,
        )
        arguments.append((this_config,))

    if args.processes > 1:
        with Pool(args.processes) as pool:
            to_write = pool.starmap(clean_gene, arguments, chunksize=1)
    else:
        to_write = [clean_gene(argument[0]) for argument in arguments]

    taxa_global = set()
    to_scrape = set()
    for _, taxa_local, _, _, _, gene_taxon_to_taxa, _, consensus in to_write:
        taxa_global.update(taxa_local)
        if generate_names:
            to_scrape.update(consensus)

    if generate_names:
        printv("Attempting to scrape names.csv", args.verbose)
        to_scrape = list(to_scrape)
        per_thread = ceil(len(to_scrape) / args.processes)
        scrape_args = [to_scrape[i:i+per_thread] for i in range(0, len(to_scrape), per_thread)]

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

    for gene, taxa_local, aa_content, nt_content, taxon_count, gene_taxon_to_taxa, path_to, _ in to_write:
        if args.gene_kick:
            if kick_gene(taxa_local, gene_kick, taxa_global):
                aa_path, nt_path = path_to
                os.remove(aa_path) # TODO Can do this better
                os.remove(nt_path)
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

            positions = None
            if args.position != {1,2,3}:
                positions = list(map(int, args.position))

            with open(output_fas, "w", encoding="UTF-8") as fp:
                for taxa, taxa_contig_sequence in taxa_sequences_global.items():
                    if positions and type_ == "nt":
                        taxa_contig_sequence = "".join(taxa_contig_sequence)
                        triplets = [taxa_contig_sequence[i:i+3] for i in range(0, len(taxa_contig_sequence), 3)]
                        taxa_contig_sequence = ["".join([triplet[i-1] for i in positions]) for triplet in triplets]
                    fp.write(">" + taxa + "\n")
                    fp.write("".join(taxa_contig_sequence) + "\n")

            with open(output_nex, "w", encoding="UTF-8") as fp:
                fp.write("#nexus\nbegin sets;\n")
                for gene in log:
                    start, end = log[gene]
                    fp.write(f"CHARSET {gene} = {start}-{end} ;\n")

                fp.write("end;\n")
    if args.count:
        out = ["Taxa,Taxon,Present,Total AA"]

        for taxon, total_taxa in total.items():
            out.append(f"{taxon_to_taxa.get(taxon, '')},{taxon},{total_taxa['p']},{total_taxa['bp']}")

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
