from __future__ import annotations

from collections.abc import Generator
from multiprocessing import Pool
from os import path as ospath
from pathlib import Path
from shutil import rmtree
from typing import Any

from pr2codon import pn2codon

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta
from phymmr_tools import find_index_pair


def do_folder(
    num_threads: int,
    input: str,
    specified_dna_table: dict,
    verbose: int,
    compress: bool,
) -> bool:
    """
    Runs pr2codon in parallel on all genes for a given input

    Args:
    ----
        num_threads (int): Number of threads to use
        input (str): Path to input folder
        specified_dna_table (dict): DNA table to use
        verbose (int): Verbosity level
        compress (bool): Whether to compress output
    Returns:
    -------
        bool: Whether all genes were processed successfully
    """

    input_path = Path(input)

    joined_align = input_path.joinpath(Path("align"))
    joined_nt_aligned = input_path.joinpath(Path("nt_aligned"))

    rmtree(joined_nt_aligned, ignore_errors=True)
    joined_nt_aligned.mkdir()

    aa_genes = sorted(
        [
            i
            for i in joined_align.iterdir()
            if i.suffix in [".fa", ".gz", ".fq", ".fastq", ".fasta"]
        ],
        key=lambda x: x.name.split(".")[0],
    )

    arguments = []
    for aa_path in aa_genes:
        aa_path = str(aa_path)
        nt_path = make_nt(aa_path)

        arguments.append(
            [
                aa_path,
                nt_path,
                Path(input, "nt_aligned", ospath.basename(nt_path)),
                specified_dna_table,
                verbose,
                compress,
            ]
        )

    printv(f"Aligning NT Files", verbose)
    with Pool(num_threads) as pool:
        result = pool.starmap(worker, arguments, chunksize=100)

    return all(result)


def make_nt(aa_path: str) -> str:
    return aa_path.replace(".aa.", ".nt.").replace("/align/", "/nt/")


def worker(
    aa_file: str,
    nt_file: str,
    out_file: Path,
    specified_dna_table: dict,
    verbose: bool,
    compress: bool,
) -> bool:
    """
    Runs pr2codon on a single gene and writes the result to a file

    Args:
    ----
        aa_file (str): Path to amino acid file
        nt_file (str): Path to nucleotide file
        out_file (Path): Path to output file
        specified_dna_table (dict): DNA table to use
        verbose (int): Verbosity level
        compress (bool): Whether to compress output
    Returns:
    -------
        bool: Whether the gene was processed successfully
    """
    gene = ospath.basename(aa_file).split(".")[0]
    printv(f"Doing: {gene}", verbose, 2)
    if not ospath.exists(aa_file):
        printv(f"ERROR CAUGHT in {gene}: AA file does not exist", verbose, 0)
        return False
    if not ospath.exists(nt_file):
        printv(f"ERROR CAUGHT in {gene}: NT file does not exist", verbose, 0)
        return False
    seqs = read_and_convert_fasta_files(gene, aa_file, nt_file, verbose)

    if seqs is False:
        return False

    outfile_stem = out_file.stem
    aligned_result = pn2codon(outfile_stem, specified_dna_table, seqs)
    aligned_lines = aligned_result.split("\n")

    records = []
    for i in range(0, len(aligned_lines) - 1, 2):
        records.append((aligned_lines[i].lstrip(">"), aligned_lines[i + 1]))

    writeFasta(str(out_file), records, compress)

    return True


def read_and_convert_fasta_files(
    gene: str,
    aa_file: str,
    nt_file: str,
    verbose: bool,
) -> dict[str, tuple[list[tuple[str, str]], list[tuple[str, str]]]]:
    """
    Grabs the aa and nt sequence pairs while also checking whether the NT file has refs or if any of the headers are missing in either files

    Args:
    ----
        gene (str): Name of gene
        aa_file (str): Path to amino acid file
        nt_file (str): Path to nucleotide file
        verbose (int): Verbosity level
    Returns:
    -------
        dict[str, tuple[list[tuple[str, str]], list[tuple[str, str]]]]: Dictionary of headers to tuples of (AA header, AA sequence) and (NT header, NT sequence)
    """
    aas = []
    nts = {}

    nt_has_refs = False
    for i, nt in enumerate(parseFasta(nt_file)):
        nt_header = nt[0]
        if i == 0 and nt_header[-1] == ".":
            nt_has_refs = True

        nts[nt_header] = nt
    # sort shouldn't affect refs, so they need a seperate list
    aa_intermediate = []

    for aa_header, aa_seq in parseFasta(aa_file, True):
        if aa_header[-1] == ".":
            aas.append((aa_header.strip(), aa_seq.strip()))
        else:
            aa_intermediate.append((aa_header.strip(), aa_seq.strip()))
    # sort aa candidates by first and last data bp
    aa_intermediate.sort(key=lambda x: find_index_pair(x[1], "-"))
    # add candidates to references
    aas.extend(aa_intermediate)

    with open(aa_file, "w") as fa:
        for header, sequence in aas:
            fa.write(">" + header + "\n")
            fa.write(sequence + "\n")

    # if no nt refs found, do not process aa refs
    aa_final = aas
    if not nt_has_refs:
        aa_final = aa_intermediate

    # Nt has header AA doesn't
    if len(nts) > len(aa_final):
        printv(
            f"ERROR CAUGHT in {gene}: There are more headers in NUC sequence FASTA ({len(nts)}) file than in PEP sequence FASTA file ({len(aa_final)})",
            verbose,
            0,
        )
        return False

    result = {}
    i = -1
    for header, sequence in aa_final:
        try:
            if header not in result:
                i += 1
            result[header] = ((header, sequence), (i, *nts[header]))
        except KeyError as e:
            printv(
                f"ERROR CAUGHT in {gene}: There is a single header in PEP sequence FASTA file that does not exist in NUC sequence FASTA file",
                verbose,
                0,
            )
            printv(f"SEQUENCE MISSING: {e}", verbose, 0)
            return False
    return result


def main(args):
    DICT_TABLES = {
        "1": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG", "TGA"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "2": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG", "AGA", "AGG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC"],
            "M": ["ATA", "ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "3": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "T": ["CTT", "CTC", "CTA", "CTG", "ACT", "ACC", "ACA", "ACG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC"],
            "M": ["ATA", "ATG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "ATT", "ATC"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "4": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "5": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC"],
            "M": ["ATA", "ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "6": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "Q": ["TAA", "TAG", "CAA", "CAG"],
            "C": ["TGT", "TGC"],
            "*": ["TGA"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["TAA", "TAG", "CAA", "CAG", "GAA", "GAG"],
        },
        "9": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC", "AAA"],
            "K": ["AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "AAA", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "10": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC", "TGA"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "11": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG", "TGA"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "12": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA"],
            "S": ["TCT", "TCC", "TCA", "TCG", "CTG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG", "TGA"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "13": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC"],
            "M": ["ATA", "ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "G": ["AGA", "AGG", "GGT", "GGC", "GGA", "GGG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "14": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG"],
            "Y": ["TAT", "TAC", "TAA"],
            "*": ["TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC", "AAA"],
            "K": ["AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "AAA", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "15": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TGA"],
            "Q": ["TAG", "CAA", "CAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["TAG", "CAA", "CAG", "GAA", "GAG"],
        },
        "16": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "TAG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TGA"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "TAG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "21": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA", "AGG"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC"],
            "M": ["ATA", "ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC", "AAA"],
            "K": ["AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "AAA", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "22": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "TAG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCG", "AGT", "AGC"],
            "*": ["TCA", "TAA", "TGA"],
            "Y": ["TAT", "TAC"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "TAG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "23": {
            "F": ["TTT", "TTC"],
            "*": ["TTA", "TAA", "TAG", "TGA"],
            "L": ["TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "24": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG", "AGG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "25": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "G": ["TGA", "GGT", "GGC", "GGA", "GGG"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "26": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG", "TGA"],
            "C": ["TGT", "TGC"],
            "W": ["TGG"],
            "A": ["CTG", "GCT", "GCC", "GCA", "GCG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "27": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "Q": ["TAA", "TAG", "CAA", "CAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["TAA", "TAG", "CAA", "CAG", "GAA", "GAG"],
        },
        "28": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "Q": ["TAA", "TAG", "CAA", "CAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["TAA", "TAG", "CAA", "CAG", "GAA", "GAG"],
        },
        "29": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC", "TAA", "TAG"],
            "C": ["TGT", "TGC"],
            "*": ["TGA"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "30": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "E": ["TAA", "TAG", "GAA", "GAG"],
            "C": ["TGT", "TGC"],
            "*": ["TGA"],
            "W": ["TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["TAA", "TAG", "GAA", "GAG", "CAA", "CAG"],
        },
        "31": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "E": ["TAA", "TAG", "GAA", "GAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["TAA", "TAG", "GAA", "GAG", "CAA", "CAG"],
        },
        "32": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TGA"],
            "W": ["TAG", "TGG"],
            "C": ["TGT", "TGC"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
        "33": {
            "F": ["TTT", "TTC"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGA"],
            "Y": ["TAT", "TAC", "TAA"],
            "*": ["TAG"],
            "C": ["TGT", "TGC"],
            "W": ["TGA", "TGG"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "H": ["CAT", "CAC"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG"],
            "I": ["ATT", "ATC", "ATA"],
            "M": ["ATG"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "N": ["AAT", "AAC"],
            "K": ["AAA", "AAG", "AGG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "X": [],
            "B": ["AAT", "AAC", "GAT", "GAC"],
            "J": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],
            "Z": ["CAA", "CAG", "GAA", "GAG"],
        },
    }

    specified_dna_table = DICT_TABLES[str(args.table)]
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    for folder in args.INPUT:
        printv(f"Processing: {ospath.basename(folder)}", args.verbose, 0)
        success = do_folder(
            args.processes,
            folder,
            specified_dna_table,
            args.verbose,
            args.compress,
        )

        if not success:
            printv(f"A fatal error has occured in {folder}.", args.verbose, 0)
            return False
        printv(f"Done! Took {time_keeper.lap():.2f}s", args.verbose, 1)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {time_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre Pal2Nal"
    raise Exception(
        msg,
    )
