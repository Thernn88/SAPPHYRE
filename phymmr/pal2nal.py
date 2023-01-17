from __future__ import annotations
from shutil import rmtree

from multiprocessing import Pool
from pathlib import Path
from typing import Any, Dict, Generator, List, Tuple
from os import path as ospath

from pro2codon import pn2codon

from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv, parseFasta, writeFasta

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


def return_aligned_paths(
    glob_paths_taxa: List[Path],
    glob_paths_genes: List[Path],
    path_aligned: Path,
    specified_dna_table: dict,
    verbose: bool,
    compress: bool,
) -> Generator[Path, Any, Any]:
    for path_nt, path_aa in zip(glob_paths_taxa, glob_paths_genes):
        if not path_nt.is_file() or not path_aa.is_file():

            continue

        nt_gene = str(path_nt.parts[-1]).split(".")[0]
        aa_gene = str(path_aa.parts[-1]).split(".")[0]

        if nt_gene != aa_gene:
            continue

        yield (
            path_aa,
            path_nt,
            path_aligned.joinpath(Path(f"{nt_gene}.nt.fa")),
            specified_dna_table,
            verbose,
            compress,
        )


def prepare_taxa_and_genes(
    input: str, specified_dna_table, verbose, compress
) -> Tuple[Generator[Tuple[Path, Path, Path], Any, Any], int]:
    input_path = Path(input)

    joined_mafft = input_path.joinpath(Path("mafft"))
    joined_nt = input_path.joinpath(Path("nt"))
    joined_nt_aligned = input_path.joinpath(Path("nt_aligned"))

    rmtree(joined_nt_aligned, ignore_errors=True)
    joined_nt_aligned.mkdir()

    glob_aa = sorted(
        [
            i
            for i in joined_mafft.iterdir()
            if i.suffix in [".fa", ".gz", ".fq", ".fastq", ".fasta"]
        ],
        key=lambda x: x.name.split(".")[0],
    )
    glob_nt = sorted(
        [
            i
            for i in joined_nt.iterdir()
            if i.suffix in [".fa", ".gz", ".fq", ".fastq", ".fasta"]
        ],
        key=lambda x: x.name.split(".")[0],
    )

    out_generator = return_aligned_paths(
        glob_nt, glob_aa, joined_nt_aligned, specified_dna_table, verbose, compress
    )

    return out_generator, len(glob_aa)


def read_and_convert_fasta_files(
    aa_file: str, nt_file: str, verbose: bool
) -> Dict[str, Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]]:

    aas = []
    nts = {}

    nt_has_refs = False
    for i, nt in enumerate(parseFasta(nt_file)):
        nt_header = nt[0]
        if i == 0:
            if nt_header[-1] == ".":
                nt_has_refs = True

        nts[nt_header] = nt

    for aa_header, aa_seq in parseFasta(aa_file):
        if not nt_has_refs and aa_header[-1] == ".":
            continue

        aa_header, aa_seq = aa_header.strip(), aa_seq.strip()

        aas.append((aa_header, aa_seq))

    result = {}
    i = -1
    for header, sequence in aas:
        try:
            if header not in result:
                i += 1
            result[header] = ((header, sequence), (i, *nts[header]))
        except KeyError as e:
            printv(
                "ERROR CAUGHT: There is a single header in PEP sequence FASTA file that does not exist in NUC sequence FASTA file",
                verbose,
                0,
            )
            printv(f"SEQUENCE MISSING: {e}", verbose, 0)
            return False
    return result


def worker(
    aa_file: str,
    nt_file: str,
    out_file: Path,
    specified_dna_table: dict,
    verbose: bool,
    compress: bool,
) -> bool:
    gene = ospath.basename(aa_file).split(".")[0]
    printv(f"Doing: {gene}", verbose, 2)
    seqs = read_and_convert_fasta_files(aa_file, nt_file, verbose)

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


def run_batch_threaded(
    num_threads: int, ls: List[List[List[Tuple[Tuple[Path, Path, Path], Dict, bool]]]]
):
    with Pool(num_threads) as pool:
        return all(pool.starmap(worker, ls, chunksize=100))


def main(args):
    specified_dna_table = DICT_TABLES[str(args.table)]
    time_keeper = TimeKeeper(KeeperMode.DIRECT)
    for folder in args.INPUT:
        printv(f"Processing: {ospath.basename(folder)}", args.verbose, 0)
        this_taxa_jobs, _ = prepare_taxa_and_genes(
            folder, specified_dna_table, args.verbose, args.compress
        )

        success = run_batch_threaded(num_threads=args.processes, ls=this_taxa_jobs)

        if not success:
            printv("A fatal error has occured.", args.verbose, 0)
            return False
        printv(f"Done! Took {time_keeper.lap():.2f}s", args.verbose, 1)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {time_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception("Cannot be called directly, please use the module:\nphymmr Pal2Nal")
