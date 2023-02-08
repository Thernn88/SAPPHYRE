# -*- coding: utf-8 -*-
# © 2022 GPLv3+ PhyMMR Team
import gzip
import os
import pyfastx
from threading import Thread
from queue import Queue
from typing import Generator
# from Bio.SeqIO.QualityIO import FastqGeneralIterator


class ConcurrentLogger(Thread):
    def __init__(self, inq: Queue):
        super().__init__(daemon=True)
        self.inq = inq

    def run(self):
        while True:
            message, verbosity, reqv = self.inq.get()
            printv(message, verbosity, reqv)
            self.inq.task_done()

    def __call__(self, msg: str, verbosity: int, reqv=1):
        self.inq.put((msg, verbosity, reqv))


def printv(msg, verbosity, reqverb=1) -> None:
    if verbosity >= reqverb:
        print(msg)


def gettempdir():
    if os.path.exists("/run/shm"):
        return "/run/shm"
    if os.path.exists("/dev/shm"):
        return "/dev/shm"
    return None

def parseFasta(path: str) -> Generator[tuple[str, str], None, None]:
    """
    Iterate over a Fasta file returning sequence records as string tuples.
    Designed in order to handle .gz and .fasta files with potential interleave.
    """
    fa = pyfastx.Fastx(path, uppercase=True)
    for header, seq in fa:
        yield header, seq


def writeFasta(path: str, records: tuple[str, str], compress=False):
    """
    Writes sequence records to a Fasta format file.
    """
    func = open
    if compress:
        path += ".gz"
        func = gzip.open

    with func(path, "wb") as fp:
        for header, sequence in records:
            fp.write(f">{header}\n{sequence}\n".encode())


def write2Line2Fasta(path: str, records: list[str], compress=False):
    func = open
    if compress:
        path += ".gz"
        func = gzip.open

    with func(path, "wb") as fp:
        for i in range(0, len(records), 2):
            fp.write(f"{records[i]}\n{records[i+1]}\n".encode())
