# -*- coding: utf-8 -*-
# © 2022 GPLv3+ PhyMMR Team
import gzip
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from pathlib import Path
from threading import Thread
from queue import Queue
from typing import Generator

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
    elif os.path.exists("/dev/shm"):
        return "/dev/shm"
    return None


def get_records(fp, type: str) -> Generator[tuple[str, str], None, None]:
    """
    Iterates over every line of a file and returns each sequence record.
    Forces sequences to be uppercase.
    """
    current_header = None
    if type == "fasta":
        for line in fp:
            if line.startswith(b">"):
                if current_header:
                    yield (
                        current_header[1:].decode(),
                        b"".join(this_sequence)
                        .replace(b" ", b"")
                        .replace(b"\r", b"")
                        .upper()
                        .decode(),
                    )
                this_sequence = []
                current_header = line.rstrip()
                continue

            this_sequence.append(line.rstrip())

        yield (
            current_header[1:].decode(),
            b"".join(this_sequence)
            .replace(b" ", b"")
            .replace(b"\r", b"")
            .upper()
            .decode(),
        )
    ### FIXME: This breaks if the name lines don't include length
    ### currently replaced with FastqGeneralIterator call in parseFasta
    elif type == "fastq":
        generator = FastqGeneralIterator(fp)
        for header, seq, description in generator:
            yield header, seq
    #     for line in fp:
    #         if line.startswith(b"@") and b"length=" in line:
    #             sequence = next(fp).rstrip()
    #             yield (
    #                 line[1:].rstrip().decode(),
    #                 sequence.replace(b" ", b"").replace(b"\r", b"").upper().decode(),
    #             )


def parseFasta(path: str) -> Generator[tuple[str, str], None, None]:
    """
    Iterate over a Fasta file returning sequence records as string tuples.
    Designed in order to handle .gz and .fasta files with potential interleave.
    """
    func = open
    suffixes = Path(path).suffixes
    if ".gz" in suffixes:
        func = gzip.open
    file_type = "fastq" if ".fastq" in suffixes or ".fq" in suffixes else "fasta"
    mode = 'rb'
    if file_type == 'fastq':
        mode = 'rt'
    return get_records(func(path, mode), file_type)


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
