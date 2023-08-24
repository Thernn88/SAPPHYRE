# © 2022 GPLv3+ Sapphyre Team
import gzip
import os
from collections.abc import Generator
from queue import Queue
from threading import Thread

import pyfastx
from needletail import parse_fastx_file

class ConcurrentLogger(Thread):
    def __init__(self, inq: Queue) -> None:
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


def parseFasta(
    path: str,
    has_interleave=False,
) -> Generator[tuple[str, str], None, None]:
    """Iterate over a Fasta file returning sequence records as string tuples.
    Designed in order to handle .gz and .fasta files with potential interleave.
    """
    if has_interleave:
        fa = pyfastx.Fastx(
            str(path),
            uppercase=True,
        )
        for header, seq in fa:
            yield header, seq
    elif str(path).rsplit(".", maxsplit=1)[-1] in {"fastq", "fq", "gz"}:
        fa = parse_fastx_file(str(path))
        for entry in fa:  # Deinterleave
            yield entry.id, entry.seq
    else:
        with open(path) as fp:
            header = None
            for line in fp:
                if line.startswith(">"):
                    header = line.lstrip(">").strip()
                else:
                    if line.strip():
                        yield header, line.strip()
                        header = None


def writeFasta(path: str, records: tuple[str, str], compress=False):
    """Writes sequence records to a Fasta format file."""
    func = open
    path = str(path)
    if compress:
        if not path.endswith(".gz"):
            path += ".gz"
        func = gzip.open
    elif path.endswith(".gz"):
        path = path.rstrip(".gz")

    with func(path, "wb") as fp:
        data = [f">{header}\n{sequence}\n" for header, sequence in records]
        fp.write("".join(data).encode())

def write2Line2Fasta(path: str, records: list[str], compress=False):
    func = open
    if compress:
        path += ".gz"
        func = gzip.open

    with func(path, "wb") as fp:
        for i in range(0, len(records), 2):
            fp.write(f"{records[i]}\n{records[i+1]}\n".encode())
