# © 2022 GPLv3+ Sapphyre Team
from isal import igzip as isal_gzip
import os
from collections.abc import Generator
from queue import Queue
from threading import Thread

import sapphyre_tools
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
    if has_interleave or str(path).rsplit(".", maxsplit=1)[-1] in {"fastq", "fq", "gz"}:
        # fa = pyfastx.Fastx(
        #     str(path),
        #     uppercase=True,
        # )
        # for tup in fa:
        #     yield tup[0], tup[1]
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
        func = isal_gzip.open
    elif path.endswith(".gz"):
        path = path.rstrip(".gz")

    with func(path, "wb") as fp:
        data = [f">{header}\n{sequence}\n" for header, sequence in records]
        fp.write("".join(data).encode())


def write2Line2Fasta(path: str, records: list[str], compress=False):
    func = open
    if compress:
        path += ".gz"
        func = isal_gzip.open

    with func(path, "wb") as fp:
        for i in range(0, len(records), 2):
            fp.write(f"{records[i]}\n{records[i+1]}\n".encode())


def find_gap_regions(consensus: str, gap="-", min_length=6) -> str:
    start = None
    indices = []
    for i, letter in enumerate(consensus):
        if letter == gap:
            if start is None:
                start = i
                continue
        else:
            if start:
                if i - start + 1> min_length:
                    indices.extend(range(start, i))
                start = None
    if start:
        if len(consensus) - start + 1 > min_length:
            indices.extend(range(start, len(consensus)))
    return indices


def cull_columns(sequences: list[str], consensus_threshold: 0.6, min_length=6) -> list[str]:
    unmasked = sapphyre_tools.dumb_consensus(sequences, consensus_threshold, 0)
    mask = set(find_gap_regions(unmasked, gap="-", min_length=min_length))
    output = []
    for seq in sequences:
        output.append("".join(seq[i] for i in range(len(unmasked)) if i not in mask))
    return output