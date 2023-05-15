import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from subprocess import PIPE, Popen

import requests
from bs4 import BeautifulSoup

from .utils import printv


def download_parallel(arguments):
    command, srr_acession, path_to_download, verbose = arguments
    printv(f"Download {srr_acession} to {path_to_download}...", verbose)

    with Popen(
        f"{command} {srr_acession} -O {path_to_download}", shell=True, stdout=PIPE,
    ) as p:
        try:
            print(p.stdout.read().decode())
        except UnicodeDecodeError:
            print(
                "ErrUnicode decoding error, UTF-8 charset does not contain the bytecode for gotten character",
            )
            sys.exit(1)


def main(args):
    cmd = "fastq-dump --gzip"
    if args.bin:
        cmd = Path(args.bin, cmd)

    csvfile = Path(args.INPUT)
    path_to_download = Path(os.getcwd(), "input", csvfile.name.removesuffix(".csv"))
    os.makedirs(path_to_download, exist_ok=True)

    with open(csvfile, encoding="utf-8") as fp:
        csv_read = csv.reader(fp, delimiter=",", quotechar='"')
        arguments = []
        for i, fields in enumerate(csv_read):
            out_fields = ['"{}"'.format(str(i).replace("ï»¿", "")) for i in fields]
            if fields[0] == "Experiment Accession":
                out_fields.append('"SRR Acession"')

            elif fields[0] != "":
                acession = fields[0]
                print(f"Searching for runs in SRA: {acession}")

                url = f"https://www.ncbi.nlm.nih.gov/sra/{acession}[accn]"
                # TODO handle network errors
                req = requests.get(url, timeout=500)
                soup = BeautifulSoup(req.content, "html.parser")
                for srr_acession in (
                    a.contents[0]
                    for a in soup.find_all("a", href=True)
                    if a["href"].startswith("//trace.ncbi.nlm.nih.gov/Traces?run")
                ):
                    print(f"Attempting to download: {srr_acession}")
                    out_fields.append(f'"{srr_acession}"')

                    # TODO: verify download is successful
                    arguments.append(
                        (cmd, srr_acession, path_to_download, args.verbose),
                    )

        with ThreadPoolExecutor(args.processes) as pool:
            pool.map(download_parallel, arguments, chunksize=1)
    return True
