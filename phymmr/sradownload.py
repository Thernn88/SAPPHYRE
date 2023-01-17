import csv
import os
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import Pool
from pathlib import Path
from subprocess import Popen, PIPE

import requests
from bs4 import BeautifulSoup

from phymmr.utils import printv


def download_parallel(arguments):
    command, srr_acession, path_to_download, verbose = arguments
    printv(f"Download {srr_acession} to {path_to_download}...", verbose)
    p = Popen(
        f"{command} {srr_acession} -O {path_to_download}", shell=True, stdout=PIPE
    )
    print(p.stdout.read().decode())


def main(args):
    cmd = "fastq-dump --gzip"
    if args.bin:
        cmd = Path(args.bin, cmd)

    csvfile = Path(args.INPUT)
    # target_folder = csvfile.name.removesuffix(".csv")
    path_to_download = Path(os.getcwd(), csvfile.name.removesuffix(".csv"))

    with open(csvfile, mode="r", encoding="utf-8") as fp:
        csv_read = csv.reader(fp, delimiter=",", quotechar='"')
        arguments = []
        for i, fields in enumerate(csv_read):
            out_fields = ['"{}"'.format(str(i).replace("ï»¿", "")) for i in fields]
            if fields[0] == "Experiment Accession":
                out_fields.append('"SRR Acession"')

            elif fields[0] != "":
                acession = fields[0]
                print(f"Searching for runs in SRA: {acession}")

                url = "https://www.ncbi.nlm.nih.gov/sra/{}[accn]".format(acession)
                # TODO handle network errors
                req = requests.get(url)
                soup = BeautifulSoup(req.content, "html.parser")
                for srr_acession in (
                    a.contents[0]
                    for a in soup.find_all("a", href=True)
                    if a["href"].startswith("//trace.ncbi.nlm.nih.gov/Traces?run")
                ):
                    print(f"Attempting to download: {srr_acession}")
                    out_fields.append(f'"{srr_acession}"')

                    # TODO: verify download is successful
                    # expected_directory = Path(path_to_download, f'{srr_acession}.fastq')
                    arguments.append(
                        (cmd, srr_acession, path_to_download, args.verbose)
                    )

        with ThreadPoolExecutor(args.processes) as pool:
            pool.map(download_parallel, arguments, chunksize=1)
    return True
