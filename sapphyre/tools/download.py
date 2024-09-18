import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from subprocess import PIPE, Popen

import openpyxl
import requests
from bs4 import BeautifulSoup

from ..utils import printv


def download_parallel_srr(arguments):
    command, srr_acession, path_to_download, verbose = arguments
    printv(f"Attempting to download: {srr_acession}", verbose)
    printv(f"Download {srr_acession} to {path_to_download}...", verbose)

    with Popen(
        f"{command} {srr_acession} -O {path_to_download}",
        shell=True,
        stdout=PIPE,
    ) as p:
        try:
            print(p.stdout.read().decode())
        except UnicodeDecodeError:
            print(
                "ErrUnicode decoding error, UTF-8 charset does not contain the bytecode for gotten character",
            )
            sys.exit(1)


def download_parallel_wgs(arguments):
    download_link, prefix, path_to_download, verbose = arguments
    printv(f"Attempting to download: {prefix}", verbose)

    file_name = download_link.split("/")[-1]

    printv(f"Download {file_name} to {path_to_download}...", verbose)
    download = requests.get(download_link)

    path = os.path.join(path_to_download, file_name)

    open(path, "wb").write(download.content)


def search_and_download(arguments):
    cmd, acession, path_to_download, verbose, wgs_flag = arguments
    if wgs_flag:
        url = f"https://www.ncbi.nlm.nih.gov/Traces/wgs/{acession}"
    else:
        url = f"https://www.ncbi.nlm.nih.gov/sra/{acession}[accn]"

    req = requests.get(url)
    soup = BeautifulSoup(req.content, "html.parser")
    
    if wgs_flag:
        em = soup.find("em", text="FASTA:")
        if not em:
            print(f"Failed to find FASTA: {acession}")
            return
        container = em.find_parent()
        for a in container.find_all("a", href=True):
            if "sra-download.ncbi.nlm.nih.gov" in a["href"]:
                download_parallel_wgs((a["href"], acession, path_to_download, verbose))
    else:
        for srr_acession in (
            a.contents[0]
            for a in soup.find_all("a", href=True)
            if a["href"].startswith("//trace.ncbi.nlm.nih.gov/Traces?run")
        ):
            download_parallel_srr((cmd, srr_acession, path_to_download, verbose))


def main(args):
    cmd = "fastq-dump --split-3 --gzip"
    if args.bin:
        cmd = Path(args.bin, cmd)

    csvfile = Path(args.INPUT)
    this_suffix = csvfile.suffix

    path_to_download = Path(
        os.getcwd(), "input", csvfile.name.removesuffix(this_suffix)
    )
    os.makedirs(path_to_download, exist_ok=True)
    
    arguments = []
    wgs_flag = args.wgs
    if this_suffix == ".csv":
        with open(csvfile, encoding="utf-8") as fp:
            csv_read = csv.reader(fp, delimiter=",", quotechar='"')
            for i, fields in enumerate(csv_read):
                if i == 0:
                    continue
                acession = fields[0]
                
                arguments.append((cmd, acession, path_to_download, args.verbose, wgs_flag))
    elif "xls" in this_suffix:
        workbook = openpyxl.load_workbook(csvfile)
        sheet = workbook.active
        for row in sheet.iter_rows(values_only=True):
            fields = list(row)
            acession = fields[0]
            arguments.append((cmd, acession, path_to_download, args.verbose, wgs_flag))

    with ThreadPoolExecutor(args.processes) as pool:
        pool.map(search_and_download, arguments, chunksize=1)

    return True