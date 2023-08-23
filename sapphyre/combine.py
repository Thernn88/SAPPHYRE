from __future__ import annotations
from math import ceil

from multiprocessing.pool import Pool
from pathlib import Path

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta


def prepend(inputs, directory):
    return [Path(directory, i) for i in inputs]


def parse_gene(path, already_grabbed_references):
    this_out = []
    for header, sequence in parseFasta(str(path)):
        if header.endswith("."):
            fields = header.split("|")
            if fields[2] not in already_grabbed_references:
                this_out.append((header, sequence.replace("-", "")))
                already_grabbed_references.add(fields[2])
        else:
            this_out.append((header, sequence.replace("-", "")))
    return {path: (this_out, already_grabbed_references)}


def write_gene(path, gene_sequences, compress):
    for gene, sequences in gene_sequences:
        gene_out = Path(path, gene)
        writeFasta(gene_out, sequences, compress)


def main(args):
    main_keeper = TimeKeeper(KeeperMode.DIRECT)
    if args.DIRECTORY:
        inputs = prepend(args.INPUT, args.DIRECTORY)
    else:
        inputs = [Path(i) for i in args.INPUT]

    aa_out = {}
    nt_out = {}
    grabbed_aa_references = {}
    grabbed_nt_references = {}

    with Pool(args.processes) as pool:
        for item in inputs:
            printv(f"Merging directory: {item}", args.verbose, 0)
            for taxa in item.iterdir():
                aa_path = Path(taxa, "aa_merged")
                nt_path = Path(taxa, "nt_merged")
                if not aa_path.exists() or not nt_path.exists():
                    printv(
                        f"WARNING: Either {aa_path} or {nt_path} doesn't exists. Abort",
                        args.verbose,
                        0,
                    )
                    continue

                arguments = []
                for aa_gene in aa_path.iterdir():
                    if aa_gene.name not in aa_out:
                        grabbed_aa_references[aa_gene.name] = set()

                    arguments.append(
                        (
                            aa_gene,
                            grabbed_aa_references[aa_gene.name],
                        ),
                    )

                aa_result = pool.starmap(parse_gene, arguments, chunksize=8)

                for result in aa_result:
                    for path, out in result.items():
                        out, already_grabbed = out
                        path = path.name
                        aa_out.setdefault(path, []).extend(out)
                        grabbed_aa_references[path] = already_grabbed

                arguments = []
                for nt_gene in nt_path.iterdir():
                    if nt_gene.name not in nt_out:
                        grabbed_nt_references[nt_gene.name] = set()
                    arguments.append(
                        (
                            nt_gene,
                            grabbed_nt_references[nt_gene.name],
                        ),
                    )

                nt_result = pool.starmap(parse_gene, arguments, chunksize=8)

                for result in nt_result:
                    for path, out in result.items():
                        out, already_grabbed = out
                        path = path.name
                        nt_out.setdefault(path, []).extend(out)
                        grabbed_nt_references[path] = already_grabbed

    aa_out_path = Path(args.output_directory, "aa")
    nt_out_path = Path(args.output_directory, "nt")
    aa_out_path.mkdir(parents=True, exist_ok=True)
    nt_out_path.mkdir(parents=True, exist_ok=True)

    aa_sequences = list(aa_out.items())
    nt_sequences = list(nt_out.items())

    if aa_sequences:

        per_thread = ceil(len(aa_sequences) / args.processes)

        aa_arguments = [(aa_out_path, aa_sequences[i:i + per_thread], args.compress) for i in range(0, len(aa_sequences), per_thread)]

        with Pool(args.processes) as pool:
            pool.starmap(write_gene, aa_arguments)
        del aa_arguments

        nt_arguments = [(nt_out_path, nt_sequences[i:i + per_thread], args.compress) for i in range(0, len(nt_sequences), per_thread)]

        with Pool(args.processes) as pool:
            pool.starmap(write_gene, nt_arguments)
        del nt_arguments

    printv(f"Finished took {main_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True
