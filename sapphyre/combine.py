from __future__ import annotations
from pathlib import Path
from multiprocessing.pool import Pool

from .utils import printv, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode


def prepend(inputs, directory):
    return [Path(directory, i) for i in inputs]


def parse_gene(path, already_grabbed_references):
    this_out = []
    for header, sequence in parseFasta(str(path)):
        if header.endswith("."):
            fields = header.split("|")
            if fields[2] not in already_grabbed_references:
                this_out.append((header, sequence))
                already_grabbed_references.add(fields[2])
        else:
            this_out.append((header, sequence))
    return {path: (this_out, already_grabbed_references)}


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
                    )
                )

            with Pool(args.processes) as pool:
                aa_result = pool.starmap(parse_gene, arguments, chunksize=1)

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
                    )
                )

            with Pool(args.processes) as pool:
                nt_result = pool.starmap(parse_gene, arguments, chunksize=1)

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

    for gene, aa_sequence in aa_out.items():
        gene_out = Path(aa_out_path, gene)
        writeFasta(gene_out, aa_sequence, args.compress)

    for gene, nt_sequence in nt_out.items():
        gene_out = Path(nt_out_path, gene)
        writeFasta(gene_out, nt_sequence, args.compress)

    printv(f"Finished took {main_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nsapphyre MergeGenes"
    )
