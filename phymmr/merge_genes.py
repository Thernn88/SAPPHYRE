from __future__ import annotations
from time import time
from pathlib import Path
from multiprocessing.pool import Pool

from .utils import printv
from .timekeeper import TimeKeeper, KeeperMode

def prepend(inputs, directory):
    return [Path(directory, i) for i in inputs]

def parse_gene(path, grab_references = False):
    this_out = []
    with path.open() as fp:
        if grab_references:
            for line in fp:
                if line != "\n":
                    this_out.append(line.strip())
        else:
            for line in fp:
                if ">" == line[0] and line[-1] != '.':
                    this_out.append(line.strip())
                    sequence = next(fp)
                    this_out.append(sequence.strip())
    return {path: this_out}

def main(args):
    main_keeper = TimeKeeper(KeeperMode.DIRECT)
    if args.DIRECTORY:
        inputs = prepend(args.INPUT, args.DIRECTORY)
    else:
        inputs = [Path(i) for i in args.INPUT]

    aa_out = {}
    nt_out = {}

    for item in inputs:
        printv(f"Merging directory: {item}", args.verbose, 0)
        for taxa in item.iterdir():
            aa_path = Path(taxa, "aa_merged")
            nt_path = Path(taxa, "nt_merged")
            if not aa_path.exists() or not nt_path.exists():
                printv(f"WARNING: Either {aa_path} or {nt_path} doesn't exists. Abort", args.verbose, 0)
                continue

            arguments = []
            for aa_gene in aa_path.iterdir():
                add_references = aa_gene not in aa_out
                arguments.append((aa_gene, add_references,))

            with Pool(args.processes) as pool:
                aa_result = pool.starmap(parse_gene, arguments, chunksize=1)

            for result in aa_result:
                aa_out = aa_out | result

            arguments = []
            for nt_gene in nt_path.iterdir():
                add_references = nt_gene not in nt_out
                arguments.append((nt_gene, add_references,))

            with Pool(args.processes) as pool:
                nt_result = pool.starmap(parse_gene, arguments, chunksize=1)

            for result in nt_result:
                nt_out = nt_out | result
        
    aa_out_path = Path(args.output_directory, "aa")
    nt_out_path = Path(args.output_directory, "nt")
    aa_out_path.mkdir(parents=True, exist_ok=True)
    nt_out_path.mkdir(parents=True, exist_ok=True)

    for gene in aa_out:
        gene_out = Path(aa_out_path, gene.name)
        gene_out.write_text("\n".join(aa_out[gene]))

    for gene in nt_out:
        gene_out = Path(nt_out_path, gene.name)
        gene_out.write_text("\n".join(nt_out[gene]))
    printv(f"Finished took {main_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr MergeGenes"
    )
