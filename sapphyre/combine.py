from __future__ import annotations

import os
from collections import defaultdict
from math import ceil
from multiprocessing.pool import Pool
from pathlib import Path

from .timekeeper import KeeperMode, TimeKeeper
from .utils import parseFasta, printv, writeFasta


def delete_all_files(path: str) -> None:
    fastas = [os.path.join(path, x) for x in os.listdir(path) if ".fa" in x]
    for fasta in fastas:
        os.remove(fasta)


def prepend(inputs, directory):
    return [Path(directory, i) for i in inputs]


def parse_gene(path):
    this_out = []
    this_refs = {}
    already_grabbed_references = set()
    present_taxon = set()
    taxon_to_target = defaultdict(list)
    for header, sequence in parseFasta(str(path)):
        if header.endswith("."):
            fields = header.split("|")
            if fields[2] not in already_grabbed_references:
                taxon_to_target[fields[1]].append(fields[2])
                already_grabbed_references.add(fields[2])
                this_refs[fields[2]] = header, sequence.replace("-", "")
        else:
            present_taxon.add(header.split("|")[1])
            this_out.append((header, sequence.replace("-", "")))
    return {path: (this_out, present_taxon, taxon_to_target, this_refs)}


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
    grabbed_aa_references = defaultdict(set)
    grabbed_nt_references = defaultdict(set)

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

                aa_result = pool.map(parse_gene, aa_path.iterdir(), chunksize=8)

                for result in aa_result:
                    for path, out_tuple in result.items():
                        (
                            out_candidates,
                            present_taxons,
                            taxon_to_target,
                            ref_sequences,
                        ) = out_tuple
                        path = path.name
                        already_grabbed = grabbed_aa_references[path]
                        references = []
                        for taxon in present_taxons:
                            for target in taxon_to_target[taxon]:
                                if target not in already_grabbed:
                                    references.append(ref_sequences[target])
                                    already_grabbed.add(target)

                        aa_out.setdefault(path, []).extend(references)
                        aa_out.setdefault(path, []).extend(out_candidates)
                        grabbed_aa_references[path].update(already_grabbed)

                arguments = []
                for nt_gene in nt_path.iterdir():
                    arguments.append(
                        (nt_gene,),
                    )

                nt_result = pool.starmap(parse_gene, arguments, chunksize=8)

                for result in nt_result:
                    for path, out_tuple in result.items():
                        (
                            out_candidates,
                            present_taxons,
                            taxon_to_target,
                            ref_sequences,
                        ) = out_tuple
                        path = path.name
                        already_grabbed = grabbed_nt_references[path]
                        references = []
                        for taxon in present_taxons:
                            for target in taxon_to_target[taxon]:
                                if target not in already_grabbed:
                                    references.append(ref_sequences[target])
                                    already_grabbed.add(target)

                        nt_out.setdefault(path, []).extend(references)
                        nt_out.setdefault(path, []).extend(out_candidates)
                        grabbed_nt_references[path].update(already_grabbed)

    aa_out_path = Path(args.output_directory, "aa")
    nt_out_path = Path(args.output_directory, "nt")
    aa_out_path.mkdir(parents=True, exist_ok=True)
    nt_out_path.mkdir(parents=True, exist_ok=True)
    # clear output directories to make reconcile safe
    delete_all_files(aa_out_path)
    delete_all_files(nt_out_path)

    aa_sequences = list(aa_out.items())
    nt_sequences = list(nt_out.items())

    if aa_sequences:
        per_thread = ceil(len(aa_sequences) / args.processes)

        aa_arguments = [
            (aa_out_path, aa_sequences[i : i + per_thread], args.compress)
            for i in range(0, len(aa_sequences), per_thread)
        ]

        with Pool(args.processes) as pool:
            pool.starmap(write_gene, aa_arguments)
        del aa_arguments

        nt_arguments = [
            (nt_out_path, nt_sequences[i : i + per_thread], args.compress)
            for i in range(0, len(nt_sequences), per_thread)
        ]

        with Pool(args.processes) as pool:
            pool.starmap(write_gene, nt_arguments)
        del nt_arguments

    printv(f"Finished took {main_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True
