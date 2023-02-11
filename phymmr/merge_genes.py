from __future__ import annotations
import json
import wrap_rocks
from pathlib import Path
from multiprocessing.pool import Pool

from .utils import printv, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode


def prepend(inputs, directory):
    return [Path(directory, i) for i in inputs]


def parse_gene(path):
    this_out = []
    for header, sequence in parseFasta(str(path)):
        if header[-1] != ".":
            this_out.append((header, sequence))
    return {path: this_out}


def get_references(gene, orthoset_db):
    core_seqs = json.loads(orthoset_db.get(f"getcore:{gene}"))

    aa_refs = []
    for taxon, target, seq in sorted(core_seqs["aa"]):
        aa_refs.append((f"{gene}|{taxon}|{target}|.",seq))

    nt_refs = []
    for taxon, target, seq in sorted(core_seqs["nt"]):
        nt_refs.append((f"{gene}|{taxon}|{target}|.",seq))

    return aa_refs, nt_refs

def main(args):
    main_keeper = TimeKeeper(KeeperMode.DIRECT)
    if args.DIRECTORY:
        inputs = prepend(args.INPUT, args.DIRECTORY)
    else:
        inputs = [Path(i) for i in args.INPUT]

    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    orthoset_db_path = Path(orthosets_dir, orthoset, "rocksdb")
    if not orthoset_db_path.exists():
        printv(
            f"WARNING: {orthoset_db_path} doesn't exists. Abort",
            args.verbose,
            0,
        )
        return False
    orthoset_db = wrap_rocks.RocksDB(str(orthoset_db_path))

    aa_out = {}
    nt_out = {}

    nt_ref_pass_through = {}

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
                    aa_refs, nt_refs = get_references(aa_gene.name.split('.')[0], orthoset_db)
                    aa_out[aa_gene.name] = aa_refs
                    nt_ref_pass_through[aa_gene.name.split('.')[0]] = nt_refs
                
                arguments.append(
                    (
                        aa_gene,
                    )
                )

            with Pool(args.processes) as pool:
                aa_result = pool.starmap(parse_gene, arguments, chunksize=1)

            for result in aa_result:
                for path, out in result.items():
                    path = path.name
                    aa_out.setdefault(path, []).extend(out)

            arguments = []
            for nt_gene in nt_path.iterdir():
                if nt_gene.name not in nt_out:
                    this_gene = nt_gene.name.split('.')[0]
                    if this_gene not in nt_ref_pass_through:
                        _, nt_refs = get_references(this_gene, orthoset_db)
                    else:
                        nt_refs = nt_ref_pass_through[this_gene]
                    nt_out[nt_gene.name] = nt_refs

                arguments.append(
                    (
                        nt_gene,
                    )
                )

            with Pool(args.processes) as pool:
                nt_result = pool.starmap(parse_gene, arguments, chunksize=1)

            for result in nt_result:
                for path, out in result.items():
                    path = path.name
                    nt_out.setdefault(path, []).extend(out)

    aa_out_path = Path(args.output_directory, "aa")
    nt_out_path = Path(args.output_directory, "nt")
    aa_out_path.mkdir(parents=True, exist_ok=True)
    nt_out_path.mkdir(parents=True, exist_ok=True)

    for gene in aa_out:
        gene_out = Path(aa_out_path, gene)
        writeFasta(gene_out, aa_out[gene], args.compress)

    for gene in nt_out:
        gene_out = Path(nt_out_path, gene)
        writeFasta(gene_out, nt_out[gene], args.compress)
    printv(f"Finished took {main_keeper.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr MergeGenes"
    )