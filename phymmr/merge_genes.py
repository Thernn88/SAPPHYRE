from __future__ import annotations

from pathlib import Path


def prepend(inputs, directory):
    return [Path(directory, i) for i in inputs]


def main(args):
    if args.DIRECTORY:
        inputs = prepend(args.INPUT, args.DIRECTORY)
    else:
        inputs = args.INPUT

    aa_out = {}
    nt_out = {}

    for item in inputs:
        print("Doing Dir:", item)
        for taxa in item.iterdir():
            aa_path = Path(taxa, "aa_merged")
            nt_path = Path(taxa, "nt_merged")
            if not aa_path.exists() or not nt_path.exists():
                if args.verbose >= 1:
                    print(
                        f"Either {aa_path} or {nt_path} doesn't exists. Ignoring"
                    )
                continue
            for aa_gene in aa_path.iterdir():
                with aa_gene.open() as fp:
                    if aa_gene not in aa_out:
                        aa_out[aa_gene] = []
                        for line in fp:
                            if line != "\n":
                                aa_out[aa_gene].append(line.strip())
                    else:
                        header = None
                        for line in fp:
                            if header:
                                aa_out[aa_gene].append(header.strip())
                                aa_out[aa_gene].append(line.strip())
                                header = None
                            else:
                                if ">" in line and line[-1] != '.':
                                    header = line

                for nt_gene in nt_path.iterdir():
                    with nt_gene.open() as fp:
                        if nt_gene not in nt_out:
                            nt_out[nt_gene] = []
                            for line in fp:
                                if line != "\n":
                                    nt_out[nt_gene].append(line.strip())
                        else:
                            header = None
                            for line in fp:
                                if header:
                                    nt_out[nt_gene].append(header.strip())
                                    nt_out[nt_gene].append(line.strip())
                                    header = None
                                else:
                                    if ">" in line and line[-1] != '.':
                                        header = line

    aa_out_path = Path(args.OUPUT_DIRECTORY, "aa")
    nt_out_path = Path(args.OUPUT_DIRECTORY, "nt")
    aa_out_path.mkdir(parents=True, exist_ok=True)
    nt_out_path.mkdir(parents=True, exist_ok=True)

    for gene in aa_out:
        gene_out = Path(aa_out_path, gene.name)
        gene_out.write_text("\n".join(aa_out[gene]))

    for gene in nt_out:
        gene_out = Path(nt_out_path, gene.name)
        gene_out.write_text("\n".join(nt_out[gene]))

    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr MergeGenes"
    )
