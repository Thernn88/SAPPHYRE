import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', nargs='+', action='append',
    help='Path to input')
parser.add_argument('-o', '--output', type=str, default='MergedGenes',
    help='Merged output.')
args = parser.parse_args()

aa_out = {}
nt_out = {}

unpacked = []
for outer_directory in args.input:#args.input:
    for inner_directory in outer_directory:
        unpacked.append(inner_directory)

for item in unpacked:
    taxas = os.listdir(item)

    print('Doing Dir:',item)
    for taxa in taxas:
        aa_path = os.path.join(item,taxa,'aa_merged')
        nt_path = os.path.join(item,taxa,'nt_merged')
        if os.path.exists(aa_path) and os.path.exists(nt_path):
            for aa_gene in os.listdir(aa_path):
                this_gene_path = os.path.join(aa_path, aa_gene)

                with open(this_gene_path) as f:
                    if aa_gene not in aa_out:
                        aa_out[aa_gene] = []
                        for line in f:
                            if line != '\n':
                                aa_out[aa_gene].append(line.strip())
                    else:
                        header = None
                        for line in f:
                            if header:
                                aa_out[aa_gene].append(header.strip())
                                aa_out[aa_gene].append(line.strip())
                                header = None
                            else:
                                if '>' in line and line.count('|') != 2:
                                    header = line
            
            for nt_gene in os.listdir(nt_path):
                this_gene_path = os.path.join(nt_path, nt_gene)

                with open(this_gene_path) as f:
                    if nt_gene not in nt_out:
                        nt_out[nt_gene] = []
                        for line in f:
                            if line != '\n':
                                nt_out[nt_gene].append(line.strip())
                    else:
                        header = None
                        for line in f:
                            if header:
                                nt_out[nt_gene].append(header.strip())
                                nt_out[nt_gene].append(line.strip())
                                header = None
                            else:
                                if '>' in line and line.count('|') != 2:
                                    header = line

aa_out_path = os.path.join(args.output,'aa')
nt_out_path = os.path.join(args.output,'nt')

if not os.path.exists(args.output): os.mkdir(args.output)
if not os.path.exists(aa_out_path): os.mkdir(aa_out_path)
if not os.path.exists(nt_out_path): os.mkdir(nt_out_path)

for gene in aa_out:
    gene_out = os.path.join(aa_out_path,gene)
    open(gene_out,'w').write('\n'.join(aa_out[gene]))

for gene in nt_out:
    gene_out = os.path.join(nt_out_path,gene)
    open(gene_out,'w').write('\n'.join(nt_out[gene]))