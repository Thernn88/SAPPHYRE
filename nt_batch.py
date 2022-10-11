import argparse
import os
from multiprocessing.pool import ThreadPool
from threading import Lock


def run_command(arg_tuple: tuple) -> None:
    string, gene, taxa, parent, aa_path, nt_path, lock = arg_tuple

    align_references(aa_path, nt_path)
    with lock:
        print(gene)
    command = string.format(gene,taxa,parent)
    os.system(command)

def align_references(aa_path: str, nt_path: str) -> None:
    order = []
    with open(aa_path) as aa_file:
        for line in aa_file:
            if line != '':
                if '>' in line:
                    fields = line.split('|')
                    if len(fields) > 3:
                        break

                    else:
                        order.append(fields[1]) # Save order of reference taxa name
    
    nt_references = {}
    nt_out_lines = []

    with open(nt_path) as nt_file:
        content = nt_file.read()
        lines = content.split('\n')

        for i in range(0,len(lines),2):
            if lines[i] != '':
                header = lines[i]
                seq = lines[i+1]

                fields = header.split('|')
                if len(fields) > 3:
                    nt_out_lines.append(header)
                    nt_out_lines.append(seq)
                
                else:
                    nt_references[fields[1]] = (header,seq)
    
    to_add = []

    for taxa_name in order:
        header, seq = nt_references[taxa_name]
        to_add.append(header)
        to_add.append(seq)

    nt_out_lines = to_add + nt_out_lines

    open(nt_path, 'w').write('\n'.join(nt_out_lines))

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, default='Parent',
                    help='Parent input path.')
parser.add_argument('-p', '--processes', type=int, default=0,
                    help='Number of threads used to call processes.')
args = parser.parse_args()

for taxa in os.listdir(args.input):
    print('Doing taxa {}'.format(taxa))
    taxa_path = os.path.join(args.input,taxa)
    nt_aligned_path = os.path.join(taxa_path,'nt_aligned')
    mafft_path = os.path.join(taxa_path,'mafft')

    if os.path.exists(mafft_path):
        if not os.path.exists(nt_aligned_path): os.mkdir(nt_aligned_path)

        genes = [gene.split('.')[0] for gene in os.listdir(mafft_path) if '.aa' in gene]

        command = 'perl pal2nal.mod.pl -output fasta {2}/{1}/mafft/{0}.aa.fa {2}/{1}/nt/{0}.nt.fa > {2}/{1}/nt_aligned/{0}.nt.fa'

        if args.processes:
            arguments = list()
            lock = Lock()
            for gene in genes:

                aa_path = os.path.join(args.input,taxa,'mafft',gene+'.aa.fa')
                nt_path = os.path.join(args.input,taxa,'nt',gene+'.nt.fa')

                arguments.append((command, gene, taxa, args.input, aa_path, nt_path, lock))
            with ThreadPool(args.processes) as pool:
                pool.map(run_command, arguments, chunksize=1)
        else:
            for gene in genes:
                print(gene)

                aa_path = os.path.join(args.input,taxa,'mafft',gene+'.aa.fa')
                nt_path = os.path.join(args.input,taxa,'nt',gene+'.nt.fa')

                align_references(aa_path, nt_path)
                os.system(command.format(gene,taxa,args.input))
    else:
        print("Can't find mafft folder for taxa {}".format(taxa))