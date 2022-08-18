import argparse
import os
from multiprocessing.pool import ThreadPool
from threading import Lock


def run_command(arg_tuple: tuple) -> None:
    string, gene, taxa, parent, lock = arg_tuple
    with lock:
        print(gene)
    command = string.format(gene,taxa,parent)
    os.system(command)


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
                arguments.append((command, gene, taxa, args.input, lock))
            with ThreadPool(args.processes) as pool:
                pool.map(run_command, arguments, chunksize=1)
        else:
            for gene in genes:
                print(gene)
                os.system(command.format(gene,taxa,args.input))
    else:
        print("Can't find mafft folder for taxa {}".format(taxa))