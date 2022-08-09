import argparse
import os
from multiprocessing.pool import ThreadPool
from threading import Lock


def run_command(arg_tuple: tuple) -> None:
    string, gene, lock = arg_tuple
    with lock:
        print(gene)
    command = string.format(gene)
    os.system(command)


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--processes', type=int, default=0,
                    help='Number of threads used to call processes.')
args = parser.parse_args()
try:
    os.mkdir('nt_aligned')
except:
    pass

genes = [gene.split('.')[0] for gene in os.listdir('aa') if '.aa' in gene]

command = 'perl pal2nal.mod.pl -output fasta aa/{0}.aa.fa nt/{0}.nt.fa > nt_aligned/{0}.nt.fa'

if args.processes:
    arguments = list()
    lock = Lock()
    for gene in genes:
        arguments.append((command, gene, lock))
    with ThreadPool(args.processes) as pool:
        pool.map(run_command, arguments, chunksize=1)
else:
    for gene in genes:
        print(gene)
        os.system(command.format(gene))