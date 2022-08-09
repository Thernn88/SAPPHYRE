import argparse
import os
from multiprocessing.pool import ThreadPool
from threading import Lock
from time import time

start = time()

def deinterleave(fasta_lines: list) -> str:
    result = []
    this_out = list()
    for line in fasta_lines:
        if line != '':
            if line[0] == '>':
                if this_out:
                    result.append(''.join(this_out))
                result.append(line)
                this_out = list()
            else:
                this_out.append(line.strip())
    if this_out != list():
        result.append(''.join(this_out))
    return result

def delete_no_pipe_headers(content: str) -> str:
    lines = content.split('\n')
    lines = deinterleave(lines)

    result = []

    for i in range(0,len(lines),2):
        header = lines[i]
        seq = lines[i+1]

        if '|' in header:
            result.append(header)
            result.append(seq)

    result = '\n'.join(result)

    return result

def run_command(arg_tuple: tuple) -> None:
    string, taxa, mafft_path, gene, lock = arg_tuple
    with lock:
        print(gene)
    command = string.format(gene, taxa)
    os.system(command)

    file = gene+'.aa.fa'
    path = os.path.join(mafft_path,file)

    content = open(path).read()

    out_content = delete_no_pipe_headers(content)

    open(path,'w').write(out_content)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, default='Taxa',
                    help='Path to taxa for input')
parser.add_argument('-p', '--processes', type=int, default=0,
                    help='Number of threads used to call processes.')
args = parser.parse_args()

mafft_path = os.path.join(args.input,'clean','mafft')
aa_path = os.path.join(args.input,'clean','aa')
if not os.path.exists(mafft_path): os.mkdir(mafft_path)

genes = [gene.split('.')[0] for gene in os.listdir(aa_path) if '.aa' in gene]

#command = 'mafft --anysymbol --auto --thread -1 --addfragments aa/{0}.aa.fa aln/{0}.aln.fa > mafft/{0}.aa.fa'
command = 'mafft-linsi --anysymbol --addfragments {1}/clean/aa/{0}.aa.fa --thread -1 aln/{0}.aln.fa > {1}/clean/mafft/{0}.aa.fa'


if args.processes:
    arguments = list()
    lock = Lock()
    for gene in genes:
        arguments.append((command, args.input, mafft_path, gene, lock))
    with ThreadPool(args.processes) as pool:
        pool.map(run_command, arguments, chunksize=1)
else:
    for gene in genes:
        print(gene)
        os.system(command.format(gene, args.input))

print('Took {:.2f}s'.format(time()-start))

