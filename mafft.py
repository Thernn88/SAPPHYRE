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
    string, taxa, mafft_path, parent, gene, lock = arg_tuple
    with lock:
        print(gene)
    command = string.format(gene, taxa, parent)
    os.system(command)

    file = gene+'.aa.fa'
    path = os.path.join(mafft_path,file)

    content = open(path).read()

    out_content = delete_no_pipe_headers(content)

    open(path,'w').write(out_content)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, default='Parent',
                    help='Path to parent folder for input')
parser.add_argument('-p', '--processes', type=int, default=0,
                    help='Number of threads used to call processes.')
parser.add_argument(
    "-oi",
    "--orthoset_input",
    type=str,
    default="orthosets",
    help="Path to directory of Orthosets folder",
)
parser.add_argument(
    "-o",
    "--orthoset",
    type=str,
    required=True,
    help="Orthoset",
)
args = parser.parse_args()

mafft_folder = 'mafft'
aa_folder = 'aa'
aln_folder = 'aln'

aln_path = os.path.join(args.orthoset_input, args.orthoset, aln_folder)

for taxa in os.listdir(args.input):
    print('Doing taxa {}'.format(taxa))
    mafft_path = os.path.join(args.input, taxa, mafft_folder)
    aa_path = os.path.join(args.input, taxa, aa_folder)
    if os.path.exists(aa_path):
        if not os.path.exists(mafft_path): os.mkdir(mafft_path)

        genes = [gene.split('.')[0] for gene in os.listdir(aa_path) if '.aa' in gene]

        #command = 'mafft --anysymbol --auto --thread -1 --addfragments {2}/{1}/aa/{0}.aa.fa --thread -1 '+aln_path+'/{0}.aln.fa > {2}/{1}/mafft/{0}.aa.fa'
        command = 'mafft-linsi --anysymbol --addfragments {2}/{1}/aa/{0}.aa.fa --thread -1 '+aln_path+'/{0}.aln.fa > {2}/{1}/mafft/{0}.aa.fa'


        if args.processes:
            arguments = list()
            lock = Lock()
            for gene in genes:
                arguments.append((command, taxa, mafft_path, args.input, gene, lock))
            with ThreadPool(args.processes) as pool:
                pool.map(run_command, arguments, chunksize=1)
        else:
            for gene in genes:
                print(gene)
                os.system(command.format(gene, taxa, args.input))

                file = gene+'.aa.fa'
                path = os.path.join(mafft_path,file)

                content = open(path).read()

                out_content = delete_no_pipe_headers(content)

                open(path,'w').write(out_content)

        print('Took {:.2f}s'.format(time()-start))
    else:
        print("Can't find aa folder for taxa {}".format(taxa))
