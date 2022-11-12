import csv 
import os
import requests
import argparse
from time import sleep
from multiprocessing.pool import Pool

def download_parallel(command, srr_acession, path_to_download, expected_directory, out_fields):
    os.system(f"{command} {srr_acession} -O {path_to_download}")
    return os.path.exists(expected_directory), out_fields

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default='Coleoptera.csv', help="CSV File Input.")
    parser.add_argument('-b', '--bin', default='C:/Users/kevin/Downloads/sratoolkit/bin', help="Path to SRA Toolkit.")
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Number of threads used to call processes.",
    )
    args = parser.parse_args()

    path_to_bin = args.bin
    base = os.path.basename(args.input).split('.')[0]

    finPath = args.input.split('.')
    finPath[0] = finPath[0]+'_Downloaded'
    finPath = '.'.join(finPath)
    failPath = args.input.split('.')
    failPath[0] = failPath[0]+'_Fails'
    failPath = '.'.join(failPath)

    open(finPath,'w',encoding='utf-8').write('')
    open(failPath,'w',encoding='utf-8').write('')

    with open(args.input, encoding='utf-8') as csvfile:
        csv_read = csv.reader(csvfile, delimiter=',', quotechar='"')
        arguments = []
        for i,fields in enumerate(csv_read):
            out_fields = ['"'+str(i).replace('ï»¿','')+'"' for i in fields]
            if fields[0] == 'Experiment Accession':
                out_fields.append('"SRR Acession"')

            elif fields[0] != '':
                acession = fields[0]
                print(f"Searching for runs in SRA: {acession}")

                url = 'https://www.ncbi.nlm.nih.gov/sra/{}[accn]'.format(acession)

                req = requests.get(url)
                
                link_elements = req.text.split('trace.ncbi.nlm.nih.gov/Traces?run=')
                
                for link in link_elements[1:]:
                    srr_acession = link.split('">')[0]
                    print(f"Attempting to download: {srr_acession}")
                    out_fields.append('"'+srr_acession+'"')

                    command = os.path.join(path_to_bin,'fastq-dump --gzip')
                    path_to_download = os.path.join(os.getcwd(),base)
                    expected_directory = os.path.join(path_to_download,srr_acession+'.fastq')
                    arguments.append((command,srr_acession,path_to_download,expected_directory,out_fields))

        with Pool(args.processes) as pool:
            downloaded_files = pool.starmap(download_parallel, arguments, chunksize=1)