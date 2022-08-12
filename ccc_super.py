import argparse
import os
import shutil
from sys import argv


def folder_check(folder: str) -> None:
    if not os.path.exists(folder):
        os.mkdir(folder)

def get_filenames(folder: str) -> list:
    result = [item for item in os.listdir(folder) if item[0] not in ['.','$']]
    return result

def sequence_is_reference(header):
    """
    Returns True if last character of header is .
    """
    return header[-1] == '.'

def clone_and_clean_files(aa_fastas: list, out_dir_clean_aa: str, nt_fastas: str, out_dir_clean_nt: str) -> None:
    rev_f_tags = {}
    anti_dupe_tags = {}
    anti_dupe_counts = {}

    aa_content = {}

    for fasta in aa_fastas:
        this_gene = os.path.basename(fasta)
        if this_gene not in aa_content:
            aa_content[this_gene] = []

        with open(fasta) as f:
            for line in f:
                if line != '\n':
                    aa_content[this_gene].append(line.strip())

    for fasta in aa_content:
        aa_new_file_path = os.path.join(out_dir_clean_aa,fasta)
        gene = fasta.split('.')[0]
        if gene not in rev_f_tags:
            rev_f_tags[gene] = {}
            anti_dupe_tags[gene] = {}

        lines = aa_content[fasta]
        
        headers = [line.split('|')[2] for line in lines if line[0] == '>']
        output = open(aa_new_file_path, 'w')  # open this once instead of every loop iteration
        already_written = set()
        for i in range(0, len(lines), 2):
            header = lines[i]
            sequence = lines[i+1]
            
            header_gene, taxaId, node, coords, frame, taxa = header.split('|')
            
            if sequence_is_reference(header):
                new_header = ''.join([header_gene, '|', taxaId, '|', node]) # its faster to change a list
                if new_header in already_written:
                    continue
                already_written.add(new_header)
            else:
                if not taxaId[-1].isnumeric():
                    taxaId = taxaId[:-1]

                new_header = [header_gene, '|', taxa, '|', taxaId, '|', node]

                if headers.count(node) > 1:
                    if node not in rev_f_tags[gene]: rev_f_tags[gene][node] = {}

                    if 'revcomp' in frame:
                        rev_f_tags[gene][node][i] = '_R'
                        new_header.append('_R')
                    else:
                        rev_f_tags[gene][node][i] = '_F'
                        
                        new_header.append('_F')

                    if not node in anti_dupe_tags[gene]: anti_dupe_tags[gene][node] = {}

                    if node in anti_dupe_counts:
                        anti_dupe_counts[node] += 1
                        count = anti_dupe_counts[node]
                    else:
                        anti_dupe_counts[node] = 1
                        count = 1

                    header_addition = '_{}'.format(count)
                    anti_dupe_tags[gene][node][i] = header_addition

                    new_header.append(header_addition)
                new_header = ''.join(new_header)

            output.write(new_header)

            output.write('\n')
            output.write(sequence)
            output.write('\n')

        output.close()

    nt_content = {}

    for fasta in nt_fastas:
        this_gene = os.path.basename(fasta)
        if this_gene not in nt_content:
            nt_content[this_gene] = []

        with open(fasta) as f:
            for line in f:
                if line != '\n':
                    nt_content[this_gene].append(line.strip())

    for fasta in nt_content:
        gene = fasta.split('.')[0]

        nt_new_file_path = os.path.join(out_dir_clean_nt,fasta)

        lines = nt_content[fasta]

        headers = [line.split('|')[2] for line in lines if '>' in line]
        output = open(nt_new_file_path,'w')
        already_written = set()
        for i in range(0,len(lines),2):
            header = lines[i]
            sequence = lines[i+1]
            
            header_gene,taxaId,node,coords,frame,taxa = header.split('|')

            if sequence_is_reference(header):
                new_header = ''.join([header_gene, '|', taxaId, '|', node]) # Format reference header
                if new_header in already_written:
                    continue
                already_written.add(new_header)
            else:
                if not taxaId[-1].isnumeric():
                    taxaId = taxaId[:-1]
                    
                new_header = [header_gene, '|', taxa, '|', taxaId, '|', node] # Format candidate header

                if headers.count(node) > 1:
                    new_header.append(rev_f_tags[gene][node][i])
                    new_header.append(anti_dupe_tags[gene][node][i])
                new_header = ''.join(new_header)
            output.write(new_header)
            output.write('\n')
            output.write(sequence)
            output.write('\n')

        output.close()

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, help="Super folder")
    parser.add_argument('-o','--output', required=True,
                        help="Location of cloned directory")
    args = parser.parse_args()

    folder_check(args.output)

    taxas = [taxa for taxa in os.listdir(args.input) if os.path.isdir(os.path.join(args.input, taxa))]
    taxa_groups = dict()

    for taxa in taxas:
        taxa_name = taxa.split('.')[0]
        if taxa_name[-1].isnumeric() and taxa_name[-2] == '_':
            taxa_name = '_'.join(taxa_name.split('_')[:-1])

        if taxa_name not in taxa_groups:
            taxa_groups[taxa_name] = []
        
        taxa_groups[taxa_name].append(taxa)

    for taxa in taxa_groups:
        to_combine_taxa = taxa_groups[taxa]

        aa_fastas = []
        nt_fasta = []

        for combine_taxa in to_combine_taxa:
            aa_in = os.path.join(os.path.join(args.input, combine_taxa), 'aa')
            aa_fastas.extend([os.path.join(aa_in, aa) for aa in get_filenames(aa_in)])

            nt_in = os.path.join(os.path.join(args.input, combine_taxa), 'nt')
            nt_fasta.extend([os.path.join(nt_in, nt) for nt in get_filenames(nt_in)])
        
        clone = os.path.join(args.output, taxa)
        folder_check(clone)

        clean_directory = os.path.join(clone, 'clean')

        folder_check(clean_directory)
        
        aa_out = os.path.join(clone, 'aa')
        nt_out = os.path.join(clone, 'nt')
        aa_clean_out = os.path.join(clean_directory, 'aa')
        nt_clean_out = os.path.join(clean_directory, 'nt')

        folder_check(aa_clean_out)
        folder_check(nt_clean_out)

        clone_and_clean_files(aa_fastas, aa_clean_out, nt_fasta, nt_clean_out)
        
if __name__ == '__main__':
    main(argv)
