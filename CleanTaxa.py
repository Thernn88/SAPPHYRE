import argparse
import os
from sys import argv
from time import time


def sequence_is_reference(header):
    """
    Returns True if last character of header is .
    """
    return header[-1] == '.'


def clean_taxa(path_to_taxa: str) -> None:
    clean_path = os.path.join(path_to_taxa,'clean')
    aa_clean_path = os.path.join(clean_path,'aa')
    nt_clean_path = os.path.join(clean_path,'nt')

    if not os.path.exists(clean_path): os.mkdir(clean_path)
    if not os.path.exists(aa_clean_path): os.mkdir(aa_clean_path)
    if not os.path.exists(nt_clean_path): os.mkdir(nt_clean_path)

    print('Doing: AA')

    aa_path = os.path.join(path_to_taxa,'aa')

    # These dictionaries ensure any changes made to headers in AA genes are exactly mirrored to NT genes
    rev_f_tags = {}
    anti_dupe_tags = {}
    anti_dupe_counts = {}

    for aa_file in os.listdir(aa_path):
        if '.fa' in aa_file:
            gene = aa_file.split('.')[0]
            if gene not in rev_f_tags:
                rev_f_tags[gene] = {}
                anti_dupe_tags[gene] = {}
                
            aa_file_path = os.path.join(aa_path,aa_file)
            aa_new_file_path = os.path.join(aa_clean_path,aa_file)

            #################################################################################################
            # file.read() is fast, but it also reads the entire file into memory at once.                   #
            # If you have a 20gb file, it will read all 20gb into memory.                                   #
            # That will kill my desktop with a single thread, and severely limit MayMay's multi-threading.  #
            # This change makes python read one line at a time, so only one line is in memory.              #
            # The split might technically make it double the file size for that line, I'm not certain.      #
            #################################################################################################

            # content = open(aa_file_path).read()
            # lines = content.split('\n')
            #
            # while '' in lines: lines.remove('') # Remove any empty lines
            lines = list()
            with open(aa_file_path) as f:
                for line in f:
                    if line != '\n':
                        lines.append(line.strip())

            new_lines = []
            header = None

            # headers = [line.split('|')[2] for line in lines if '>' in line]
            headers = [line.split('|')[2] for line in lines if line[0] == '>']
            output = open(aa_new_file_path, 'w')  # open this once instead of every loop iteration
            for i in range(0, len(lines), 2):
                header = lines[i]
                sequence = lines[i+1]
                
                header_gene, taxaId, node, coords, frame, taxa = header.split('|')
                
                if sequence_is_reference(header):
                    # new_header = header_gene+'|'+taxaId+'|'+node # Format reference header
                    new_header = [header_gene, '|', taxaId, '|', node]  # its faster to change a list

                else:
                    if not taxaId[-1].isnumeric():
                        taxaId = taxaId[:-1]
                        
                    # new_header = header_gene+'|'+taxa+'|'+taxaId+'|'+node # Format candidate header
                    new_header = [header_gene, '|', taxa, '|', taxaId, '|', node]

                    if headers.count(node) > 1:
                        if node not in rev_f_tags[gene]: rev_f_tags[gene][node] = {}
                        
                        if 'revcomp' in frame:
                            # Ensure change made in this gene, for this node, for this exact sequence, is made in the NT file
                            rev_f_tags[gene][node][i] = '_R'

                            #####################################################################################
                            # Strings are actually immutable. Once one is created, you can't change it.         #
                            # The string += operation is convenient, but it is a lie. You can't add to a string.#
                            # Python makes a new one every time and hides it from you to be nice.               #
                            # Every time you use +=, the run time is len(string) + len(addition), not 1.        #
                            #####################################################################################

                            # new_header += '_R'  # this is actually O(N)
                            new_header.append('_R')  # this is basically O(1), we don't care about amortization right now
                        else:
                            # Ensure change made in this gene, for this node, for this exact sequence, is made in the NT file
                            rev_f_tags[gene][node][i] = '_F'
                            
                            # new_header += '_F'
                            new_header.append('_F')

                        if not node in anti_dupe_tags[gene]: anti_dupe_tags[gene][node] = {}

                        if node in anti_dupe_counts:
                            anti_dupe_counts[node] += 1
                            count = anti_dupe_counts[node]
                        else:
                            anti_dupe_counts[node] = 1
                            count = 1

                        # Ensure change made in this gene, for this node, for this exact sequence, is made in the NT file
                        header_addition = '_{}'.format(count)
                        anti_dupe_tags[gene][node][i] = header_addition

                        # new_header += header_addition
                        new_header.append(header_addition)

                # new_lines.append(new_header)
                # new_lines.append(sequence)
                output.write(new_header)
                output.write('\n')
                output.write(sequence)
                output.write('\n')
            # open(aa_new_file_path,'w').write('\n'.join(new_lines))
            output.close()

    print('Doing: NT')
    nt_path = os.path.join(path_to_taxa,'nt')
    for nt_file in os.listdir(nt_path):
        if '.fa' in nt_file:
            gene = nt_file.split('.')[0]

            nt_file_path = os.path.join(nt_path,nt_file)
            nt_new_file_path = os.path.join(nt_clean_path,nt_file)

            # content = open(nt_file_path).read()
            # lines = content.split('\n')
            lines = list()
            with open(nt_file_path) as f:
                for line in f:
                    if line != '\n':
                        lines.append(line.strip())

            # while '' in lines: lines.remove('') # Remove any empty lines

            new_lines = []
            
            headers = [line.split('|')[2] for line in lines if '>' in line]
            output = open(nt_new_file_path,'w')
            for i in range(0,len(lines),2):
                header = lines[i]
                sequence = lines[i+1]
                
                header_gene,taxaId,node,coords,frame,taxa = header.split('|')
                
                if sequence_is_reference(header):
                    # new_header = header_gene+'|'+taxaId+'|'+node # Format reference header
                    new_header = [header_gene, '|', taxaId, '|', node] # Format reference header

                else:
                    if not taxaId[-1].isnumeric():
                        taxaId = taxaId[:-1]
                        
                    new_header = header_gene+'|'+taxa+'|'+taxaId+'|'+node # Format candidate header
                    new_header = [header_gene, '|', taxa, '|', taxaId, '|', node] # Format candidate header

                    if headers.count(node) > 1:
                        new_header += rev_f_tags[gene][node][i] + anti_dupe_tags[gene][node][i]
                output.write(new_header)
                output.write('\n')
                output.write(sequence)
                output.write('\n')
                # new_lines.append(new_header)
                # new_lines.append(sequence)
            output.close()
            # open(nt_new_file_path,'w').write('\n'.join(new_lines))


def main(argv):
    start_time = time()
    parser = argparse.ArgumentParser("Remove duplicates across multiple taxa.")
    parser.add_argument('-i','--input',type=str, default='parent',
        help='Path to parent folder of taxa')
    args = parser.parse_args()

    this_input = args.input
    parent_folder = os.listdir(this_input)

    for taxa in parent_folder:
        print('Doing Taxa:', taxa)

        taxa_path = os.path.join(this_input, taxa)

        clean_taxa(taxa_path)
    
    print('Done took {:.2f}s'.format( time() - start_time ))


if __name__ == '__main__':
    main(argv)