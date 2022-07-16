import argparse
import os
parser = argparse.ArgumentParser("Containment check post alignment.")
parser.add_argument('-b','--base',type=str, default='in1',
    help='Path to base')
parser.add_argument('-i','--input', nargs='+', action='append',
    help='Path to input')
parser.add_argument('-o','--output',type=str, default='SuperSet',
    help='Path to output')
args = parser.parse_args()

merge_aa = {}
merge_nt = {}

done_aa = {}
done_nt = {}

unpacked = [args.base]
for outer in args.input:
    for inner in outer:
        unpacked.append(inner)


for item in unpacked:
    if item == args.base:
        taxas = [''] #Skip sub-dir check
    else:
        taxas = os.listdir(item)

    print('Doing Dir:',item)
    for taxa in taxas:
        if taxa == '':
            folder_check = True
            base_path = item+'/'
        else:
            folder_check = os.path.isdir(item+'/'+taxa)
            base_path = item+'/'+taxa+'/'
        if folder_check:
            print('Doing Taxa:',taxa)
            print('Doing: AA')

            for aa_file in os.listdir(base_path+'aa'):
                if '.fa' in aa_file:
                    gene = aa_file.split('.')[0]
                    lines = open(base_path+'aa/'+aa_file).read().split('\n')

                    if gene not in merge_aa.keys():
                        merge_aa[gene] = lines
                        for line in lines:
                            if '>' in line:
                                done_aa[line] = False
                    else:
                        header = None
                        out_lines = []
                        for line in lines:
                            if header:
                                try:
                                    done_aa[header] == False
                                except KeyError:
                                    done_aa[header] = False
                                    out_lines.append(header)
                                    out_lines.append(line)
                                header = None
                            else:
                                if line.count('|') == 3:
                                    header = line

                        merge_aa[gene] += out_lines
                        

            print('Doing: NT')
            for nt_file in os.listdir(base_path+'nt'):
                if '.fa' in nt_file:
                    gene = nt_file.split('.')[0]
                    lines = open(base_path+'nt/'+nt_file).read().split('\n')

                    if gene not in merge_nt.keys():
                        merge_nt[gene] = lines
                        for line in lines:
                            if '>' in line:
                                done_nt[line] = False
                    else:
                        header = None
                        out_lines = []
                        for line in lines:
                            if header:
                                try:
                                    done_nt[header] == False
                                except:
                                    done_nt[header] = False
                                    out_lines.append(header)
                                    out_lines.append(line)
                                header = None
                            else:
                                if line.count('|') == 3:
                                    header = line

                        merge_nt[gene] += out_lines

try:
    os.mkdir(args.output)
except:
    pass
try:
    os.mkdir(args.output+'/aa')
except:
    pass
try:
    os.mkdir(args.output+'/nt')
except:
    pass

for gene in merge_aa.keys():
    open(args.output+'/aa/'+gene+'.aa.fas','w').write('\n'.join(merge_aa[gene]))

for gene in merge_nt.keys():
    open(args.output+'/nt/'+gene+'.nt.fas','w').write('\n'.join(merge_nt[gene]))
