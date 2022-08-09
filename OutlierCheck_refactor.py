import argparse
from copy import deepcopy
import os
from statistics import mean
from multiprocessing.pool import Pool
import numpy as np
import blosum_distance as bd
from itertools import combinations
from time import time

class Record():

    def __init__(self, head, seq, raw_seq=None):
        self.id = head
        self.sequence = seq
        if raw_seq == None:
            self.raw = seq
        else:
            self.raw = raw_seq
    
    def __hash__(self):
        return hash(self.id+self.sequence)

    def __str__(self):
        return self.sequence


def taxa_sort(lines: list) -> list:
    """
    Iterates over a list of candidates and creates Records. Sorts by
    taxa, then makes a fasta list. Returns the list.
    """
    records = list()
    for i in range(0, len(lines), 2):
        specimen = Record(lines[i], lines[i+1])
        records.append(specimen)
    records.sort(key=lambda x: (x.id.split('|')[2],x.id.split('|')[3]))
    output = list()
    for record in records:
        output.append(record.id)
        output.append(record.sequence)
    return output


def original_sort(headers, lines) -> list:
    """
    Returns candidate sequences to their original order.
    """
    output = list()
    record_dict = dict()
    for i in range(0, len(lines), 2):
        record_dict[lines[i]] = lines[i+1]
    for header in headers:
        sequence = record_dict.get(header, False)
        if sequence:
            output.append(header)
            output.append(sequence)
    return output


def folder_check(path: str) -> None:
    if not os.path.exists(path):
        os.mkdir(path)
    aa_folder = os.path.join(path, 'aa')
    if not os.path.exists(aa_folder):
        os.mkdir(aa_folder)
    nt_folder = os.path.join(path, 'nt')
    if not os.path.exists(nt_folder):
        os.mkdir(nt_folder)
    logs_folder = os.path.join(path, 'logs')
    if not os.path.exists(logs_folder):
        os.mkdir(logs_folder)


def is_reference_header(header: str) -> bool:
    """
    Counts the | pipe characters in a string. If two are found, returns true.
    Otherwise returns false.
    """
    result = header.count('|') == 2
    return result


def get_headers(lines: list) -> list:
    """
    Returns a list of every other line in the provided argument. Used to get
    header names from a list of sequences.
    """
    result = list()
    for i in range(0, len(lines), 2):
        result.append(lines[i])
    return result


def split_sequences(lines: list, excluded: set) -> tuple:
    """
    Reads over a fasta record in the given list and returns a tuple of two smaller lists.
    The first returned list is the reference sequences found, the second returned list
    is the candidate sequences found.
    """
    bad_names = {'bombyx_mori', 'danaus_plexippus'}
    references = list()
    candidates = list()
    for i in range(0, len(lines), 2):
        header = lines[i]
        sequence = lines[i+1]
        if is_reference_header(header):
            if header.split('|')[1].lower() in bad_names:
                excluded.add(header.strip())
            references.append(header.strip())
            references.append(sequence.strip())
        else:
            candidates.append(header.strip())
            candidates.append(sequence.strip())
    return references, candidates


def make_indices(sequence: str, gap_character='-') -> tuple:
    """
    Finds the index of the first and last non-gap bp in a sequence.
    Returns the start value and the end values + 1 as a tuple.
    """
    start = None
    end = None
    for i,character in enumerate(sequence):
        if character != gap_character:
            start = i
            break
    for i in range(len(sequence)-1, -1, -1):
        if sequence[i] != gap_character:
            end = i+1
            break
    if start == None or end == None:
        raise ValueError()
    return start, end


def sequence_has_data(sequence: str) -> bool:
    """
    Returns True if the string contains a non-gap character.
    Otherwise, returns False.
    """
    result = False
    for character in sequence:
        if character != '-':
            result = True
            break
    return result


def constrain_data_lines(lines: list, start: int, end: int) -> tuple:
    """
    Given a start and end value, iterates over the list of sequences and
    trims the non-header lines to given values. No return, mutates the original data.
    """
    full = list()
    heads = list()
    for i in range(0, len(lines), 2):
        newline = lines[i+1][start:end]
        if sequence_has_data(newline):
            full.append(lines[i])
            full.append(newline)
            heads.append(lines[i])
    return (full, heads)


def convert_to_record_objects(lines: list) -> list:
    """
    Given a list of stings from a fasta file, returns a list of Sequence objects
    from the biopython module. This allows us to make a MultipleSequenceAlignment
    object later.
    """
    result = list()
    for i in range(0, len(lines), 2):
        record_object = Record(lines[i], lines[i+1])
        result.append(record_object)
    return result


def find_index_groups(references: list, candidates: list) -> tuple:
    """
    Iterate over a list of candidate fastas as lines of text and finds their start
    and stop indices. Makes a tuple out of the pairs, then uses the
    tuple as a key in two dictionaries. One dictionary stores lists of
    candidates with identical indices, and the other dictionary stores
    the ref set after constraining to those indices.
    """
    candidate_dict = dict()
    for i in range(0, len(candidates), 2):
        header = candidates[i]
        sequence = candidates[i+1]
        raw_seq = sequence
        index_tuple = make_indices(sequence)
        start, stop = index_tuple
        lines = [candidates[i], candidates[i+1]]
        lines, _ = constrain_data_lines(lines, start, stop)
        cand_seq = Record(lines[0],lines[1],raw_seq)
        made_already = candidate_dict.get(index_tuple, False)
        if not made_already:
            seq_set = set()
            seq_set.add(cand_seq)
            candidate_dict[index_tuple] = seq_set
        else:
            made_already.add(cand_seq)
            candidate_dict[index_tuple] = made_already
    # after processing candidates, make appropriate ref sets
    reference_dict = dict()
    raw_ref_dict = dict()
    for key in candidate_dict:
        start, stop = key
        ref_lines = deepcopy(references)
        raw_ref_dict[key] = ref_lines
        ref_lines, ref_headers = constrain_data_lines(ref_lines, start, stop)
        reference_dict[key] = ref_lines
    return reference_dict, candidate_dict


def make_ref_mean(matrix: list, ignore_zeros=False) -> float:
    """
    Iterates over a distance matrix and calculates the mean value of all found
    distances. Returns the value as a float. If ignore_zeros is enabled, ignores
    any distance value of zero.
    """
    sum = 0
    zeros_found = 0
    total_number = 0
    for row in matrix:
        for column in row:
            sum += column
            total_number += 1
            if column == 0:
                zeros_found += 1
    if ignore_zeros:
        total_number -= zeros_found
    mean = sum / total_number
    return mean


def candidate_pairwise_calls(candidate: Record, refs: list) -> list:
    """
    Calls calc._pairwise on a candidate and each ref sequence in the list.
    Returns the distances as list. Used to avoid recalculating the ref distances
    for any given candidate index.
    """
    result = list()
    for ref in refs:
        result.append(bd.blosum62_distance(str(candidate), str(ref)))
    result.append(0.0)
    return result


def compare_means(references: list, candidates: list, threshold: float,excluded_headers: set, keep_refs: bool, sort: str) -> tuple:
    """
    For each candidate record, finds the index of the first non-gap bp and makes
    matching cuts in the reference sequences. Afterwards finds the mean of the trimmed
    data.
    """
    regulars = list()
    outliers = list()
    if keep_refs:
        for line in references:
            regulars.append(line)
    ref_dict, candidates_dict = find_index_groups(references, candidates)
    to_add_later =list()
    for index_pair in ref_dict:
        # start, stop = index_pair
        current_refs = ref_dict[index_pair]
        candidates_at_index = candidates_dict[index_pair]
        # first we have to calculate the reference distances to make the ref mean
        convert_to_record_objects(current_refs)
        ref_alignments = [seq for seq in convert_to_record_objects(current_refs) if seq.id not in excluded_headers]
        ref_distances = list()
        for seq1, seq2 in combinations(ref_alignments, 2):
            ref1 = str(seq1)
            ref2 = str(seq2)
            ref_distances.append(bd.blosum62_distance(ref1, ref2))
        # First quartile (Q1)
        try:
            Q1 = np.percentile(ref_distances, 25, method = 'midpoint')
        except IndexError:
            Q1 = 0.0
        # Third quartile (Q3)
        try:
            Q3 = np.percentile(ref_distances, 75, method = 'midpoint')
        except IndexError:
            Q3 = 0.0
        # Interquartile range (IQR)
        IQR = Q3 - Q1
        upper_bound = Q3 + (threshold * IQR) + .02
        intermediate_list = list()
        for candidate in candidates_at_index:
            candidate_distances = candidate_pairwise_calls(candidate, ref_alignments)
            mean_distance = mean(candidate_distances)
            header = candidate.id
            raw_sequence = candidate.raw
            grade = 'Fail'
            if mean_distance <= upper_bound:
                if sort == 'original':
                    to_add_later.append(header)
                    to_add_later.append(raw_sequence)
                elif sort == 'cluster':
                    intermediate_list.append(header)
                    intermediate_list.append(raw_sequence)
                grade = 'Pass'
            outliers.append((header, mean_distance, upper_bound, grade, IQR))
        if sort == 'cluster':
            intermediate_list = taxa_sort(intermediate_list)
            to_add_later.extend(intermediate_list)
    return regulars, to_add_later, outliers


def make_nt_name(path: str) -> str:
    folder, name = os.path.split(path)
    fields = name.split('.')
    name, ext = fields[0], fields[-1]
    name = name + '.nt.'
    top_level = '/'.join(folder.split('/')[:-1])
    nt_folder = os.path.join(top_level, 'nt')
    result = os.path.join(nt_folder, name)
    return result, name, top_level


def deinterleave(fasta_lines: list) -> list:
    result = []
    this_out = list()
    for line in fasta_lines:
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


def delete_empty_columns(raw_fed_sequences: list) -> list:
    """
    Iterates over each sequence and deletes columns
    that consist of 100% dashes.
    """
    result = []
    sequences = []
    raw_sequences = [i.replace('\n','') for i in raw_fed_sequences]

    while '' in raw_sequences:
        raw_sequences.remove('')

    for i in range(0, len(raw_sequences), 2):
        sequences.append(raw_sequences[i+1])

    positions_to_keep = []
    if sequences != []:
        for i in range(len(sequences[0])):
            for sequence in sequences:
                if sequence[i] != '-':
                    positions_to_keep.append(i)
                    break
        for i in range(0, len(raw_sequences), 2):
            result.append(raw_sequences[i]+'\n')
            
            sequence = [raw_sequences[i+1][x] for x in positions_to_keep]
            sequence.append('\n')
            sequence = ''.join(sequence)
        
            result.append(sequence)

    return result


def remove_excluded_sequences(lines: list, excluded: set) -> list:
    """
    Given a list of fasta lines and a set of headers to be excluded from output,
    returns a list of all valid headers and sequences. Use before the delete_column
    call in the nt portion.
    """
    output = list()
    for i in range(0, len(lines), 2):
        if lines[i].strip() not in excluded:
            output.append(lines[i])
            output.append(lines[i+1])
    return output


def main_process(args_input, args_output, args_threshold, args_references, sort: str):
    if not args_references:
        keep_refs = True
    else:
        keep_refs = False

    file_input = args_input
    filename = os.path.basename(file_input)
    name = filename.split('.')[0]

    threshold = args_threshold/100
    aa_output = os.path.join(args_output, 'aa')
    aa_output = os.path.join(aa_output, filename)

    outliers_csv = os.path.join(args_output+'/logs', 'outliers_'+name+'.csv')
    outliers_csv = open(outliers_csv, 'w+')

    lines = list()
    with open(file_input) as f:
        lines = f.readlines()
        lines = deinterleave(lines)
        # print(file_input)
    to_be_excluded = set()
    reference_sequences, candidate_sequences = split_sequences(lines, to_be_excluded)
    candidate_headers = [header for header in candidate_sequences if header[0] == '>']
    raw_regulars, to_add, outliers = compare_means(reference_sequences, candidate_sequences, threshold, to_be_excluded, keep_refs, sort)
    if sort == 'original':
        to_add = original_sort(candidate_headers, to_add)

    for line in to_add:
        raw_regulars.append(line)


    regulars = delete_empty_columns(raw_regulars)

    if len(to_add) > 0:  # If candidate added to fasta
        aa_output = open(aa_output,'w+')
        aa_output.writelines(regulars)
        to_be_excluded = set()
        for outlier in outliers:
            header, distance, ref_dist, grade, iqr = outlier
            if grade == 'Fail':
                to_be_excluded.add(header)

            header = header[1:]
            result = [header, str(distance), str(ref_dist),str(iqr), grade]
            outliers_csv.write(','.join(result)+'\n')
        
        nt_input_path, nt_filename, top_level_folder = make_nt_name(args_input)
        nt_output_path = os.path.join(args_output, 'nt')
        allowed_extensions = {'fa', 'fas', 'fasta'}
        for ext in allowed_extensions:
            nt_ext_name = nt_filename+ext
            nt_trial_path = nt_input_path+ext
            if os.path.exists(nt_trial_path):
                nt_input_path = nt_trial_path
                nt_output_path = os.path.join(nt_output_path, nt_ext_name)
                break
        nt_output_handle = open(nt_output_path,'w+')
        with open(nt_input_path) as f:
            lines = f.readlines()
            de_lines = deinterleave(lines)
            non_empty_lines = remove_excluded_sequences(de_lines, to_be_excluded)
            non_empty_lines = delete_empty_columns(non_empty_lines)
        for i in range(0, len(non_empty_lines), 2):
            nt_output_handle.write(non_empty_lines[i])
            nt_output_handle.write(non_empty_lines[i+1])
        nt_output_handle.close()


def run_command(arg_tuple: tuple) -> None:
    input, output, threshold, references_args, sort = arg_tuple
    main_process(input, output, threshold, references_args, sort)


if __name__ == '__main__':
    start = time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-aa', '--aa_input', default='aa', help="Source of AA files")
    parser.add_argument('-o', '--output', default='outlier', help="Output folder")
    parser.add_argument('-p', '--processes', type=int, default=0,
                        help='Number of threads used to call processes.')
    parser.add_argument('-t', '--threshold', type=int, default=50,
                        help='Greater than reference mean to be counted as an outlier. Default is 2x.')
    parser.add_argument('--no-references', action='store_true',
                        help='Disable output of reference sequences')
    parser.add_argument('-s', '--sort', choices=['cluster', 'original'], default='original',
                        help="Sort candidate output by cluster and taxa, or preserver original order.")
    args = parser.parse_args()
    allowed_extensions = {'fa', 'fas', 'fasta'}
    file_inputs = [args.aa_input+'/'+gene for gene in os.listdir(args.aa_input) if '.aa' in gene and gene.split('.')[-1] in allowed_extensions]

    folder_check(args.output)

    if args.processes:
        arguments = list()
        for gene in file_inputs:
            arguments.append((gene,args.output,args.threshold,args.no_references, args.sort))

        with Pool(args.processes) as pool:
            pool.map(run_command, arguments, chunksize=1)
    else:
        for gene in file_inputs:
            print(gene)
            main_process(gene,args.output,args.threshold,args.no_references, args.sort)

    log_folder_path = os.path.join(args.output, 'logs')
    global_csv_path = os.path.join(log_folder_path, 'outliers_global.csv')

    logs = [x for x in os.listdir(log_folder_path) if 'outliers_' in x and 'global' not in x]
    with open(global_csv_path, 'w') as global_csv:
        global_csv.write('Gene,Header,Mean_Dist,Ref_Mean,IQR\n')
        for log in logs:
            log_file_path = os.path.join(log_folder_path, log)
            with open(log_file_path) as log_f:
                for line in log_f:
                    if line.strip().split(',')[-1] == 'Fail':
                        global_csv.write(line)
                        if line[-1] != '\n':
                            global_csv.write('\n')
    time_taken = time()
    time_taken = time_taken-start

    print('Finished in {} seconds'.format(round(time_taken)))
    #open('runtime.txt','w').writelines('Finished in {} seconds'.format(round(time_taken)))
