"""
Merges all sequences per taxa into single sequence in each gene

PyLint 9.61/10
"""
import argparse
import os
from multiprocessing.pool import Pool
from time import time
from typing import Union
import wrap_rocks
import json
from Bio.Seq import Seq

def make_seq_dict(sequences: list, data_region: tuple) -> dict:
    """
    Creates a dictionary of each sequence present at each coordinate of a list of sequences
    """

    seq_dict = {}

    for i in range(data_region[0], data_region[1] + 1): # exclusive, so add 1
        seq_dict[i] = []

    for sequence in sequences:
        for cursor in range(sequence[0], sequence[1] + 1):  # exclusive, so add 1
            seq_dict[cursor].append((sequence[2], sequence[3]))
    return seq_dict


def most_common_element_with_count(iterable) -> tuple:
    """
    Returns the most common element and corresponding
    count within an iterable
    """
    counts = {}
    winner = ("dummy", -1)
    for element in iterable:
        count = counts.get(element, 0) + 1
        counts[element] = count
        if count > winner[1]:
            winner = (element, count)
    return winner


def make_nt_name(file_name: str) -> str:
    """
    Universal rename for aa file format to nt.
    """
    file = file_name.split(".")
    file[file.index("aa")] = "nt"
    result = ".".join(file)
    return result


def is_reference(header: str) -> bool:
    """
    Returns True if the reference header identifier is present in the header
    """

    return header[-1] == '.'


def parse_fasta(text_input: str) -> list:
    """
    Returns references from raw fasta text input.
    """
    lines = text_input.split("\n")

    references = []
    candidates = []

    while "" in lines:
        lines.remove("") # Remove blank lines

    end_of_references = False
    for i in range(0, len(lines), 2):
        header = lines[i]
        sequence = lines[i + 1]

        if end_of_references is False:
            if is_reference(header):
                references.append((header, sequence))
            else:
                end_of_references = True

        if end_of_references is True:
            candidates.append((header, sequence))

    return references, candidates

def get_start_end(sequence: str) -> tuple:
    """
    Returns index of first and last none dash character in sequence.
    """
    start = None
    end = None
    for i, character in enumerate(sequence):
        if character != "-":
            start = i
            break
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != "-":
            end = i
            break
    return (start, end)

def expand_region(original: tuple, expansion: tuple) -> tuple:
    """
    Expands two (start, end) tuples to cover the entire region
    """

    start = original[0]
    if expansion[0] < start:
        start = expansion[0]

    end = original[1]
    if expansion[1] > end:
        end = expansion[1]

    return (start, end)


def disperse_into_overlap_groups(taxa_pair: list) -> dict:
    """
    Splits list of (header,sequence) into overlap based groups.

    Returns (overlap region, sequences in overlap region)
    """

    result = []
    current_group = []

    current_region = None

    for start, end, header, sequence in taxa_pair:
        this_object = (start, end, header, sequence)

        if current_region is None or find_overlap((start, end), current_region) is None:
            if current_group:
                result.append((current_region, current_group))

            current_region = (start, end)
            current_group = [this_object]
        else:
            current_group.append(this_object)
            current_region = expand_region(current_region, (start, end))


    if current_group:
        result.append((current_region, current_group))

    return result


def find_overlap(tuple_a: tuple, tuple_b: tuple) -> Union[tuple, None]:
    """
    Takes two start/end pairs and returns the overlap.
    """
    start = max(tuple_a[0], tuple_b[0])
    end = min(tuple_a[1], tuple_b[1])
    if end - start < 0:
        return None
    return start, end


def calculate_split(sequence_a: str, sequence_b: str, comparison_sequence: str) -> int:
    """
    Iterates over each position in the overlap range of sequence A and sequence B and
    creates a frankenstein sequence of sequence A + Sequence B joined at each
    position in the overlap.

    Final split position = pos in overlap with highest score.

    Score is determined by the amount of characters that are the same between each
    position in the frankenstein sequence and the comparison sequence.
    """

    pair_a = get_start_end(sequence_a)
    pair_b = get_start_end(sequence_b)

    overlap_start, overlap_end = find_overlap(pair_a, pair_b)

    sequence_a_overlap = sequence_a[overlap_start : overlap_end + 1]
    sequence_b_overlap = sequence_b[overlap_start : overlap_end + 1]
    comparison_overlap = comparison_sequence[overlap_start : overlap_end + 1]

    highest_score = 0
    base_score = 0
    highest_scoring_pos = 0

    for i, character in enumerate(sequence_b_overlap):
        if character == comparison_overlap[i]:
            base_score += 1
    highest_score = base_score

    for i, character_a in enumerate(sequence_a_overlap):
        if sequence_b_overlap[i] == comparison_overlap[i]:
            base_score -= 1
        if character_a == comparison_overlap[i]:
            base_score += 1
        if base_score >= highest_score:
            highest_score = base_score
            highest_scoring_pos = i
    return highest_scoring_pos + overlap_start


def directory_check(target_output_path) -> str:
    """
    Creates necessary directories for merge output.
    """

    aa_merged = os.path.join(target_output_path, "aa_merged")
    nt_merged = os.path.join(target_output_path, "nt_merged")

    if not os.path.exists(aa_merged):
        os.mkdir(aa_merged)

    if not os.path.exists(nt_merged):
        os.mkdir(nt_merged)

    if os.path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif os.path.exists("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = os.path.join(target_output_path, "tmp")
        if not os.path.exists(tmp_path):
            os.mkdir(tmp_path)
    
    return tmp_path


def main(
    gene,
    output_dir,
    aa_path,
    nt_path,
    dupe_tmp_file,
    fallback_taxa,
    debug,
    majority,
    minimum_mr_amount,
) -> None:
    """
    Merge main loop. Opens fasta file, parses sequences and merges based on taxa
    """
    print(gene)
    already_calculated_splits = {}

    

    with open(dupe_tmp_file) as dupe_tmp_in:
        dupe_counts = json.load(dupe_tmp_in)

    for protein in ["aa", "nt"]:
        if protein == "aa":
            with open(aa_path, encoding="UTF-8") as aa_in:
                fasta_text = aa_in.read()
        else:
            with open(nt_path, encoding="UTF-8") as nt_in:
                fasta_text = nt_in.read()
            gene = make_nt_name(gene)

        references, candidates = parse_fasta(fasta_text)

        gene_out = []
        comparison_sequences = {}
        for header, sequence in references:
            gene_out.append(header)
            gene_out.append(sequence)
            this_taxa = header.split("|")[1]
            comparison_sequences[this_taxa] = sequence

        sequenecs_to_merge = []
        for header, sequence in candidates:
            start, end = get_start_end(sequence)

            this_object = (start, end, header, sequence)
            sequenecs_to_merge.append(this_object)

        sequenecs_to_merge.sort(key = lambda x : x [0])

        overlap_groups = disperse_into_overlap_groups(sequenecs_to_merge)

        for overlap_region, this_sequences in overlap_groups:

            # Use the header of the sequence that starts first as the base
            base_header = this_sequences[0][2]

            this_gene, this_taxa, this_taxa_id, node, coords, frame = base_header.split("|")
            base_header = "|".join([this_gene, this_taxa, this_taxa_id, node])

            # Create a list of each header and sequence that is present in this overlap region
            consists_of = []
            
            for sequence in this_sequences:
                # If the sequence had a dupe removed during dataset preperation dedupe
                # add 1 to it's MR count.
                count = 1 if sequence not in dupe_counts else dupe_counts[sequence]
                
                consists_of.append((sequence[2], sequence[3], count))

            # Gets last pipe of each component of a merge
            stitch = "&&".join(
                [sequence[2].split("|")[3] for sequence in this_sequences[1:]]
            )

            
            if stitch != "":
                final_header = f"{base_header}&&{stitch}"
            else:  # If only single component aka no merge occurs don't change header
                final_header = base_header

            # The overlap region is the first and last index of data
            data_start, data_end = overlap_region

            # As the fasta is aligned the trailing end of all sequences can be assumed to be the same as the length of the last sequence
            trailing_end = len(this_sequences[-1][3])

            # Add gap characters until the first data index
            new_merge = ["-"] * data_start

            # Create a hashmap of the headers present at every position in the range of
            # the overlap region
            current_point_seqs = make_seq_dict(this_sequences, overlap_region)


            for cursor in range(data_start, data_end+1):
                # Grab the sequences present at this current index
                sequences_at_current_point = current_point_seqs[cursor]
                amount_of_sequences = len(sequences_at_current_point)

                if amount_of_sequences == 1:
                    # If there's only one sequence just add the character present at that point
                    new_merge.append(sequences_at_current_point[0][1][cursor])

                elif amount_of_sequences > 1:
                    # If there is more than one sequence at this current index
                    splits = amount_of_sequences - 1

                    # Grab most occuring taxon
                    taxons_of_split = [
                        header.split("|")[1]
                        for (header, _) in sequences_at_current_point
                    ]

                    most_occuring = most_common_element_with_count(taxons_of_split)
                    if most_occuring[1] == 1:  # No taxa occurs more than once
                        comparison_taxa = fallback_taxa
                    else:
                        comparison_taxa = most_occuring[0]

                    # Grab the reference sequence for the mode taxon
                    comparison_sequence = comparison_sequences.get(
                        comparison_taxa, comparison_sequences[fallback_taxa]
                    )

                    # Set next character to the first possible sequence.

                    next_character = sequences_at_current_point[0][1][
                        cursor
                    ] 

                    # For every 
                    for split_count in range(splits):
                        header_a, sequence_a = sequences_at_current_point[
                            split_count
                        ]
                        header_b, sequence_b = sequences_at_current_point[
                            split_count + 1
                        ]

                        split_key = header_a + header_b

                        if protein == "aa":
                            if split_key in already_calculated_splits:
                                split_position = already_calculated_splits[
                                    split_key
                                ]
                            else:
                                split_position = calculate_split(
                                    sequence_a, sequence_b, comparison_sequence
                                )
                                already_calculated_splits[
                                    split_key
                                ] = split_position

                        elif protein == "nt":
                            split_position = (
                                already_calculated_splits[split_key] * 3
                            )

                        if cursor >= split_position:
                            # Iterate over each split if cursor is past calculated split
                            # position add from sequence B. We only want to add from one
                            # sequence out of every possible split so we calculate which
                            # sequence to add from here then add the character to the
                            # final merge in the next line.
                            next_character = sequence_b[cursor]

                    new_merge.append(next_character)

            new_merge.extend(["-"] * (trailing_end - data_end))

            # Doing MR Outlier Check
            #
            # Each position of the new merge is checked with each position of overlapping
            # candidates at that point. If majority (args.majority percent) of characters
            # are the same but differ from the character chosen in the merge it will replace
            # it with the majority rules letter.

            dna_codons = {
                "GCT": "A",
                "GCC": "A",
                "GCA": "A",
                "GCG": "A",
                "TGT": "C",
                "TGC": "C",
                "GAT": "D",
                "GAC": "D",
                "GAA": "E",
                "GAG": "E",
                "TTT": "F",
                "TTC": "F",
                "GGT": "G",
                "GGC": "G",
                "GGA": "G",
                "GGG": "G",
                "CAT": "H",
                "CAC": "H",
                "ATA": "I",
                "ATT": "I",
                "ATC": "I",
                "AAA": "K",
                "AAG": "K",
                "TTA": "L",
                "TTG": "L",
                "CTT": "L",
                "CTC": "L",
                "CTA": "L",
                "CTG": "L",
                "ATG": "M",
                "AAT": "N",
                "AAC": "N",
                "CCT": "P",
                "CCC": "P",
                "CCA": "P",
                "CCG": "P",
                "CAA": "Q",
                "CAG": "Q",
                "CGT": "R",
                "CGC": "R",
                "CGA": "R",
                "CGG": "R",
                "AGA": "R",
                "AGG": "R",
                "TCT": "S",
                "TCC": "S",
                "TCA": "S",
                "TCG": "S",
                "AGT": "S",
                "AGC": "S",
                "ACT": "T",
                "ACC": "T",
                "ACA": "T",
                "ACG": "T",
                "GTT": "V",
                "GTC": "V",
                "GTA": "V",
                "GTG": "V",
                "TGG": "W",
                "TAT": "Y",
                "TAC": "Y",
                "TAA": "*",
                "TAG": "*",
                "TGA": "*",
                "---": "---",
            }

            if protein == "aa":
                if debug:
                    majority_assignments = ['-' * data_start]  # Used for debug log
                start_ends = {}
                for i in range(data_start, data_end+1):
                    char = new_merge[i]
                    candidate_characters = {}
                    total_characters = 0
                    mode = -999
                    mode_char = None
                    for header, sequence, count in consists_of:
                        if header not in start_ends:
                            start_ends[header] = get_start_end(sequence)
                        start, end = start_ends[header]

                        if start <= i <= end:
                            total_characters += count
                            this_char = sequence[i]
                            this_count = (candidate_characters[this_char] + 1) if this_char in candidate_characters else count
                            candidate_characters[this_char] = this_count
                            if this_count > mode:
                                mode_char = this_char
                                mode = this_count

                    mr_won = False
                    if total_characters >= minimum_mr_amount:
                        perc_appearing = mode / total_characters
                        if perc_appearing >= majority:
                            if mode_char != char:
                                new_merge[i] = mode_char
                                if debug:
                                    majority_assignments.append(mode_char)
                                mr_won = True
                    if not mr_won:
                        if debug:
                            majority_assignments.append("-")
                if debug:
                    majority_assignments.append("-" * (trailing_end - data_end))
            elif protein == "nt":
                length = len(new_merge)
                new_merge = ["".join(new_merge[i : i + 3]) for i in range(0, length, 3)]
                
                start_ends = {}
                
                triplet_data_start = int(data_start / 3)
                triplet_data_end = int(data_end / 3)

                if debug:
                    majority_assignments = ["---" * triplet_data_start]  # Used for debug log

                for raw_i in range(triplet_data_start, triplet_data_end+1):
                    i = (raw_i * 3) 
                    char = new_merge[raw_i]
                    
                    candidate_characters = []
                    total_characters = 0
                    

                    for header, sequence, count in consists_of:
                        if header not in start_ends:
                            start_ends[header] = get_start_end(sequence)
                        start, end = start_ends[header]

                        if start <= i <= end:
                            this_char = sequence[i : i + 3]
                            total_characters += count

                            candidate_characters.append(this_char)        

                    mr_won = False

                    if total_characters >= minimum_mr_amount:
                        # Translate all NT triplets into AA
                        translated_characters = [
                            dna_codons[j] for j in candidate_characters
                        ]

                        # Figure out what AA occurs the most
                        (
                            most_occuring_translated_char,
                            translated_char_count,
                        ) = most_common_element_with_count(translated_characters)

                        #Check to see if its majority
                        perc_appearing = translated_char_count / total_characters

                        if perc_appearing >= majority:
                            # Grab all the NT that map to that AA that triggered majority
                            candidate_chars_mapping_to_same_dna = [
                                j
                                for j in candidate_characters
                                if dna_codons[j] == most_occuring_translated_char
                            ]

                            # Grab the most occuring NT that maps to that AA from the current seq
                            (
                                mode_cand_raw_character,
                                _,
                            ) = most_common_element_with_count(candidate_chars_mapping_to_same_dna)
                            if mode_cand_raw_character != char:

                                new_merge[raw_i] = mode_cand_raw_character
                                if debug:
                                    majority_assignments.append(mode_cand_raw_character)

                                mr_won = True
                    if not mr_won:
                        if debug:
                            majority_assignments.append("---")
                if debug:
                    majority_assignments.append("---" * (length - data_end))  # Used for debug log

            new_merge = Seq("").join(new_merge)

            gene_out.append(final_header)
            gene_out.append(str(new_merge))

            if debug:
                # If debug enabled add each component under final merge
                gene_out.append(">" + this_taxa_id + "|MajorityRulesAssigned")
                gene_out.append("".join(majority_assignments))
                for header, sequence, _ in consists_of:
                    node, _, frame = header.split('|')[3:]
                    gene_out.append(f">{node}|{frame}")
                    gene_out.append(sequence)

        if protein == "aa":
            output_path = os.path.join(output_dir, "aa_merged", gene)
        else:
            output_path = os.path.join(output_dir, "nt_merged", gene)

        with open(output_path, "w", encoding="UTF-8") as fp:
            fp.write("\n".join(gene_out))


def run_command(arg_tuple) -> None:
    """
    Calls the main() function parallel in each thread
    """
    (
        gene,
        output_dir,
        aa_path,
        nt_path,
        comparison,
        debug, 
        dupe_tmp_file,
        majority,
        majority_count,
    ) = arg_tuple
    main(
        gene, output_dir, aa_path, nt_path, comparison, debug, dupe_tmp_file, majority, majority_count
    )


if __name__ == "__main__":
    start_time = time()
    parser = argparse.ArgumentParser("Merges all sequences per taxa in genes.")
    parser.add_argument("-i", "--input", default="Parent", help="Path to parent input")
    parser.add_argument(
        "-aa",
        "--aa_input",
        type=str,
        default="aa",
        help="Path to directory of AA folder",
    )
    parser.add_argument(
        "-nt",
        "--nt_input",
        type=str,
        default="nt",
        help="Path to directory of NT folder",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled displays each component of merged headers.",
    )
    parser.add_argument(
        "-c",
        "--comparison",
        type=str,
        default="Drosophila_melanogaster",
        help="""Fallback Comparison Taxa. Sequence in which
        Sequence A and Sequence B is compared to in split calculation.""",
    )
    parser.add_argument(
        "-m",
        "--majority",
        type=float,
        default=0.66,
        help="Percentage for majority ruling.",
    )
    parser.add_argument(
        "-mc",
        "--majority_count",
        type=int,
        default=4,
        help="Percentage for majority ruling.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=0,
        help="Number of threads used to call processes.",
    )

    args = parser.parse_args()    

    for taxa in os.listdir(args.input):
        print(f"Doing taxa, {taxa}")
        taxa_path = os.path.join(args.input, taxa)
        input_path = os.path.join(taxa_path, "outlier")
        aa_input = os.path.join(input_path, args.aa_input)
        nt_input = os.path.join(input_path, args.nt_input)

        if os.path.exists(aa_input):
            tmp_dir = directory_check(taxa_path)

            dupe_tmp_file = os.path.join(tmp_dir, "DupeSeqs.tmp")

            rocks_db_path = os.path.join(taxa_path, "rocksdb")
            rocksdb_db = wrap_rocks.RocksDB(rocks_db_path)

            with open(dupe_tmp_file, 'w') as dupe_tmp_out:
                dupe_tmp_out.write(rocksdb_db.get('getall:dupes'))

            file_inputs = []
            for item in os.listdir(aa_input):
                if ".fa" in item:
                    file_inputs.append(item)

            if args.processes:
                arguments = []
                for target_gene in file_inputs:
                    target_aa_path = os.path.join(aa_input, target_gene)
                    target_nt_path = os.path.join(nt_input, make_nt_name(target_gene))
                    arguments.append(
                        (
                            target_gene,
                            taxa_path,
                            target_aa_path,
                            target_nt_path,
                            dupe_tmp_file,
                            args.comparison,
                            args.debug,
                            args.majority,
                            args.majority_count,
                        )
                    )
                with Pool(args.processes) as pool:
                    pool.map(run_command, arguments, chunksize=1)
            else:
                for target_gene in file_inputs:
                    target_aa_path = os.path.join(aa_input, target_gene)
                    target_nt_path = os.path.join(nt_input, make_nt_name(target_gene))
                    main(
                        target_gene,
                        taxa_path,
                        target_aa_path,
                        target_nt_path,
                        dupe_tmp_file,
                        args.comparison,
                        args.debug,
                        args.majority,
                        args.majority_count,
                    )

            timed = round(time() - start_time)
            print(f"Finished in {timed} seconds")

            if os.path.exists(dupe_tmp_file):
                os.remove(dupe_tmp_file)
        else:
            print(f"Can't find aa folder for taxa, {taxa}")
