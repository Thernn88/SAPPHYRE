"""
Merges all sequences per taxa into single sequence in each gene

PyLint 9.61/10
"""
import argparse
import os
from multiprocessing.pool import Pool
from time import time
from typing import Union


def make_seq_dict(sequences: list, trailing_end: int) -> dict:
    """
    Creates a dictionary of each sequence present at each coordinate of a list of sequences
    """

    seq_dict = {}

    for i in range(trailing_end):
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


def format_taxa(raw_taxa: str) -> str:
    """
    Formats rerun taxas to origianl taxa.
    """

    if "R" in raw_taxa.split("_")[-1]:

        after_r = raw_taxa.split("_")[-1].replace("R", "")
        character_after_r = after_r[0]
        if character_after_r.isdigit():
            print("_R Triggered. Comment line 87 to disable this")
            # If _R# in header name strip
            formatted_taxa = "_".join(raw_taxa.split("_")[:-1])
    else:
        taxa_min = formatted_taxa.count("_") > 2 or formatted_taxa.count("_") == 1
        if formatted_taxa.split("_")[-1].isdigit() and taxa_min:
            # If character after last underscore is digit and header consists of
            # 1 or more than 2 headers

            # Strip last underscore
            formatted_taxa = "_".join(formatted_taxa.split("_")[:-1])
            print("_# Triggered. Comment line 93 to disable this")

    return formatted_taxa


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


def disperse_into_individual_taxa(sequence_pair: list) -> dict:
    """
    Splits list of (header,sequence) into taxa based dict.
    """

    result = {}

    for header, sequence in sequence_pair:
        original_taxa = header.split("|")[2]
        taxa_formatted = format_taxa(original_taxa)

        if taxa_formatted not in result:
            result[taxa_formatted] = []

        start, end = get_start_end(sequence)

        this_object = (start, end, header, sequence)

        result[taxa_formatted].append(this_object)

    return result


def expand_region(original: tuple, expansion: tuple) -> tuple:
    """
    Expands two (start, end) tuples to cover the entire region
    """

    start = original[0]
    if expansion[0] > start:
        start = expansion[0]

    end = original[1]
    if expansion[1] > end:
        end = expansion[1]

    return (start, end)


def disperse_into_overlap_groups(taxa_pair: list) -> dict:
    """
    Splits list of (header,sequence) into overlap based groups.
    """

    result = []
    current_group = []

    current_region = None

    for start, end, header, sequence in taxa_pair:
        this_object = (start, end, header, sequence)

        if current_region is None or find_overlap((start, end), current_region) is None:
            if current_group:
                result.append(current_group)

            current_group = [this_object]
        else:
            current_group.append(this_object)

        if current_region is not None:
            current_region = expand_region(current_region, (start, end))
        else:
            current_region = (start, end)

    if current_group:
        result.append(current_group)

    return result


def find_overlap(tuple_a: tuple, tuple_b: tuple) -> Union[tuple, None]:
    """
    Takes two start/end pairs and returns the overlap.
    """
    start = max(tuple_a[0], tuple_b[0])
    end = min(tuple_a[1], tuple_b[1])
    if end - start < 0:
        return None, None
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


def directory_check(target_output_path) -> None:
    """
    Creates necessary directories for merge output.
    """

    aa_merged = os.path.join(target_output_path, "aa_merged")
    nt_merged = os.path.join(target_output_path, "nt_merged")

    if not os.path.exists(aa_merged):
        os.mkdir(aa_merged)

    if not os.path.exists(nt_merged):
        os.mkdir(nt_merged)


def main(
    gene,
    output_dir,
    aa_path,
    nt_path,
    fallback_taxa,
    debug,
    majority,
    minimum_mr_amount,
) -> None:
    """
    Merge main loop. Opens fasta file, parses sequences and merges based on taxa
    """
    already_calculated_splits = {}

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

        overlap_regions_to_sequences = disperse_into_individual_taxa(candidates)

        for this_overlap_sequences in overlap_regions_to_sequences.values():
            this_overlap_sequences.sort(key=lambda x: x[0])

            overlap_groups = disperse_into_overlap_groups(this_overlap_sequences)

            for this_sequences in overlap_groups:
                base_header = this_sequences[0][2]

                this_gene, this_taxa, this_taxa_id, node = base_header.split("|")
                
                # Formats base header taxa to match every component's formatted taxa
                formatted_taxa_id = format_taxa(this_taxa_id)
                base_header = "|".join([this_gene, taxa, formatted_taxa_id, node])

                consists_of = [
                    (sequence[2], sequence[3]) for sequence in this_sequences
                ]

                # Gets last pipe of each component of a merge
                stitch = "&&".join(
                    [sequence[2].split("|")[-1] for sequence in this_sequences[1:]]
                )

                if stitch != "":
                    final_header = base_header + "&&" + stitch
                else:  # No merge occurs
                    final_header = base_header

                trailing_end = len(this_sequences[-1][3])

                new_merge = []

                current_point_seqs = make_seq_dict(this_sequences, trailing_end)

                for cursor in range(trailing_end):

                    sequences_at_current_point = current_point_seqs[cursor]
                    if len(sequences_at_current_point) == 0:
                        new_merge.append("-")
                    elif len(sequences_at_current_point) == 1:
                        new_merge.append(sequences_at_current_point[0][1][cursor])
                    elif len(sequences_at_current_point) > 1:

                        splits = int(len(sequences_at_current_point) / 2)
                        taxas_of_split = [
                            header.split("|")[1]
                            for (header, sequence) in sequences_at_current_point
                        ]

                        most_occuring = most_common_element_with_count(taxas_of_split)
                        if most_occuring[1] == 1:  # No taxa occurs more than once
                            comparison_taxa = fallback_taxa
                        else:
                            comparison_taxa = most_occuring[0]

                        comparison_sequence = comparison_sequences.get(
                            comparison_taxa, comparison_sequences[fallback_taxa]
                        )

                        this_splits_headers = []
                        this_splits_sequences = []

                        for split_count in range(splits):
                            header_a, sequence_a = sequences_at_current_point[
                                split_count
                            ]
                            header_b, sequence_b = sequences_at_current_point[
                                split_count + 1
                            ]

                            split_key = header_a + header_b

                            this_splits_headers.append((header_a, header_b))
                            this_splits_sequences.append((sequence_a, sequence_b))

                        next_character = this_splits_sequences[0][0][
                            cursor
                        ]  # Add from sequence A if cursor is not past any split position

                        for i, pair in enumerate(this_splits_headers):
                            header_a, header_b = pair
                            sequence_a, sequence_b = this_splits_sequences[i]

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
                    majority_assignments = []  # Used for debug log
                    start_ends = {}
                    for i, char in enumerate(new_merge):
                        candidate_characters = []
                        candidate_sequence = []
                        for header, sequence in consists_of:
                            if header not in start_ends:
                                start_ends[header] = get_start_end(sequence)
                            start, end = start_ends[header]

                            if start <= i <= end:
                                candidate_characters.append(sequence[i])
                                candidate_sequence.append(header)

                        mr_won = False
                        if len(candidate_characters) >= minimum_mr_amount:
                            (
                                mode_cand_character,
                                mode_count,
                            ) = most_common_element_with_count(candidate_characters)
                            perc_appearing = mode_count / len(candidate_characters)
                            if perc_appearing >= majority:
                                if mode_cand_character != char:
                                    new_merge[i] = mode_cand_character
                                    majority_assignments.append(mode_cand_character)
                                    mr_won = True
                        if not mr_won:
                            majority_assignments.append("-")
                    new_merge = "".join(new_merge)

                elif protein == "nt":
                    candidate_sequence = {}
                    new_merge = "".join(new_merge)
                    length = len(new_merge)
                    new_merge = [new_merge[i : i + 3] for i in range(0, length, 3)]
                    majority_assignments = []  # Used for debug log
                    start_ends = {}
                    for i, char in enumerate(new_merge):
                        adjusted_i = range(0, length, 3)[i]
                        candidate_characters = []
                        candidate_sequence = []
                        for header, sequence in consists_of:
                            if header not in start_ends:
                                start_ends[header] = get_start_end(sequence)
                            start, end = start_ends[header]

                            if start <= adjusted_i <= end:
                                candidate_characters.append(
                                    sequence[adjusted_i : adjusted_i + 3]
                                )
                                candidate_sequence.append(header)

                        mr_won = False

                        if len(candidate_characters) >= minimum_mr_amount:

                            translated_characters = [
                                dna_codons[i] for i in candidate_characters
                            ]

                            (
                                most_occuring_translated_char,
                                translated_char_count,
                            ) = most_common_element_with_count(translated_characters)

                            perc_appearing = translated_char_count / len(
                                translated_characters
                            )
                            candidate_characters = [
                                i
                                for i in candidate_characters
                                if dna_codons[i] == most_occuring_translated_char
                            ]

                            (
                                mode_cand_raw_character,
                                _,
                            ) = most_common_element_with_count(candidate_characters)

                            if perc_appearing >= majority:
                                if mode_cand_raw_character != char:

                                    new_merge[i] = mode_cand_raw_character
                                    majority_assignments.append(mode_cand_raw_character)

                                    mr_won = True
                        if not mr_won:
                            majority_assignments.append("---")
                    new_merge = "".join(new_merge)

                gene_out.append(final_header)
                gene_out.append("".join(new_merge))

                if debug:
                    # If debug enabled add each component under final merge
                    gene_out.append(">" + this_taxa_id + "|MajorityRulesAssigned")
                    gene_out.append("".join(majority_assignments))
                    for header, sequence in consists_of:
                        this_taxa_id = header.split("|")[-2]
                        last_pipe = this_taxa_id + "|" + header.split("|")[-1]
                        gene_out.append(">" + last_pipe)
                        gene_out.append(sequence)

        if protein == "aa":
            output_path = os.path.join(output_dir, "aa_merged", gene)
        else:
            output_path = os.path.join(output_dir, "nt_merged", gene)

        with open(output_path, "w", encoding="UTF-8") as output_file:
            output_file.write("\n".join(gene_out))


def run_command(arg_tuple: tuple) -> None:
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
        majority,
        majority_count,
    ) = arg_tuple
    main(
        gene, output_dir, aa_path, nt_path, comparison, debug, majority, majority_count
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
            directory_check(taxa_path)

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
                    print(target_gene)
                    target_aa_path = os.path.join(aa_input, target_gene)
                    target_nt_path = os.path.join(nt_input, make_nt_name(target_gene))
                    main(
                        target_gene,
                        taxa_path,
                        target_aa_path,
                        target_nt_path,
                        args.comparison,
                        args.debug,
                        args.majority,
                        args.majority_count,
                    )

            timed = round(time() - start_time)
            print(f"Finished in {timed} seconds")
        else:
            print(f"Can't find aa folder for taxa, {taxa}")
