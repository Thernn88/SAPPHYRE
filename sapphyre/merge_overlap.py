"""
Merges all sequences per taxa into single sequence in each gene

PyLint 9.61/10
"""
from __future__ import annotations

import json
import os
from multiprocessing.pool import Pool
from pathlib import Path
from shutil import rmtree
from typing import Union, Literal

import wrap_rocks
from Bio.Seq import Seq

from .utils import printv, parseFasta, writeFasta
from .timekeeper import TimeKeeper, KeeperMode

TMP_PATH = None
if os.path.exists("/run/shm"):
    TMP_PATH = "/run/shm"
elif os.path.exists("/dev/shm"):
    TMP_PATH = "/dev/shm"

DNA_CODONS = {
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


def make_nt_name(x):
    return str(x).replace(".aa.", ".nt.")


def make_seq_dict(sequences: list, data_region: tuple) -> dict:
    """
    Creates a dictionary of each sequence present at each coordinate of a list of sequences
    """
    seq_dict = {
        i: [] for i in range(data_region[0], data_region[1] + 1)
    }  # exclusive, so add 1

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


def parse_fasta(path: str) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """
    Returns references from raw fasta text input.
    """
    references: list[tuple[str, str]] = []
    candidates: list[tuple[str, str]] = []
    end_of_references = False
    try:
        for header, sequence in parseFasta(path):
            if end_of_references is False:
                # the reference header identifier is present in the header
                if header[-1] == ".":
                    references.append((header, sequence))
                else:
                    end_of_references = True

            if end_of_references is True:
                candidates.append((header, sequence))
    except ValueError as e:
        print(f"Fatal error in {os.path.basename(path)}: {e}")

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
    return start, end


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

    return start, end


def disperse_into_overlap_groups(taxa_pair: list) -> list[tuple]:
    """
    Splits list of (header,sequence) into overlap based groups.

    Returns (overlap region, sequences in overlap region)
    """
    result = []
    current_group = []
    current_region = None

    for start, end, header, sequence in taxa_pair:
        if "NODE" in header.split("|")[-1]:
            node = header.split("|")[-1]
            frame = ""
        else:
            node, frame = header.split("|")[-2:]
        this_object = (start, end, header, sequence, node, frame)

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


def find_overlap(
    tuple_a: tuple, tuple_b: tuple, allowed_deviation: int = 1
) -> Union[tuple, None]:
    """
    Takes two start/end pairs and returns the overlap.
    """
    start = max(tuple_a[0], tuple_b[0])
    end = min(tuple_a[1], tuple_b[1])
    if end - start < -allowed_deviation:
        return None
    return start, end


def calculate_split(sequence_a: str, sequence_b: str, comparison_sequence: str) -> int:
    """
    Iterates over each position in the overlap range of sequence A and sequence B and
    creates a frankenstein sequence of sequence A + Sequence B joined at each
    position in the overlap.

    Final split position = pos in overlap with the highest score.

    Score is determined by the amount of characters that are the same between each
    position in the frankenstein sequence and the comparison sequence.
    """
    pair_a = get_start_end(sequence_a)
    pair_b = get_start_end(sequence_b)

    overlap_start, overlap_end = find_overlap(pair_a, pair_b, 0)

    sequence_a_overlap = sequence_a[overlap_start : overlap_end + 1]
    sequence_b_overlap = sequence_b[overlap_start : overlap_end + 1]
    comparison_overlap = comparison_sequence[overlap_start : overlap_end + 1]

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
    aa_merged_path = os.path.join(target_output_path, "aa_merged")
    nt_merged_path = os.path.join(target_output_path, "nt_merged")
    rmtree(aa_merged_path, ignore_errors=True)
    rmtree(nt_merged_path, ignore_errors=True)
    os.mkdir(aa_merged_path)
    os.mkdir(nt_merged_path)
    if not TMP_PATH:
        tmp_path = os.path.join(target_output_path, "tmp")
        os.makedirs(tmp_path, exist_ok=True)
        return tmp_path
    return TMP_PATH


def grab_merge_start_end(taxa_pair: list) -> list[tuple]:
    """
    Grabs start and end of merge sequences.
    """

    merge_start = min(seq[0] for seq in taxa_pair)
    merge_end = max(seq[1] for seq in taxa_pair)
    sequences = []
    for start, end, header, sequence in taxa_pair:
        if "NODE" in header.split("|")[-1]:
            node = header.split("|")[-1]
            frame = ""
        else:
            node, frame = header.split("|")[-2:]
        this_object = (start, end, header, sequence, node, frame)
        sequences.append(this_object)

    return [((merge_start, merge_end), sequences)]


def do_protein(
    protein: Literal["aa", "nt"],
    path,
    output_dir: Path,
    prep_dupe_counts,
    rep_dupe_counts,
    ref_stats,
    already_calculated_splits,
    gene,
    majority,
    minimum_mr_amount,
    ignore_overlap_chunks,
    debug=None,
):
    references, candidates = parse_fasta(path)

    gene_out = []

    def get_ref(header: str) -> str:
        return header.split("|")[1]

    # Grab all the reference sequences
    comparison_sequences = {}
    for header, sequence in references:
        gene_out.append((header, sequence))
        if protein == "aa":
            taxon = header.split("|")[1]
            comparison_sequences[taxon] = sequence

    # Grab all the candidate sequences and sort into taxa id based groups
    taxa_groups = {}

    def get_taxa(header: str) -> str:
        return header.split("|")[2]

    for header, sequence in candidates:
        start, end = get_start_end(sequence)
        this_object = (start, end, header, sequence)
        taxa = get_taxa(header)
        taxa_groups.setdefault(taxa, []).append(this_object)

    for _, sequences_to_merge in taxa_groups.items():
        # Sort by start position
        sequences_to_merge.sort(key=lambda x: x[0])

        # Disperse sequences into clusters of overlap
        if ignore_overlap_chunks:
            overlap_groups = grab_merge_start_end(sequences_to_merge)
        else:
            overlap_groups = disperse_into_overlap_groups(sequences_to_merge)

        for overlap_region, this_sequences in overlap_groups:
            # Use the header of the sequence that starts first as the base
            base_header = this_sequences[0][2]

            try:
                this_gene, _, this_taxa_id, node, frame = base_header.split("|")
            except ValueError:
                this_gene, _, this_taxa_id, node = base_header.split("|")

            # Create a list of each header and sequence that is present in this overlap region
            consists_of = []

            for sequence in this_sequences:
                # if protein == "aa":
                count = prep_dupe_counts.get(sequence[4], 1) + sum(
                    [
                        prep_dupe_counts.get(header, 1)
                        for header in rep_dupe_counts.get(sequence[4], [])
                    ]
                )
                # else:
                # count = dupe_counts.get(sequence[4], 0) + 1
                consists_of.append((sequence[2], sequence[3], count))

            # Gets last pipe of each component of a merge
            stitch = "&&".join(
                [sequence[2].split("|")[3] for sequence in this_sequences[1:]]
            )

            taxons = [i[2].split("|")[1] for i in this_sequences]
            taxons_filtered = [i for i in taxons if i in ref_stats]
            if not taxons_filtered:
                this_taxa = taxons[0]
            else:
                most_occuring = most_common_element_with_count(taxons_filtered)
                if len(taxons_filtered) > 1 and most_occuring[1] > 1:
                    this_taxa = most_occuring[0]
                elif len(taxons_filtered) == 1:
                    this_taxa = taxons_filtered[0]
                else:
                    if ref_stats:
                        # Grab most common
                        for reference in ref_stats:
                            if reference in taxons_filtered:
                                this_taxa = reference
                                break
                    else:
                        # No ref stats, grab first
                        this_taxa = most_occuring[0]

            if ignore_overlap_chunks:
                base_header = "|".join(
                    [this_gene, this_taxa, this_taxa_id, "contig_sequence"]
                )
                final_header = base_header
            else:
                base_header = "|".join([this_gene, this_taxa, this_taxa_id, node])

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

            # Create a hashmap to store the split quality taxons
            quality_taxons_here = {}

            for cursor in range(data_start, data_end + 1):
                sequences_at_current_point = current_point_seqs[cursor]
                amount_of_seqs_at_cursor = len(sequences_at_current_point)

                if amount_of_seqs_at_cursor == 0:
                    new_merge.append("-")  # Gap between sequences
                elif amount_of_seqs_at_cursor == 1:
                    new_merge.append(sequences_at_current_point[0][1][cursor])
                elif amount_of_seqs_at_cursor > 1:
                    # If there is more than one sequence at this current index
                    splits = amount_of_seqs_at_cursor - 1

                    if protein == "aa":
                        # Grab most occurring taxon
                        if comparison_sequences:
                            taxons_of_split = [
                                get_ref(header)
                                for (header, _) in sequences_at_current_point
                                if get_ref(header) in comparison_sequences
                            ]

                        comparison_taxa = None
                        most_occuring = most_common_element_with_count(taxons_of_split)
                        if len(taxons_of_split) > 1 and most_occuring[1] > 1:
                            comparison_taxa = most_occuring[0]
                        elif len(taxons_of_split) == 1:
                            comparison_taxa = taxons_of_split[0]
                        else:
                            if ref_stats:
                                for reference in ref_stats:
                                    if not reference in comparison_sequences:
                                        continue
                                    if taxons_of_split:
                                        if reference in taxons_of_split:
                                            comparison_taxa = reference
                                            break
                                    else:
                                        comparison_taxa = reference
                                        break
                            else:
                                if taxons_of_split:
                                    comparison_taxa = most_occuring[0] if most_occuring[1] != -1 else taxons_of_split[0]
                                else:
                                    headers_here = "".join([header for (header, _) in sequences_at_current_point])
                                    key = f"{data_start}{data_end}{headers_here}"
                                    if key in quality_taxons_here:
                                        comparison_taxa = quality_taxons_here[key]
                                    else:
                                        quality_taxons = []
                                        for taxon, sequence in comparison_sequences.items():
                                            data_region = sequence[data_start: data_end+1]
                                            data = len(data_region) - data_region.count("-")
                                            quality_taxons.append((data, taxon))
                                    
                                        comparison_taxa = max(quality_taxons, key = lambda x: x[0])[1]
                                        quality_taxons_here[key] = comparison_taxa
                            
                        # Grab the reference sequence for the mode taxon
                        if comparison_taxa:
                            comparison_sequence = comparison_sequences[comparison_taxa]

                    next_character = sequences_at_current_point[0][1][cursor]

                    # Iterate over each split if cursor is past calculated split
                    # position add from sequence B. We only want to add from one
                    # sequence out of every possible split, so we calculate which
                    # sequence to add from here then add the character to the
                    # final merge in the next line.
                    for split_count in range(splits - 1, -1, -1):
                        header_a, sequence_a = sequences_at_current_point[split_count]
                        header_b, sequence_b = sequences_at_current_point[
                            split_count + 1
                        ]

                        split_key = header_a + header_b

                        if protein == "aa":
                            if split_key in already_calculated_splits:
                                split_position = already_calculated_splits[split_key]
                            else:
                                split_position = calculate_split(
                                    sequence_a, sequence_b, comparison_sequence
                                )
                                already_calculated_splits[split_key] = split_position

                        elif protein == "nt":

                            split_position = already_calculated_splits[split_key] * 3

                        if cursor >= split_position:
                            next_character = sequence_b[cursor]
                            break

                    new_merge.append(next_character)

            new_merge.extend(["-"] * (trailing_end - data_end - 1))

            # Doing MR Outlier Check
            #
            # Each position of the new merge is checked with each position of overlapping
            # candidates at that point. If majority (args.majority percent) of characters
            # are the same but differ from the character chosen in the merge it will replace
            # it with the majority rules letter.

            if protein == "aa":
                if debug:
                    majority_assignments = ["-"] * len(new_merge)  # Used for debug log

                start_ends = {}
                for i in range(data_start, data_end + 1):
                    candidate_characters = {}
                    total_characters = 0
                    mode = -999
                    mode_char = None

                    char = new_merge[i]

                    for header, sequence, count in consists_of:
                        if header not in start_ends:
                            start_ends[header] = get_start_end(sequence)
                        start, end = start_ends[header]

                        # If sequence has data at this position
                        if start <= i <= end:
                            total_characters += count
                            this_char = sequence[i]

                            char_total = candidate_characters.get(this_char, 0) + count
                            candidate_characters[this_char] = char_total
                            if char_total > mode:
                                mode_char = this_char
                                mode = char_total

                    if (
                        total_characters >= minimum_mr_amount
                        and mode / total_characters >= majority
                        and mode_char != char
                    ):
                        new_merge[i] = mode_char
                        if debug:
                            majority_assignments[i] = "X" if mode_char == "-" else mode_char

            elif protein == "nt":
                length = len(new_merge)
                new_merge = ["".join(new_merge[i : i + 3]) for i in range(0, length, 3)]

                start_ends = {}

                triplet_data_start = int(data_start / 3)
                triplet_data_end = int(data_end / 3)

                if debug:
                    majority_assignments = ["---"] * len(
                        new_merge
                    )  # Used for debug log

                for raw_i in range(triplet_data_start, triplet_data_end + 1):
                    i = raw_i * 3
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

                            candidate_characters.extend([this_char] * count)

                    if total_characters >= minimum_mr_amount:
                        # Translate all NT triplets into AA
                        translated_characters = [
                            DNA_CODONS[j] for j in candidate_characters
                        ]

                        # Figure out what AA occurs the most
                        (
                            most_occuring_translated_char,
                            translated_char_count,
                        ) = most_common_element_with_count(translated_characters)

                        # Check to see if its majority
                        perc_appearing = translated_char_count / total_characters

                        if perc_appearing >= majority:
                            # Grab all the NT that map to that AA that triggered majority
                            candidate_chars_mapping_to_same_dna = [
                                j
                                for j in candidate_characters
                                if DNA_CODONS[j] == most_occuring_translated_char
                            ]

                            # Grab the most occurring NT that maps to that AA from the current seq
                            (
                                mode_cand_raw_character,
                                _,
                            ) = most_common_element_with_count(
                                candidate_chars_mapping_to_same_dna
                            )
                            if mode_cand_raw_character != char:
                                new_merge[raw_i] = mode_cand_raw_character
                                if debug:
                                    majority_assignments[
                                        raw_i
                                    ] = "XXX" if mode_cand_raw_character == "---" else mode_cand_raw_character

            new_merge = str(Seq("").join(new_merge))
            gene_out.append((final_header, new_merge))
            if debug:
                # If debug enabled add each component under final merge
                gene_out.append(
                    (
                        f"{this_taxa_id}|MajorityRulesAssigned",
                        "".join(majority_assignments),
                    )
                )
                for header, sequence, count in consists_of:
                    if "NODE" in header.split("|")[-1]:
                        node = header.split("|")[-1]
                        frame = ""
                    else:
                        node, frame = header.split("|")[-2:]
                    gene_out.append((f"{node}|{frame}|{count}", sequence))

    if protein == "nt" and gene_out:  # Remove empty columns
        to_keep = set()
        for i in range(len(gene_out[-1][1])):
            keep_this = False
            for record in gene_out:  # Every second element will be a sequence
                if record[1][i] != "-":
                    keep_this = True
                    break
            if keep_this:
                to_keep.add(i)

        for i, record in enumerate(gene_out):
            gene_out[i] = (
                record[0],
                "".join([let for i, let in enumerate(record[1]) if i in to_keep]),
            )

    output_path = os.path.join(output_dir, f"{protein}_merged", gene.rstrip(".gz"))

    return output_path, gene_out


def do_gene(
    gene,
    output_dir: Path,
    aa_path,
    nt_path,  # this one
    prep_dupe_counts,
    rep_dupe_counts,
    ref_stats,
    debug,
    majority,
    minimum_mr_amount,
    verbosity,
    ignore_overlap_chunks,
    compress,
) -> None:
    """
    Merge main loop. Opens fasta file, parses sequences and merges based on taxa
    """
    already_calculated_splits = {}
    printv(f"Doing: {gene}", verbosity, 2)

    aa_path, aa_data = do_protein(
        "aa",
        aa_path,
        output_dir,
        prep_dupe_counts,
        rep_dupe_counts,
        ref_stats,
        already_calculated_splits,
        gene,
        majority,
        minimum_mr_amount,
        ignore_overlap_chunks,
        debug=debug,
    )

    # print(already_calculated_splits)

    nt_path, nt_data = do_protein(
        "nt",
        nt_path,
        output_dir,
        prep_dupe_counts,
        rep_dupe_counts,
        ref_stats,
        already_calculated_splits,
        make_nt_name(gene),
        majority,
        minimum_mr_amount,
        ignore_overlap_chunks,
        debug=debug,
    )
    if nt_data:
        writeFasta(aa_path, aa_data, compress)
        writeFasta(nt_path, nt_data, compress)


def run_command(arg_tuple: tuple) -> None:
    """
    Calls the do_gene() function parallel in each thread
    """
    do_gene(*arg_tuple)


def do_folder(folder: Path, args):
    folder_time = TimeKeeper(KeeperMode.DIRECT)

    printv(f"Processing: {os.path.basename(folder)}", args.verbose, 0)
    input_path = Path(folder)
    aa_input = Path(input_path, args.aa_input)
    nt_input = Path(input_path, args.nt_input)

    if not os.path.exists(aa_input):
        print(aa_input)
        printv(f"WARNING: Can't find aa folder for taxa, {folder}", args.verbose, 0)
        return

    tmp_dir = directory_check(folder)
    dupe_tmp_file = Path(tmp_dir, "DupeSeqs.tmp")
    rocks_db_path = Path(folder, "rocksdb", "sequences", "nt")
    if rocks_db_path.exists():
        rocksdb_db = wrap_rocks.RocksDB(str(rocks_db_path))
        prepare_dupe_counts = json.loads(rocksdb_db.get("getall:gene_dupes"))
        reporter_dupe_counts = json.loads(rocksdb_db.get("getall:reporter_dupes"))
        ref_stats = rocksdb_db.get("getall:valid_refs").split(',')
    else:
        prepare_dupe_counts = {}
        reporter_dupe_counts = {}
        ref_stats = []

    target_genes = []
    for item in Path(aa_input).iterdir():
        if item.suffix in [".fa", ".gz", ".fq", ".fastq", ".fasta"]:
            target_genes.append(item.name)

    target_genes.sort(key=lambda x: Path(aa_input, x).stat().st_size, reverse=True)

    if args.processes > 1:
        arguments = []
        for target_gene in target_genes:
            prep_dupes_in_this_gene = prepare_dupe_counts.get(
                target_gene.split(".")[0], {}
            )
            rep_dupes_in_this_gene = reporter_dupe_counts.get(
                target_gene.split(".")[0], {}
            )
            target_aa_path = Path(aa_input, target_gene)
            target_nt_path = Path(nt_input, make_nt_name(target_gene))
            arguments.append(
                (
                    target_gene,
                    folder,
                    target_aa_path,
                    target_nt_path,
                    prep_dupes_in_this_gene,
                    rep_dupes_in_this_gene,
                    ref_stats,
                    args.debug,
                    args.majority,
                    args.majority_count,
                    args.verbose,
                    args.ignore_overlap_chunks,
                    args.compress,
                )
            )
        with Pool(args.processes) as pool:
            pool.map(run_command, arguments, chunksize=1)
    else:
        for target_gene in target_genes:
            prep_dupes_in_this_gene = prepare_dupe_counts.get(
                target_gene.split(".")[0], {}
            )
            rep_dupes_in_this_gene = reporter_dupe_counts.get(
                target_gene.split(".")[0], {}
            )
            target_aa_path = os.path.join(aa_input, target_gene)
            target_nt_path = os.path.join(nt_input, make_nt_name(target_gene))
            do_gene(
                target_gene,
                folder,
                target_aa_path,
                target_nt_path,
                prep_dupes_in_this_gene,
                rep_dupes_in_this_gene,
                ref_stats,
                args.debug,
                args.majority,
                args.majority_count,
                args.verbose,
                args.ignore_overlap_chunks,
                args.compress,
            )
    printv(f"Done! Took {folder_time.differential():.2f}s", args.verbose)

    if os.path.exists(dupe_tmp_file):
        os.remove(dupe_tmp_file)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for folder in args.INPUT:
        do_folder(Path(folder), args)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nsapphyre MergeOverlap"
    )