import hashlib
import os
import sqlite3
import wrap_rocks
from Bio.Seq import Seq
from time import time
import argparse
import math
from multiprocessing.pool import Pool
import shutil
import uuid
import json

T_global_start = time()


class NodeRecord:
    def __init__(self, header):
        self.header = header
        self.score = -math.inf
        self.cdna_start = None
        self.cdna_end = None
        self.cdna_sequence = ""
        self.aa_start = None
        self.aa_end = None
        self.aa_sequence = ""

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.score == other.score

    def __gt__(self, other):
        return self.score > other.score

    def __ge__(self, other):
        return self.score >= other.score

    def __le__(self, other):
        return self.score <= other.score


def clear_output_path(input):
    """
    Clears protein output paths
    """

    for protein_path in ["aa", "nt"]:
        protein_output = os.path.join(input, protein_path)
        for item in os.listdir(protein_output):
            item_path = os.path.join(protein_output, item)
            if ".fa" in item:
                os.remove(item_path)


def get_set_id(orthoset_db_con, orthoset):
    """
    Retrieves orthoset id from orthoset db
    """
    orthoset_id = None

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute("SELECT * FROM orthograph_set_details;")

    for row in rows:
        id, name, description = row

        if name == orthoset:
            orthoset_id = id

    if orthoset_id == None:
        raise Exception("Orthoset {} id cant be retrieved".format(orthoset))
    else:
        return orthoset_id


def get_taxa_in_set(set_id, orthoset_db_con):
    reference_taxa = []

    query = f'''SELECT DISTINCT {orthoset_set_details}.name, {orthoset_taxa}.name
        FROM {orthoset_seqpairs}
        INNER JOIN {orthoset_taxa}
            ON {orthoset_seqpairs}.taxid = {orthoset_taxa}.id
        INNER JOIN {orthoset_orthologs}
            ON {orthoset_orthologs}.sequence_pair = {orthoset_seqpairs}.id 
        INNER JOIN {orthoset_set_details}
            ON {orthoset_orthologs}.setid = {orthoset_set_details}.id
        WHERE {orthoset_set_details}.id = "{set_id}"'''

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute(query)

    reference_taxa = [row[1] for row in rows]

    return reference_taxa

def get_scores_list(score_threshold, min_length, num_threads):
    batches = rocksdb_db.get("hmmbatch:all")
    batches = batches.split(',')

    score_based_results = {}
    ufr_out = [["Gene", "Hash", "Header", "Score", "Start", "End"]]

    for batch_i in batches:
        batch_rows = rocksdb_db.get(f"hmmbatch:{batch_i}")
        batch_rows = json.loads(batch_rows)
        for this_row in batch_rows:
            orthoid = this_row["gene"]

            header = this_row["header"].strip()
            this_row["header"] = header

            env_start = this_row["env_start"]
            env_end = this_row["env_end"]

            this_row["hmmhit"] = ""

            length = env_end - env_start

            hmm_score = this_row["score"]
            hmm_start = this_row["hmm_start"]
            hmm_end = this_row["hmm_end"]

            if float(hmm_score) > score_threshold:
                if length > min_length:

                    components = header.split("|")
                    if "[" in components[-1]:
                        out_header = "|".join(components[:-1]) + " " + components[-1]
                    else:
                        out_header = header

                    ufr_out.append(
                        [
                            orthoid,
                            out_header,
                            str(hmm_score),
                            str(hmm_start),
                            str(hmm_end),
                        ]
                    )

                    if hmm_score not in score_based_results:
                        score_based_results[hmm_score] = []
                    score_based_results[hmm_score].append(this_row)


    ufr_out = sorted(ufr_out, key=lambda x: (x[0], x[1], x[3], x[4]))

    return score_based_results, ufr_out


def get_blastresults_for_hmmsearch_id(hmmsearch_id):
    key = "blastfor:{}".format(hmmsearch_id)

    db_entry = rocksdb_db.get(key)

    # TODO: CHECK IF ITS POSSIBLE FOR NO BLAST RESULTS OR IF THIS IS A BLAST BUG - THIS IS POSSIBLE
    if db_entry == None:
        return []
    else:
        result = json.loads(db_entry)
        return result


def is_reciprocal_match(
    blast_results, reference_taxa
):
    reftaxon_count = {ref_taxa: 0 for ref_taxa in reference_taxa}

    for match in blast_results:
        reftaxon = match["reftaxon"]

        if reftaxon in reftaxon_count:
            reftaxon_count[reftaxon] = 1

            if strict_search_mode:
                total_count = sum(
                    [reftaxon_count[reftaxon] for reftaxon in reference_taxa]
                )
                if total_count == len(reference_taxa):
                    return match
            else:
                return match
    return None


def calculate_length_of_ali(result):
    return abs(result["ali_end"] - result["ali_start"]) + 1


def transcript_not_long_enough(result, minimum_transcript_length):
    length = calculate_length_of_ali(result)

    if length < minimum_transcript_length:
        return True
    else:
        return False


def coordinate_overlap(other_hits, new_hit):
    for hit in other_hits:
        if "ali_start" in hit:
            starts_before_ends_within = (
                new_hit["ali_start"] < hit["ali_start"]
                and new_hit["ali_end"] > hit["ali_start"]
            )
            starts_before_ends_after = (
                new_hit["ali_start"] < hit["ali_start"]
                and new_hit["ali_end"] > hit["ali_end"]
            )
            starts_within_ends_within = (
                new_hit["ali_start"] > hit["ali_start"]
                and new_hit["ali_end"] < hit["ali_end"]
            )

            if (
                starts_before_ends_within
                or starts_before_ends_after
                or starts_within_ends_within
            ):
                return True
    return False


def get_reference_sequence(hit_id, orthoset_db_con):
    orthoset_db_cur = orthoset_db_con.cursor()

    query = f'''SELECT
            {orthoset_aaseqs}.{db_col_sequence}, 
            {orthoset_taxa}.{db_col_name}
        FROM {orthoset_aaseqs}
        INNER JOIN {orthoset_taxa}
            ON {orthoset_aaseqs}.{db_col_taxid} = {orthoset_taxa}.id
        WHERE {orthoset_aaseqs}.{db_col_id} = "{hit_id}"'''

    rows = orthoset_db_cur.execute(query)

    for row in rows:
        return (row[0], row[1])


def reverse_complement(nt_seq):
    seq = Seq(nt_seq)
    return str(seq.reverse_complement())


def get_nucleotide_transcript_for(header):
    base_header = get_baseheader(header).strip()
    hash_of_header = hashlib.sha256(base_header.encode()).hexdigest()

    row_data = rocksdb_db.get(hash_of_header)
    this_header, sequence = row_data.split("\n")

    if "revcomp" in header:
        return base_header, reverse_complement(sequence)
    else:
        return base_header, sequence


def crop_to_hmm_alignment(seq, header, hit):
    start = hit["ali_start"] - 1  # Adjust for zero based number
    start = start * 3
    end = hit["ali_end"] * 3

    return seq[start:end]


### EXONERATE FUNCTIONS


def parse_results(result_file_content):

    cdna_sequence = ""
    aa_sequence = ""
    cdna_start = None
    cdna_end = None
    aa_start = None
    aa_end = None

    lines = result_file_content.split("\n")
    current_seq = ""
    for line in lines:
        if ">" in line:
            header, start, end = line.split(" ")
            current_seq = header.replace(">", "")
            if current_seq == "cdna":
                cdna_start = start + 1
                cdna_end = end
            elif current_seq == "aa":
                aa_start = start + 1
                aa_end = end
        else:
            if current_seq == "cdna":
                cdna_sequence += line
            elif current_seq == "aa":
                aa_sequence += line

    if cdna_sequence == "":
        cdna_sequence = None
    if aa_sequence == "":
        aa_sequence = None

    return (cdna_sequence, cdna_start, cdna_end, aa_sequence, aa_start, aa_end)


def fastaify(headers, sequences, input_path, tmp_path):
    name = str(uuid.uuid4()) + ".tmp"  # Super unique name
    path = os.path.join(tmp_path, name)

    this_out = []
    for i, seq in enumerate(sequences):
        this_out.append(">" + headers[i])
        this_out.append(seq)

    open(path, "w").write("\n".join(this_out))

    return path


def translate_cdna(cdna_seq, nt_table):
    if cdna_seq == None:
        return None

    if len(cdna_seq) % 3 != 0:
        print("WARNING: NT Sequence length is not divisable by 3")

    split_seq = [cdna_seq[x : x + 3] for x in range(0, len(cdna_seq), 3)]
    translation = "".join([nt_table[i] for i in split_seq])

    return translation


def parse_nodes(lines):
    """
    vulgar => cdna => aa => vulgar
    0 => 1 => 2 => 0
    """
    nodes = dict()
    this_node = None
    state = 0
    for line in lines:
        if "vulgar:" in line:  # if vulgar line, start new record
            current_header = line.split(" . ")[1].split(" ")[0]
            this_node = NodeRecord(current_header)
            this_node.score = int(line.split()[9])
            previously_seen = nodes.get(current_header, False)
            if previously_seen:
                nodes[current_header] = max(this_node, previously_seen)
            else:
                nodes[current_header] = this_node
        elif ">cdna" in line:  # if cdna, set cdna values
            fields = line.split()
            this_node.cdna_start = int(fields[1])
            this_node.cdna_end = int(fields[2])
            state = 1
        elif ">aa" in line:  # if aa in line, set aa values
            fields = line.split()
            this_node.aa_start = int(fields[1])
            this_node.aa_end = int(fields[2])
            state = 2
        elif state == 1:  # if cdna data line
            this_node.cdna_sequence += line
        elif state == 2:  # if aa data line
            this_node.aa_sequence += line
    return nodes


def parse_multi_results(result_file_content):
    result = []
    lines = result_file_content.split("\n")
    nodes = parse_nodes(lines)
    for key in nodes:
        node = nodes[key]
        result.append(
            (
                node.header,
                node.cdna_sequence,
                node.cdna_start,
                node.cdna_end,
                node.aa_sequence,
                node.aa_start,
                node.aa_end,
            )
        )
    return result


def get_multi_orf(query, targets, input_path, score_threshold, tmp_path, extend=False):
    exonerate_ryo = ">cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas"
    genetic_code = 1
    exonerate_model = "protein2genome"
    exhaustive = ""

    headers = [i["header"].strip().replace(" ", "|") for i in targets]
    if extend:
        sequences = [i["est_sequence_complete"] for i in targets]
    else:
        sequences = [i["est_sequence_hmm_region"] for i in targets]

    queryfile = fastaify(["query"], [query], input_path, tmp_path)
    targetfile = fastaify(headers, sequences, input_path, tmp_path)

    outfile = os.path.join(tmp_path, str(uuid.uuid4()) + ".exonerateout")

    exonerate_cmd = f"exonerate --score {score_threshold} --ryo '{exonerate_ryo}' --subopt 0 --geneticcode {genetic_code} --model '{exonerate_model}' --querytype 'protein' --targettype 'dna' --verbose 0 --showalignment 'no' --showvulgar 'yes' --query '{queryfile}' --target '{targetfile}' > {outfile}"

    os.system(exonerate_cmd)

    result_content = open(outfile).read()

    results = parse_multi_results(result_content)

    os.remove(queryfile)
    os.remove(targetfile)
    os.remove(outfile)

    return results


def extended_orf_contains_original_orf(hit):
    if (
        hit["extended_orf_cdna_start"] > hit["orf_cdna_end_on_transcript"]
        or hit["extended_orf_cdna_end"] < hit["orf_cdna_start_on_transcript"]
    ):
        return False
    else:
        return True


def get_overlap_length(candidate):
    if (
        candidate["extended_orf_aa_start_on_transcript"] <= candidate["ali_end"]
        and candidate["extended_orf_aa_end_on_transcript"] >= candidate["ali_start"]
    ):
        overlap_start = (
            candidate["extended_orf_aa_start_on_transcript"]
            if candidate["extended_orf_aa_start_on_transcript"] > candidate["ali_start"]
            else candidate["ali_start"]
        )
        overlap_end = (
            candidate["extended_orf_aa_end_on_transcript"]
            if candidate["extended_orf_aa_end_on_transcript"] > candidate["ali_end"]
            else candidate["ali_end"]
        )

        overlap_length = overlap_end - overlap_start
        return overlap_length
    else:
        return 0


def overlap_by_orf(candidate):
    orf_length = abs(
        candidate["extended_orf_aa_end"] - candidate["extended_orf_aa_start"]
    )
    overlap_length = get_overlap_length(candidate)

    return overlap_length / orf_length


def remove_extended_orf(hit):
    hit.pop("extended_orf_aa_sequence")
    hit.pop("extended_orf_aa_start")
    hit.pop("extended_orf_aa_end")
    hit.pop("extended_orf_cdna_sequence")
    hit.pop("extended_orf_cdna_start")
    hit.pop("extended_orf_cdna_end")

    return hit


### FIN


def get_ortholog_group(orthoset_id, orthoid, orthoset_db_con):
    query = f"""SELECT 
            {orthoset_taxa}.{db_col_name},
            {orthoset_aaseqs}.{db_col_header},
            {orthoset_aaseqs}.{db_col_sequence}
        FROM {orthoset_aaseqs}
        INNER JOIN {orthoset_seqpairs}
            ON {orthoset_aaseqs}.{db_col_id} = {orthoset_seqpairs}.{db_col_aaseq}
        INNER JOIN {orthoset_orthologs}
            ON {orthoset_seqpairs}.{db_col_id} = {orthoset_orthologs}.{db_col_seqpair}
        INNER JOIN {orthoset_taxa}
            ON {orthoset_aaseqs}.{db_col_taxid} = {orthoset_taxa}.{db_col_id}
        AND   {orthoset_orthologs}.{db_col_setid} = ?
        AND   {orthoset_orthologs}.{db_col_orthoid} = ?
        ORDER BY {orthoset_taxa}.{db_col_name}"""

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute(
        query,
        (
            orthoset_id,
            orthoid,
        ),
    )

    return rows

def format_candidate_header(gene, taxa_name, taxa_id, sequence_id, coords, frame):
    return header_seperator.join([gene, taxa_name, taxa_id, sequence_id, coords, frame])

def format_reference_header(gene, taxa_name, taxa_id, identifier = '.'):
    return header_seperator.join([gene, taxa_name, taxa_id, identifier])

def print_core_sequences(orthoid, core_sequences):
    core_sequences = sorted(core_sequences)

    result = []
    for core in core_sequences:
        header = format_reference_header(orthoid, core[0], core[1])

        sequence = core[2]

        result.append(">" + header + "\n")
        result.append(sequence + "\n")

    return result


def get_rf(header):
    raw_header = get_baseheader(header)
    header_part = get_translate(header)
    if "translate" in header_part or "revcomp" in header_part:
        frame = header_part

    return raw_header, frame


def print_unmerged_sequences(
    hits, orthoid, type, minimum_seq_data_length, species_name, kicks
):
    result = []
    kicks_result = set()
    for i, hit in enumerate(hits):
        this_hdr, rf = get_rf(hit["header"])

        start = (
            hit["extended_orf_aa_start_on_transcript"]
            if "extended_orf_aa_start_on_transcript" in hit
            else hit["orf_aa_start_on_transcript"]
        )
        end = (
            hit["extended_orf_aa_end_on_transcript"]
            if "extended_orf_aa_end_on_transcript" in hit
            else hit["orf_aa_end_on_transcript"]
        )

        header = format_candidate_header(orthoid, hit["reftaxon"], species_name, this_hdr, f'{round(start)}-{round(end)}', rf)

        if type == "nt":
            seq = (
                hit["extended_orf_cdna_sequence"]
                if "extended_orf_cdna_sequence" in hit
                else hit["orf_cdna_sequence"]
            )
        elif type == "aa":
            seq = (
                hit["extended_orf_aa_sequence"]
                if "extended_orf_aa_sequence" in hit
                else hit["orf_aa_sequence"]
            )

        if type == "aa":
            if len(seq) - seq.count("-") > minimum_seq_data_length:
                result.append(">" + header + "\n")
                result.append(seq + "\n")
            else:
                kicks_result.add(i)
        elif type == "nt":
            if i not in kicks:
                result.append(">" + header + "\n")
                result.append(seq + "\n")

    return kicks_result, result


def get_ortholog_group_nucleotide(orthoset_id, orthoid, orthoset_db_con):
    query = f"""SELECT 
            {orthoset_taxa}.{db_col_name},
            {orthoset_ntseqs}.{db_col_header},
            {orthoset_ntseqs}.{db_col_sequence}
        FROM {orthoset_ntseqs}
        INNER JOIN {orthoset_seqpairs}
            ON {orthoset_ntseqs}.{db_col_id} = {orthoset_seqpairs}.{db_col_ntseq}
        INNER JOIN {orthoset_orthologs}
            ON {orthoset_seqpairs}.{db_col_id} = {orthoset_orthologs}.{db_col_seqpair}
        INNER JOIN {orthoset_taxa}
            ON {orthoset_ntseqs}.{db_col_taxid} = {orthoset_taxa}.{db_col_id}
        AND   {orthoset_orthologs}.{db_col_setid} = ?
        AND   {orthoset_orthologs}.{db_col_orthoid} = ?
        ORDER BY {orthoset_taxa}.{db_col_name}"""

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute(
        query,
        (
            orthoset_id,
            orthoid,
        ),
    )

    return rows


def get_difference(scoreA, scoreB):
    """
    Returns decimal difference of two scores
    """

    try:
        if scoreA / scoreB > 1:
            return scoreA / scoreB

        elif scoreB / scoreA > 1:
            return scoreB / scoreA

        elif scoreA == scoreB:
            return 1

    except ZeroDivisionError:
        return 0


def get_baseheader(header):
    """
    Returns header content before first whitespace.
    """
    baseheader = header.split("|")[0].strip()
    return baseheader


def get_translate(header):
    """
    Returns header content before first whitespace.
    """
    translate = header.split("|")[1]
    return translate


def run_exonerate(arg_tuple):
    (
        gene,
        list_of_hits,
        orthoset_db_path,
        input_path,
        min_score,
        orthoset_id,
        aa_out_path,
        min_length,
        species_name,
        nt_out_path,
        tmp_path,
        exonerate_verbose,
    ) = arg_tuple
    exonerate_gene_multi(
        gene,
        list_of_hits,
        orthoset_db_path,
        input_path,
        min_score,
        orthoset_id,
        aa_out_path,
        min_length,
        species_name,
        nt_out_path,
        tmp_path,
        exonerate_verbose,
    )


def exonerate_gene_multi(
    orthoid,
    list_of_hits,
    orthoset_db_path,
    input_path,
    min_score,
    orthoset_id,
    aa_out_path,
    min_length,
    species_name,
    nt_out_path,
    tmp_path,
    exonerate_verbose,
):
    T_gene_start = time()

    summary_out_path = os.path.join(tmp_path, "{}-summary.txt".format(orthoid))
    summary_out_log = []

    nt_table = {
        "TTT": "F",
        "TCT": "S",
        "TAT": "Y",
        "TGT": "C",
        "TTC": "F",
        "TCC": "S",
        "TAC": "Y",
        "TGC": "C",
        "TTA": "L",
        "TCA": "S",
        "TAA": "*",
        "TGA": "*",
        "TTG": "L",
        "TCG": "S",
        "TAG": "*",
        "TGG": "W",
        "CTT": "L",
        "CCT": "P",
        "CAT": "H",
        "CGT": "R",
        "CTC": "L",
        "CCC": "P",
        "CAC": "H",
        "CGC": "R",
        "CTA": "L",
        "CCA": "P",
        "CAA": "Q",
        "CGA": "R",
        "CTG": "L",
        "CCG": "P",
        "CAG": "Q",
        "CGG": "R",
        "ATT": "I",
        "ACT": "T",
        "AAT": "N",
        "AGT": "S",
        "ATC": "I",
        "ACC": "T",
        "AAC": "N",
        "AGC": "S",
        "ATA": "I",
        "ACA": "T",
        "AAA": "K",
        "AGA": "R",
        "ATG": "M",
        "ACG": "T",
        "AAG": "K",
        "AGG": "R",
        "GTT": "V",
        "GCT": "A",
        "GAT": "D",
        "GGT": "G",
        "GTC": "V",
        "GCC": "A",
        "GAC": "D",
        "GGC": "G",
        "GTA": "V",
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",
        "GTG": "V",
        "GCG": "A",
        "GAG": "E",
        "GGG": "G",
    }

    orthoset_db_con = sqlite3.connect(orthoset_db_path)

    if exonerate_verbose:
        print("Exonerating and doing output for: ", orthoid)
    reftaxon_related_transcripts = {}
    reftaxon_to_proteome_sequence = {}
    for hit in list_of_hits:
        proteome_sequence, this_reftaxon = get_reference_sequence(
            hit["target"], orthoset_db_con
        )

        reftaxon_to_proteome_sequence[this_reftaxon] = proteome_sequence

        hit["reftaxon"] = this_reftaxon
        hit["proteome_sequence"] = proteome_sequence

        est_header, est_sequence_complete = get_nucleotide_transcript_for(hit["header"])
        est_sequence_hmm_region = crop_to_hmm_alignment(
            est_sequence_complete, est_header, hit
        )

        hit["est_header"] = est_header
        hit["est_sequence_complete"] = est_sequence_complete
        hit["est_sequence_hmm_region"] = est_sequence_hmm_region

        if this_reftaxon not in reftaxon_related_transcripts:
            reftaxon_related_transcripts[this_reftaxon] = []

        reftaxon_related_transcripts[this_reftaxon].append(hit)

    output_sequences = []

    total_results = 0

    for taxon_hit in reftaxon_related_transcripts:
        hits = reftaxon_related_transcripts[taxon_hit]
        query = reftaxon_to_proteome_sequence[taxon_hit]

        results = get_multi_orf(query, hits, input_path, min_score, tmp_path)

        total_results += len(hits)

        if extend_orf:
            extended_results = get_multi_orf(
                query, hits, input_path, min_score, tmp_path, extend=True
            )

        for hit in hits:
            matching_alignment = [
                i for i in results if i[0].strip() == hit["header"].strip()
            ]

            if len(matching_alignment) == 0:
                (
                    orf_aa_sequence,
                    orf_cdna_sequence,
                    orf_cdna_start,
                    orf_cdna_end,
                    orf_aa_start,
                    orf_aa_end,
                ) = (None, None, None, None, None, None)
            else:
                (
                    current_header,
                    orf_cdna_sequence,
                    orf_cdna_start,
                    orf_cdna_end,
                    orf_aa_sequence,
                    orf_aa_start,
                    orf_aa_end,
                ) = matching_alignment[0]

                orf_aa_sequence = translate_cdna(orf_cdna_sequence, nt_table)

                hit["orf_aa_sequence"] = orf_aa_sequence
                hit["orf_cdna_sequence"] = orf_cdna_sequence
                hit["orf_cdna_start"] = orf_cdna_start
                hit["orf_cdna_end"] = orf_cdna_end
                hit["orf_aa_start"] = orf_aa_start
                hit["orf_aa_end"] = orf_aa_end

                if orf_cdna_sequence == None:
                    continue

                else:
                    hit["orf_cdna_start_on_transcript"] = (
                        hit["orf_cdna_start"] + (hit["ali_start"] * 3) - 3
                    )
                    hit["orf_cdna_end_on_transcript"] = (
                        hit["orf_cdna_end"] + (hit["ali_start"] * 3) - 3
                    )
                    hit["orf_aa_start_on_transcript"] = (
                        hit["orf_cdna_start"] + (hit["ali_start"] * 3) - 3
                    ) / 3
                    hit["orf_aa_end_on_transcript"] = (
                        hit["orf_cdna_end"] + (hit["ali_start"] * 3) - 3
                    ) / 3

                    if extend_orf:
                        matching_extended_alignment = [
                            i for i in extended_results if i[0] == hit["header"]
                        ]
                        if matching_extended_alignment != []:
                            (
                                current_header,
                                extended_orf_cdna_sequence,
                                extended_orf_cdna_start,
                                extended_orf_cdna_end,
                                extended_orf_aa_sequence,
                                extended_orf_aa_start,
                                extended_orf_aa_end,
                            ) = matching_extended_alignment[0]

                            extended_orf_aa_start_on_transcript = (
                                extended_orf_cdna_start - 1
                            ) / 3 + 1
                            # assert extended_orf_aa_start_on_transcript.is_integer()
                            extended_orf_aa_end_on_transcript = (
                                extended_orf_cdna_end - 1
                            ) / 3 + 1
                            # assert extended_orf_aa_end_on_transcript.is_integer()
                            extended_orf_aa_sequence = translate_cdna(
                                extended_orf_cdna_sequence, nt_table
                            )

                            hit["extended_orf_aa_sequence"] = extended_orf_aa_sequence
                            hit[
                                "extended_orf_cdna_sequence"
                            ] = extended_orf_cdna_sequence
                            hit["extended_orf_cdna_start"] = extended_orf_cdna_start
                            hit["extended_orf_cdna_end"] = extended_orf_cdna_end
                            hit["extended_orf_aa_start"] = extended_orf_aa_start
                            hit["extended_orf_aa_end"] = extended_orf_aa_end
                            hit["extended_orf_aa_start"] = extended_orf_aa_start
                            hit["extended_orf_aa_end"] = extended_orf_aa_end

                            hit[
                                "extended_orf_aa_start_on_transcript"
                            ] = extended_orf_aa_start_on_transcript
                            hit[
                                "extended_orf_aa_end_on_transcript"
                            ] = extended_orf_aa_end_on_transcript

                            if extended_orf_contains_original_orf(hit):
                                orf_overlap = overlap_by_orf(hit)
                                orf_overlap = int(orf_overlap * 1000000)
                                if orf_overlap < orf_overlap_minimum:
                                    hit = remove_extended_orf(hit)

                            else:
                                hit = remove_extended_orf(hit)
                        else:
                            print("Failed to extend orf on {}".format(hit["header"]))

                    output_sequences.append(hit)

    if len(output_sequences) > 0:
        output_sequences = sorted(output_sequences, key=lambda d: d["hmm_start"])

        core_sequences = get_ortholog_group(orthoset_id, orthoid, orthoset_db_con)

        this_aa_out = []
        this_aa_path = os.path.join(aa_out_path, orthoid + ".aa.fa")
        this_aa_out.extend(print_core_sequences(orthoid, core_sequences))
        kicks, output = print_unmerged_sequences(
            output_sequences, orthoid, "aa", min_length, species_name, kicks=set()
        )
        
        if output:
            this_aa_out.extend(output)
            open(this_aa_path, "w").writelines(this_aa_out)

            this_out = [orthoid]
            for hit in output_sequences:
                length = (
                    len(hit["extended_orf_aa_sequence"])
                    if "extended_orf_aa_sequence" in hit
                    else len(hit["orf_aa_sequence"])
                )
                this_out.append("{}[{} aa]".format(hit["header"], length))

            summary_out_log.append("\t".join(this_out))

            core_sequences_nt = get_ortholog_group_nucleotide(
                orthoset_id, orthoid, orthoset_db_con
            )

            this_nt_out = []
            this_nt_path = os.path.join(nt_out_path, orthoid + ".nt.fa")

            this_nt_out.extend(print_core_sequences(orthoid, core_sequences_nt))
            na_kicks, output = print_unmerged_sequences(
                output_sequences, orthoid, "nt", min_length, species_name, kicks=kicks
            )
            this_nt_out.extend(output)

            orthoid_summary_out = []
            for hit in output_sequences:
                if "orf_aa_sequence" in hit:
                    if hit["orf_aa_sequence"] != None:
                        header = hit["header"]
                        length = (
                            len(hit["extended_orf_aa_sequence"])
                            if "extended_orf_aa_sequence" in hit
                            else len(hit["orf_aa_sequence"])
                        )
                        orthoid_summary_out.append((header, length))

            open(this_nt_path, "w").writelines(this_nt_out)

    open(summary_out_path, "w").write("\n".join(summary_out_log))

    if exonerate_verbose:
        print(
            "{} took {:.2f}s. Had {} sequences".format(
                orthoid, time() - T_gene_start, len(output_sequences)
            )
        )

def reciprocal_search(
    hmmresults,
    list_of_wanted_orthoids,
    reference_taxa,
    score,
    reciprocal_verbose,
):
    if reciprocal_verbose:
        T_reciprocal_start = time()
        print("Ensuring reciprocal hit for hmmresults in {}".format(score))

    results = []
    this_fails = []
    for result in hmmresults:
        orthoid = result["gene"]

        if list_of_wanted_orthoids != []:
            if orthoid not in list_of_wanted_orthoids:
                continue

        result_hmmsearch_id = result["hmm_id"]

        blast_results = get_blastresults_for_hmmsearch_id(result_hmmsearch_id)

        this_match = is_reciprocal_match(
            blast_results, reference_taxa
        )

        if this_match == None:
            this_fails.append(
                [
                    result["gene"],
                    result["hmmhit"],
                    result["header"],
                    str(result["score"]),
                    str(result["hmm_start"]),
                    str(result["hmm_end"]),
                    "Reciprocal mismatch",
                ]
            )
        else:
            this_match.update(result)  # Persist hmmresult data
            this_match["gene"] = orthoid
            this_match["score"] = score
            results.append(this_match)

    if reciprocal_verbose:
        print(
            "Checked reciprocal hits for {}. Took {:.2f}s.".format(
                score, time() - T_reciprocal_start
            )
        )

    return {"Results": results, "Kicks": this_fails}


####
db_col_name = "name"
db_col_target = "target"
db_col_score = "score"
db_col_evalue = "evalue"
db_col_start = "start"
db_col_end = "end"
db_col_env_start = "env_start"
db_col_env_end = "env_end"
db_col_ali_start = "ali_start"
db_col_ali_end = "ali_end"
db_col_hmm_start = "hmm_start"
db_col_hmm_end = "hmm_end"
db_col_header = "header"
db_col_hmmsearch_id = "hmmsearch_id"
db_col_id = "id"
db_col_digest = "digest"
db_col_orthoid = "ortholog_gene_id"
db_col_taxid = "taxid"
db_col_query = "query"
db_col_setid = "setid"
db_col_sequence = "sequence"
db_col_aaseq = "aa_seq"
db_col_seqpair = "sequence_pair"
db_col_ntseq = "nt_seq"

orthoset_set_details = "orthograph_set_details"
orthoset_taxa = "orthograph_taxa"
orthoset_seqpairs = "orthograph_sequence_pairs"
orthoset_orthologs = "orthograph_orthologs"
orthoset_aaseqs = "orthograph_aaseqs"
orthoset_ntseqs = "orthograph_ntseqs"

####
# Misc Settings
####

strict_search_mode = False
orthoid_list_file = None
frameshift_correction = True
extend_orf = True
orf_overlap_minimum = 0.15
orf_overlap_minimum = int(orf_overlap_minimum * 1000000)
clear_output = False

header_seperator = "|"

####

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="PhyMMR/Acroceridae/SRR6453524.fa",
        help="Path to directory of Input folder",
    )
    parser.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="PhyMMR/orthosets",
        help="Path to directory of Orthosets folder",
    )
    parser.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    parser.add_argument(
        "-ml", "--min_length", type=int, default=30, help="Minimum Transcript Length"
    )
    parser.add_argument(
        "-ms", "--min_score", type=float, default=40, help="Minimum Hit Domain Score"
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Number of threads used to call processes.",
    )
    parser.add_argument("-v", "--verbose", type=int, default=2, help="Verbose debug.")
    parser.add_argument("-d", "--debug", type=int, default=0, help="Verbose debug.")

    args = parser.parse_args()

    debug = args.debug != 0

    num_threads = args.processes

    ####
    # Filter settings
    ####

    min_length = args.min_length
    min_score = args.min_score
    verbose = range(0, args.verbose + 1)

    ####

    input_path = args.input
    species_name = os.path.basename(input_path).split(".")[0]

    print("Doing {}.".format(species_name))

    if 1 in verbose:
        T_init_db = time()

    if os.path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif os.path.exists("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = os.path.join(input_path, "tmp")
        if not os.path.exists(tmp_path):
            os.mkdir(tmp_path)

    if orthoid_list_file:
        list_of_wanted_orthoids = open(orthoid_list_file).read().split("\n")
        wanted_orthoid_only = True
    else:
        list_of_wanted_orthoids = []
        wanted_orthoid_only = False

    orthoset_path = args.orthoset_input
    orthoset = args.orthoset

    orthoset_db_path = os.path.join(orthoset_path, orthoset + ".sqlite")
    orthoset_db_con = sqlite3.connect(orthoset_db_path)

    cache_size = 16000000

    orthoset_db_con.execute("PRAGMA journal_mode = OFF;")
    orthoset_db_con.execute("PRAGMA synchronous = 0;")
    orthoset_db_con.execute(f"PRAGMA cache_size = {cache_size};")
    orthoset_db_con.execute("PRAGMA locking_mode = EXCLUSIVE;")
    orthoset_db_con.execute("PRAGMA temp_store = MEMORY;")

    orthoset_id = get_set_id(orthoset_db_con, orthoset)

    if clear_output:
        clear_output_path(input_path)

    aa_out = "aa"
    nt_out = "nt"

    aa_out_path = os.path.join(input_path, aa_out)
    nt_out_path = os.path.join(input_path, nt_out)

    if os.path.exists(aa_out_path):
        shutil.rmtree(aa_out_path)
    os.mkdir(aa_out_path)

    if os.path.exists(nt_out_path):
        shutil.rmtree(nt_out_path)
    os.mkdir(nt_out_path)

    rocks_db_path = os.path.join(input_path, "rocksdb")

    if 2 in verbose:
        T_reference_taxa = time()
        print(
            "Initialized databases. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing reference taxa in set.".format(
                time() - T_global_start, time() - T_init_db
            )
        )

    reference_taxa = get_taxa_in_set(orthoset_id, orthoset_db_con)

    if 2 in verbose:
        T_hmmresults = time()
        print(
            "Got referenca taxa in set. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing hmmsearch hits".format(
                time() - T_global_start, time() - T_reference_taxa
            )
        )

    rocksdb_db = wrap_rocks.RocksDB(rocks_db_path)
    score_based_results, ufr_rows = get_scores_list(
        min_score, min_length, num_threads
    )
    # gene_based_results,header_based_results,ufr_rows = get_scores_list(min_score,taxa_db_path,orthoset_db_path,min_length,orthoset_id)

    ufr_path = os.path.join(input_path, "unfiltered-hits.csv")

    ufr_out = ["Gene,Header,Score,Start,End\n"]
    for row in ufr_rows:
        ufr_out.append(",".join(row) + "\n")

    if debug:
        open(ufr_path, "w").writelines(ufr_out)

    ####################################
    if 2 in verbose:
        print(
            "Got hmmresults. Elapsed time {:.2f}s. Took {:.2f}s.".format(
                time() - T_global_start, time() - T_hmmresults
            )
        )

    if 1 in verbose:
        T_reciprocal_search = time()
        if args.verbose != 1:
            print(
                "Retrieved data from DB. Elapsed time {:.2f}s. Took {:.2f}s. Doing reciprocal check.".format(
                    time() - T_global_start, time() - T_hmmresults
                )
            )

    scores = list(score_based_results.keys())
    scores.sort(reverse=True)  # Ascending
    transcripts_mapped_to = {}

    reciprocal_verbose = 4 in verbose
    brh_path = os.path.join(input_path, "best-reciprocal-hits.txt")
    filtered_sequences_log = []

    arguments = list()
    for score in scores:
        hmmresults = score_based_results[score]
        arguments.append(
            (
                hmmresults,
                list_of_wanted_orthoids,
                reference_taxa,
                score,
                reciprocal_verbose,
            )
        )
    with Pool(num_threads) as pool:
        reciprocal_data = pool.starmap(reciprocal_search, arguments, chunksize=1)

    brh_file_out = []

    for data in reciprocal_data:
        filtered_sequences_log.extend(data["Kicks"])
        for this_match in data["Results"]:
            orthoid = this_match["gene"]

            if transcript_not_long_enough(this_match, min_length):
                continue

            if orthoid not in transcripts_mapped_to:
                transcripts_mapped_to[orthoid] = []

            transcripts_mapped_to[orthoid].append(this_match)

            brh_file_out.append(
                "\t".join(
                    [
                        this_match["gene"],
                        this_match["header"],
                        str(this_match["ali_start"]),
                        str(this_match["ali_end"]),
                        str(this_match["score"]),
                        str(this_match["hmm_evalue"]),
                        str(this_match["hmm_start"]),
                        str(this_match["hmm_end"]),
                    ]
                )
                + "\n"
            )

    brh_file_out.sort()
    open(brh_path, "w").writelines(brh_file_out)

    if 2 in verbose:
        T_internal_search = time()
        print(
            "Reciprocal check done, found {} reciprocal hits. Elapsed time {:.2f}s. Took {:.2f}s. Exonerating genes.".format(
                len(brh_file_out), time() - T_global_start, time() - T_reciprocal_search
            )
        )

    # Clear summary files
    for item in os.listdir(tmp_path):
        item_path = os.path.join(tmp_path, item)
        if "-summary.txt" in item:
            os.remove(item_path)

    exonerate_verbose = 3 in verbose
    summary_file_out = []
    T_exonerate_genes = time()

    # Disperse into reftaxons
    if num_threads == 1:
        for orthoid in transcripts_mapped_to:
            list_of_hits = transcripts_mapped_to[orthoid]
            exonerate_gene_multi(
                orthoid,
                list_of_hits,
                orthoset_db_path,
                input_path,
                min_score,
                orthoset_id,
                aa_out_path,
                min_length,
                species_name,
                nt_out_path,
                tmp_path,
                exonerate_verbose,
            )

    else:
        arguments = list()
        for orthoid in transcripts_mapped_to:
            list_of_hits = transcripts_mapped_to[orthoid]
            arguments.append(
                (
                    orthoid,
                    list_of_hits,
                    orthoset_db_path,
                    input_path,
                    min_score,
                    orthoset_id,
                    aa_out_path,
                    min_length,
                    species_name,
                    nt_out_path,
                    tmp_path,
                    exonerate_verbose,
                )
            )
        with Pool(num_threads) as pool:
            pool.map(run_exonerate, arguments, chunksize=1)

    summary_file_out = []
    for item in os.listdir(tmp_path):
        item_path = os.path.join(tmp_path, item)
        if "-summary.txt" in item:
            item_content = open(item_path).read()
            summary_file_out.append(item_content)
            os.remove(item_path)

    summary_file_path = os.path.join(input_path, "summary.txt")
    open(summary_file_path, "w").write("\n".join(summary_file_out))

    if 1 in verbose:
        print(
            "Done. Final time {:.2f}s. Exonerate took {:.2f}s.".format(
                time() - T_global_start, time() - T_exonerate_genes
            )
        )

    if args.verbose == 0:
        print("Done took {:.2f}s.".format(time() - T_global_start))
