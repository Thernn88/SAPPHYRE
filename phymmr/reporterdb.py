from __future__ import annotations
import json
import math
import os
import shutil
import sqlite3
import uuid
from collections import namedtuple
from multiprocessing.pool import Pool
from time import time
from typing import List

import phymmr_tools
import wrap_rocks
import xxhash
from Bio.Seq import Seq

T_global_start = time()


class Hit:
    __slots__ = (
        "hmm_id",
        "header",
        "gene",
        "score",
        "hmm_start",
        "hmm_end",
        "ali_start",
        "ali_end",
        "est_sequence_hmm_region",
        "est_sequence_complete",
        "extended_orf_cdna_start",
        "orf_cdna_end_on_transcript",
        "extended_orf_cdna_end",
        "orf_cdna_start_on_transcript",
        "extended_orf_aa_start_on_transcript",
        "extended_orf_aa_end_on_transcript",
        "extended_orf_aa_end",
        "extended_orf_aa_start",
        "orf_aa_start_on_transcript",
        "orf_aa_end_on_transcript",
        "extended_orf_cdna_sequence",
        "orf_cdna_sequence",
        "extended_orf_aa_sequence",
        "orf_aa_sequence",
        "reftaxon",
        "proteome_sequence",
        "est_header",
        "est_sequence_complete",
        "est_sequence_hmm_region",
        "orf_cdna_start",
        "orf_cdna_end",
        "orf_aa_start",
        "orf_aa_end",
        "target",
    )

    def __init__(self, as_json):
        self.hmm_id = as_json["hmm_id"]
        self.header = as_json["header"]
        self.gene = as_json["gene"]
        self.score = as_json["score"]
        self.hmm_start = as_json["hmm_start"]
        self.hmm_end = as_json["hmm_end"]
        self.ali_start = as_json["ali_start"]
        self.ali_end = as_json["ali_end"]

    def remove_extended_orf(self):
        self.extended_orf_aa_sequence = None
        self.extended_orf_aa_start = None
        self.extended_orf_aa_end = None
        self.extended_orf_cdna_sequence = None
        self.extended_orf_cdna_start = None
        self.extended_orf_cdna_end = None

    def add_extended_orf(self, exonerate_record):
        (
            _,
            self.extended_orf_cdna_sequence,
            self.extended_orf_cdna_start,
            self.extended_orf_cdna_end,
            self.extended_orf_aa_sequence,
            self.extended_orf_aa_start,
            self.extended_orf_aa_end,
        ) = exonerate_record

        self.extended_orf_aa_start_on_transcript = (
            (self.extended_orf_cdna_start - 1) / 3
        ) + 1

        self.extended_orf_aa_end_on_transcript = (
            (self.extended_orf_cdna_end - 1) / 3
        ) + 1

        self.extended_orf_aa_sequence = translate_cdna(self.extended_orf_cdna_sequence)

    def add_orf(self, exonerate_record):
        (
            _,
            self.orf_cdna_sequence,
            self.orf_cdna_start,
            self.orf_cdna_end,
            self.orf_aa_sequence,
            self.orf_aa_start,
            self.orf_aa_end,
        ) = exonerate_record

        self.orf_aa_sequence = translate_cdna(self.orf_cdna_sequence)

        self.orf_cdna_start_on_transcript = (
            self.orf_cdna_start + (self.ali_start * 3) - 3
        )
        self.orf_cdna_end_on_transcript = self.orf_cdna_end + (self.ali_start * 3) - 3
        self.orf_aa_start_on_transcript = (
            self.orf_cdna_start + (self.ali_start * 3) - 3
        ) / 3
        self.orf_aa_end_on_transcript = (
            self.orf_cdna_end + (self.ali_start * 3) - 3
        ) / 3


class NodeRecord:
    def __init__(self, header, is_extension):
        self.header = header
        self.score = -math.inf
        self.cdna_start = None
        self.cdna_end = None
        self.cdna_sequence = ""
        self.aa_start = None
        self.aa_end = None
        self.aa_sequence = ""
        self.is_extension = is_extension

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

def get_set_id(orthoset_db_con, orthoset):
    """
    Retrieves orthoset id from orthoset's name.

    Args:
        orthoset_db_con: SQLite db pointer
        orthoset (str): name of the thing.

    Returns:
        (int) id
    """
    orthoset_id = None

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute(
        'SELECT id FROM orthograph_set_details WHERE name = "{}";'.format(orthoset)
    )

    if not rows:
        raise Exception("Orthoset {} id cant be retrieved".format(orthoset))

    # Return first result
    return next(rows)[0]


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

    return [row[1] for row in rows]


def get_scores_list(score_threshold, min_length, rocks_hits_db, debug):
    batches = rocks_hits_db.get("hmmbatch:all")
    batches = batches.split(",")

    score_based_results = {}
    ufr_out = [["Gene", "Hash", "Header", "Score", "Start", "End"]]

    for batch_i in batches:
        batch_rows = rocks_hits_db.get(f"hmmbatch:{batch_i}")
        batch_rows = json.loads(batch_rows)
        for this_row in batch_rows:
            length = this_row["env_end"] - this_row["env_start"]
            hmm_score = this_row["score"]

            if hmm_score > score_threshold and length > args.min_length:
                if debug:
                    hmm_start = this_row["hmm_start"]
                    hmm_end = this_row["hmm_end"]
                    components = this_row["header"].split("|")
                    if "[" in components[-1]:
                        out_header = "|".join(components[:-1]) + " " + components[-1]
                    else:
                        out_header = this_row["header"]

                    ufr_out.append(
                        [
                            this_row["gene"],
                            out_header,
                            str(hmm_score),
                            str(hmm_start),
                            str(hmm_end),
                        ]
                    )

                this_hit = Hit(this_row)

                score_based_results.setdefault(hmm_score, [])
                score_based_results[hmm_score].append(this_hit)

    if debug:
        ufr_out = sorted(ufr_out, key=lambda x: (x[0], x[1], x[3], x[4]))

    return score_based_results, ufr_out


def get_blastresults_for_hmmsearch_id(hmmsearch_id, rocks_hits_db):
    key = "blastfor:{}".format(hmmsearch_id)
    db_entry = rocks_hits_db.get(key)

    if not db_entry:
        return []

    result = json.loads(db_entry)
    return result


def get_reftaxon_name(hit_id, orthoset_db_path):
    query = f"""SELECT {orthoset_taxa}.{db_col_name}
FROM {orthoset_taxa}
    INNER JOIN {orthoset_aaseqs} ON {orthoset_taxa}.{db_col_id} = {orthoset_aaseqs}.{db_col_taxid}
WHERE {orthoset_aaseqs}.{db_col_id} = ?"""

    orthoset_db_con = sqlite3.connect(orthoset_db_path)
    orthoset_db_cur = orthoset_db_con.cursor()

    rows = orthoset_db_cur.execute(query, (hit_id,))

    # Return first result
    return next(rows)[0]


def calculate_length_of_ali(
    result,
):  # XXX: could be merged with transcript_not_long_enough
    return abs(result.ali_end - result.ali_start) + 1


def transcript_not_long_enough(result, minimum_transcript_length):
    length = calculate_length_of_ali(result)
    return length < minimum_transcript_length


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

    return next(rows)


def reverse_complement(nt_seq):
    return phymmr_tools.bio_revcomp(nt_seq)


def get_nucleotide_transcript_for(header, rocks_sequence_db):
    base_header = get_baseheader(header).strip()
    hash_of_header = xxhash.xxh64_hexdigest(base_header)

    row_data = rocks_sequence_db.get(hash_of_header)
    _, sequence = row_data.split("\n")

    if "revcomp" in header:
        return base_header, reverse_complement(sequence)
    return base_header, sequence


def crop_to_hmm_alignment(seq, header, hit):
    start = hit.ali_start - 1  # Adjust for zero based number
    start = start * 3
    end = hit.ali_end * 3

    return seq[start:end]


### EXONERATE FUNCTIONS


def fastaify(headers, sequences, tmp_path):
    name = str(uuid.uuid4()) + ".tmp"  # Super unique name
    path = os.path.join(tmp_path, name)

    with open(path, "w") as fp:
        for header, seq in zip(headers, sequences):
            fp.write(">" + header + "\n" + seq + "\n")

    return path


def translate_cdna(cdna_seq):
    if cdna_seq == None:
        return None

    if len(cdna_seq) % 3 != 0:
        print("WARNING: NT Sequence length is not divisable by 3")

    return str(Seq(cdna_seq).translate())


def parse_nodes(lines):
    """
    vulgar => cdna => aa => vulgar
    0 => 1 => 2 => 0
    """
    nodes = dict()
    this_node = None
    state = 0
    for line in lines:
        line = line.strip()
        if "vulgar:" in line:  # if vulgar line, start new record
            raw_header = line.split(" . ")[1].split(" ")[0]
            is_extension = False
            if "extended_" in raw_header:
                current_header = raw_header.replace("extended_", "")
                is_extension = True
            else:
                current_header = raw_header

            this_node = NodeRecord(current_header, is_extension)
            this_node.score = int(line.split()[9])
            previously_seen = nodes.get(raw_header, False)
            if previously_seen:
                nodes[raw_header] = max(this_node, previously_seen)
            else:
                nodes[raw_header] = this_node
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


def parse_multi_results(handle):
    extended_result = []
    result = []
    nodes = parse_nodes(handle)
    for key in nodes:
        node = nodes[key]
        if node.is_extension:
            extended_result.append(
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
        else:
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
    return extended_result, result


def get_multi_orf(query, targets, score_threshold, tmp_path, include_extended=False):
    exonerate_ryo = ">cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas"
    genetic_code = 1
    exonerate_model = "protein2genome"

    headers = [i.header for i in targets]
    if include_extended:
        headers.extend(["extended_" + i.header for i in targets])

    sequences = [i.est_sequence_hmm_region for i in targets]
    if include_extended:
        sequences.extend([i.est_sequence_complete for i in targets])

    queryfile = fastaify(["query"], [query], tmp_path)
    targetfile = fastaify(headers, sequences, tmp_path)

    outfile = os.path.join(tmp_path, str(uuid.uuid4()) + ".exonerateout")

    exonerate_cmd = f"exonerate --score {score_threshold} --ryo '{exonerate_ryo}' --subopt 0 --geneticcode {genetic_code} --model '{exonerate_model}' --querytype 'protein' --targettype 'dna' --verbose 0 --showalignment 'no' --showvulgar 'yes' --query '{queryfile}' --target '{targetfile}' > {outfile}"

    os.system(exonerate_cmd)

    with open(outfile) as fp:
        extended_results, results = parse_multi_results(fp)

    os.remove(queryfile)
    os.remove(targetfile)
    os.remove(outfile)

    return extended_results, results


def extended_orf_contains_original_orf(hit):
    return (
        hit.extended_orf_cdna_start > hit.orf_cdna_end_on_transcript
        or hit.extended_orf_cdna_end < hit.orf_cdna_start_on_transcript
    )


def get_overlap_length(candidate):
    if (
        candidate.extended_orf_aa_start_on_transcript <= candidate.ali_end
        and candidate.extended_orf_aa_end_on_transcript >= candidate.ali_start
    ):
        overlap_start = (
            candidate.extended_orf_aa_start_on_transcript
            if candidate.extended_orf_aa_start_on_transcript > candidate.ali_start
            else candidate.ali_start
        )
        overlap_end = (
            candidate.extended_orf_aa_end_on_transcript
            if candidate.extended_orf_aa_end_on_transcript > candidate.ali_end
            else candidate.ali_end
        )

        overlap_length = overlap_end - overlap_start
        return overlap_length
    return 0


def overlap_by_orf(candidate):
    orf_length = abs(candidate.extended_orf_aa_end - candidate.extended_orf_aa_start)
    overlap_length = get_overlap_length(candidate)

    return overlap_length / orf_length


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


def format_candidate_header(gene, taxa_name, taxa_id, sequence_id, frame):
    return header_seperator.join([gene, taxa_name, taxa_id, sequence_id, frame])


def format_reference_header(gene, taxa_name, taxa_id, identifier="."):
    return header_seperator.join([gene, taxa_name, taxa_id, identifier])


def print_core_sequences(orthoid, core_sequences):
    core_sequences = sorted(core_sequences)

    result = []
    for core in core_sequences:
        header = format_reference_header(orthoid, core[0], core[1])
        result.append(">" + header + "\n")
        result.append(core[2] + "\n")  # append the sequence

    return result


def get_rf(header):
    raw_header = get_baseheader(header)
    header_part = get_translate(header)
    if "translate" in header_part or "revcomp" in header_part:
        frame = header_part

    return raw_header, frame


def print_unmerged_sequences(
    hits, orthoid, type, minimum_seq_data_length, taxa_id, kicks
):
    result = []
    kicks_result = set()
    exact_hit_mapped_already = set()
    for i, hit in enumerate(hits):
        base_header, rf = get_rf(hit.header)

        header = format_candidate_header(
            orthoid,
            hit.reftaxon,
            taxa_id,
            base_header,
            rf,
        )

        if type == "nt":
            seq = (
                hit.extended_orf_cdna_sequence
                if hit.extended_orf_cdna_sequence is not None
                else hit.orf_cdna_sequence
            )

            #Deinterleave
            seq = seq.replace('\n','')

            if i not in kicks:
                result.append(">" + header + "\n")
                result.append(seq + "\n")
        elif type == "aa":
            seq = (
                hit.extended_orf_aa_sequence
                if hit.extended_orf_aa_sequence is not None
                else hit.orf_aa_sequence
            )

            seq = seq.replace('\n','')
            unique_hit = base_header+seq

            if not unique_hit in exact_hit_mapped_already and len(seq) - seq.count("-") > minimum_seq_data_length:
                result.append(">" + header + "\n")
                result.append(seq + "\n")
                exact_hit_mapped_already.add(unique_hit)
            else:
                kicks_result.add(i)
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
    Returns header content before first pipe.
    """
    baseheader = header.split("|")[0].strip()
    return baseheader


def get_translate(header):
    """
    Returns header content of second pipe.
    """
    translate = header.split("|")[1]
    return translate


def get_match(header, results):
    for result in results:
        if result[0] == header:
            return result
    return None


ExonerateArgs = namedtuple("ExonerateArgs",
    [
        "orthoid",
        "list_of_hits",
        "orthoset_db_path",
        "min_score",
        "orthoset_id",
        "aa_out_path",
        "min_length",
        "taxa_id",
        "nt_out_path",
        "tmp_path",
        "verbose",
        "rdb_kw", # kwargs: {'rocks_sequence_db': ..., 'rocks_hits_db': ...}
    ]
)


def run_exonerate(arg_tuple: ExonerateArgs):
    exonerate_gene_multi(arg_tuple)


def exonerate_gene_multi(eargs: ExonerateArgs):
    if eargs.verbose >= 2:
        print("Exonerating and doing output for: ", eargs.orthoid)
        T_gene_start = time()

    orthoset_db_con = sqlite3.connect(eargs.orthoset_db_path)
    reftaxon_related_transcripts = {}
    reftaxon_to_proteome_sequence = {}
    for hit in eargs.list_of_hits:
        this_reftaxon = hit.reftaxon

        est_header, est_sequence_complete = get_nucleotide_transcript_for(
            hit.header, rdb_kw.rocks_sequence_db
        )
        est_sequence_hmm_region = crop_to_hmm_alignment(
            est_sequence_complete, est_header, hit
        )

        hit.est_header = est_header
        hit.est_sequence_complete = est_sequence_complete
        hit.est_sequence_hmm_region = est_sequence_hmm_region

        if this_reftaxon not in reftaxon_related_transcripts:
            reftaxon_to_proteome_sequence[this_reftaxon] = hit.proteome_sequence
            reftaxon_related_transcripts[this_reftaxon] = []

        reftaxon_related_transcripts[this_reftaxon].append(hit)

    output_sequences = []
    total_results = 0
    for taxon_hit, hits in reftaxon_related_transcripts.items():
        total_results += len(hits)
        query = reftaxon_to_proteome_sequence[taxon_hit]
        extended_results, results = get_multi_orf(
            query, hits, eargs.min_score, eargs.tmp_path, include_extended=extend_orf
        )

        for hit in hits:
            matching_alignment = get_match(hit.header, results)
            if matching_alignment:
                hit.add_orf(matching_alignment)

                if hit.orf_cdna_sequence:
                    if extend_orf:
                        matching_extended_alignment = get_match(
                            hit.header, extended_results
                        )

                        if matching_extended_alignment:
                            hit.add_extended_orf(matching_extended_alignment)

                            if extended_orf_contains_original_orf(hit):
                                orf_overlap = overlap_by_orf(hit)
                                if orf_overlap < orf_overlap_minimum:
                                    hit.remove_extended_orf()
                            else:
                                hit.remove_extended_orf()
                        else:
                            print("Failed to extend orf on {}".format(hit.header))
                    output_sequences.append(hit)
    if len(output_sequences) > 0:
        output_sequences = sorted(output_sequences, key=lambda d: d.hmm_start)
        core_sequences = get_ortholog_group(
            eargs.orthoset_id, eargs.orthoid, orthoset_db_con
        )
        this_aa_out = []
        this_aa_path = os.path.join(eargs.aa_out_path, eargs.orthoid + ".aa.fa")
        this_aa_out.extend(print_core_sequences(eargs.orthoid, core_sequences))
        kicks, output = print_unmerged_sequences(
            output_sequences,
            eargs.orthoid,
            "aa",
            eargs.min_length,
            eargs.taxa_id,
            kicks=set()
        )

        if output:
            this_aa_out.extend(output)
            with open(this_aa_path, "w") as fp:
                fp.writelines(this_aa_out)

            core_sequences_nt = get_ortholog_group_nucleotide(
                eargs.orthoset_id, eargs.orthoid, orthoset_db_con
            )

            this_nt_out = []
            this_nt_path = os.path.join(eargs.nt_out_path, eargs.orthoid + ".nt.fa")

            this_nt_out.extend(print_core_sequences(eargs.orthoid, core_sequences_nt))
            na_kicks, output = print_unmerged_sequences(
                output_sequences,
                eargs.orthoid,
                "nt",
                eargs.min_length,
                eargs.taxa_id,
                kicks=kicks
            )
            this_nt_out.extend(output)

            with open(this_nt_path, "w") as fp:
                fp.writelines(this_nt_out)

    if eargs.verbose >= 2:
        print(
            "{} took {:.2f}s. Had {} sequences".format(
                eargs.orthoid, time() - T_gene_start, len(output_sequences)
            )
        )

def is_reciprocal_match(blast_results, reference_taxa: List[str]):
    reftaxon_count = {ref_taxa: 0 for ref_taxa in reference_taxa}
    ref_taxon_to_target = {}

    for result in blast_results:
        ref_taxon = result["reftaxon"]

        if ref_taxon in reftaxon_count:
            if not strict_search_mode:
                return ref_taxon, result["ref_sequence"]
            if ref_taxon not in ref_taxon_to_target:
                ref_taxon_to_target[ref_taxon] = (ref_taxon, result["ref_sequence"])
            reftaxon_count[ref_taxon] += 1  # only need the one
            if all(reftaxon_count.values()):  # Everything's been counted
                return ref_taxon_to_target[
                    max(reftaxon_count)
                ]  # Grab most hit reftaxon
    return None, None

def reciprocal_search(
    hmmresults, list_of_wanted_orthoids, reference_taxa, score, args, rocks_hits_db
):
    if args.verbose >= 3:
        T_reciprocal_start = time()
        print("Ensuring reciprocal hit for hmmresults in {}".format(score))

    results = []
    for result in hmmresults:
        orthoid = result.gene

        if list_of_wanted_orthoids and orthoid not in list_of_wanted_orthoids:
            continue

        result_hmmsearch_id = result.hmm_id
        blast_results = get_blastresults_for_hmmsearch_id(
            result_hmmsearch_id, rocks_hits_db
        )
        this_match_reftaxon, this_match_ref_sequence = is_reciprocal_match(
            blast_results, reference_taxa
        )

        if this_match_reftaxon:
            result.proteome_sequence = this_match_ref_sequence
            result.reftaxon = this_match_reftaxon
            results.append(result)

    if args.verbose >= 3:
        print(
            "Checked reciprocal hits for {}. Took {:.2f}s.".format(
                score, time() - T_reciprocal_start
            )
        )
    return results


def do_taxa(path, taxa_id, argsn **kwargs):
    print("Doing {}.".format(taxa_id))
    if args.verbose >= 1:
        T_init_db = time()

    num_threads = args.processes
    if not isinstance(num_threads, int) or num_threads < 1:
        num_threads = 1

    if os.path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif os.path.exists("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = os.path.join(path, "tmp")
        os.makedirs(tmp_path, exist_ok=True)

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

    aa_out = "aa"
    nt_out = "nt"

    aa_out_path = os.path.join(path, aa_out)
    nt_out_path = os.path.join(path, nt_out)

    if clear_output:
        if os.path.exists(aa_out_path):
            shutil.rmtree(aa_out_path)
        if os.path.exists(nt_out_path):
            shutil.rmtree(nt_out_path)

    os.makedirs(aa_out_path, exist_ok=True)
    os.makedirs(nt_out_path, exist_ok=True)

    if args.verbose >= 1:
        T_reference_taxa = time()
        print(
            "Initialized databases. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing reference taxa in set.".format(
                time() - T_global_start, time() - T_init_db
            )
        )

    reference_taxa = get_taxa_in_set(orthoset_id, orthoset_db_con)

    if args.verbose >= 1:
        T_hmmresults = time()
        print(
            "Got reference taxa in set. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing hmmresults".format(
                time() - T_global_start, time() - T_reference_taxa
            )
        )

    score_based_results, ufr_rows = get_scores_list(
        args.min_score, args.min_length, kwargs.rocks_hits_db, args.debug
    )

    if args.debug:
        ufr_path = os.path.join(path, "unfiltered-hits.csv")

        ufr_out = ["Gene,Header,Score,Start,End\n"]
        for row in ufr_rows:
            ufr_out.append(",".join(row) + "\n")

        open(ufr_path, "w").writelines(ufr_out)

    ####################################
    if args.verbose >= 1:
        print(
            "Got hmmresults. Elapsed time {:.2f}s. Took {:.2f}s.".format(
                time() - T_global_start, time() - T_hmmresults
            )
        )
        print(
            "Retrieved data from DB. Elapsed time {:.2f}s. Took {:.2f}s. Doing reciprocal check.".format(
                time() - T_global_start, time() - T_hmmresults
            )
        )
        T_reciprocal_search = time()

    scores = list(score_based_results.keys())
    scores.sort(reverse=True)  # Ascending
    transcripts_mapped_to = {}

    arguments = []
    for score in scores:
        hmmresults = score_based_results[score]
        arguments.append(
            (
                hmmresults,
                list_of_wanted_orthoids,
                reference_taxa,
                score,
                args,
                kwargs.rocks_hits_db,
            )
        )
    with Pool(num_threads) as pool:
        reciprocal_data = pool.starmap(reciprocal_search, arguments, chunksize=1)

    brh_count = 0

    for data in reciprocal_data:
        brh_count += len(data)
        for this_match in data:
            orthoid = this_match.gene

            if transcript_not_long_enough(this_match, args.min_length):
                continue

            if orthoid not in transcripts_mapped_to:
                transcripts_mapped_to[orthoid] = []

            transcripts_mapped_to[orthoid].append(this_match)

    if args.verbose >= 1:
        T_internal_search = time()
        print(
            "Reciprocal check done, found {} reciprocal hits. Elapsed time {:.2f}s. Took {:.2f}s. Exonerating genes.".format(
                brh_count, time() - T_global_start, time() - T_reciprocal_search
            )
        )

    T_exonerate_genes = time()

    if num_threads > 1:
        arguments = []
        func = arguments.append
    else:
        func = exonerate_gene_multi

    # this sorting the list so that the ones with the most hits are first
    for orthoids in sorted(
        transcripts_mapped_to,
        key=lambda k: len(transcripts_mapped_to[k]),
        reverse=True
    ):
        func(
            ExonerateArgs(
                orthoids,
                transcripts_mapped_to[orthoids],
                orthoset_db_path,
                args.min_score,
                orthoset_id,
                aa_out_path,
                args.min_length,
                taxa_id,
                nt_out_path,
                tmp_path,
                args.verbose,
                kwargs
            )
        )

    if num_threads > 1:
        with Pool(num_threads) as pool:
            pool.map(run_exonerate, arguments, chunksize=1)

    if args.verbose >= 1:
        print(
            "Done. Final time {:.2f}s. Exonerate took {:.2f}s.".format(
                time() - T_global_start, time() - T_exonerate_genes
            )
        )

    if args.verbose == 0:
        print("Done took {:.2f}s.".format(time() - T_global_start))


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

# TODO Make these argparse variables
strict_search_mode = False
orthoid_list_file = None
extend_orf = True
orf_overlap_minimum = 0.15
clear_output = True

header_seperator = "|"

####

def main(args):
    if not all(os.path.exists(i) for i in args.INPUT):
        print("ERROR: All folders passed as argument must exists.")
        return False
    for input_path in args.INPUT:
        print(f"### Processing path '{input_path}'.")
        rocks_db_path = os.path.join(input_path, "rocksdb")
        rocks_sequence_db = wrap_rocks.RocksDB(os.path.join(rocks_db_path, "sequences"))
        rocks_hits_db = wrap_rocks.RocksDB(os.path.join(rocks_db_path, "hits"))
        do_taxa(
            path=input_path,
            taxa_id=os.path.basename(input_path).split(".")[0],
            args=args,
            rocks_sequence_db=rocks_sequence_db,
            rocks_hits_db=rocks_hits_db,
        )
    return True


if __name__ == '__main__':
    main()
