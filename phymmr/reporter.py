from __future__ import annotations

import json
import math
import os
import shutil
import sys
import uuid
from collections import namedtuple
from multiprocessing.pool import Pool
from time import time
from typing import List, Optional
from tempfile import TemporaryDirectory, NamedTemporaryFile
from concurrent.futures import ThreadPoolExecutor, wait

import phymmr_tools
import xxhash
from Bio.Seq import Seq

from . import rocky
from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv, gettempdir

MainArgs = namedtuple(
    "MainArgs",
    [
        'verbose',
        'processes',
        'debug',
        'INPUT',
        'orthoset_input',
        'orthoset',
        'min_length',
        'min_score',
    ]
)

class Result:
    __slots__ = (
        "hmm_id",
        "gene",
        "ref_taxon",
        )

    def __init__(self, as_json) -> None:
        self.hmm_id = as_json["hmmId"]
        self.gene = as_json["gene"]
        self.ref_taxon = as_json["refTaxon"]

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
        "est_header",
        "est_sequence_hmm_region",
        "est_sequence_complete",
        "est_sequence_complete",
        "est_sequence_hmm_region",
        "reftaxon",
        "first_closest",
        "second_closest",
        "mapped_to",
        "first_alignment",
        "first_extended_alignment",
        "second_alignment",
        "second_extended_alignment"
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
        self.second_closest = None
        self.mapped_to = None
        self.first_alignment = None
        self.first_extended_alignment = None
        self.second_alignment = None
        self.second_extended_alignment = None

    def add_orf(self, exonerate_record): 
        exonerate_record.orf_cdna_start_on_transcript = (
            exonerate_record.orf_cdna_start + (self.ali_start * 3) - 3
        )
        exonerate_record.orf_cdna_end_on_transcript = exonerate_record.orf_cdna_end + (self.ali_start * 3) - 3
        exonerate_record.orf_aa_start_on_transcript = (
            exonerate_record.orf_cdna_start + (self.ali_start * 3) - 3
        ) / 3
        exonerate_record.orf_aa_end_on_transcript = (
            exonerate_record.orf_cdna_end + (self.ali_start * 3) - 3
        ) / 3
        self.first_alignment = exonerate_record

class NodeRecord:
    def __init__(self, header, is_extension):
        self.header = header
        self.score = -math.inf
        self.orf_cdna_start = None
        self.orf_cdna_end = None
        self.orf_cdna_sequence = ""
        self.orf_aa_start = None
        self.orf_aa_end = None
        self.orf_aa_sequence = ""
        self.is_extension = is_extension
        self.orf_cdna_end_on_transcript = None
        self.orf_aa_start_on_transcript = None
        self.extended_orf_cdna_sequence = None
        self.extended_orf_aa_sequence = None
        
    def check_extend(self):
        if self.is_extension:
            self.extended_orf_cdna_sequence = self.orf_cdna_sequence
            self.extended_orf_cdna_start = self.orf_cdna_start
            self.extended_orf_cdna_end = self.orf_cdna_end
            self.extended_orf_aa_sequence = self.orf_aa_sequence
            self.extended_orf_aa_start = self.orf_aa_start
            self.extended_orf_aa_end = self.orf_aa_end

        self.extended_orf_aa_start_on_transcript = (
            (self.extended_orf_cdna_start - 1) / 3
        ) + 1

        self.extended_orf_aa_end_on_transcript = (
            (self.extended_orf_cdna_end - 1) / 3
        ) + 1

        self.extended_orf_aa_sequence = translate_cdna(self.extended_orf_cdna_sequence)

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

def get_hmmresults(score_threshold, min_length, rocks_hits_db, debug):
    batches = rocks_hits_db.get("hmmbatch:all")
    batches = batches.split(",")

    gene_based_results = {}
    ufr_out = [["Gene", "Hash", "Header", "Score", "Start", "End"]]
    for batch_i in batches:
        batch_rows = rocks_hits_db.get(f"hmmbatch:{batch_i}")
        batch_rows = json.loads(batch_rows)
        for this_row in batch_rows:
            length = this_row["env_end"] - this_row["env_start"]
            hmm_score = this_row["score"]
            gene = this_row["gene"]
            if hmm_score > score_threshold and length >= min_length:
                if debug:
                    hmm_start = this_row["hmm_start"]
                    hmm_end = this_row["hmm_end"]
                    components = this_row["header"].split("|")
                    if "[" in components[-1]:
                        out_header = "|".join(components[:-1]) + " " + components[-1]
                    else:
                        out_header = this_row["header"]

                    if debug:
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

                gene_based_results.setdefault(gene, [])
                gene_based_results[gene].append(this_hit)

    if debug:
        ufr_out = sorted(ufr_out, key=lambda x: (x[0], x[1], x[3], x[4]))
    return gene_based_results, ufr_out


def get_blastresults(rocks_hits_db):
    batches = rocks_hits_db.get("blastbatch:all")
    batches = batches.split(",")

    blast_results = {}

    for batch_i in batches:
        batch_rows = rocks_hits_db.get(f"blastbatch:{batch_i}")
        batch_rows = json.loads(batch_rows)
        for result in batch_rows:
            gene = result["gene"]
            hmm_id = result["hmmId"]
            this_result = Result(result)
            blast_results.setdefault(gene, {})
            blast_results[gene].setdefault(hmm_id, [])
            blast_results[gene][hmm_id].append(this_result)

    return blast_results


def get_reference_data(rocks_hits_db):
    raw_data = rocks_hits_db.get("getall:refseqs")

    processed = json.loads(raw_data)

    return processed

def reverse_complement(nt_seq):
    return phymmr_tools.bio_revcomp(nt_seq)


def get_nucleotide_transcript_for(header):
    base_header = header.split("|")[0]
    hash_of_header = xxhash.xxh64_hexdigest(base_header)

    row_data = rocky.get_rock("rocks_nt_db").get(hash_of_header)
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

def translate_cdna(cdna_seq):
    if not cdna_seq:
        return None

    if len(cdna_seq) % 3 != 0:
        printv("WARNING: NT Sequence length is not divisable by 3", 0)

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
            this_node.orf_cdna_start = int(fields[1])
            this_node.orf_cdna_end = int(fields[2])
            state = 1
        elif ">aa" in line:  # if aa in line, set aa values
            fields = line.split()
            this_node.orf_aa_start = int(fields[1])
            this_node.orf_aa_end = int(fields[2])
            state = 2
        elif state == 1:  # if cdna data line
            this_node.orf_cdna_sequence += line
        elif state == 2:  # if aa data line
            this_node.orf_aa_sequence += line
    return nodes


def parse_multi_results(handle):
    extended_result = []
    result = []
    handle.seek(0)
    nodes = parse_nodes(handle)
    for node in nodes.values():
        if node.is_extension:
            node.check_extend()
            extended_result.append(
                (
                    node
                )
            )
        else:
            node.orf_aa_sequence = translate_cdna(node.orf_cdna_sequence)
            result.append(
                (
                    node
                )
            )
    return extended_result, result

def get_multi_orf(query, targets, score_threshold, include_extended):
    exonerate_ryo = ">cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas"
    genetic_code = 1
    exonerate_model = "protein2genome"

    sequences = [(i.header, i.est_sequence_hmm_region) for i in targets]
    if include_extended:
        sequences.extend([("extended_" + i.header, i.est_sequence_complete) for i in targets])

    with TemporaryDirectory(dir=gettempdir()) as tmpdir, NamedTemporaryFile(dir=tmpdir, mode="w+") as tmpquery, NamedTemporaryFile(dir=tmpdir, mode="w+") as tmptarget, NamedTemporaryFile(dir=tmpdir, mode="r") as tmpout:
        tmptarget.write("".join([f">{header}\n{sequence}\n" for header, sequence in sequences]))
        tmptarget.flush() # Flush the internal buffer so it can be read by exonerate

        tmpquery.write(f">query\n{query}\n")
        tmpquery.flush() # Flush the internal buffer so it can be read by exonerate

        exonerate_cmd = f"exonerate --score {score_threshold} --ryo '{exonerate_ryo}' --subopt 0 --geneticcode {genetic_code} --model '{exonerate_model}' --querytype 'protein' --targettype 'dna' --verbose 0 --showalignment 'no' --showvulgar 'yes' --query '{tmpquery.name}' --target '{tmptarget.name}' > {tmpout.name}"

        os.system(exonerate_cmd)

        extended_results, results = parse_multi_results(tmpout)

    return extended_results, results


def extended_orf_contains_original_orf(extended_alignment, alignment):
    return (
        extended_alignment.extended_orf_cdna_start > alignment.orf_cdna_end_on_transcript
        or extended_alignment.extended_orf_cdna_end < alignment.orf_cdna_start_on_transcript
    )


def get_overlap_length(candidate, alignment):
    if (
        alignment.extended_orf_aa_start_on_transcript <= candidate.ali_end
        and alignment.extended_orf_aa_end_on_transcript >= candidate.ali_start
    ):
        overlap_start = (
            alignment.extended_orf_aa_start_on_transcript
            if alignment.extended_orf_aa_start_on_transcript > candidate.ali_start
            else candidate.ali_start
        )
        overlap_end = (
            alignment.extended_orf_aa_end_on_transcript
            if alignment.extended_orf_aa_end_on_transcript > candidate.ali_end
            else candidate.ali_end
        )

        overlap_length = overlap_end - overlap_start
        return overlap_length
    return 0


def overlap_by_orf(candidate, alignment):
    orf_length = abs(alignment.extended_orf_aa_end - alignment.extended_orf_aa_start)
    overlap_length = get_overlap_length(candidate, alignment)

    return overlap_length / orf_length


### FIN


def get_ortholog_group(orthoid, orthoset_db):
    core_seqs = json.loads(orthoset_db.get(f"getcore:{orthoid}"))
    return core_seqs["aa"], core_seqs["nt"]

def format_candidate_header(gene, taxa_name, taxa_id, sequence_id, frame):
    return header_seperator.join([gene, taxa_name, taxa_id, sequence_id, frame])


def format_reference_header(gene, taxa_name, taxa_id, identifier="."):
    return header_seperator.join([gene, taxa_name, taxa_id, identifier])


def print_core_sequences(orthoid, core_sequences):
    core_sequences = sorted(core_sequences)

    result = []
    for core in core_sequences:
        header = format_reference_header(orthoid, core[0], core[1])
        result.append(f">{header}\n{core[2]}\n")

    return result

def print_unmerged_sequences(
    hits, orthoid, taxa_id
):
    aa_result = []
    nt_result = []
    header_maps_to_where = {}
    header_mapped_x_times = {}
    base_header_mapped_already = {}
    exact_hit_mapped_already = set()
    for hit in hits:
        base_header, reference_frame = hit.header.split('|')

        header = format_candidate_header(
            orthoid,
            hit.reftaxon,
            taxa_id,
            base_header,
            reference_frame,
        )

        alignment = hit.first_alignment if hit.mapped_to == hit.first_closest else hit.second_alignment
        extended_alignment = hit.first_extended_alignment if hit.mapped_to == hit.first_closest else hit.second_extended_alignment
        nt_seq = (
            extended_alignment.extended_orf_cdna_sequence
            if extended_alignment and extended_alignment.extended_orf_cdna_sequence is not None
            else alignment.orf_cdna_sequence
        )

        # De-interleave
        nt_seq = nt_seq.replace('\n', '')

        aa_seq = (
            extended_alignment.extended_orf_aa_sequence
            if extended_alignment and extended_alignment.extended_orf_aa_sequence is not None
            else alignment.orf_aa_sequence
        )

        aa_seq = aa_seq.replace('\n', '')
        unique_hit = base_header+aa_seq

        if (
            unique_hit not in exact_hit_mapped_already
        ):
            if base_header in base_header_mapped_already:
                already_mapped_header, already_mapped_sequence = base_header_mapped_already[base_header]

                if len(aa_seq) > len(already_mapped_sequence):
                    if already_mapped_sequence in aa_seq:
                        aa_result[header_maps_to_where[already_mapped_header]] = f">{header}\n{aa_seq}\n" 
                        nt_result[header_maps_to_where[already_mapped_header]] = f">{header}\n{nt_seq}\n" 
                        continue
                else:
                    if aa_seq in already_mapped_sequence:
                        continue

                if base_header in header_mapped_x_times:
                    # Make header unique
                    old_header = base_header
                    header = format_candidate_header(
                        orthoid,
                        hit.reftaxon,
                        taxa_id,
                        base_header+f"_{header_mapped_x_times[old_header]}",
                        reference_frame,
                    )

                    header_mapped_x_times[base_header] += 1
            else:
                base_header_mapped_already[base_header] = header, aa_seq

            header_maps_to_where[header] = len(aa_result) # Save the index of the sequence output
            aa_result.append(f">{header}\n{aa_seq}\n" )
            nt_result.append(f">{header}\n{nt_seq}\n" )

            header_mapped_x_times.setdefault(base_header, 1)
            exact_hit_mapped_already.add(unique_hit)

    return aa_result, nt_result


def get_difference(score_a, score_b):
    """
    Returns decimal difference of two scores
    """

    try:
        if score_a / score_b > 1:
            return score_a / score_b
        if score_b / score_a > 1:
            return score_b / score_a
        if score_a == score_b:
            return 1
    except ZeroDivisionError:
        return 0

def parse_results(results):
    res = {}
    for result in results:
        res[result.header] = result
    
    return res

ExonerateArgs = namedtuple(
    "ExonerateArgs",
    [
        "orthoid",
        "list_of_hits",
        "min_score",
        "aa_out_path",
        "min_length",
        "taxa_id",
        "nt_out_path",
        "tmp_path",
        "verbose",
        "reference_sequences",
    ]
)

def exonerate_gene_multi(eargs: ExonerateArgs):

    t_gene_start = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Exonerating and doing output for: {eargs.orthoid}", eargs.verbose, 2)

    reftaxon_related_transcripts = {i: [] for i in eargs.reference_sequences.keys()}
    for hit in eargs.list_of_hits:
        this_reftaxon = hit.reftaxon

        est_header, est_sequence_complete = get_nucleotide_transcript_for(
            hit.header
        )
        est_sequence_hmm_region = crop_to_hmm_alignment(
            est_sequence_complete, est_header, hit
        )

        hit.est_header = est_header
        hit.est_sequence_complete = est_sequence_complete
        hit.est_sequence_hmm_region = est_sequence_hmm_region

        reftaxon_related_transcripts[this_reftaxon].append(hit)

        # If it doesn't align to closest ref fall back to this ref and try again
        if hit.second_closest is not None:
            reftaxon_related_transcripts[hit.second_closest].append(hit)

    output_sequences = []
    total_results = 0
    for taxon_hit, hits in reftaxon_related_transcripts.items():
        hits_to_exonerate = [i for i in hits if i.mapped_to is None]
        if len(hits_to_exonerate) == 0:
            continue

        total_results += len(hits)
        query = eargs.reference_sequences[taxon_hit]
        extended_results, results = get_multi_orf(
            query, hits_to_exonerate, eargs.min_score, include_extended=extend_orf
        )

        extended_results = parse_results(extended_results)
        results = parse_results(results)

        for hit in hits_to_exonerate:
            # We still want to see the first results before rerun
            if taxon_hit == hit.second_closest and hit.mapped_to is None:
                hit.second_alignment = results.get(hit.header, None)
                hit.second_extended_alignment = extended_results.get(hit.header, None)
                continue

            if hit.header in results:
                matching_alignment = results.get(hit.header, None)
                if matching_alignment.orf_cdna_sequence:
                    hit.add_orf(matching_alignment)
                    if extend_orf:
                        if hit.header in extended_results:
                            hit.first_extended_alignment = extended_results.get(hit.header, None)

                            if extended_orf_contains_original_orf(hit.first_extended_alignment, matching_alignment):
                                orf_overlap = overlap_by_orf(hit, hit.first_extended_alignment)
                                if orf_overlap >= orf_overlap_minimum:
                                    continue
                            
                            # Does not contain original orf or does not overlap enough.
                            hit.first_extended_alignment = None

                    aa_seq = (
                            hit.first_extended_alignment.extended_orf_aa_sequence
                            if hit.first_extended_alignment and hit.first_extended_alignment.extended_orf_aa_sequence is not None
                            else hit.first_alignment.orf_aa_sequence
                        )

                    if len(aa_seq) >= eargs.min_length:
                        hit.mapped_to = hit.reftaxon
                        output_sequences.append(hit)
                    else:
                        if hit.second_alignment is not None: # If first run doesn't pass check if rerun does
                            matching_alignment = hit.second_alignment
                            if matching_alignment:
                                hit.add_orf(matching_alignment)

                                if matching_alignment.orf_cdna_sequence:
                                    if extend_orf:
                                        if hit.second_extended_alignment:
                                            if extended_orf_contains_original_orf(hit.second_extended_alignment, matching_alignment):
                                                orf_overlap = overlap_by_orf(hit, hit.second_extended_alignment)
                                                if orf_overlap >= orf_overlap_minimum:
                                                    continue
                                            
                                            # Does not contain original orf or does not overlap enough.
                                            hit.second_extended_alignment = None

                                    aa_seq = (
                                            hit.second_extended_alignment.extended_orf_aa_sequence
                                            if hit.second_extended_alignment and hit.second_extended_alignment.extended_orf_aa_sequence is not None
                                            else hit.second_alignment.orf_aa_sequence
                                        )
                                        
                                    if len(aa_seq) >= eargs.min_length:
                                        hit.reftaxon = hit.second_closest
                                        hit.mapped_to = hit.reftaxon
                                        output_sequences.append(hit)
    if len(output_sequences) > 0:
        output_sequences = sorted(output_sequences, key=lambda d: d.hmm_start)
        core_sequences, core_sequences_nt = get_ortholog_group(eargs.orthoid, rocky.get_rock("rocks_orthoset_db"))
        this_aa_path = os.path.join(eargs.aa_out_path, eargs.orthoid + ".aa.fa")
        aa_output, nt_output = print_unmerged_sequences(
            output_sequences,
            eargs.orthoid,
            eargs.taxa_id,
        )
        
        if aa_output:
            with open(this_aa_path, "w") as fp:
                fp.writelines(print_core_sequences(eargs.orthoid, core_sequences))
                fp.writelines(aa_output)

            this_nt_path = os.path.join(eargs.nt_out_path, eargs.orthoid + ".nt.fa")

            with open(this_nt_path, "w") as fp:
                fp.writelines(print_core_sequences(eargs.orthoid, core_sequences_nt))
                fp.writelines(nt_output)

    printv(f"{eargs.orthoid} took {t_gene_start.differential():.2f}s. Had {len(output_sequences)} sequences", eargs.verbose, 2)
    return len(output_sequences)

def reciprocal_search(
    arg_tuple
):
    hmmresults, gene_blast_results, reference_taxa, gene, verbose = arg_tuple
    t_reciprocal_start = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Ensuring reciprocal hit for hmmresults in {gene}", verbose, 2)
    results = []

    for hit in hmmresults:
        hmm_id = str(hit.hmm_id)

        if hmm_id not in gene_blast_results:
            continue

        blast_results = gene_blast_results[hmm_id]
        reftaxon_count = {ref_taxa: 0 for ref_taxa in reference_taxa}
        this_match_reftaxon = None
        this_second_match = None

        for result in blast_results:
            ref_taxon = result.ref_taxon

            if ref_taxon in reftaxon_count:
                if not this_match_reftaxon:
                    this_match_reftaxon = ref_taxon
                elif not this_second_match and ref_taxon != this_match_reftaxon:
                    this_second_match = ref_taxon
                
                if not strict_search_mode and this_match_reftaxon and this_second_match:
                    break
                elif all(reftaxon_count.values()):
                    break

                reftaxon_count[ref_taxon] = 1  # only need the one
        
        if not strict_search_mode and this_match_reftaxon or all(reftaxon_count.values()):
            hit.reftaxon = this_match_reftaxon
            hit.first_closest = this_match_reftaxon
            if this_second_match:
                hit.second_closest = this_second_match
            results.append(hit)

    printv(f"Checked reciprocal hits for {gene}. Took {t_reciprocal_start.differential():.2f}s.", verbose, 2)
    return results


def do_taxa(path, taxa_id, args):
    printv(f"Processing: {taxa_id}", args.verbose, 0)
    time_keeper = TimeKeeper(KeeperMode.DIRECT)

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

    if orthoid_list_file:  # TODO orthoid_list_file must be passed as argument
        with open(orthoid_list_file) as fp:
            list_of_wanted_orthoids = fp.read().split("\n")
    else:
        list_of_wanted_orthoids = []

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

    printv(f"Initialized databases. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reference taxa in set.", args.verbose)

    reference_taxa = json.loads(rocky.get_rock("rocks_orthoset_db").get("getall:taxainset"))

    printv(f"Got reference taxa in set. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing hmmresults.", args.verbose)

    gene_based_results, ufr_rows = get_hmmresults(
        args.min_score, args.min_length, rocky.get_rock("rocks_hits_db"), args.debug
    )

    if args.debug:
        ufr_path = os.path.join(path, "unfiltered-hits.csv")
        ufr_out = ["Gene,Header,Score,Start,End\n"]
        for row in ufr_rows:
            ufr_out.append(",".join(row) + "\n")
        with open(ufr_path, "w") as fp:
            fp.writelines(ufr_out)

    printv(f"Got hmmresults. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing Blast results.", args.verbose)

    blast_resultss = get_blastresults(rocky.get_rock("rocks_hits_db"))

    printv(f"Got Blast results. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Doing reciprocal check.", args.verbose)

    genes = list(gene_based_results.keys())
    genes.sort(key = lambda x: len(x), reverse=True)  # Ascending
    transcripts_mapped_to = {}

    argmnts = []
    for gene, hmmresults in gene_based_results.items():
        if list_of_wanted_orthoids and gene not in list_of_wanted_orthoids:
            continue
        argmnts.append(
            (
                hmmresults,
                blast_resultss.get(gene, []),
                reference_taxa,
                gene,
                args.verbose,
            )
        )
    with ThreadPoolExecutor(num_threads) as pool:
        reciprocal_data = pool.map(reciprocal_search, argmnts)
    # with Pool(num_threads) as pool:
    #     reciprocal_data = pool.starmap(reciprocal_search, argmnts, chunksize=1)

    brh_count = 0
    for data in reciprocal_data:
        brh_count += len(data)
        for this_match in data:
            orthoid = this_match.gene
            if orthoid not in transcripts_mapped_to:
                transcripts_mapped_to[orthoid] = []

            transcripts_mapped_to[orthoid].append(this_match)

    printv(f"Reciprocal check done, found {brh_count} reciprocal hits. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing references sequences.", args.verbose)

    gene_reference_data = get_reference_data(rocky.get_rock("rocks_orthoset_db"))

    printv(f"Got reference data. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Exonerating genes.", args.verbose)

    arguments: list[Optional[ExonerateArgs]] = []
    for orthoid in sorted(
        transcripts_mapped_to,
        key=lambda k: len(transcripts_mapped_to[k]),
        reverse=True
    ):
        arguments.append(
            (ExonerateArgs(
                orthoid,
                transcripts_mapped_to[orthoid],
                args.min_score,
                aa_out_path,
                args.min_length,
                taxa_id,
                nt_out_path,
                tmp_path,
                args.verbose,
                gene_reference_data[orthoid]
            ),)
        )

    # this sorting the list so that the ones with the most hits are first
    if num_threads > 1:
        if num_threads > 1:
            with Pool(num_threads) as pool:
                recovered = pool.starmap(exonerate_gene_multi, arguments, chunksize=1)
    
    else:
        recovered = [exonerate_gene_multi(i[0]) for i in arguments]

    printv(f"Done! Took {time_keeper.differential():.2f}s overall. Exonerate took {time_keeper.lap():.2f}s and found {sum(recovered)} sequences.", args.verbose)

# TODO Make these argparse variables
strict_search_mode = False
orthoid_list_file = None
extend_orf = True
orf_overlap_minimum = 0.15
clear_output = True

header_seperator = "|"


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exist.", args.verbose, 0)
        return False
    rocky.create_pointer("rocks_orthoset_db", os.path.join(args.orthoset_input, args.orthoset, "rocksdb"))
    for input_path in args.INPUT:
        rocks_db_path = os.path.join(input_path, "rocksdb")
        rocky.create_pointer("rocks_nt_db", os.path.join(rocks_db_path, "sequences", "nt"))
        rocky.create_pointer("rocks_hits_db", os.path.join(rocks_db_path, "hits"))
        do_taxa(
            path=input_path,
            taxa_id=os.path.basename(input_path).split(".")[0],
            args=args,
        )
        rocky.close_pointer("rocks_nt_db")
        rocky.close_pointer("rocks_hits_db")
    rocky.close_pointer("rocks_orthoset_db")
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr Reporter"
    )
