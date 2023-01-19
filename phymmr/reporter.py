from __future__ import annotations

import json
import math
import os
import shutil
from collections import namedtuple
from multiprocessing.pool import Pool
from typing import Optional
from tempfile import TemporaryDirectory, NamedTemporaryFile

import phymmr_tools
import xxhash
from Bio.Seq import Seq

from . import rocky
from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv, gettempdir, writeFasta

MainArgs = namedtuple(
    "MainArgs",
    [
        "verbose",
        "processes",
        "debug",
        "INPUT",
        "orthoset_input",
        "orthoset",
        "min_length",
        "min_score",
        "compress",
    ],
)


class Hit:
    __slots__ = (
        "header",
        "gene",
        "score",
        "ali_start",
        "ali_end",
        "s_ali_start",
        "s_ali_end",
        "f_est_sequence_trimmed",
        "s_est_sequence_trimmed",
        "est_sequence_complete",
        "reftaxon",
        "f_ref_taxon",
        "s_ref_taxon",
        "second_alignment",
        "second_extended_alignment",
        "first_alignment",
        "first_extended_alignment",
        "mapped_to",
        "attempted_first",
    )

    def __init__(self, first_hit, second_hit):
        self.header = first_hit["header"]
        self.gene = first_hit["gene"]
        self.ali_start = int(first_hit["ali_start"])
        self.ali_end = int(first_hit["ali_end"])

        self.f_ref_taxon = first_hit["ref_taxon"]
        if second_hit:
            self.s_ref_taxon = second_hit["ref_taxon"]
            self.s_est_sequence_trimmed = second_hit["trim_seq"]
            self.s_ali_start = int(second_hit["ali_start"])
            self.s_ali_end = int(second_hit["ali_end"])
        else:
            self.s_ref_taxon = None
            self.s_est_sequence_trimmed = None
            self.s_ali_start = None
            self.s_ali_end = None

        self.reftaxon = first_hit["ref_taxon"]
        self.est_sequence_complete = first_hit["full_seq"]
        self.f_est_sequence_trimmed = first_hit["trim_seq"]
        self.second_alignment = None
        self.second_extended_alignment = None
        self.first_alignment = None
        self.first_extended_alignment = None
        self.mapped_to = None
        self.attempted_first = False

    def add_orf(self, exonerate_record):
        exonerate_record.orf_cdna_start_on_transcript = (
            exonerate_record.orf_cdna_start + self.ali_start - 3
        )
        exonerate_record.orf_cdna_end_on_transcript = (
            exonerate_record.orf_cdna_end + self.ali_start - 3
        )
        exonerate_record.orf_aa_start_on_transcript = (
            exonerate_record.orf_cdna_start + self.ali_start - 3
        ) / 3
        exonerate_record.orf_aa_end_on_transcript = (
            exonerate_record.orf_cdna_end + self.ali_start - 3
        ) / 3
        if not self.first_alignment:
            self.first_alignment = exonerate_record
        else:
            self.second_alignment = exonerate_record


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
        self.extended_orf_cdna_start_on_transcript = None
        self.extended_orf_cdna_end_on_transcript = None

    def check_extend(self):
        if self.is_extension:
            self.extended_orf_cdna_sequence = self.orf_cdna_sequence
            self.extended_orf_cdna_start = self.orf_cdna_start
            self.extended_orf_cdna_end = self.orf_cdna_end
            self.extended_orf_aa_sequence = self.orf_aa_sequence
            self.extended_orf_aa_start = self.orf_aa_start
            self.extended_orf_aa_end = self.orf_aa_end

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


def get_diamondhits(rocks_hits_db, list_of_wanted_orthoids):
    gene_based_results = {}
    for gene in rocks_hits_db.get("getall:presentgenes").split(","):
        if list_of_wanted_orthoids and gene not in list_of_wanted_orthoids:
            continue

        gene_based_results[gene] = [
            Hit(this_data["f"], this_data["s"])
            for this_data in json.loads(rocks_hits_db.get(f"gethits:{gene}"))
        ]

    return gene_based_results


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
        return reverse_complement(sequence)
    return sequence


def crop_to_alignment(seq, hit):
    start = hit.ali_start - 1  # Adjust for zero based number
    end = hit.ali_end

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
            extended_result.append((node))
        else:
            node.orf_aa_sequence = translate_cdna(node.orf_cdna_sequence)
            result.append((node))
    return extended_result, result


def get_multi_orf(query, targets, score_threshold, include_extended, taxon):
    exonerate_ryo = ">cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas"
    genetic_code = 1
    exonerate_model = "protein2genome"

    sequences = []
    for i in targets:
        this_sequence = (
            i.f_est_sequence_trimmed
            if taxon == i.f_ref_taxon
            else i.s_est_sequence_trimmed
        )
        sequences.append((i.header, this_sequence))
    if include_extended:
        sequences.extend(
            [("extended_" + i.header, i.est_sequence_complete) for i in targets]
        )

    with TemporaryDirectory(dir=gettempdir()) as tmpdir, NamedTemporaryFile(
        dir=tmpdir, mode="w+"
    ) as tmpquery, NamedTemporaryFile(
        dir=tmpdir, mode="w+"
    ) as tmptarget, NamedTemporaryFile(
        dir=tmpdir, mode="r"
    ) as tmpout:
        tmptarget.write(
            "".join([f">{header}\n{sequence}\n" for header, sequence in sequences])
        )
        tmptarget.flush()  # Flush the internal buffer so it can be read by exonerate

        tmpquery.write(f">query\n{query}\n")
        tmpquery.flush()  # Flush the internal buffer so it can be read by exonerate

        exonerate_cmd = f"exonerate --score {score_threshold} --ryo '{exonerate_ryo}' --subopt 0 --geneticcode {genetic_code} --model '{exonerate_model}' --querytype 'protein' --targettype 'dna' --verbose 0 --showalignment 'no' --showvulgar 'yes' --query '{tmpquery.name}' --target '{tmptarget.name}' > {tmpout.name}"

        os.system(exonerate_cmd)

        extended_results, results = parse_multi_results(tmpout)

    return extended_results, results


def extended_orf_contains_original_orf(extended_alignment, alignment):
    return (
        extended_alignment.extended_orf_cdna_start
        >= alignment.orf_cdna_end_on_transcript
        or extended_alignment.extended_orf_cdna_end
        >= alignment.orf_cdna_start_on_transcript
    )


def get_overlap_length(candidate, alignment):
    if (
        alignment.extended_orf_aa_start <= candidate.ali_end
        and alignment.extended_orf_aa_end >= candidate.ali_start
    ):
        overlap_start = (
            alignment.extended_orf_aa_start
            if alignment.extended_orf_aa_start > candidate.ali_start
            else candidate.ali_start
        )
        overlap_end = (
            alignment.extended_orf_aa_end
            if alignment.extended_orf_aa_end > candidate.ali_end
            else candidate.ali_end
        )

        overlap_length = overlap_end - overlap_start
        return overlap_length
    return 0


def overlap_by_orf(candidate, alignment):
    orf_length = abs(alignment.extended_orf_aa_end - alignment.extended_orf_aa_start)
    overlap_length = min(alignment.extended_orf_cdna_end, candidate.ali_end) - max(
        alignment.extended_orf_cdna_start, candidate.ali_start
    )
    if overlap_length < 0:
        overlap_length = 0

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
    result = []
    for core in sorted(core_sequences):
        header = format_reference_header(orthoid, core[0], core[1])
        result.append((header, core[2]))

    return result


def print_unmerged_sequences(hits, orthoid, taxa_id, minimum_length):
    aa_result = []
    nt_result = []
    header_maps_to_where = {}
    header_mapped_x_times = {}
    base_header_mapped_already = {}
    exact_hit_mapped_already = set()
    for hit in hits:

        base_header, reference_frame = hit.header.split("|")

        header = format_candidate_header(
            orthoid,
            hit.reftaxon,
            taxa_id,
            base_header,
            reference_frame,
        )

        nt_seq = hit.s_est_sequence_trimmed if len(hit.f_est_sequence_trimmed) < minimum_length else hit.f_est_sequence_trimmed
        if len(nt_seq) < minimum_length:
            continue
        aa_seq = translate_cdna(nt_seq)

        unique_hit = base_header + aa_seq
        if unique_hit not in exact_hit_mapped_already:

            if base_header in base_header_mapped_already:
                (
                    already_mapped_header,
                    already_mapped_sequence,
                ) = base_header_mapped_already[base_header]

                if len(aa_seq) > len(already_mapped_sequence):
                    if already_mapped_sequence in aa_seq:
                        aa_result[header_maps_to_where[already_mapped_header]] = (
                            header,
                            aa_seq,
                        )
                        nt_result[header_maps_to_where[already_mapped_header]] = (
                            header,
                            nt_seq,
                        )
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
                        base_header + f"_{header_mapped_x_times[old_header]}",
                        reference_frame,
                    )

                    header_mapped_x_times[base_header] += 1
            else:
                base_header_mapped_already[base_header] = header, aa_seq

            header_maps_to_where[header] = len(
                aa_result
            )  # Save the index of the sequence output
            aa_result.append((header, aa_seq))
            nt_result.append((header, nt_seq))

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
        if result.header not in res:
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
        "compress",
    ],
)


def exonerate_gene_multi(eargs: ExonerateArgs):

    t_gene_start = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Exonerating and doing output for: {eargs.orthoid}", eargs.verbose, 2)

    output_sequences = eargs.list_of_hits

    if len(output_sequences) > 0:
        core_sequences, core_sequences_nt = get_ortholog_group(
            eargs.orthoid, rocky.get_rock("rocks_orthoset_db")
        )
        this_aa_path = os.path.join(eargs.aa_out_path, eargs.orthoid + ".aa.fa")
        aa_output, nt_output = print_unmerged_sequences(
            output_sequences,
            eargs.orthoid,
            eargs.taxa_id,
            eargs.min_length
        )

        if aa_output:
            aa_core_sequences = print_core_sequences(eargs.orthoid, core_sequences)
            writeFasta(this_aa_path, aa_core_sequences + aa_output, eargs.compress)

            this_nt_path = os.path.join(eargs.nt_out_path, eargs.orthoid + ".nt.fa")

            nt_core_sequences = print_core_sequences(eargs.orthoid, core_sequences_nt)
            writeFasta(this_nt_path, nt_core_sequences + nt_output, eargs.compress)

    printv(
        f"{eargs.orthoid} took {t_gene_start.differential():.2f}s. Had {len(output_sequences)} sequences",
        eargs.verbose,
        2,
    )
    return len(output_sequences)


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

    printv(
        f"Initialized databases. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reciprocal diamond hits.",
        args.verbose,
    )

    transcripts_mapped_to = get_diamondhits(
        rocky.get_rock("rocks_hits_db"), list_of_wanted_orthoids
    )

    printv(
        f"Got hits. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Grabbing reference data.",
        args.verbose,
    )

    gene_reference_data = get_reference_data(rocky.get_rock("rocks_orthoset_db"))

    printv(
        f"Got reference data. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Exonerating genes.",
        args.verbose,
    )

    arguments: list[Optional[ExonerateArgs]] = []
    for orthoid in sorted(
        transcripts_mapped_to, key=lambda k: len(transcripts_mapped_to[k]), reverse=True
    ):
        arguments.append(
            (
                ExonerateArgs(
                    orthoid,
                    transcripts_mapped_to[orthoid],
                    args.min_score,
                    aa_out_path,
                    args.min_length,
                    taxa_id,
                    nt_out_path,
                    tmp_path,
                    args.verbose,
                    gene_reference_data[orthoid],
                    args.compress,
                ),
            )
        )

    # this sorting the list so that the ones with the most hits are first
    if num_threads > 1:
        if num_threads > 1:
            with Pool(num_threads) as pool:
                recovered = pool.starmap(exonerate_gene_multi, arguments, chunksize=1)

    else:
        recovered = [exonerate_gene_multi(i[0]) for i in arguments]

    printv(
        f"Done! Took {time_keeper.differential():.2f}s overall. Exonerate took {time_keeper.lap():.2f}s and found {sum(recovered)} sequences.",
        args.verbose,
    )


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
    rocky.create_pointer(
        "rocks_orthoset_db", os.path.join(args.orthoset_input, args.orthoset, "rocksdb")
    )
    for input_path in args.INPUT:
        rocks_db_path = os.path.join(input_path, "rocksdb")
        rocky.create_pointer(
            "rocks_nt_db", os.path.join(rocks_db_path, "sequences", "nt")
        )
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
