from __future__ import annotations

import json
import os
import shutil
from collections import namedtuple
from multiprocessing.pool import Pool
from typing import Optional

import phymmr_tools
from Bio.Seq import Seq

from . import rocky
from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv, writeFasta

MainArgs = namedtuple(
    "MainArgs",
    [
        "verbose",
        "processes",
        "debug",
        "INPUT",
        "orthoset_input",
        "orthoset",
        "compress",
    ],
)


class Hit:
    __slots__ = (
        "header",
        "gene",
        "ali_start",
        "ali_end",
        "est_sequence",
        "ref_taxon",
    )

    def __init__(self, first_hit, gene):
        self.header = first_hit["header"]
        self.gene = gene
        self.ali_start = int(first_hit["ali_start"])
        self.ali_end = int(first_hit["ali_end"])
        self.ref_taxon = first_hit["ref_taxon"]

        self.est_sequence = first_hit["seq"]

    def trim_to_coords(self):
        self.est_sequence = self.est_sequence[self.ali_start - 1 : self.ali_end]
        if "revcomp" in self.header:
            self.est_sequence = phymmr_tools.bio_revcomp(self.est_sequence)


def get_diamondhits(rocks_hits_db, list_of_wanted_orthoids):
    gene_based_results = {}
    for gene in rocks_hits_db.get("getall:presentgenes").split(","):
        if list_of_wanted_orthoids and gene not in list_of_wanted_orthoids:
            continue

        gene_based_results[gene] = [
            Hit(this_data, gene)
            for this_data in json.loads(rocks_hits_db.get(f"gethits:{gene}"))
        ]

    return gene_based_results


def get_reference_data(rocks_hits_db):
    return json.loads(rocks_hits_db.get("getall:refseqs"))

def get_gene_variants(rocks_hits_db):
    return json.loads(rocks_hits_db.get("getall:target_variants"))

def get_toprefs(rocks_nt_db):
    return rocks_nt_db.get("getall:valid_refs").split(",")

def translate_cdna(cdna_seq):
    if not cdna_seq:
        return None

    if len(cdna_seq) % 3 != 0:
        printv("WARNING: NT Sequence length is not divisable by 3", 0)

    return str(Seq(cdna_seq).translate())


def get_ortholog_group(orthoid, orthoset_db):
    core_seqs = json.loads(orthoset_db.get(f"getcore:{orthoid}"))
    return core_seqs["aa"], core_seqs["nt"]


def format_candidate_header(gene, taxa_name, taxa_id, sequence_id, frame):
    return header_seperator.join([gene, taxa_name, taxa_id, sequence_id, frame])


def format_reference_header(gene, taxa_name, taxa_id, identifier="."):
    return header_seperator.join([gene, taxa_name, taxa_id, identifier])


def print_core_sequences(orthoid, core_sequences, target_taxon, top_refs):
    result = []
    for core in sorted(core_sequences):
        if target_taxon:
            if not core[1] in target_taxon:
                continue
        else:
            if not core[0] in top_refs:
                continue

        header = format_reference_header(orthoid, core[0], core[1])
        result.append((header, core[2]))

    return result


def print_unmerged_sequences(hits, orthoid, taxa_id):
    aa_result = []
    nt_result = []
    header_maps_to_where = {}
    header_mapped_x_times = {}
    base_header_mapped_already = {}
    seq_mapped_already = {}
    aa_reinsert_dupes = {}
    aa_rd_set = {}
    exact_hit_mapped_already = set()
    dupes = {}
    aa_dupe_count = 0
    for hit in hits:

        base_header, reference_frame = hit.header.split("|")

        header = format_candidate_header(
            orthoid,
            hit.ref_taxon,
            taxa_id,
            base_header,
            reference_frame,
        )
        nt_seq = hit.est_sequence
        aa_seq = translate_cdna(nt_seq)

        unique_hit = base_header + aa_seq

        if nt_seq in seq_mapped_already:
            mapped_to = seq_mapped_already[nt_seq]
            dupes.setdefault(mapped_to, []).append(base_header)
            continue
        seq_mapped_already[nt_seq] = base_header

        if aa_seq in aa_rd_set:
            aa_dupe_count += 1
            aa_reinsert_dupes.setdefault(aa_rd_set[aa_seq], []).append(header)
            continue
        else:
            aa_rd_set[aa_seq] = header

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
                        hit.ref_taxon,
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

    return dupes, aa_result, nt_result, aa_reinsert_dupes, aa_dupe_count


OutputArgs = namedtuple(
    "OutputArgs",
    [
        "gene",
        "list_of_hits",
        "aa_out_path",
        "taxa_id",
        "nt_out_path",
        "verbose",
        "reference_sequences",
        "compress",
        "target_taxon",
        "top_refs",

    ],
)


def trim_and_write(oargs: OutputArgs):
    t_gene_start = TimeKeeper(KeeperMode.DIRECT)
    printv(f"Doing output for: {oargs.gene}", oargs.verbose, 2)

    for sequence in oargs.list_of_hits:
        sequence.trim_to_coords()

    core_sequences, core_sequences_nt = get_ortholog_group(
        oargs.gene, rocky.get_rock("rocks_orthoset_db")
    )

    this_aa_path = os.path.join(oargs.aa_out_path, oargs.gene + ".aa.fa")
    this_gene_dupes, aa_output, nt_output, aa_reinsert_dupes, aa_dupe_count = print_unmerged_sequences(
        oargs.list_of_hits,
        oargs.gene,
        oargs.taxa_id,
    )

    if aa_output:
        aa_core_sequences = print_core_sequences(oargs.gene, core_sequences, oargs.target_taxon, oargs.top_refs)
        writeFasta(this_aa_path, aa_core_sequences + aa_output, oargs.compress)

        this_nt_path = os.path.join(oargs.nt_out_path, oargs.gene + ".nt.fa")

        nt_core_sequences = print_core_sequences(oargs.gene, core_sequences_nt, oargs.target_taxon, oargs.top_refs)
        writeFasta(this_nt_path, nt_core_sequences + nt_output, oargs.compress)

    printv(
        f"{oargs.gene} took {t_gene_start.differential():.2f}s. Had {len(aa_output)} sequences",
        oargs.verbose,
        2,
    )
    return oargs.gene, this_gene_dupes, len(aa_output), aa_reinsert_dupes, aa_dupe_count


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
    target_taxon = get_gene_variants(rocky.get_rock("rocks_hits_db"))
    top_refs = get_toprefs(rocky.get_rock("rocks_nt_db"))

    printv(
        f"Got reference data. Elapsed time {time_keeper.differential():.2f}s. Took {time_keeper.lap():.2f}s. Exonerating genes.",
        args.verbose,
    )

    arguments: list[Optional[OutputArgs]] = []
    for orthoid in sorted(
        transcripts_mapped_to, key=lambda k: len(transcripts_mapped_to[k]), reverse=True
    ):
        arguments.append(
            (
                OutputArgs(
                    orthoid,
                    transcripts_mapped_to[orthoid],
                    aa_out_path,
                    taxa_id,
                    nt_out_path,
                    args.verbose,
                    gene_reference_data[orthoid],
                    args.compress,
                    set(target_taxon.get(orthoid, [])),
                    top_refs,
                ),
            )
        )

    # this sorting the list so that the ones with the most hits are first
    if num_threads > 1:
        if num_threads > 1:
            with Pool(num_threads) as pool:
                recovered = pool.starmap(trim_and_write, arguments, chunksize=1)

    else:
        recovered = [trim_and_write(i[0]) for i in arguments]

    final_count = 0
    final_aa_dupes_count = 0
    this_gene_based_dupes = {}
    this_gene_inserts = {}
    for gene, nt_dupes, amount, aa_dupes, aa_dupes_amount in recovered:
        final_count += amount
        final_aa_dupes_count += aa_dupes_amount
        this_gene_based_dupes[gene] = nt_dupes
        this_gene_inserts[gene] = aa_dupes

    key = "getall:reporter_dupes"
    data = json.dumps(this_gene_based_dupes)
    rocky.get_rock("rocks_nt_db").put(key, data)

    key = "getall:insert_dupes"
    data = json.dumps(this_gene_inserts)
    rocky.get_rock("rocks_nt_db").put(key, data)

    printv(
        f"Done! Took {time_keeper.differential():.2f}s overall. Exonerate took {time_keeper.lap():.2f}s and found {final_count} sequences but kicked {final_aa_dupes_count} aa dupes.",
        args.verbose,
    )


# TODO Make these argparse variables
orthoid_list_file = None
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
