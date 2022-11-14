from __future__ import annotations

import itertools
import json
import math
import os
import subprocess
from multiprocessing.pool import Pool
from time import time

import wrap_rocks
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .utils import printv
from .timekeeper import TimeKeeper, KeeperMode

class Hit:
    __slots__ = (
        "header",
        "base_header",
        "gene",
        "evalue",
        "score",
        "hmm_start",
        "hmm_end",
        "ali_start",
        "ali_end",
        "env_start",
        "env_end",
        "uuid",
        "hmm_sequence",
        "hmm_id",
    )

    def __init__(
        self,
        target,
        query,
        evalue,
        score,
        hmm_start,
        hmm_end,
        ali_start,
        ali_end,
        env_start,
        env_end,
    ):
        self.header = str(target).replace(">", "").replace(" ", "|")
        self.base_header = self.header.split('|')[0].strip()
        self.gene = str(query)
        self.evalue = float(evalue)
        self.score = float(score)
        self.hmm_start = int(hmm_start)
        self.hmm_end = int(hmm_end)
        self.ali_start = int(ali_start)
        self.ali_end = int(ali_end)
        self.env_start = int(env_start)
        self.env_end = int(env_end)
        self.uuid = None
        self.hmm_sequence = None
        self.hmm_id = None

    def __str__(self):
        return "\t".join(
            [
                self.header,
                self.gene,
                self.evalue,
                self.score,
                self.hmm_start,
                self.hmm_end,
                self.ali_start,
                self.ali_end,
                self.env_start,
                self.env_end,
            ]
        )

    def __repr__(self) -> str:
        return f"{self.header} {self.score}"

    def list_values(self, species_id):
        """
        Given a species_id, returns a list of values suitable for the
        hmmsearch table insertion.
        """
        log_evalue = -999
        e_float = float(self.evalue)
        if e_float != 0:
            log_evalue = math.log(e_float)
        return [
            str(species_id),
            self.gene,
            self.header,
            self.score,
            str(e_float),
            str(log_evalue),
            self.hmm_start,
            self.hmm_end,
            self.ali_start,
            self.ali_end,
            self.env_start,
            self.env_end,
        ]

    def to_json(self):
        return {
            "header": self.header,
            "gene": self.gene,
            "score": self.score,
            "hmm_evalue": self.evalue,
            "hmm_start": self.hmm_start,
            "hmm_end": self.hmm_end,
            "ali_start": self.ali_start,
            "ali_end": self.ali_end,
            "env_start": self.env_start,
            "env_end": self.env_end,
            "uuid": self.uuid,
            "hmm_id": self.hmm_id,
            "hmm_sequence": self.hmm_sequence
        }


def get_difference(scoreA, scoreB):
    """
    Returns decimal difference of two scores
    """
    try:
        if scoreA / scoreB > 1:
            return scoreA / scoreB
        if scoreB / scoreA > 1:
            return scoreB / scoreA
        if scoreA == scoreB:
            return 1
    except ZeroDivisionError:
        return 0


def get_overlap(a_start, a_end, b_start, b_end):
    overlap_end = min(a_end, b_end)
    overlap_start = max(a_start, b_start)
    amount = (overlap_end - overlap_start) + 1  # inclusive

    return 0 if amount < 0 else amount


def internal_filter_gene(this_gene_hits, gene, min_overlap_internal, score_diff_internal, debug, requires_sort):
    if requires_sort:
        this_gene_hits.sort(key=lambda hit: hit.score, reverse=True)

    filtered_sequences_log = []

    for i, hit_a in enumerate(this_gene_hits):
        if not hit_a:
            continue
        for j in range(len(this_gene_hits) - 1, i, -1):
            hit_b = this_gene_hits[j]
            if hit_b:
                if hit_a.base_header != hit_b.base_header:
                    if ((hit_a.score / hit_b.score) if hit_b.score != 0 else 0) < score_diff_internal:
                        break

                    amount_of_overlap = get_overlap(hit_a.hmm_start, hit_a.hmm_end, hit_b.hmm_start, hit_b.hmm_end)
                    distance = (hit_b.hmm_end - hit_b.hmm_start) + 1  # Inclusive
                    percentage_of_overlap = amount_of_overlap / distance

                    if percentage_of_overlap >= min_overlap_internal:
                        this_gene_hits[j] = None
                        if debug:
                            filtered_sequences_log.append(
                                [
                                    hit_b.gene, hit_b.header, str(hit_b.score), str(hit_b.hmm_start),
                                    str(hit_b.hmm_end), 'Internal Overlapped with Lowest Score', hit_a.gene,
                                    hit_a.header, str(hit_a.score), str(hit_a.hmm_start), str(hit_a.hmm_end)
                                ]
                            )

    this_out_data = {'Passes': [i for i in this_gene_hits if i is not None], 'Log': filtered_sequences_log,
                     'Gene': gene}
    return this_out_data


def multi_filter_dupes(
    this_hits,
    debug,
    min_overlap_multi,
    score_diff_multi,
):
    kick_happend = True
    filtered_sequences_log = []

    this_hits.sort(key=lambda data: (data.score, data.gene), reverse=True)

    while kick_happend:
        kick_happend = False

        master = this_hits[0]
        # candidates = [i for i in this_hits[1:] if i is not None]
        candidates = [i for i in this_hits[1:]]

        master_env_start = master.env_start
        master_env_end = master.env_end

        for i, candidate in enumerate(candidates, 1):
            if candidate is None:  # guard statement
                continue
            if candidate.gene == master.gene:  # From same gene = Pseudomaster
                # Remove
                if debug:
                    filtered_sequences_log.append(
                        [
                            candidate.gene,
                            candidate.header,
                            str(candidate.score),
                            str(candidate.hmm_start),
                            str(candidate.hmm_end),
                            "Pseudomaster",
                            master.gene,
                            master.header,
                            str(master.score),
                            str(master.hmm_start),
                            str(master.hmm_end),
                        ]
                    )

                this_hits[i] = None
                candidates[i - 1] = None
                # Extend master range
                kick_happend = True
                if candidate.env_start < master_env_start:
                    master_env_start = candidate.env_start
                if candidate.env_end > master_env_start:
                    master_env_end = candidate.env_end
            else:
                break

        miniscule_score = False

        for i, candidate in enumerate(candidates, 1):
            if candidate:
                distance = (master_env_end - master_env_start) + 1  # Inclusive
                amount_of_overlap = get_overlap(master_env_start, master_env_end, candidate.env_start,
                                                candidate.env_end)
                percentage_of_overlap = amount_of_overlap / distance

                if percentage_of_overlap >= min_overlap_multi:
                    score_difference = get_difference(master.score, candidate.score)
                    if score_difference >= score_diff_multi:
                        kick_happend = True
                        this_hits[i] = None
                        if debug:
                            filtered_sequences_log.append(
                                [
                                    candidate.gene,
                                    candidate.header,
                                    str(candidate.score),
                                    str(candidate.env_start),
                                    str(candidate.env_end),
                                    "Multi Overlapped with Lowest Score",
                                    master.gene,
                                    master.header,
                                    str(master.score),
                                    str(master.env_start),
                                    str(master.env_end),
                                ]
                            )
                    else:
                        miniscule_score = True
                        break

        if (
            miniscule_score
        ):  # Remove all overlapping candidates if it's score is a miniscule difference of the masters
            for i, candidate in enumerate(candidates, 1):
                if candidate:
                    distance = (master_env_end - master_env_start) + 1  # Inclusive
                    amount_of_overlap = get_overlap(master_env_start, master_env_end, candidate.env_start,
                                                    candidate.env_end)
                    percentage_of_overlap = amount_of_overlap / distance
                    if percentage_of_overlap >= min_overlap_multi:
                        kick_happend = True
                        this_hits[i] = None
                        if debug:
                            filtered_sequences_log.append(
                                [
                                    candidate.gene,
                                    candidate.header,
                                    str(candidate.score),
                                    str(candidate.hmm_start),
                                    str(candidate.hmm_end),
                                    "Multi Overlapped with Miniscule Score",
                                    master.gene,
                                    master.header,
                                    str(master.score),
                                    str(master.hmm_start),
                                    str(master.hmm_end),
                                ]
                            )
    multi_data = {
        "Log": filtered_sequences_log,
        "Remaining": [i for i in this_hits if i is not None],
    }
    return multi_data


def internal_multi_filter(flagged_headers, this_gene_hits, minimum_overlap_multi_internal, debug, gene):
    this_gene_hits.sort(key=lambda hit: hit.score, reverse=True)

    bh_based_results = {}
    for i, hit in enumerate(this_gene_hits):  # Iterate once and make hashmap
        if not hit.base_header in bh_based_results:
            bh_based_results[hit.base_header] = [i]
        else:
            bh_based_results[hit.base_header].append(i)

    filtered_sequences_log = []

    for b_header in flagged_headers:
        if b_header in bh_based_results:  # Not kicked during multi
            bh_hits = bh_based_results[b_header]

            for hit_a_index, hit_b_index in itertools.combinations(bh_hits, 2):
                hit_a = this_gene_hits[hit_a_index]
                if not hit_a:
                    continue

                hit_b = this_gene_hits[hit_b_index]
                if not hit_b:
                    continue

                overlap_amount = get_overlap(hit_a.ali_start, hit_a.ali_end, hit_b.ali_start, hit_b.ali_end)
                distance = (hit_b.ali_end - hit_b.ali_start) + 1

                overlap_percent = overlap_amount / distance

                if overlap_percent > minimum_overlap_multi_internal:

                    this_gene_hits[hit_b_index] = None
                    if debug:
                        filtered_sequences_log.append(
                            [hit_b.gene, hit_b.header, str(hit_b.score), str(hit_b.ali_start), str(hit_b.ali_end),
                             'Internal Multi Overlapped with Lowest Score', hit_a.gene, hit_a.header, str(hit_a.score),
                             str(hit_a.ali_start), str(hit_a.ali_end)])

    this_out_data = {'Passes': [i for i in this_gene_hits if i is not None], 'Log': filtered_sequences_log,
                     'Gene': gene}

    return this_out_data


def get_hmm_name(hmmfile: str) -> str:
    """
    The name is always the second line of the hmm file. This function reads the second
    line of the file and splits the line. The second field is the name.
    Returns the name as a string.
    """
    with open(hmmfile) as hf:
        hf.readline()
        name_line = hf.readline()
        name = name_line.split()[1].strip()
        return name


def make_exclusion_list(path: str) -> set:
    """
    Use on a line delimited text file to make the exclusion list.
    Reads the file at the given path line by line and stores all the values
    as strings in a set. Returns the set.
    :param path: A path string to the excluded text file
    :return: A set containing the names of all excluded orthologs
    """
    excluded = set()
    with open(path) as infile:
        for line in infile:
            excluded.add(line.strip)
    return excluded


def make_temp_protfile(
    ortholog: str,
    temp: str,
    sequences_db_conn
) -> str:
    """
    Looks up the digest and sequence in the orthodb est table, then
    writes to the protfile. Returns path to protfile as a str.
    Creates side effects.
    """
    prot_name = ortholog + "_prot.tmp"
    prot_path = os.path.join(temp, prot_name)

    recipe = sequences_db_conn.get("getall:prot").split(',')

    with open(prot_path, "w") as prot_file_handle:
        for component in recipe:
            prot_file_handle.write(sequences_db_conn.get(f"getprot:{component}"))

    return prot_path


def parse_domtbl_fields(fields: list) -> Hit:
    hit = Hit(
        fields[0],
        fields[3],
        fields[12],
        fields[13],
        fields[15],
        fields[16],
        fields[17],
        fields[18],
        fields[19],
        fields[20],
    )
    return hit


def domtbl_dupe_check(domtbl_path: str) -> None:
    """Reads the given domtbl file and removes duplicate data rows.
    Overwrites old file with new version."""
    data_lines_seen = {}
    output = []
    with open(domtbl_path, 'r+') as f:
        for line in f:
            if line[0] == "#":
                output.append(line)
                continue
            seen_before = data_lines_seen.get(line, False)
            if not seen_before:
                output.append(line)
                data_lines_seen[line] = True
        f.seek(0)
        f.writelines(output)
        f.truncate()


def get_hits_from_domtbl(domtbl_path: str, score, evalue) -> list:
    domtbl_dupe_check(domtbl_path)
    hits = []
    with open(domtbl_path) as f:
        for line in f:
            if line and line[0] == "#":
                continue
            fields = line.strip().split()

            if float(fields[13]) >= score:
                hits.append(parse_domtbl_fields(fields))
    return hits


def empty_domtbl_file(path: str) -> bool:
    """
    Aborted runs result in a domtbl file with 0 bytes. If such a file is
    found, return True. Otherwise returns False.
    If an empty domtbl already exists, always overwrite.
    """
    return os.path.getsize(path) == 0


def search_prot(
    prot_file: str,
    domtbl_path: str,
    hmm_file: str,
    evalue: float,
    score: float,
    ovw: bool,
    verbose: bool,
    prog="hmmsearch",
    threads=1,
) -> list:
    if os.path.exists(domtbl_path) and not ovw:
        # always overwrite the empty domtbl files
        if not empty_domtbl_file(domtbl_path):
            return get_hits_from_domtbl(domtbl_path, score, evalue)
    if score:
        option = "-T"
        threshold = str(score)
    else:
        option = "-E"
        threshold = str(evalue)
    # $searchprog --domtblout $outfile $threshold_option --cpu $num_threads $hmmfile $protfile
    command = [
        prog,
        "--domtblout",
        domtbl_path,
        option,
        threshold,
        "--cpu",
        str(threads),
        hmm_file,
        prot_file,
    ]

    p = subprocess.run(command, stdout=subprocess.PIPE)
    if p.returncode != 0:  # non-zero return code means an error
        printv(f"{domtbl_path}:hmmsearch error code {p.returncode}", verbose, 0)
    else:
        printv(f"Searched {os.path.basename(domtbl_path)}", verbose, 2)
    return get_hits_from_domtbl(domtbl_path, score, evalue)


def hmm_search(
    hmm_file: str, domtbl_dir: str, evalue, score, prot: str, ovw: bool, verbose: bool
) -> None:
    """
    Reimplements hmmsearch loop in lines 468 to 538 in orthograph analyzer.
    """
    # hmmfile contains dir in path at this point
    hmm_name = get_hmm_name(hmm_file)
    domtbl_path = os.path.join(domtbl_dir, hmm_name + ".hmm.domtbl")

    hits = search_prot(prot, domtbl_path, hmm_file, evalue, score, ovw, verbose)
    return hits


def run_process(args, input_path: str) -> None:
    """
    Runs the hmmsearch process on an individual input path.
    Allows the use of nargs in main().
    """
    MAX_HMM_BATCH_SIZE = args.max_hmm_batch_size
    debug = args.debug != 0
    score_diff_multi = args.score_diff_multi
    min_overlap_multi = args.min_overlap_multi
    score_diff_internal = args.score_diff_internal
    min_overlap_internal = args.min_overlap_internal

    time_keeper = TimeKeeper(KeeperMode.DIRECT)

    num_threads = args.processes
    if os.path.exists("/run/shm"):
        in_ram = "/run/shm"
    else:
        in_ram = "/dev/shm"
    # temp_dir = os.path.join(input_path, "tmp")
    temp_dir = os.path.join(in_ram, "tmp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # get ortholog
    ortholog = os.path.basename(input_path).split(".")[0]

    # if database not set, make expected path to database
    db_path = os.path.join(input_path, "rocksdb")
    sequences_path = os.path.join(db_path, "sequences")
    hits_path = os.path.join(db_path, "hits")

    sequences_db_conn = wrap_rocks.RocksDB(sequences_path)

    # path to domtbl directory
    domtbl_dir = os.path.join(input_path, "hmmsearch")

    os.makedirs(domtbl_dir, exist_ok=True)

    # path to .hmm file directory
    orthoset_dir = os.path.join(args.orthoset_input, args.orthoset)
    hmm_dir = os.path.join(orthoset_dir, "hmms")

    # make a set of excluded orthologs if argument was included
    excluded = set()
    if args.excluded_list:
        excluded = make_exclusion_list(args.excluded_list)

    # make a set of wanted orthologs if argument was included
    wanted = set()
    if args.wanted_list:
        wanted = make_exclusion_list(args.wanted_list)

    # make a list of valid ortholog names, excluded hidden files
    printv("Finding hmm files", args.verbose)
    hmm_list = [
        hmm for hmm in os.listdir(hmm_dir) if ".hmm" in hmm and hmm not in excluded
    ]
    if wanted:
        hmm_list = [hmm for hmm in hmm_list if hmm.split(".")[0] in wanted]

    hmm_list.sort()

    # rejoin file names with directory path
    hmm_list = [os.path.join(hmm_dir, hmm) for hmm in hmm_list]

    # make protfile for hmm_search later
    printv("Creating protfile", args.verbose)
    time_keeper.lap()

    protfile = make_temp_protfile(
        ortholog, temp_dir, sequences_db_conn
    )

    printv(f"Wrote prot file in {time_keeper.lap():.2f}s", args.verbose)


    hmm_results = []  # list of lists containing Hit objects

    if num_threads > 1:
        arg_tuples = []
        for hmm in hmm_list:
            arg_tuples.append(
                (hmm, domtbl_dir, args.evalue, args.score, protfile, args.overwrite, args.verbose,)
            )
        with Pool(num_threads) as search_pool:
            hmm_results = search_pool.starmap(hmm_search, arg_tuples)
    else:
        for hmm in hmm_list:
            hmm_results.append(
                hmm_search(hmm, domtbl_dir, args.evalue, args.score, protfile, args.overwrite, args.verbose))

    printv(f"Search time: {time_keeper.lap():.2f}", args.verbose)

    filter_timer = TimeKeeper(KeeperMode.DIRECT)

    f_duplicates = {}
    gene_based_results = {}
    header_based_results = {}
    count = 0

    # Disperse hits into genes
    for hit_group in hmm_results:
        for hit in hit_group:
            count += 1
            this_gene = hit.gene
            if this_gene not in gene_based_results:
                gene_based_results[this_gene] = []

            hit.uuid = hit.header + f"_hit_{count}"

            gene_based_results[this_gene].append(hit)

    printv(f"Filtering {count} hits", args.verbose)

    printv("Doing multi-gene dupe filter. Searching for duplicates", args.verbose)
    required_internal_multi_genes = {}
    rimg_set = set()

    for gene in gene_based_results:
        this_gene_baseheaders = set()
        requires_internal_multi_filter = {}
        this_requires = False

        for hit in gene_based_results[gene]:
            if "revcomp" in hit.header:
                base_header = hit.base_header
                raw_length = int(base_header.split("_length_")[1]) / 3
                length = math.floor(raw_length)

                new_env_start = length - int(hit.env_end)
                new_env_end = length - int(hit.env_start)
                hit.env_start = new_env_start
                hit.env_end = new_env_end

            if hit.header not in f_duplicates:
                f_duplicates[hit.header] = []
            if hit.header not in header_based_results:
                header_based_results[hit.header] = []
            if hit.base_header not in this_gene_baseheaders:
                this_gene_baseheaders.add(hit.base_header)
            else:
                requires_internal_multi_filter[hit.base_header] = True  # As a dict to avoid dupe headers
                this_requires = True

            f_duplicates[hit.header].append(hit)
            header_based_results[hit.header].append(hit)

        if this_requires:
            required_internal_multi_genes[gene] = list(requires_internal_multi_filter.keys())
            rimg_set.add(gene)  # Used for internal sort

    headers = list(f_duplicates.keys())
    possible_dupes = 0
    for header in headers:
        if len(f_duplicates[header]) > 1:
            unique_genes = list(
                dict.fromkeys([i.gene for i in f_duplicates[header]])
            )
            if len(unique_genes) <= 1:  # But the same gene
                f_duplicates.pop(header)
            else:
                possible_dupes += len(f_duplicates[header])
        else:
            f_duplicates.pop(header)

    for header_left in f_duplicates:
        header_based_results[header_left] = []

    printv(f"Found {possible_dupes} potential dupes. Search took {filter_timer.lap():.2f}s. Filtering dupes", args.verbose)

    if debug:
        filtered_sequences_log = []

    if num_threads == 1:
        for header in f_duplicates:
            this_hits = f_duplicates[header]
            data = multi_filter_dupes(
                this_hits,
                debug,
                min_overlap_multi,
                score_diff_multi,
            )

            if debug:
                filtered_sequences_log.extend(data["Log"])
            header_based_results[data["Remaining"][0].header] = data["Remaining"]

    else:
        arguments = []
        for header in f_duplicates:
            this_hits = f_duplicates[header]
            arguments.append(
                (
                    this_hits,
                    debug,
                    min_overlap_multi,
                    score_diff_multi,
                )
            )

        with Pool(num_threads) as pool:
            multi_data = pool.starmap(multi_filter_dupes, arguments, chunksize=1)

        for data in multi_data:
            if debug:
                filtered_sequences_log.extend(data["Log"])
            header_based_results[data["Remaining"][0].header] = data["Remaining"]

    printv(f"Done! Took {filter_timer.lap():.2f}s", args.verbose)
    transcripts_mapped_to = {}

    for header in header_based_results:
        for match in header_based_results[header]:
            if match.gene not in transcripts_mapped_to:
                transcripts_mapped_to[match.gene] = []
            transcripts_mapped_to[match.gene].append(match)

    printv(f"Looking for internal duplicates over {len(required_internal_multi_genes)} flagged genes", args.verbose)

    if num_threads == 1:
        internal_multi_data = []
        for gene in required_internal_multi_genes:
            if gene in transcripts_mapped_to:
                this_gene_transcripts = transcripts_mapped_to[gene]
                check_headers = required_internal_multi_genes[gene]
                internal_multi_data.append(  #
                    internal_multi_filter(
                        check_headers,
                        this_gene_transcripts,
                        args.minimum_overlap_internal_multi,
                        debug,
                        gene,
                    )
                )

    else:
        arguments = []
        for gene in required_internal_multi_genes:
            if gene in transcripts_mapped_to:
                this_gene_transcripts = transcripts_mapped_to[gene]
                check_headers = required_internal_multi_genes[gene]
                arguments.append(
                    (
                        check_headers,
                        this_gene_transcripts,
                        args.minimum_overlap_internal_multi,
                        debug,
                        gene,
                    )
                )
        with Pool(num_threads) as pool:
            internal_multi_data = pool.starmap(internal_multi_filter, arguments, chunksize=1)

    for data in internal_multi_data:
        gene = data["Gene"]

        transcripts_mapped_to[gene] = data["Passes"]
        if debug:
            filtered_sequences_log.extend(data["Log"])

    printv(f"Done! Took {filter_timer.lap():.2f}s", args.verbose)
    printv("Doing internal overlap filter", args.verbose)

    total_hits = 0
    if num_threads == 1:
        internal_data = []
        for gene in transcripts_mapped_to:
            this_gene_transcripts = transcripts_mapped_to[gene]
            sort = gene not in rimg_set
            internal_data.append(
                internal_filter_gene(
                    this_gene_transcripts,
                    gene,
                    min_overlap_internal,
                    score_diff_internal,
                    debug,
                    sort
                )
            )
    else:
        arguments = []
        for gene in transcripts_mapped_to:
            this_gene_transcripts = transcripts_mapped_to[gene]
            sort = gene not in rimg_set
            arguments.append(
                (
                    this_gene_transcripts,
                    gene,
                    min_overlap_internal,
                    score_diff_internal,
                    debug,
                    sort
                )
            )
        with Pool(num_threads) as pool:
            internal_data = pool.starmap(internal_filter_gene, arguments, chunksize=1)

    transcripts_mapped_to = {}
    for data in internal_data:
        gene = data["Gene"]

        if gene not in transcripts_mapped_to:
            transcripts_mapped_to[gene] = []

        transcripts_mapped_to[gene] = data["Passes"]
        total_hits += len(data["Passes"])
        if debug:
            filtered_sequences_log.extend(data["Log"])

    if debug:
        filtered_sequences_log_out = []
        for line in filtered_sequences_log:
            filtered_sequences_log_out.append(",".join(line) + "\n")

        filtered_sequences_log_path = os.path.join(input_path, "filtered-hits.csv")
        with open(filtered_sequences_log_path, "w") as fp:
            fp.writelines(filtered_sequences_log_out)

    printv(f"Done! Took {filter_timer.lap():.2f}s", args.verbose)
    printv(f"Filtering took {filter_timer.differential():.2f}s overall", args.verbose)

    # Grab the seq dict
    sequence_dict = {}
    header = None
    time_keeper.lap()
    with open(protfile) as prot_file_handle:
        fasta_file = SimpleFastaParser(prot_file_handle)
        for header, sequence in fasta_file:
            sequence_dict[header] = sequence

    if os.path.exists(protfile):
        os.remove(protfile)

    printv(f"Read prot file in {time_keeper.lap():.2f}s", args.verbose)

    # Seperate hmm hits into seperate levels
    global_hmm_obj_recipe = []
    current_batch = []
    current_hit_count = 1
    hit_id = 0
    batch_i = 0

    hits_db_conn = wrap_rocks.RocksDB(hits_path)

    dupe_counts = json.loads(sequences_db_conn.get("getall:dupes"))
    dupes_per_gene = {}

    for gene in transcripts_mapped_to:
        dupe_count_divvy = {}

        for hit in transcripts_mapped_to[gene]:
            if hit.header in sequence_dict:
                # Make dupe count gene based
                if hit.base_header in dupe_counts:  # NT Dupe
                    dupe_count_divvy[hit.base_header] = dupe_counts[hit.base_header]
                if hit.header in dupe_counts:  # AA Dupe
                    dupe_count_divvy[hit.header] = dupe_counts[hit.header]

                hit_id += 1
                hit.hmm_sequence = "".join(sequence_dict[hit.header])
                hit.hmm_id = hit_id
                if current_hit_count >= MAX_HMM_BATCH_SIZE:
                    data = json.dumps(current_batch)
                    del current_batch

                    key = f'hmmbatch:{batch_i}'

                    global_hmm_obj_recipe.append(str(batch_i))

                    batch_i += 1

                    hits_db_conn.put(key, data)

                    del data
                    current_batch = []
                    current_hit_count = 1

                current_batch.append(hit.to_json())
                current_hit_count += 1

        if len(dupe_count_divvy) > 1:
            dupes_per_gene[gene] = dupe_count_divvy
    key = "getall:gene_dupes"
    data = json.dumps(dupes_per_gene)
    sequences_db_conn.put(key, data)

    del sequence_dict

    if current_batch:
        data = json.dumps(current_batch)
        key = f'hmmbatch:{batch_i}'
        global_hmm_obj_recipe.append(str(batch_i))
        hits_db_conn.put(key, data)

    del current_batch

    # insert key to grab all hmm objects into the database
    key = 'hmmbatch:all'
    data = ','.join(global_hmm_obj_recipe)
    hits_db_conn.put(key, data)

    batches = len(global_hmm_obj_recipe)
    kicks = count - total_hits
    printv(f"Inserted {total_hits} hits over {batches} batch(es) in {time_keeper.lap():.2f} seconds. Kicked {kicks} hits during filtering", args.verbose)
    printv(f"Done! Took {time_keeper.differential():.2f}s overall", args.verbose, 0)


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
    for input_path in args.INPUT:
        printv(f"Processing: {os.path.basename(input_path)}", args.verbose)
        run_process(args, input_path)
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr Hmmsearch"
    )
