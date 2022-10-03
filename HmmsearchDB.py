import argparse
import math
from multiprocessing.pool import Pool
import os
import subprocess
from sys import argv
from time import time
import wrap_rocks
import json
from tqdm import tqdm


class Hit:
    __slots__ = (
        "target",
        "query",
        "evalue",
        "score",
        "hmm_start",
        "hmm_end",
        "ali_start",
        "ali_end",
        "env_start",
        "env_end",
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
        self.target = target
        self.query = query
        self.evalue = evalue
        self.score = score
        self.hmm_start = hmm_start
        self.hmm_end = hmm_end
        self.ali_start = ali_start
        self.ali_end = ali_end
        self.env_start = env_start
        self.env_end = env_end

    def __str__(self):
        return "\t".join(
            [
                self.target,
                self.query,
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

    def to_json(self):
        return {
            "header": self.target,
            "gene": self.query,
            "score": self.score,
            "hmm_evalue": self.evalue,
            "hmm_start": self.hmm_start,
            "hmm_end": self.hmm_end,
            "ali_start": self.ali_start,
            "ali_end": self.ali_end,
            "env_start": self.env_start,
            "env_end": self.env_end,
        }

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
            self.query,
            self.target,
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

def get_baseheader(header):
    """
    Returns header content before first whitespace.
    """
    baseheader = header.split("|")[0].strip()
    return baseheader

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

def get_overlap(a_start, a_end, b_start, b_end):
    overlap_end = min(a_end, b_end)
    overlap_start = max(a_start, b_start)
    amount = (overlap_end - overlap_start) + 1 #inclusive 

    return 0 if amount < 0 else amount

def internal_filter_gene(this_gene_hits, gene, min_overlap_internal, score_diff_internal, filter_verbose):
    this_gene_hits.sort(key=lambda hit: hit['score'], reverse=True)

    filtered_sequences_log = []

    for i,hit_a in enumerate(this_gene_hits):
        if not hit_a:
            continue
        for j in range(len(this_gene_hits)-1, i, -1):
            hit_b = this_gene_hits[j]
            if hit_b:
                if get_baseheader(hit_a['header']) != get_baseheader(hit_b['header']):
                    if get_difference(hit_a['score'], hit_b['score']) < score_diff_internal:
                        break
                    
                    amount_of_overlap = get_overlap(hit_a['hmm_start'], hit_a['hmm_end'], hit_b['hmm_start'], hit_b['hmm_end'])
                    distance = (hit_b['hmm_end'] - hit_b['hmm_start']) + 1 # Inclusive
                    percentage_of_overlap = amount_of_overlap / distance

                    if (percentage_of_overlap >= min_overlap_internal):
                        this_gene_hits[j] = None
                        if filter_verbose: 
                            filtered_sequences_log.append([hit_b['gene'],hit_b['header'],str(hit_b['score']),str(hit_b['hmm_start']),str(hit_b['hmm_end']),'Internal Overlapped with Lowest Score',hit_a['gene'],hit_a['header'],str(hit_a['score']),str(hit_a['hmm_start']),str(hit_a['hmm_end'])])

    this_out_data = {'Passes':[i for i in this_gene_hits if i is not None], 'Log':filtered_sequences_log, 'gene':gene}

    return this_out_data

def run_internal_filter(
    this_gene_transcripts,
    orthoid,
    min_overlap_internal,
    score_diff_internal,
    filter_verbose,
):
    return internal_filter_gene(
        this_gene_transcripts,
        orthoid,
        min_overlap_internal,
        score_diff_internal,
        filter_verbose,
    )

def run_multi_filter(
    this_hits,
    filter_verbose,
    min_overlap_multi,
    score_diff_multi,
):
    return multi_filter_dupes(
        this_hits,
        filter_verbose,
        min_overlap_multi,
        score_diff_multi,
    )

def multi_filter_dupes(
    this_hits,
    filter_verbose,
    min_overlap_multi,
    score_diff_multi,
):
    kick_happend = True
    filtered_sequences_log = []
    this_kicks = []

    this_hits.sort(key=lambda data: (data["score"], data["gene"]), reverse = True)

    while kick_happend:
        kick_happend = False

        master = this_hits[0]
        candidates = this_hits[1:]

        master_env_start = master["env_start"]
        master_env_end = master["env_end"]

        for candidate in candidates:
            if candidate["gene"] == master["gene"]:  # From same gene = Pseudomaster
                # Remove
                this_hits.remove(candidate)
                this_kicks.append(candidate)
                if filter_verbose:
                    filtered_sequences_log.append(
                        [
                            candidate["gene"],
                            candidate["header"],
                            str(candidate["score"]),
                            str(candidate["hmm_start"]),
                            str(candidate["hmm_end"]),
                            "Pseudomaster",
                            master["gene"],
                            master["header"],
                            str(master["score"]),
                            str(master["hmm_start"]),
                            str(master["hmm_end"]),
                        ]
                    )
                candidates.remove(candidate)
                # Extend master range
                kick_happend = True
                if candidate["env_start"] < master_env_start:
                    master_env_start = candidate["env_start"]
                if candidate["env_end"] > master_env_start:
                    master_env_end = candidate["env_end"]
            else:
                break

        miniscule_score = False

        for candidate in candidates:
            distance = (master_env_end - master_env_start) + 1 # Inclusive
            amount_of_overlap = get_overlap(master_env_start, master_env_end, candidate["env_start"], candidate["env_end"])
            percentage_of_overlap = amount_of_overlap / distance

            if percentage_of_overlap >= min_overlap_multi:
                score_difference = get_difference(master["score"], candidate["score"])
                if score_difference >= score_diff_multi:
                    kick_happend = True
                    this_hits.remove(candidate)
                    this_kicks.append(candidate)
                    if filter_verbose:
                        filtered_sequences_log.append(
                            [
                                candidate["gene"],
                                candidate["header"],
                                str(candidate["score"]),
                                str(candidate["env_start"]),
                                str(candidate["env_end"]),
                                "Multi Overlapped with Lowest Score",
                                master["gene"],
                                master["header"],
                                str(master["score"]),
                                str(master["env_start"]),
                                str(master["env_end"]),
                            ]
                        )
                else:
                    miniscule_score = True
                    break

        if (
            miniscule_score
        ):  # Remove all overlapping candidates if it's score is a miniscule difference of the masters
            for candidate in candidates:
                distance = (master_env_end - master_env_start) + 1 # Inclusive
                amount_of_overlap = get_overlap(master_env_start, master_env_end, candidate["env_start"], candidate["env_end"])
                percentage_of_overlap = amount_of_overlap / distance
                if percentage_of_overlap >= min_overlap_multi:
                    kick_happend = True
                    this_hits.remove(candidate)
                    this_kicks.append(candidate)
                    if filter_verbose:
                        filtered_sequences_log.append(
                            [
                                candidate["gene"],
                                candidate["header"],
                                str(candidate["score"]),
                                str(candidate["hmm_start"]),
                                str(candidate["hmm_end"]),
                                "Multi Overlapped with Miniscule Score",
                                master["gene"],
                                master["header"],
                                str(master["score"]),
                                str(master["hmm_start"]),
                                str(master["hmm_end"]),
                            ]
                        )

    multi_data = {
        "Kicks": this_kicks,
        "Log": filtered_sequences_log,
        "Remaining": this_hits,
    }
    return multi_data

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
    num_threads: int,
    ortholog: str,
    temp: str,
    force_prot: bool,
    specimen_type=1,
    table="orthograph_ests",
    translate_program = "fastatranslate",
    genetic_code = 1
) -> str:
    """
    Looks up the digest and sequence in the orthodb est table, then
    writes to the protfile. Returns path to protfile as a str.
    Creates side effects.
    """
    prot_name = ortholog + "_prot.tmp"
    prot_path = os.path.join(temp, prot_name)
    # return existing path if protfile exists and we aren't forcing a new one
    if not force_prot and os.path.exists(prot_path):
        return prot_path

    prepared_name = ortholog + "_prep.tmp"
    prepared_path = os.path.join(temp, prepared_name)

    open(prot_path,"w")

    start = time()

    recipe = db_conn.get("getall:prepared").split(',')

    with open(prepared_path, "w+") as prot_file_handle:
        for component in recipe:
            prot_file_handle.write(db_conn.get(component))

    os.system(
                    f"{translate_program} --geneticcode {genetic_code} '{prepared_path}' > '{prot_path}'"
                )
        
    print("Wrote & translated prot file in {:.2f}s.".format(time() - start))
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
    data_lines_seen = dict()
    output = list()
    with open(domtbl_path) as f:
        for line in f:
            if line[0] == "#":
                output.append(line)
                continue
            seen_before = data_lines_seen.get(line, False)
            if not seen_before:
                output.append(line)
                data_lines_seen[line] = True
    with open(domtbl_path, "w+") as f:
        f.writelines(output)


def get_hits_from_domtbl(domtbl_path: str, score, evalue) -> list:
    domtbl_dupe_check(domtbl_path)
    hits = list()
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
    evalue,
    score,
    ovw,
    prog="hmmsearch",
    threads=1,
) -> str:
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
        print(f"{domtbl_path}:hmmsearch error code {p.returncode}")
    else:
        print(f"{domtbl_path}")
    return get_hits_from_domtbl(domtbl_path, score, evalue)


def hmm_search(
    hmm_file: str, domtbl_dir: str, evalue, score, prot: str, ovw: bool
) -> None:
    """
    Reimplements hmmsearch loop in lines 468 to 538 in orthograph analyzer.
    """
    # hmmfile contains dir in path at this point
    hmm_name = get_hmm_name(hmm_file)
    domtbl_path = os.path.join(domtbl_dir, hmm_name + ".hmm.domtbl")

    hits = search_prot(prot, domtbl_path, hmm_file, evalue, score, ovw)
    return hits

def de_interleave_record(s):
    return s.split('\n',1)

def main(argv):
    global db_conn
    global_start = time()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default="PhyMMR/Acroceridae/SRR6453524.fa",
        help="Path to input directory.",
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
        "-ovw",
        "--overwrite",
        default=False,
        action="store_true",
        help="Remake domtbl files even if previous file exists.",
    )
    parser.add_argument(
        "-s", "--score", type=float, default=40, help="Score threshold. Defaults to 40"
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=0,
        help="Evalue threshold. Defaults to 0",
    )
    parser.add_argument(
        "--excluded-list",
        default=False,
        help="File containing names of genes to be excluded",
    )
    parser.add_argument(
        "--wanted-list", default=False, help="File containing list of wanted genes"
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="""Number of worker processes launched.
                        Defaults to 1.""",
    )
    parser.add_argument(
        "--remake-protfile",
        default=False,
        action="store_true",
        help="Force creation of a new protfile even if one already exists.",
    )
    parser.add_argument(
        "-sdm",
        "--score_diff_multi",
        type=float,
        default=1.05,
        help="Multi-gene Score Difference Adjustment",
    )
    parser.add_argument(
        "-mom",
        "--min_overlap_multi",
        type=float,
        default=0.3,
        help="Multi-gene Minimum Overlap Adjustment",
    )
    parser.add_argument(
        "-sdi",
        "--score_diff_internal",
        type=float,
        default=1.5,
        help="Internal Score Difference Adjustmen",
    )
    parser.add_argument(
        "-moi",
        "--min_overlap_internal",
        type=float,
        default=0.9,
        help="Internal Minimum Overlap Adjustment",
    )
    parser.add_argument(
        "-m",
        "--max_hmm_batch_size",
        default=250000, 
        type=int, 
        help="Max hits per hmmsearch batch in db. Default: 500 thousand.",
    )
    parser.add_argument("-d", "--debug", type=int, default=0, help="Verbose debug.")

    args = parser.parse_args()

    MAX_HMM_BATCH_SIZE = args.max_hmm_batch_size
    debug = args.debug != 0
    score_diff_multi = args.score_diff_multi
    min_overlap_multi = args.min_overlap_multi
    score_diff_internal = args.score_diff_internal
    min_overlap_internal = args.min_overlap_internal

    num_threads = args.processes
   
    if os.path.exists("/run/shm"):
        in_ram = "/run/shm"
    else:
        in_ram = "/dev/shm"
    # temp_dir = os.path.join(args.input, "tmp")
    temp_dir = os.path.join(in_ram, "tmp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # get ortholog
    ortholog = os.path.basename(args.input).split(".")[0]

    # if database not set, make expected path to database
    db_path = os.path.join(args.input, "rocksdb")
    db_conn = wrap_rocks.RocksDB(db_path)

    # path to domtbl directory
    domtbl_dir = os.path.join(args.input, "hmmsearch")

    if not os.path.exists(domtbl_dir):
        os.mkdir(domtbl_dir)

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
    print("Finding hmm files")
    hmm_list = [
        hmm for hmm in os.listdir(hmm_dir) if ".hmm" in hmm and hmm not in excluded
    ]
    if wanted:
        hmm_list = [hmm for hmm in hmm_list if hmm.split(".")[0] in wanted]

    hmm_list.sort()

    # rejoin file names with directory path
    hmm_list = [os.path.join(hmm_dir, hmm) for hmm in hmm_list]

    # make protfile for hmm_search later
    print("Checking protfile")

    protfile = make_temp_protfile(
        num_threads, ortholog, temp_dir, args.remake_protfile, specimen_type=1
    )

    # clean prot file and store a hashmap of sequences for later use in filtering

    sequence_dict = {}
    start = time()
    with open(protfile) as prot_file_handle:
        content = prot_file_handle.read().strip().replace(' ','|')

        records = content.split('>')[1:]
        records = map(de_interleave_record, records)
        for header,sequence in records:
            sequence_dict[header] = sequence

    with open(protfile, 'w') as prot_file_ovw:
        prot_file_ovw.write(content)
    
    # save the sequence dict for later.
    sequence_dict_location = os.path.join(temp_dir, "SeqDict.tmp")
    with open(sequence_dict_location, "w") as seq_dict_handle:
        json.dump(sequence_dict, seq_dict_handle)

    del sequence_dict

    print('Cleaned and read prot file in {:.2f}s'.format(time()-start))

    arg_tuples = list()
    for hmm in hmm_list:
        arg_tuples.append(
            (hmm, domtbl_dir, args.evalue, args.score, protfile, args.overwrite)
        )

    start = time()
    hmm_results = list()  # list of lists containing Hit objects
    with Pool(num_threads) as search_pool:
        hmm_results = search_pool.starmap(hmm_search, arg_tuples)
    print("Search time: {:.2f}".format(time() - start))

    if os.path.exists(protfile):
        os.remove(protfile)

    filter_start = time()

    f_duplicates = {}
    gene_based_results = {}
    header_based_results = {}
    count = 0

    for hit_group in hmm_results:
        for hit in hit_group:
            count += 1
            this_orthoid = hit.query
            if this_orthoid not in gene_based_results:
                gene_based_results[this_orthoid] = []

            this_hit_data = hit.to_json()
            this_hit_data["uuid"] = this_hit_data["header"] + f"_hit_{count}"

            gene_based_results[this_orthoid].append(this_hit_data)

    print(f"Filtering {count} hits.")

    print('Filtering multi-gene dupes')

    for orthoid in gene_based_results:
        for hit in gene_based_results[orthoid]:
            if "revcomp" in hit["header"]:
                base_header = get_baseheader(hit["header"])
                raw_length = int(base_header.split("_length_")[1]) / 3
                length = math.floor(raw_length)

                new_env_start = length - int(hit["env_end"])
                new_env_end = length - int(hit["env_start"])
                hit["env_start"] = new_env_start
                hit["env_end"] = new_env_end

            if hit["header"] not in f_duplicates:
                f_duplicates[hit["header"]] = []
            if hit["header"] not in header_based_results:
                header_based_results[hit["header"]] = []

            f_duplicates[hit["header"]].append(hit)
            header_based_results[hit["header"]].append(hit)

    headers = list(f_duplicates.keys())
    for header in headers:
        if len(f_duplicates[header]) > 1:
            unique_genes = list(
                dict.fromkeys([i["gene"] for i in f_duplicates[header]])
            )
            if len(unique_genes) <= 1:  # But the same gene
                f_duplicates.pop(header)
        else:
            f_duplicates.pop(header)

    for header_left in f_duplicates:
        header_based_results[header_left] = []

    filter_verbose = debug
    filtered_sequences_log = []

    if num_threads == 1:
        for header in f_duplicates:
            this_hits = f_duplicates[header]
            data = multi_filter_dupes(
                this_hits,
                filter_verbose,
                min_overlap_multi,
                score_diff_multi,
            )

            filtered_sequences_log.extend(data["Log"])
            header_based_results[data["Remaining"][0]["header"]] = data["Remaining"]

    else:
        arguments = list()
        for header in f_duplicates:
            this_hits = f_duplicates[header]
            arguments.append(
                (
                    this_hits,
                    filter_verbose,
                    min_overlap_multi,
                    score_diff_multi,
                )
            )

        with Pool(num_threads) as pool:
            multi_data = pool.starmap(run_multi_filter, arguments, chunksize=1)

        for data in multi_data:
            filtered_sequences_log.extend(data["Log"])

            hits_to_add = []
            for hit in data["Remaining"]:
                hits_to_add.append(hit)

            header_based_results[data["Remaining"][0]["header"]] = hits_to_add

    transcripts_mapped_to = {}


    for header in header_based_results:
        for match in header_based_results[header]:
            if match["gene"] not in transcripts_mapped_to:
                transcripts_mapped_to[match["gene"]] = []
            transcripts_mapped_to[match["gene"]].append(match)

    print('Doing internal filtering')

    total_hits = 0
    if num_threads == 1:
        internal_data = []
        for orthoid in transcripts_mapped_to:
            this_gene_transcripts = transcripts_mapped_to[orthoid]
            internal_data.append(
                internal_filter_gene(
                    this_gene_transcripts,
                    orthoid,
                    min_overlap_internal,
                    score_diff_internal,
                    filter_verbose,
                )
            )

    else:
        arguments = list()
        for orthoid in transcripts_mapped_to:
            this_gene_transcripts = transcripts_mapped_to[orthoid]
            arguments.append(
                (
                    this_gene_transcripts,
                    orthoid,
                    min_overlap_internal,
                    score_diff_internal,
                    filter_verbose,
                )
            )
        with Pool(num_threads) as pool:
            internal_data = pool.starmap(run_internal_filter, arguments, chunksize=1)


    transcripts_mapped_to = {}
    duplicates = {}
    total_kicks = 0

    for data in internal_data:
        gene = data["gene"]
        if gene not in duplicates.keys():
            duplicates[gene] = set()

        transcripts_mapped_to[gene] = []

        for this_match in data["Passes"]:
            if this_match["uuid"] in duplicates[gene]:
                filtered_sequences_log.append(
                    [
                        this_match["gene"],
                        this_match["header"],
                        str(this_match["score"]),
                        str(this_match["hmm_start"]),
                        str(this_match["hmm_end"]),
                        "Duplicate UUID",
                    ]
                )
                continue
            else:
                duplicates[gene].add(this_match["uuid"])
            transcripts_mapped_to[gene].append(this_match)

        total_hits += len(data["Passes"])
        filtered_sequences_log.extend(data["Log"])

    filtered_sequences_log_out = []
    for line in filtered_sequences_log:
        filtered_sequences_log_out.append(",".join(line) + "\n")

    if debug:
        filtered_sequences_log_path = os.path.join(args.input, "filtered-hits.csv")
        open(filtered_sequences_log_path, "w").writelines(filtered_sequences_log_out)

    print('Done! Filtering took {:.2f}s'.format(time() - filter_start))

    # Grab that seq dict
    with open(sequence_dict_location) as seq_dict_handle:
        sequence_dict = json.load(seq_dict_handle)

    os.remove(sequence_dict_location)

    end = time()

    # Seperate hmm hits into seperate levels
    global_hmm_obj_recipe = []
    current_batch = []
    current_hit_count = 1
    hit_id = 0
    batch_i= 0 
    for gene in transcripts_mapped_to:
        for hit in transcripts_mapped_to[gene]:
            hit_id += 1
            hit["hmm_sequence"] = "".join(sequence_dict[hit["header"]])
            hit["hmm_id"] = hit_id
            if current_hit_count >= MAX_HMM_BATCH_SIZE:
                data = json.dumps(current_batch)
                del current_batch

                key = f'hmmbatch:{batch_i}'

                global_hmm_obj_recipe.append(str(batch_i))

                batch_i += 1

                db_conn.put(key, data)

                del data
                current_batch = []
                current_hit_count = 1

            current_batch.append(hit)
            current_hit_count += 1

    if current_batch:
        data = json.dumps(current_batch)
        key = f'hmmbatch:{batch_i}'

        global_hmm_obj_recipe.append(str(batch_i))

        db_conn.put(key, data)

    del current_batch

    # insert key to grab all hmm objects into the database

    key = 'hmmbatch:all'
    data = ','.join(global_hmm_obj_recipe)

    db_conn.put(key, data)

    db_end = time()
    print("Inserted {} hits over {} batch(es) in {:.2f} seconds. Kicked {} hits during filtering".format(total_hits, len(global_hmm_obj_recipe), db_end - end, count-total_hits))
    print("Took {:.2f}s overall".format(time() - global_start))
    print("Cleaning temp files. Closing DB.")
    clear_command = f'rm {os.path.join(temp_dir,"*")}'
    os.system(clear_command)

if __name__ == "__main__":
    main(argv)
