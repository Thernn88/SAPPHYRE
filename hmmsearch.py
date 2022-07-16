import argparse
import inserts
import math
from multiprocessing.pool import Pool
import os
import sqlite3
import subprocess
from sys import argv
from time import time


class Hit:
    __slots__ = ('target', 'query', 'evalue', 'score', 'hmm_start', 'hmm_end',
                 'ali_start', 'ali_end', 'env_start', 'env_end')

    def __init__(self, target, query, evalue, score, hmm_start, hmm_end,
                 ali_start, ali_end, env_start, env_end):
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
        return "\t".join([self.target, self.query, self.evalue, self.score, self.hmm_start,
                          self.hmm_end, self.ali_start, self.ali_end, self.env_start, self.env_end])

    def list_values(self, species_id):
        """
        Given a species_id, returns a list of values suitable for the
        hmmsearch table insertion.
        """
        log_evalue = -999
        e_float = float(self.evalue)
        if e_float != 0:
            log_evalue = math.log(e_float)
        return [str(species_id), self.query, self.target, self.score, str(e_float), str(log_evalue), self.hmm_start,
                self.hmm_end, self.ali_start, self.ali_end, self.env_start, self.env_end]


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


def get_orthoset_id_from_database(orthoset: str, ortho_db) -> int:
    """
    Gets the orthoset id number from the database.
    Should only find one result.
    Should only run once.
    """
    try:
        print("Getting ortholog set id number.")
        db = sqlite3.connect(ortho_db)
        cursor = db.cursor()
        query = """
        SELECT id
        FROM orthograph_species_info
        WHERE name = ?;
        """
        cursor.execute(query, (orthoset,))
        result = cursor.fetchall()
        db.close()
        if len(result) > 1:
            print(f"Warning: found multiple sets with name {orthoset}.\n" +
                  "Returning id number of first result.")
        return result[0][0]
    except sqlite3.Error:
        print("Cannot get set id from ortholog database")
        raise IOError


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


def make_hmm_query(ests="orthograph_ests", hmmsearch="orthograph_hmmsearch") -> str:
    """
    Generate sql query needed to get hmm results. Mirrors query in
    Sqlite.pm Returns a str containing the query.
    """
    query = """SELECT {0}.digest, {0}.sequence, {1}.env_start, {1}.env_end, {1}.id
    FROM {0}
    INNER JOIN {1}
    ON {1}.target = {0}.digest
    WHERE {1}.query = ?
    AND {1}.taxid = ?
    """.format(ests, hmmsearch)
    return query


def get_hmm_results(database: str) -> list:
    """
    Runs a query and returns a list of results
    """
    try:
        db = sqlite3.connect(database)
        cursor = db.cursor()
        query = make_hmm_query()
        cursor.execute(query)
        record = cursor.fetchall()
        db.close()
        return record

    except sqlite3.Error:
        raise TimeoutError("Cannot query database")

def make_temp_protfile(species_id: int, database: str, temp: str, force_prot: bool, specimen_type=1,
                       table="orthograph_ests") -> str:
    """
    Looks up the digest and sequence in the orthodb est table, then
    writes to the protfile. Returns path to protfile as a str.
    Creates side effects.
    """
    prot_name = os.path.basename(database).split('.')[0] + "_prot.tmp"
    prot_path = os.path.join(temp, prot_name)
    # return existing path if protfile exists and we aren't forcing a new one
    if not force_prot and os.path.exists(prot_path):
        return prot_path
    db = sqlite3.connect(database)
    cursor = db.cursor()
    query = """
    SELECT digest, sequence
    FROM {0}
    WHERE taxid = ?
    AND type = ?
    """.format(table)
    cursor.execute(query, (species_id, specimen_type))
    prot_file_handle = open(prot_path, 'w+')
    # fetch results 1000 at a time
    while True:
        rows = cursor.fetchmany(1000)
        if not rows:
            break
        for row in rows:
            prot_file_handle.write(f">{row[0]}\n{row[1]}\n")
    # close connections and handles
    prot_file_handle.close()
    db.close()
    return prot_path


def is_valid_result_line(fields, score, evalue) -> bool:

    if float(fields[13]) >= score:
        if float(fields[12]) <= evalue:
            return True
    return False


def parse_domtbl_fields(fields: list) -> Hit:
    hit = Hit(fields[0], fields[3], fields[12], fields[13], fields[15], fields[16],
              fields[17], fields[18], fields[19], fields[20])
    return hit


def domtbl_dupe_check(domtbl_path: str) -> None:
    """Reads the given domtbl file and removes duplicate data rows.
    Overwrites old file with new version."""
    data_lines_seen = dict()
    output = list()
    with open(domtbl_path) as f:
        for line in f:
            if line[0] == '#':
                output.append(line)
                continue
            seen_before = data_lines_seen.get(line, False)
            if not seen_before:
                output.append(line)
                data_lines_seen[line] = True
    with open(domtbl_path, 'w+') as f:
        f.writelines(output)



def get_hits_from_domtbl(domtbl_path: str, score, evalue) -> list:
    domtbl_dupe_check(domtbl_path)
    hits = list()
    with open(domtbl_path) as f:
        for line in f:
            if line and line[0] == '#':
                continue
            fields = line.strip().split()
            # if is_valid_result_line(fields, score, evalue):
            hits.append(parse_domtbl_fields(fields))
    return hits


def empty_domtbl_file(path: str) -> bool:
    """
    Aborted runs result in a domtbl file with 0 bytes. If such a file is
    found, return True. Otherwise returns False.
    If an empty domtbl already exists, always overwrite.
    """
    return os.path.getsize(path) == 0


def search_prot(prot_file: str, domtbl_path: str, hmm_file: str, evalue, score, ovw, prog="hmmsearch", threads=1) -> str:
    if os.path.exists(domtbl_path) and not ovw:
        # always overwrite the empty domtbl files
        if not empty_domtbl_file(domtbl_path):
            print(f"{domtbl_path}")
            return get_hits_from_domtbl(domtbl_path, score, evalue)
    if score:
        option = "-T"
        threshold = str(score)
    else:
        option = '-E'
        threshold = str(evalue)
    # $searchprog --domtblout $outfile $threshold_option --cpu $num_threads $hmmfile $protfile
    command = [
        prog,
        "--domtblout", domtbl_path,
        option, threshold,
        "--cpu", str(threads),
        hmm_file,
        prot_file
    ]
    p = subprocess.run(command, stdout=subprocess.PIPE)
    if p.returncode != 0:  # non-zero return code means an error
        print(f'{domtbl_path}:hmmsearch error code {p.returncode}')
    else:
        print(f'{domtbl_path}')
    return get_hits_from_domtbl(domtbl_path, score, evalue)


def hmm_search(hmm_file: str, domtbl_dir: str, evalue, score, prot: str, ovw: bool) -> None:
    """
    Reimplements hmmsearch loop in lines 468 to 538 in orthograph analyzer.
    """
    # hmmfile contains dir in path at this point
    hmm_name = get_hmm_name(hmm_file)
    domtbl_path = os.path.join(domtbl_dir, hmm_name + ".hmm.domtbl")
    hits = search_prot(prot, domtbl_path, hmm_file, evalue, score, ovw)
    return hits


def make_insert_query(table="orthograph_hmmsearch") -> str:
    query = """INSERT OR IGNORE INTO {0} (
    'taxid',
    'query',
    'target',
    'score',
    'evalue',
    'log_evalue',
    'hmm_start',
    'hmm_end',
    'ali_start',
    'ali_end',
    'env_start',
    'env_end'
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """.format(table)
    return query

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default="Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa",
                        help="Path to input directory.")
    parser.add_argument('-oi', '--orthoset_input', type=str, default='Syrphidae/orthosets',
                        help='Path to directory of Orthosets folder')
    parser.add_argument('-o', '--orthoset', type=str, default='Ortholog_set_Mecopterida_v4',
                        help='Orthoset')
    parser.add_argument('-ovw', '--overwrite', default=False, action='store_true',
                        help="Remake domtbl files even if previous file exists.")
    parser.add_argument('-s', '--score', type=float, default=40,
                        help="Score threshold. Defaults to 40")
    parser.add_argument('-e', '--evalue', type=float, default=0,
                        help="Evalue threshold. Defaults to 0.00001")
    parser.add_argument('--excluded-list', default=False,
                        help='File containing names of genes to be excluded')
    parser.add_argument('--wanted-list', default=False,
                        help='File containing list of wanted genes')
    parser.add_argument('-p', '--processes', type=int, default=1,
                        help="""Number of worker processes launched.
                        Defaults to 1.""")
    parser.add_argument('--remake-protfile', default=False, action='store_true',
                        help='Force creation of a new protfile even if one already exists.')
    args = parser.parse_args()
    if os.path.exists("/run/shm"):
        in_ram = "/run/shm"
    else:
        in_ram = "/dev/shm"
    # temp_dir = os.path.join(args.input, "tmp")
    temp_dir = os.path.join(in_ram, "tmp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # if database not set, make expected path to database
    ortholog = os.path.basename(args.input).split(".")[0]
    db_path = ortholog + ".sqlite"
    db_path = os.path.join(args.input, db_path)

    # path to domtbl directory
    domtbl_dir = os.path.join(args.input, "hmmsearch")

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

    # get the ortholog set id number
    set_id = get_orthoset_id_from_database(ortholog, db_path)

    # make a list of valid ortholog names, excluded hidden files
    print("Finding hmm files")
    hmm_list = [hmm for hmm in os.listdir(hmm_dir)
                if '.hmm' in hmm and hmm not in excluded]
    if wanted:
        hmm_list = [hmm for hmm in hmm_list if hmm.split('.')[0] in wanted]

    hmm_list.sort()

    # rejoin file names with directory path
    hmm_list = [os.path.join(hmm_dir, hmm) for hmm in hmm_list]

    # make protfile for hmm_search later
    print("Checking protfile")
    protfile = make_temp_protfile(set_id, db_path, temp_dir, args.remake_protfile, specimen_type=1)

    arg_tuples = list()
    for hmm in hmm_list:
        arg_tuples.append((hmm, domtbl_dir, args.evalue, args.score, protfile, args.overwrite))
    start = time()
    hmm_results = list()  # list of lists containing Hit objects
    with Pool(args.processes) as search_pool:
        hmm_results = search_pool.starmap(hmm_search, arg_tuples)
    end = time()
    print("Search time:", end - start)

    # insert results into the database
    db_conn = sqlite3.connect(db_path)
    cursor = db_conn.cursor()
    print("Clearing table")
    # clear indices for faster inserts
    cursor.execute("BEGIN")
    cursor.execute("DROP INDEX IF EXISTS orthograph_hmmsearch_query;")
    cursor.execute("DROP INDEX IF EXISTS orthograph_hmmsearch_score;")
    cursor.execute("DROP INDEX IF EXISTS orthograph_hmmsearch_target;")
    cursor.execute("DROP INDEX IF EXISTS orthograph_hmmsearch_taxid;")
    cursor.execute("DROP INDEX IF EXISTS orthograph_hmmsearch_evalue;")
    cursor.execute("END;")
    cursor.execute("DELETE FROM orthograph_hmmsearch;")
    db_conn.commit()

    print("Inserting results into database")
    vectors = list()
    for hit_group in hmm_results:
        for hit in hit_group:
            vectors.append(hit.list_values(set_id))
    inserts.hmmsearch_inserts(db_path, vectors, True)
    cursor.execute('BEGIN')
    cursor.execute("CREATE INDEX orthograph_hmmsearch_evalue ON orthograph_hmmsearch (log_evalue);")
    cursor.execute("CREATE INDEX orthograph_hmmsearch_query ON orthograph_hmmsearch (query);")
    cursor.execute("CREATE INDEX orthograph_hmmsearch_score ON orthograph_hmmsearch (score);")
    cursor.execute("CREATE INDEX orthograph_hmmsearch_target ON orthograph_hmmsearch (target);")
    cursor.execute("CREATE INDEX orthograph_hmmsearch_taxid ON orthograph_hmmsearch (taxid);")
    cursor.execute('END')
    db_conn.commit()
    db_conn.close()

    db_end = time()
    print("Insertion time:", db_end - end)
    print("Cleaning temp files.")
    clear_command = f'rm {os.path.join(temp_dir,"*")}'
    os.system(clear_command)


if __name__ == '__main__':
    main(argv)
