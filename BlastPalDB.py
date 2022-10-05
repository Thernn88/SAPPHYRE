import os
import math
from multiprocessing.pool import Pool
import argparse
from time import time,sleep
from shutil import rmtree
import json
import wrap_rocks
import sqlite3

class Blast:
    __slots__ = (
        "query_id",
        "subject_id",
        "score",
        "evalue",
        "log_evalue",
        "blast_start",
        "blast_end",
        "hmmsearch_id",
        "reftaxon"
    )

    def __init__(self, query_id, subject_id, evalue, bit_score, q_start, q_end):
        try:
            log_evalue = str(math.log(float(evalue)))
        except: #Undefined
            log_evalue = '-999'
            
        self.query_id, self.hmmsearch_id = query_id.split('_hmmid')
        self.subject_id = int(subject_id)
        self.score = float(bit_score)
        self.evalue = float(evalue)
        self.log_evalue = float(log_evalue)
        self.blast_start = int(q_start)
        self.blast_end = int(q_end)

    def to_json(self):
        return {'target':self.subject_id,
                'score':self.score,
                'evalue':self.evalue,
                'log_evalue':self.log_evalue,
                'blast_start':self.blast_start,
                'blast_end':self.blast_end,
                'reftaxon':self.reftaxon}

def do(gene, tmp_path, this_gene_sequences, blast_path, blast_db_path, blast_minimum_score, blast_minimum_evalue, prog = 'blastp', evalue_threshold = .00001, score_threshold = 40, num_threads = 1):
    blast_file_name = f'{gene}.blast'
    this_blast_path = os.path.join(blast_path,blast_file_name)
    result = this_blast_path + ".done"

    if not os.path.exists(result) or os.path.getsize(result) == 0:
        if not os.path.exists(this_blast_path):
            print('Blasted:',gene)    
            target_tmp_path = os.path.join(tmp_path,gene+'.fa')

            with open(target_tmp_path,'w') as target_handle:
                for target, sequence in this_gene_sequences:
                    target_handle.write('>'+target+'\n'+sequence+'\n')

            if os.path.exists(this_blast_path):
                open(this_blast_path,'w').write('')

            cmd = "{prog} -outfmt '7 qseqid sseqid evalue bitscore qstart qend' -evalue '{evalue_threshold}' -threshold '{score_threshold}' -num_threads '{num_threads}' -db '{db}' -query '{queryfile}' -out '{outfile}'"
            cmd = cmd.format(prog = prog,
                             evalue_threshold = evalue_threshold,
                             score_threshold = score_threshold,
                             num_threads = 2,
                             db = blast_db_path,
                             queryfile = target_tmp_path,
                             outfile = this_blast_path)
            os.system(cmd)
            os.remove(target_tmp_path)
            
        os.rename(this_blast_path, result)

    gene_out = {}
    this_return = []

    this_content = open(result).read()
    if this_content != '':
        for line in this_content.split('\n'):
            fields = line.split('\t')
            if len(fields) == 6:
                this_blast = Blast(*fields)

                if this_blast.score >= blast_minimum_score and this_blast.evalue <= blast_minimum_evalue:
                    this_return.append(this_blast)
        
    return this_return

def get_orthologs_for_set_hashref(set_id, orthoset_db_con):
    result = set()
    query = f'''SELECT DISTINCT
        orthograph_orthologs.ortholog_gene_id,
        orthograph_aaseqs.id
        FROM orthograph_orthologs
        INNER JOIN orthograph_sequence_pairs
            ON orthograph_orthologs.sequence_pair = orthograph_sequence_pairs.id
        INNER JOIN orthograph_aaseqs
            ON orthograph_sequence_pairs.aa_seq = orthograph_aaseqs.id
        INNER JOIN orthograph_set_details
            ON orthograph_orthologs.setid = orthograph_set_details.id
        WHERE orthograph_set_details.id = "{set_id}"'''

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute(query)

    for row in rows:
        gene, id = row
        result.add(id)

    return result

def get_set_id(orthoset_db_con, orthoset):
    """
    Retrieves orthoset id from orthoset db
    """
    orthoset_id = None

    orthoset_db_cur = orthoset_db_con.cursor()
    rows = orthoset_db_cur.execute("SELECT * FROM orthograph_set_details;")

    for row in rows:
        id, name, _ = row

        if name == orthoset:
            return id

    if orthoset_id == None:
        raise Exception("Orthoset {} id cant be retrieved".format(orthoset))
        
def get_reftaxon_names(orthoset_db_path):
    result = {}
    query = f"SELECT orthograph_aaseqs.id, orthograph_taxa.name FROM orthograph_taxa INNER JOIN orthograph_aaseqs ON orthograph_taxa.id = orthograph_aaseqs.taxid"

    orthoset_db_con = sqlite3.connect(orthoset_db_path)
    orthoset_db_cur = orthoset_db_con.cursor()

    rows = orthoset_db_cur.execute(query)

    for row in rows:
        result[row[0]] = row[1]

    return result

def main():
    start = time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',type=str, default='PhyMMR/Acroceridae/SRR6453524.fa',
        help='Path to directory of Input folder')
    parser.add_argument('-oi','--orthoset_input',type=str, default='PhyMMR/orthosets',
        help='Path to directory of Orthosets folder')
    parser.add_argument('-o','--orthoset',type=str, default='Ortholog_set_Mecopterida_v4',
        help='Orthoset')    
    parser.add_argument('-bs', '--blast_minimum_score', type=float, default=40.0,
        help='Minimum score filter in blast.')
    parser.add_argument('-be', '--blast_minimum_evalue', type=float, default=0.00001,
        help='Minimum evalue filter in blast.')
    parser.add_argument('-ovw', '--overwrite', action="store_true",
        help='Overwrite existing blast results.')
    parser.add_argument('-p', '--processes', type=int, default=1,
        help='Number of threads used to call processes.')
    args = parser.parse_args()



    print('Grabbing HMM data from db.')

    input_path = args.input
    orthoset = args.orthoset
    orthosets_dir = args.orthoset_input

    orthoset_db_path = os.path.join(orthosets_dir, orthoset + ".sqlite")

    #make dirs
    blast_path = os.path.join(input_path,'blast')
    if args.overwrite:
        rmtree(blast_path)
    if not os.path.exists(blast_path):
        os.mkdir(blast_path)


    if os.path.exists("/run/shm"):
        tmp_path = "/run/shm"
    elif os.path.exists ("/dev/shm"):
        tmp_path = "/dev/shm"
    else:
        tmp_path = os.path.join(input_path,'tmp')
        if not os.path.exists(tmp_path):
            os.mkdir(tmp_path)

    blast_db_path = os.path.join(orthosets_dir,orthoset,'blast',orthoset)

    db_path = os.path.join(input_path,'rocksdb')
    db = wrap_rocks.RocksDB(db_path)

    gene_to_hits = {}

    grab_hmm_start = time()

    global_hmm_object_raw = db.get('hmmbatch:all')
    global_hmm_batches = global_hmm_object_raw.split(',')
    hit_count = 0

    for batch_i in global_hmm_batches:
        key = f'hmmbatch:{batch_i}'
        hmm_json = db.get(key)
        hmm_hits = json.loads(hmm_json)

        for hmm_object in hmm_hits:
            hit_count += 1
            hmm_id = hmm_object['hmm_id']
            header = hmm_object['header'].strip()
            gene = hmm_object['gene']

            if gene not in gene_to_hits:
                gene_to_hits[gene] = []

            gene_to_hits[gene].append((header+f'_hmmid{hmm_id}', hmm_object['hmm_sequence']))

    genes = list(gene_to_hits.keys())

    print('Grabbed HMM Data. Took: {:.2f}s. Grabbed {} hits.'.format(time()-grab_hmm_start, hit_count))

    del global_hmm_object_raw

    num_threads = args.processes

    to_write = []

    blast_start = time()

    #Run
    if num_threads == 1:
        for gene in genes:
            to_write.append(do(gene, tmp_path, gene_to_hits[gene], blast_path, blast_db_path, args.blast_minimum_score, args.blast_minimum_evalue))
    else:                                        
        arguments = list()
        for gene in genes:
            arguments.append((gene, tmp_path, gene_to_hits[gene],blast_path,blast_db_path, args.blast_minimum_score, args.blast_minimum_evalue))

        with Pool(num_threads) as pool:
            to_write = pool.starmap(do, arguments, chunksize=1)

    print('Blast done! Took {:.2f}s. Filtering wanted results'.format(time()-blast_start))
    filter_start = time()

    hmmsearch_id_results = {}

    orthoset_db_con = sqlite3.connect(orthoset_db_path)

    orthoset_id = get_set_id(orthoset_db_con, orthoset)
    aaseqs_in_orthoid = get_orthologs_for_set_hashref(orthoset_id, orthoset_db_con)

    pre_kick = 0
    total = 0
    reftaxons = get_reftaxon_names(orthoset_db_path)

    for batch in to_write:
        pre_kick += len(batch)
        for blast_result in batch:
            if blast_result.subject_id in aaseqs_in_orthoid:
                total += 1

                blast_result.reftaxon = reftaxons[blast_result.subject_id]

                if blast_result.hmmsearch_id not in hmmsearch_id_results:
                    hmmsearch_id_results[blast_result.hmmsearch_id] = [blast_result.to_json()]
                else:
                    hmmsearch_id_results[blast_result.hmmsearch_id].append(blast_result.to_json())

    print('Done! took {:.2f}s. Writing to DB.'.format(time() - filter_start))
    write_start = time()

    for hmmsearch_id, json_list in hmmsearch_id_results.items():
        data = json.dumps(json_list)
        key = 'blastfor:{}'.format(hmmsearch_id)
        db.put(key,data)

    print('Writing {} results took {:.2f}s after kicking {} results.'.format(total ,time()-write_start, pre_kick-total))
    print('Took {:.2f}s overall.'.format(time()-start))

if __name__ == "__main__":
    main()
