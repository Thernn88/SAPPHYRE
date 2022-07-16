import sqlite3
import os
import math
from multiprocessing.pool import Pool
from multiprocessing import Process
import argparse
from time import time,sleep
import sys
from shutil import rmtree

def do(gene, tmp_path, this_gene_sequences, blast_path,blast_db_path):
	blast_file_name = '{gene}.blast'.format(gene=gene)
	this_blast_path = os.path.join(blast_path,blast_file_name)

	if not os.path.exists(this_blast_path + ".done"):
		target_tmp_path = os.path.join(tmp_path,gene+'.fa')
		this_tmp_out = []

		for target,header,sequence in this_gene_sequences:
			this_tmp_out.append('>'+target+'\n')
			this_tmp_out.append(sequence+'\n')

		open(target_tmp_path,'w').writelines(this_tmp_out)

		open(this_blast_path,'w').write('')

		cmd = "blastp -outfmt '7 qseqid sseqid evalue bitscore qstart qend' -evalue '{evalue_threshold}' -threshold '{score_threshold}' -num_threads '{num_threads}' -db '{db}' -query '{queryfile}' -out '{outfile}'"
		cmd = cmd.format(evalue_threshold = .00001,
							score_threshold = 40,
							num_threads = 1,
							db = blast_db_path,
							queryfile = target_tmp_path,
							outfile = this_blast_path)
		os.system(cmd)
	
		os.remove(target_tmp_path)
		os.rename(this_blast_path, this_blast_path + ".done")
		print('Blasted:',gene)	

def run_command(arg_tuple: tuple) -> None:
    gene, tmp_path, gene_to_hashes, blast_path,blast_db_path = arg_tuple
    do(gene, tmp_path, gene_to_hashes, blast_path,blast_db_path)

def output_thread(genes,db_path,hashes_to_gene_and_id,blast_minimum_score,blast_minimum_evalue,blast_path):
	species_id = 1

	db_con = sqlite3.connect(db_path)
	db_cur = db_con.cursor()
    
	db_cur.execute("PRAGMA journal_mode = MEMORY")
	db_cur.execute('DROP TABLE orthograph_blast;')
	db_cur.execute("CREATE TABLE orthograph_blast (id INTEGER PRIMARY KEY, taxid UNSIGNED INTEGER NOT NULL, query TEXT(32) NOT NULL, target UNSIGNED INTEGER NOT NULL, score DOUBLE NOT NULL, evalue TEXT(8) NOT NULL, log_evalue DOUBLE NOT NULL, start UNSIGNED INTEGER NOT NULL, end UNSIGNED INTEGER NOT NULL, hmmsearch_id UNSIGNED INTEGER NOT NULL)")

	while len(genes) != 0:
		for gene in genes:
			blast_file_name = '{gene}.blast'.format(gene=gene)
			this_blast_path = os.path.join(blast_path,blast_file_name) + ".done"
			if os.path.exists(this_blast_path):
				this_content = open(this_blast_path).read()
				if this_content != '':
					hits = None
					this_hits = []
					for line in this_content.split('\n'):
						if '#' in line:
							if 'hits found' in line:
								hits = int(line.split(' ')[1])
						elif len(line.split('\t')) == 6:
							query_id, subject_id, evalue, bit_score, q_start, q_end = line.split('\t')

							if float(bit_score) >= blast_minimum_score and float(evalue) <= blast_minimum_evalue:
								try:
									log_evalue = str(math.log(float(evalue)))
								except: #Undefined
									log_evalue = '-999'

								gene,hmmsearch_id = hashes_to_gene_and_id[gene][query_id]

								hmmsearch_id = str(hmmsearch_id)

								query_id = query_id.split('_hmmstart')[0]
						
								sql_call = "INSERT OR IGNORE INTO orthograph_blast ('taxid','query','target','score','evalue','log_evalue','start','end','hmmsearch_id') VALUES (?,?,?,?,?,?,?,?,?)"
								db_cur.execute(sql_call,(species_id,query_id,subject_id,bit_score,evalue,log_evalue,q_start,q_end,hmmsearch_id))

					genes.remove(gene)
					print('Got output for:',gene)

	db_con.commit()
	db_cur.execute("BEGIN")
	db_cur.execute("CREATE INDEX orthograph_blast_evalue ON orthograph_blast (log_evalue)")
	db_cur.execute("CREATE INDEX orthograph_blast_hmmsearch_id ON orthograph_blast (hmmsearch_id)")
	db_cur.execute("CREATE INDEX orthograph_blast_query ON orthograph_blast (query)")
	db_cur.execute("CREATE INDEX orthograph_blast_target ON orthograph_blast (target)")
	db_cur.execute("CREATE INDEX orthograph_blast_taxid ON orthograph_blast (taxid)")
	db_cur.execute("END")
	db_con.commit()
	db_con.close()

if __name__ == '__main__':
	start = time()
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input',type=str, default='Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa',
        help='Path to directory of Input folder')
	parser.add_argument('-oi','--orthoset_input',type=str, default='Syrphidae/orthosets',
        help='Path to directory of Orthosets folder')
	parser.add_argument('-o','--orthoset',type=str, default='Ortholog_set_Mecopterida_v4',
        help='Orthoset')	
	parser.add_argument('-hs', '--hmm_minimum_score', type=float, default=40.0,
        help='Minimum score filter in hmm grab.')
	parser.add_argument('-bs', '--blast_minimum_score', type=float, default=40.0,
        help='Minimum score filter in blast.')
	parser.add_argument('-be', '--blast_minimum_evalue', type=float, default=0.00001,
        help='Minimum evalue filter in blast.')
	parser.add_argument('-ovw', '--overwrite', action="store_true",
        help='Overwrite existing blast results.')
	parser.add_argument('-p', '--processes', type=int, default=4,
        help='Number of threads used to call processes.')
	args = parser.parse_args()

	minimum_score = args.hmm_minimum_score

	input_path = args.input
	orthoset = args.orthoset
	orthosets_dir = args.orthoset_input

	orthoset_db_path = os.path.join(orthosets_dir,orthoset+'.sqlite')

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

	#Get some data from orthograph-analyzer
	db_basename = os.path.basename(input_path).split('.')[0]+'.sqlite'
	db_path = os.path.join(input_path,db_basename)

	db_con = sqlite3.connect(db_path)
	db_cur = db_con.cursor()

	gene_to_hashes = {}
	hashes_to_gene_and_id = {}
	raw_hash_to_gene = {}

	est1 = time()

	rows = db_cur.execute('SELECT digest, header, sequence FROM orthograph_ests;')
	data = {}
	genes = []

	for row in rows: 
		data[row[0]] = row[1], row[2]

	print('Grabbed EST Data. Took {:.2f}s.'.format(time()-est1))

	hmm1 = time()
	query = f'''SELECT 
					* 
				FROM 
					orthograph_hmmsearch
				WHERE
					orthograph_hmmsearch.score >= {minimum_score};
					'''
	rows = db_cur.execute(query)
	hmm_data = 0
	for row in rows:
		hmm_data += 1
		id,taxid,query,target,score,evalue,log_evalue,env_start,env_end,ali_start,ali_end,hmm_start,hmm_end = row
		header,sequence = data[target]
		if query not in gene_to_hashes:
			genes.append(query)
			gene_to_hashes[query] = []
			hashes_to_gene_and_id[query] = {}

		target_with_coords = target+'_hmmstart{}_hmmend{}'.format(str(hmm_start),str(hmm_end))

		if target_with_coords not in gene_to_hashes[query]:
			gene_to_hashes[query].append((target_with_coords,header,sequence))

		hashes_to_gene_and_id[query][target_with_coords] = query,id

	print('Grabbed HMM Data. Took: {:.2f}s. Grabbed {} rows.'.format(time()-hmm1,hmm_data))

	num_threads = args.processes

	#Run
	if num_threads == 1:
		for gene in gene_to_hashes:
			do(gene, tmp_path, gene_to_hashes[gene], blast_path, blast_db_path)
		output_thread(genes,db_path,hashes_to_gene_and_id,args.blast_minimum_score,args.blast_minimum_evalue,blast_path)
			
	else:
		p = Process(target=output_thread, args=(genes,db_path,hashes_to_gene_and_id,args.blast_minimum_score,args.blast_minimum_evalue,blast_path))
		p.start()
		
		arguments = list()
		for gene in gene_to_hashes:
			arguments.append((gene, tmp_path, gene_to_hashes[gene] ,blast_path,blast_db_path))
		with Pool(num_threads-1) as pool:
			pool.map(run_command, arguments, chunksize=1)

		p.join()

	db_con.close()
	print('Done. Took {:.2f}s.'.format(time()-start))
