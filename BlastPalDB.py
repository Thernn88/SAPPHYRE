import os
import math
from multiprocessing.pool import Pool
import argparse
from time import time,sleep
import sys
from shutil import rmtree
import json
import wrap_rocks

def do(gene, tmp_path, this_gene_sequences, blast_path, blast_db_path, blast_minimum_score, blast_minimum_evalue):
	blast_file_name = '{gene}.blast'.format(gene=gene)
	this_blast_path = os.path.join(blast_path,blast_file_name)
	result = this_blast_path + ".done"

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
		
		os.rename(this_blast_path, result)
		os.remove(target_tmp_path)

	gene_out = {}
	this_return = []

	this_content = open(result).read()
	if this_content != '':
		for line in this_content.split('\n'):
			if len(line.split('\t')) == 6:
				query_id, subject_id, evalue, bit_score, q_start, q_end = line.split('\t')

				if float(bit_score) >= blast_minimum_score and float(evalue) <= blast_minimum_evalue:
					try:
						log_evalue = str(math.log(float(evalue)))
					except: #Undefined
						log_evalue = '-999'

					gene,hmmsearch_id = hashes_to_gene_and_id[gene][query_id]

					hmmsearch_id = str(hmmsearch_id)

					query_id = query_id.split('_hmmstart')[0]

					this_out = {'target':int(subject_id),
								'score':float(bit_score),
								'evalue':float(evalue),
								'log_evalue':float(log_evalue),
								'blast_start':int(q_start),
								'blast_end':int(q_end)}

				
					if hmmsearch_id not in gene_out:
						gene_out[hmmsearch_id] = []
					gene_out[hmmsearch_id].append(this_out)

		for hmmsearch_id in gene_out:
			this_out_results = gene_out[hmmsearch_id]
			key = 'blastfor:{}'.format(hmmsearch_id)

			data = json.dumps(this_out_results)

			this_return.append((key, data, len(this_out_results)))

	print('Blasted:',gene)	
	
	return this_return

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

	print('Grabbing HMM data from db.')

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

	db_path = os.path.join(input_path,'rocksdb')
	db = wrap_rocks.RocksDB(db_path)

	gene_to_hashes = {}
	hashes_to_gene_and_id = {}

	grab_hmm_start = time()

	global_hmm_object_raw = db.get('hmmsearch:all')
	global_hmm_object = json.loads(global_hmm_object_raw)

	for hmm_id in global_hmm_object:
		key = 'hmmsearch:'+str(hmm_id)
		hmm_json = db.get(key)
		hmm_object = json.loads(hmm_json)

		header = hmm_object['header'].strip()
		sequence = hmm_object['hmm_sequence']

		gene = hmm_object['gene']
		hmm_start = hmm_object['hmm_start']
		hmm_end = hmm_object['hmm_end']
		

		target_with_coords = header+'_hmmstart{}_hmmend{}'.format(str(hmm_start),str(hmm_end))

		if gene not in gene_to_hashes:
			gene_to_hashes[gene] = []
			hashes_to_gene_and_id[gene] = {}

		hashes_to_gene_and_id[gene][target_with_coords] = (gene,hmm_id)

		gene_to_hashes[gene].append((target_with_coords,header,sequence))

	genes = list(gene_to_hashes.keys())

	print('Grabbed HMM Data. Took: {:.2f}s. Grabbed {} rows.'.format(time()-grab_hmm_start,len(global_hmm_object)))

	num_threads = args.processes
	if num_threads > 1:
		num_threads = math.floor(num_threads / 2)

	to_write = []

	#Run
	if num_threads == 1:
		for gene in genes:
			to_write.append(do(gene, tmp_path, gene_to_hashes[gene], blast_path, blast_db_path, args.blast_minimum_score, args.blast_minimum_evalue))
	else:										
		arguments = list()
		for gene in genes:
			arguments.append((gene, tmp_path, gene_to_hashes[gene],blast_path,blast_db_path, args.blast_minimum_score, args.blast_minimum_evalue))

		with Pool(num_threads) as pool:
			to_write = pool.starmap(do, arguments, chunksize=1)

	print('Writing to DB.')

	write_start = time()

	i = 0
	for batch in to_write:
		for key,data,count in batch:
			db.put(key,data)
			i += count
	
	print('Done. Took {:.2f}s. Writing {} results took {:.2f}s.'.format(time()-start,i,time()-write_start))