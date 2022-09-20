import os
import math
from multiprocessing.pool import Pool
import argparse
from time import time,sleep
from shutil import rmtree
import json
import wrap_rocks

def do(gene, tmp_path, this_gene_sequences, blast_path, blast_db_path, blast_minimum_score, blast_minimum_evalue, prog = 'blastp', evalue_threshold = .00001, score_threshold = 40, num_threads = 1):
	blast_file_name = f'{gene}.blast'
	this_blast_path = os.path.join(blast_path,blast_file_name)
	result = this_blast_path + ".done"

	if not os.path.exists(result):
		if not os.path.exists(this_blast_path):
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
							 num_threads = num_threads,
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
				query_id, subject_id, evalue, bit_score, q_start, q_end = fields

				if float(bit_score) >= blast_minimum_score and float(evalue) <= blast_minimum_evalue:
					try:
						log_evalue = str(math.log(float(evalue)))
					except: #Undefined
						log_evalue = '-999'

					query_id,hmmsearch_id = query_id.split('_hmmid')

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

def main():
	start = time()
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input',type=str, default='Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa',
		help='Path to directory of Input folder')
	parser.add_argument('-oi','--orthoset_input',type=str, default='Syrphidae/orthosets',
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

	global_hmm_object_raw = db.get('hmmsearch:all')
	global_hmm_object = json.loads(global_hmm_object_raw)

	for hmm_id in global_hmm_object:
		key = 'hmmsearch:'+str(hmm_id)
		hmm_json = db.get(key)
		hmm_object = json.loads(hmm_json)

		header = hmm_object['header'].strip()
		gene = hmm_object['gene']

		if gene not in gene_to_hits:
			gene_to_hits[gene] = []

		gene_to_hits[gene].append((header+f'_hmmid{hmm_id}', hmm_object['hmm_sequence']))

	genes = list(gene_to_hits.keys())

	print('Grabbed HMM Data. Took: {:.2f}s. Grabbed {} rows.'.format(time()-grab_hmm_start,len(global_hmm_object)))

	del global_hmm_object_raw
	del global_hmm_object

	num_threads = args.processes

	if num_threads > 1:
		num_threads = math.floor(num_threads / 2)

	to_write = []

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

	print('Writing to DB.')

	write_start = time()

	i = 0
	for batch in to_write:
		for key,data,count in batch:
			db.put(key,data)
			i += count
	
	print('Done. Took {:.2f}s. Writing {} results took {:.2f}s.'.format(time()-start,i,time()-write_start))

if __name__ == "__main__":
	main()