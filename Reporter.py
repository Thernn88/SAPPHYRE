import os
import sqlite3
from Bio.Seq import Seq
from time import time
import argparse
import math
from multiprocessing.pool import Pool
import shutil
import uuid
import json

T_global_start = time()

class NodeRecord:

	def __init__(self, header):
		self.header = header
		self.score = -math.inf
		self.cdna_start = None
		self.cdna_end = None
		self.cdna_sequence = ''
		self.aa_start = None
		self.aa_end = None
		self.aa_sequence = ''

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

def clear_output_path(input) :
	'''
	Clears protein output paths
	'''

	for protein_path in ['aa','nt']:
		protein_output = os.path.join(input,protein_path)
		for item in os.listdir(protein_output):
			item_path = os.path.join(protein_output,item)
			if '.fa' in item:
				os.remove(item_path)

def get_set_id(orthoset_db_con, orthoset) :
	'''
	Retrieves orthoset id from orthoset db
	'''
	orthoset_id = None

	orthoset_db_cur = orthoset_db_con.cursor()
	rows = orthoset_db_cur.execute('SELECT * FROM orthograph_set_details;')

	for row in rows:
		id,name,description = row

		if name == orthoset:
			orthoset_id = id

	if orthoset_id == None:
		raise Exception('Orthoset {} id cant be retrieved'.format(orthoset))
	else:
		return orthoset_id

def get_species_id(taxa_db_con, hmmsearch_table) :
	taxa_db_cur = taxa_db_con.cursor()
	rows = taxa_db_cur.execute(f'SELECT taxid FROM {hmmsearch_table};')

	species_id = None

	for row in rows:
		species_id = row[0]
		if type(species_id) == int:
			return species_id

def get_taxa_in_set(set_id, orthoset_db_con) :
	reference_taxa = []
	
	query = f'''SELECT DISTINCT {orthoset_set_details}.name, {orthoset_taxa}.name
		FROM {orthoset_seqpairs}
		INNER JOIN {orthoset_taxa}
			ON {orthoset_seqpairs}.taxid = {orthoset_taxa}.id
		INNER JOIN {orthoset_orthologs}
			ON {orthoset_orthologs}.sequence_pair = {orthoset_seqpairs}.id 
		INNER JOIN {orthoset_set_details}
			ON {orthoset_orthologs}.setid = {orthoset_set_details}.id
		WHERE {orthoset_set_details}.id = "{set_id}"'''

	orthoset_db_cur = orthoset_db_con.cursor()
	rows = orthoset_db_cur.execute(query)

	reference_taxa = [row[1] for row in rows]

	return reference_taxa

def get_orthologs_for_set_hashref(set_id, orthoset_db_con) :
	result = set()

	query = f'''SELECT DISTINCT
		{orthoset_orthologs}.ortholog_gene_id,
		{orthoset_aaseqs}.id 
		FROM {orthoset_orthologs} 
		INNER JOIN {orthoset_seqpairs} 
			ON {orthoset_orthologs}.sequence_pair = {orthoset_seqpairs}.id
		INNER JOIN {orthoset_aaseqs}
			ON {orthoset_seqpairs}.aa_seq = {orthoset_aaseqs}.id
		INNER JOIN {orthoset_set_details} 
			ON {orthoset_orthologs}.setid = {orthoset_set_details}.id
		WHERE {orthoset_set_details}.id = "{set_id}"'''

	orthoset_db_cur = orthoset_db_con.cursor()
	rows = orthoset_db_cur.execute(query)

	for row in rows:
		gene,id = row
		result.add(id)

	return result

def get_scores_list(score_threshold,taxa_db_path,orthoset_db_path,min_length,orthoset_id) :
	query = f'''SELECT DISTINCT
			{orthoset_orthologs}.{db_col_orthoid},
			{taxa_hmmsearch_table}.{db_col_id},
			{taxa_hmmsearch_table}.{db_col_target},
			{taxa_hmmsearch_table}.{db_col_score},
			{taxa_hmmsearch_table}.{db_col_evalue},
			{taxa_hmmsearch_table}.{db_col_hmm_start},
			{taxa_hmmsearch_table}.{db_col_hmm_end},
			{taxa_hmmsearch_table}.{db_col_ali_start},
			{taxa_hmmsearch_table}.{db_col_ali_end},
			{taxa_hmmsearch_table}.{db_col_env_start},
			{taxa_hmmsearch_table}.{db_col_env_end},
			{taxa_ests_table}.{db_col_header}
		FROM {taxa_hmmsearch_table}
		LEFT JOIN {taxa_ests_table}
			ON {taxa_hmmsearch_table}.{db_col_target} = {taxa_ests_table}.{db_col_digest}
		LEFT JOIN {orthoset_orthologs}
			ON {taxa_hmmsearch_table}.{db_col_query} = {orthoset_orthologs}.{db_col_orthoid}
		LEFT JOIN {taxa_species_info}
			ON {taxa_hmmsearch_table}.{db_col_taxid} = {taxa_species_info}.{db_col_id}
		LEFT JOIN {orthoset_set_details}
			ON {orthoset_orthologs}.{db_col_setid} = {orthoset_set_details}.{db_col_id}
		WHERE {taxa_ests_table}.{db_col_digest}		  IS NOT NULL
			AND {orthoset_orthologs}.{db_col_orthoid}	IS NOT NULL
			AND {taxa_species_info}.{db_col_id}			  IS NOT NULL
			AND {orthoset_set_details}.{db_col_id}	   IS NOT NULL
			AND {orthoset_set_details}.{db_col_id}	   = {orthoset_id}
			AND {taxa_hmmsearch_table}.{db_col_score} > {score_threshold}'''

	taxa_db_con = sqlite3.connect(taxa_db_path)
	taxa_db_cur = taxa_db_con.cursor()

	taxa_db_cur.execute(f'ATTACH DATABASE "{orthoset_db_path}" as "species_database";')

	rows = taxa_db_cur.execute(query)

	gene_based_result = {}
	header_based_results = {}
	
	ufr_out = [['Gene','Hash','Header','Score','Start','End']]

	for row in rows:
		orthoid,hmmsearch_id,hmmhit,hmm_score,hmm_evalue,hmm_start,hmm_end,ali_start,ali_end,env_start,env_end,header = row
		length = env_end-env_start

		if length > min_length:
			this_row = {'orthoid':orthoid, 'hmmsearch_id':hmmsearch_id, 'hmmhit':hmmhit, 'hmm_score': hmm_score, 'hmm_evalue':hmm_evalue, 'hmm_start':hmm_start, 'hmm_end':hmm_end,
					'ali_start':ali_start, 'ali_end':ali_end, 'env_start':env_start, 'env_end':env_end, 'header':header}

			ufr_out.append([orthoid,hmmhit,header,str(hmm_score),str(hmm_start),str(hmm_end)])

			if orthoid not in gene_based_result:
				gene_based_result[orthoid] = []
			gene_based_result[orthoid].append(this_row)
			
			if header not in header_based_results:
				header_based_results[header] = []
			header_based_results[header].append(this_row)
			
	ufr_out = sorted(ufr_out, key = lambda x: (x[0], x[1], x[3], x[4], x[5]))

	return gene_based_result, header_based_results, ufr_out

def get_blastresults_for_hmmsearch_ids(taxa_db_con,hmmsearch_id) :
	result = []
	taxa_db_cur = taxa_db_con.cursor()

	query = f'''SELECT DISTINCT
			{taxa_hmmsearch_table}.{db_col_id},
			{taxa_blast_table}.{db_col_target},
			{taxa_blast_table}.{db_col_score},
			{taxa_blast_table}.{db_col_evalue},
			{taxa_blast_table}.{db_col_start},
			{taxa_blast_table}.{db_col_end},
			{taxa_hmmsearch_table}.{db_col_target},
			{taxa_hmmsearch_table}.{db_col_evalue},
			{taxa_hmmsearch_table}.{db_col_env_start},
			{taxa_hmmsearch_table}.{db_col_env_end},
			{taxa_hmmsearch_table}.{db_col_ali_start},
			{taxa_hmmsearch_table}.{db_col_ali_end},
			{taxa_hmmsearch_table}.{db_col_hmm_start},
			{taxa_hmmsearch_table}.{db_col_hmm_end},
			{taxa_ests_table}.{db_col_header}
		FROM 
			{taxa_hmmsearch_table}
		LEFT JOIN {taxa_blast_table} ON {taxa_hmmsearch_table}.{db_col_id} = {taxa_blast_table}.{db_col_hmmsearch_id}
		LEFT JOIN {taxa_ests_table} ON {taxa_ests_table}.{db_col_digest} = {taxa_hmmsearch_table}.{db_col_target}
		WHERE 
			{taxa_hmmsearch_table}.{db_col_id} IS NOT NULL
			AND {taxa_blast_table}.{db_col_hmmsearch_id} IS NOT NULL
			AND {taxa_hmmsearch_table}.{db_col_target} IS NOT NULL
			AND {taxa_hmmsearch_table}.{db_col_id} = ?
		ORDER BY {taxa_blast_table}.{db_col_score} DESC'''

	rows = taxa_db_cur.execute(query,(hmmsearch_id,))
	for row in rows:
		hmmsearch_id,blast_hit,blast_score,blast_evalue,blast_start,blast_end,hmmhit,hmm_evalue,env_start,env_end,ali_start,ali_end,hmm_start,hmm_end,header = row
		this_row = {'blast_hit':blast_hit,'blast_score':blast_score,'blast_evalue':blast_evalue,'blast_start':blast_start,'blast_end':blast_end,'hmmhit':hmmhit,'hmm_evalue':hmm_evalue,
					'env_start':env_start,'env_end':env_end,'ali_start':ali_start,'ali_end':ali_end,'hmm_start':hmm_start,'hmm_end':hmm_end,'header':header}

		result.append(this_row)
	
	return result

def get_reftaxon_name(hit_id, taxa_db_path, orthoset_db_path) :
	query = f"SELECT {orthoset_taxa}.{db_col_name} FROM {orthoset_taxa} INNER JOIN {orthoset_aaseqs} ON {orthoset_taxa}.{db_col_id} = {orthoset_aaseqs}.{db_col_taxid} WHERE {orthoset_aaseqs}.{db_col_id} = ?"

	taxa_db_con = sqlite3.connect(taxa_db_path)
	taxa_db_cur = taxa_db_con.cursor()

	taxa_db_cur.execute(f'ATTACH DATABASE "{orthoset_db_path}" as "species_database";')

	rows = taxa_db_cur.execute(query, (hit_id,))

	for row in rows:
		return row[0]

def is_reciprocal_match(orthoid, blast_results, taxa_db_path, orthoset_db_path, reference_taxa, aaseqs_in_orthoid = None):
	reftaxon_count = {ref_taxa:0 for ref_taxa in reference_taxa}

	for i,result in enumerate(blast_results):
		if aaseqs_in_orthoid != None:
			if not result['blast_hit'] in aaseqs_in_orthoid:
				continue

		reftaxon = get_reftaxon_name(result['blast_hit'], taxa_db_path, orthoset_db_path);
		match = result

		if reftaxon in reftaxon_count:
			reftaxon_count[reftaxon] = 1

			if strict_search_mode:
				total_count = sum([reftaxon_count[reftaxon] for reftaxon in reference_taxa])
				if total_count == len(reference_taxa):
					return match
			else:
				return match

	return None

def calculate_length_of_ali(result) :
	return abs(result['ali_end'] - result['ali_start']) + 1

def transcript_not_long_enough(result, minimum_transcript_length) :
	length = calculate_length_of_ali(result)

	if length < minimum_transcript_length:
		return True
	else:
		return False

def coordinate_overlap(other_hits, new_hit):
	for hit in other_hits:
		if 'ali_start' in hit:
			starts_before_ends_within = new_hit['ali_start'] < hit['ali_start'] and new_hit['ali_end'] > hit['ali_start']
			starts_before_ends_after = new_hit['ali_start'] < hit['ali_start'] and new_hit['ali_end'] > hit['ali_end']
			starts_within_ends_within = new_hit['ali_start'] > hit['ali_start'] and new_hit['ali_end'] < hit['ali_end']

			if starts_before_ends_within or starts_before_ends_after or starts_within_ends_within:
				return True
	return False

def get_reference_sequence(hit_id, orthoset_db_con) :
	orthoset_db_cur = orthoset_db_con.cursor()
	
	query = f'''SELECT
			{orthoset_aaseqs}.{db_col_sequence}, 
			{orthoset_taxa}.{db_col_name}
		FROM {orthoset_aaseqs}
		INNER JOIN {orthoset_taxa}
			ON {orthoset_aaseqs}.{db_col_taxid} = {orthoset_taxa}.id
		WHERE {orthoset_aaseqs}.{db_col_id} = "{hit_id}"'''

	rows = orthoset_db_cur.execute(query)

	for row in rows:
		return (row[0],row[1])

def reverse_complement(nt_seq) :
	seq = Seq(nt_seq)
	return str(seq.reverse_complement())

def get_nucleotide_transcript_for(header, taxa_db_con):
	formatted_header = header.split(' ')[0]
	query = f'''SELECT {db_col_sequence}
		FROM {taxa_ests_table}
		WHERE {db_col_header} = ?'''

	taxa_db_cur = taxa_db_con.cursor()

	rows = taxa_db_cur.execute(query,(formatted_header,))

	revcomp = 'revcomp' in header
	

	for row in rows:
		if revcomp:
			return formatted_header,reverse_complement(row[0])
		else:
			return formatted_header,row[0]

def crop_to_hmm_alignment(seq, header, hit):
	start = hit['ali_start']-1 #Adjust for zero based number
	start = start * 3
	end = hit['ali_end'] * 3


	return seq[start:end]

def get_transcript_for(hmmhit, taxa_db_con) :
	query = f'''SELECT {db_col_sequence}
		FROM {taxa_ests_table}
		WHERE {db_col_digest} = ?'''
	
	taxa_db_cur = taxa_db_con.cursor()

	rows = taxa_db_cur.execute(query,(hmmhit,))

	for row in rows:
		return row[0]

### EXONERATE FUNCTIONS

def parse_results(result_file_content) :
	
	cdna_sequence = ''
	aa_sequence = ''
	cdna_start = None
	cdna_end = None
	aa_start = None
	aa_end = None

	lines = result_file_content.split('\n')
	current_seq = ''
	for line in lines:
		if '>' in line:
			header,start,end = line.split(' ')
			current_seq = header.replace('>','')
			if current_seq == 'cdna':
				cdna_start = start + 1
				cdna_end = end
			elif current_seq == 'aa':
				aa_start = start + 1
				aa_end = end
		else:
			if current_seq == 'cdna':
				cdna_sequence += line
			elif current_seq == 'aa':
				aa_sequence += line

	if cdna_sequence == '':
		cdna_sequence = None
	if aa_sequence == '':
		aa_sequence = None

	return (cdna_sequence,cdna_start,cdna_end,aa_sequence,aa_start,aa_end)

def fastaify(headers, sequences, input_path, tmp_path) :
	name = str(uuid.uuid4())+'.tmp' #Super unique name
	path = os.path.join(tmp_path,name)
	
	this_out = []
	for i,seq in enumerate(sequences):
		this_out.append('>'+headers[i])
		this_out.append(seq)

	open(path,'w').write('\n'.join(this_out))

	return path

def translate_cdna(cdna_seq, nt_table):
	if cdna_seq == None:
		return None

	if len(cdna_seq) % 3 != 0:
		print('WARNING: NT Sequence length is not divisable by 3')
	
	split_seq = [cdna_seq[x:x+3] for x in range(0, len(cdna_seq), 3)]
	translation = ''.join([nt_table[i] for i in split_seq])

	return translation

def parse_nodes(lines) :
	"""
	vulgar => cdna => aa => vulgar
	0 => 1 => 2 => 0
	"""
	nodes = dict()
	this_node = None
	state = 0
	for line in lines:
		if "vulgar:" in line:  # if vulgar line, start new record
			current_header = line.split(' . ')[1].split(' ')[0].replace('|',' ')
			this_node = NodeRecord(current_header)
			this_node.score = int(line.split()[9])
			previously_seen = nodes.get(current_header, False)
			if previously_seen:
				nodes[current_header] = max(this_node, previously_seen)
			else:
				nodes[current_header] = this_node
		elif ">cdna" in line:  # if cdna, set cdna values
			fields = line.split()
			this_node.cdna_start = int(fields[1])
			this_node.cdna_end = int(fields[2])
			state = 1
		elif ">aa" in line:  # if aa in line, set aa values
			fields = line.split()
			this_node.aa_start = int(fields[1])
			this_node.aa_end = int(fields[2])
			state = 2
		elif state == 1:  # if cdna data line
			this_node.cdna_sequence += line
		elif state == 2:  # if aa data line
			this_node.aa_sequence += line
	return nodes

def parse_multi_results(result_file_content):
	result = []
	lines = result_file_content.split('\n')
	nodes = parse_nodes(lines)
	for key in nodes:
		node = nodes[key]
		result.append((node.header, node.cdna_sequence, node.cdna_start, node.cdna_end, node.aa_sequence, node.aa_start, node.aa_end))
	return result

def get_multi_orf(query, targets, input_path, score_threshold, tmp_path, extend = False):
	exonerate_ryo = '>cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas'
	genetic_code = 1
	exonerate_model = 'protein2genome'
	exhaustive = ''

	headers = [i['header'].replace(' ','|') for i in targets]
	if extend:
		sequences = [i['est_sequence_complete'] for i in targets]
	else:
		sequences = [i['est_sequence_hmm_region'] for i in targets]

	queryfile = fastaify(['query'], [query], input_path, tmp_path)
	targetfile = fastaify(headers, sequences, input_path, tmp_path)
	
	outfile = os.path.join(tmp_path,str(uuid.uuid4())+'.exonerateout')
	
	exonerate_cmd = f"exonerate --score {score_threshold} --ryo '{exonerate_ryo}' --subopt 0 --geneticcode {genetic_code} --model '{exonerate_model}' --querytype 'protein' --targettype 'dna' --verbose 0 --showalignment 'no' --showvulgar 'yes' --query '{queryfile}' --target '{targetfile}' > {outfile}"

	os.system(exonerate_cmd)

	result_content = open(outfile).read()

	results = parse_multi_results(result_content)

	os.remove(queryfile)
	os.remove(targetfile)
	os.remove(outfile)

	return results

def extended_orf_contains_original_orf(hit):
	if hit['extended_orf_cdna_start'] > hit['orf_cdna_end_on_transcript'] or hit['extended_orf_cdna_end'] < hit['orf_cdna_start_on_transcript']:
		return False
	else:
		return True

def get_overlap_length(candidate):
	if candidate['extended_orf_aa_start_on_transcript'] <= candidate['ali_end'] and candidate['extended_orf_aa_end_on_transcript'] >= candidate['ali_start']:
		overlap_start = candidate['extended_orf_aa_start_on_transcript'] if candidate['extended_orf_aa_start_on_transcript'] > candidate['ali_start'] else candidate['ali_start']
		overlap_end = candidate['extended_orf_aa_end_on_transcript'] if candidate['extended_orf_aa_end_on_transcript'] > candidate['ali_end'] else candidate['ali_end']

		overlap_length = overlap_end - overlap_start
		return overlap_length
	else:
		return 0

def overlap_by_orf(candidate):
	orf_length = abs(candidate['extended_orf_aa_end'] - candidate['extended_orf_aa_start'])
	overlap_length = get_overlap_length(candidate)

	return overlap_length / orf_length

def remove_extended_orf(hit):
	hit.pop('extended_orf_aa_sequence')
	hit.pop('extended_orf_aa_start')
	hit.pop('extended_orf_aa_end')
	hit.pop('extended_orf_cdna_sequence')
	hit.pop('extended_orf_cdna_start')
	hit.pop('extended_orf_cdna_end')

	return hit

### FIN

def get_ortholog_group(orthoset_id, orthoid, orthoset_db_con):
	query = f'''SELECT 
			{orthoset_taxa}.{db_col_name},
			{orthoset_aaseqs}.{db_col_header},
			{orthoset_aaseqs}.{db_col_sequence}
		FROM {orthoset_aaseqs}
		INNER JOIN {orthoset_seqpairs}
			ON {orthoset_aaseqs}.{db_col_id} = {orthoset_seqpairs}.{db_col_aaseq}
		INNER JOIN {orthoset_orthologs}
			ON {orthoset_seqpairs}.{db_col_id} = {orthoset_orthologs}.{db_col_seqpair}
		INNER JOIN {orthoset_taxa}
			ON {orthoset_aaseqs}.{db_col_taxid} = {orthoset_taxa}.{db_col_id}
		AND   {orthoset_orthologs}.{db_col_setid} = ?
		AND   {orthoset_orthologs}.{db_col_orthoid} = ?
		ORDER BY {orthoset_taxa}.{db_col_name}'''

	orthoset_db_cur = orthoset_db_con.cursor()
	rows = orthoset_db_cur.execute(query,(orthoset_id,orthoid,))

	return rows

def format_transcript_header(hdr,start,end,rf,reftax):
	return header_seperator.join([hdr,f'{start}-{end}'.replace('.0',''),rf,reftax])

def format_header(orthoid,reftaxon,taxon,header):
	return header_seperator.join([orthoid,taxon,header])

def print_core_sequences(orthoid, core_sequences):
	core_sequences = sorted(core_sequences)

	result = []
	for core in core_sequences:
		transcript_header = format_transcript_header(core[1],1,len(core[2]),'.','.')
		header = format_header(orthoid,'.',core[0],transcript_header)
		sequence = core[2]

		result.append('>'+header+'\n')
		result.append(sequence+'\n')
	
	return result

def get_rf(header):
	raw_header = header.split(' ')[0]
	header_part = header.split(' ')[1]
	if 'translate' in header_part or 'revcomp' in header_part:
		frame = header_part
	
	return raw_header,frame

def print_unmerged_sequences(hits, orthoid, type, minimum_seq_data_length, species_name, kicks) :
	result = []
	kicks_result = set()
	for i,hit in enumerate(hits):
		this_hdr, tmp_rf = get_rf(hit['header'])

		rf = tmp_rf if type != 'nt' else '.'

		start = hit['extended_orf_aa_start_on_transcript'] if 'extended_orf_aa_start_on_transcript' in hit else hit['orf_aa_start_on_transcript']
		end = hit['extended_orf_aa_end_on_transcript'] if 'extended_orf_aa_end_on_transcript' in hit else hit['orf_aa_end_on_transcript']

		transcript_hdr = format_transcript_header(this_hdr,round(start),round(end),rf,hit['reftaxon'])

		hdr = format_header(orthoid, hit['reftaxon'],species_name,transcript_hdr)

		if type == 'nt':
			seq = hit['extended_orf_cdna_sequence'] if 'extended_orf_cdna_sequence' in hit else hit['orf_cdna_sequence']
		elif type == 'aa':
			seq = hit['extended_orf_aa_sequence'] if 'extended_orf_aa_sequence' in hit else hit['orf_aa_sequence']

		constant_header = '|'.join(hdr.split('|')[:4])

		if type == 'aa':
			if len(seq) - seq.count('-') > minimum_seq_data_length:
				result.append('>'+hdr+'\n')
				result.append(seq+'\n')
			else:
				kicks_result.add(i)
		elif type == 'nt':
			if i not in kicks:
				result.append('>'+hdr+'\n')
				result.append(seq+'\n')
		
	return kicks_result,result

def get_ortholog_group_nucleotide(orthoset_id, orthoid, orthoset_db_con):
	query = f'''SELECT 
			{orthoset_taxa}.{db_col_name},
			{orthoset_ntseqs}.{db_col_header},
			{orthoset_ntseqs}.{db_col_sequence}
		FROM {orthoset_ntseqs}
		INNER JOIN {orthoset_seqpairs}
			ON {orthoset_ntseqs}.{db_col_id} = {orthoset_seqpairs}.{db_col_ntseq}
		INNER JOIN {orthoset_orthologs}
			ON {orthoset_seqpairs}.{db_col_id} = {orthoset_orthologs}.{db_col_seqpair}
		INNER JOIN {orthoset_taxa}
			ON {orthoset_ntseqs}.{db_col_taxid} = {orthoset_taxa}.{db_col_id}
		AND   {orthoset_orthologs}.{db_col_setid} = ?
		AND   {orthoset_orthologs}.{db_col_orthoid} = ?
		ORDER BY {orthoset_taxa}.{db_col_name}'''

	orthoset_db_cur = orthoset_db_con.cursor()
	rows = orthoset_db_cur.execute(query,(orthoset_id,orthoid,))

	return rows

def get_difference(scoreA, scoreB) :
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

def get_baseheader(header) :
	"""
	Returns header content before first whitespace.
	"""
	baseheader = header.split(' ')[0]
	return baseheader

def run_exonerate(arg_tuple):
    gene,list_of_hits,taxa_db_path,orthoset_db_path,input_path,min_score,orthoset_id,aa_out_path,min_length,species_name,nt_out_path,tmp_path,exonerate_verbose = arg_tuple
    exonerate_gene_multi(gene,list_of_hits,taxa_db_path,orthoset_db_path,input_path,min_score,orthoset_id,aa_out_path,min_length,species_name,nt_out_path,tmp_path,exonerate_verbose)

def exonerate_gene_multi(orthoid,list_of_hits,taxa_db_path,orthoset_db_path,input_path,min_score,orthoset_id,aa_out_path,min_length,species_name,nt_out_path,tmp_path,exonerate_verbose):
	T_gene_start = time()

	summary_out_path = os.path.join(tmp_path,'{}-summary.txt'.format(orthoid))
	summary_out_log = []

	nt_table =  {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
			'TTA': 'L',
			'TCA': 'S', 'TAA': '*', 'TGA': '*', 'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L',
			'CCT': 'P',
			'CAT': 'H', 'CGT': 'R', 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L', 'CCA': 'P',
			'CAA': 'Q',
			'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N',
			'AGT': 'S',
			'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
			'ATG': 'M',
			'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 'GTC': 'V',
			'GCC': 'A',
			'GAC': 'D', 'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GTG': 'V', 'GCG': 'A',
			'GAG': 'E',
			'GGG': 'G'}

	taxa_db_con = sqlite3.connect(taxa_db_path)
	orthoset_db_con = sqlite3.connect(orthoset_db_path)

	if exonerate_verbose:
		print('Exonerating and doing output for: ',orthoid)
	reftaxon_related_transcripts = {} 
	reftaxon_to_proteome_sequence = {}
	for hit in list_of_hits:
		proteome_sequence,this_reftaxon = get_reference_sequence(hit['blast_hit'], orthoset_db_con)

		reftaxon_to_proteome_sequence[this_reftaxon] = proteome_sequence

		hit['reftaxon'] = this_reftaxon
		hit['proteome_sequence'] = proteome_sequence

		est_header,est_sequence_complete = get_nucleotide_transcript_for(hit['header'],taxa_db_con)
		est_sequence_hmm_region = crop_to_hmm_alignment(est_sequence_complete,est_header,hit)

		hit['est_header'] = est_header
		hit['est_sequence_complete'] = est_sequence_complete
		hit['est_sequence_hmm_region'] = est_sequence_hmm_region

		if this_reftaxon not in reftaxon_related_transcripts:
			reftaxon_related_transcripts[this_reftaxon] = []
		
		reftaxon_related_transcripts[this_reftaxon].append(hit)

	output_sequences = []

	total_results = 0

	for taxon_hit in reftaxon_related_transcripts:
		hits = reftaxon_related_transcripts[taxon_hit]
		query = reftaxon_to_proteome_sequence[taxon_hit]
		
		results = get_multi_orf(query, hits, input_path, min_score, tmp_path)

		total_results += len(hits)

		if extend_orf:
			extended_results = get_multi_orf(query, hits, input_path, min_score, tmp_path, extend = True)

		for hit in hits:
			matching_alignment = [i for i in results if i[0] == hit['header']]

			if len(matching_alignment) == 0:
				orf_aa_sequence,orf_cdna_sequence,orf_cdna_start,orf_cdna_end,orf_aa_start,orf_aa_end = None,None,None,None,None,None
			else:
				current_header,orf_cdna_sequence,orf_cdna_start,orf_cdna_end,orf_aa_sequence,orf_aa_start,orf_aa_end = matching_alignment[0]

				orf_aa_sequence = translate_cdna(orf_cdna_sequence, nt_table)

				hit['orf_aa_sequence'] = orf_aa_sequence
				hit['orf_cdna_sequence'] = orf_cdna_sequence
				hit['orf_cdna_start'] = orf_cdna_start
				hit['orf_cdna_end'] = orf_cdna_end
				hit['orf_aa_start'] = orf_aa_start
				hit['orf_aa_end'] = orf_aa_end

				if orf_cdna_sequence == None:
					continue

				else:
					hit['orf_cdna_start_on_transcript'] = hit['orf_cdna_start'] + ( hit['ali_start'] * 3 ) - 3
					hit['orf_cdna_end_on_transcript']   = hit['orf_cdna_end'] + ( hit['ali_start'] * 3 ) - 3
					hit['orf_aa_start_on_transcript']   = ( hit['orf_cdna_start'] + ( hit['ali_start'] * 3 ) - 3 ) / 3
					hit['orf_aa_end_on_transcript']	 = ( hit['orf_cdna_end'] + ( hit['ali_start'] * 3 ) - 3 ) / 3

					if extend_orf:
						matching_extended_alignment = [i for i in extended_results if i[0] == hit['header']]
						if matching_extended_alignment != []:
							current_header,extended_orf_cdna_sequence,extended_orf_cdna_start,extended_orf_cdna_end,extended_orf_aa_sequence,extended_orf_aa_start,extended_orf_aa_end = matching_extended_alignment[0]

							extended_orf_aa_start_on_transcript = ( extended_orf_cdna_start - 1 ) / 3 + 1
							extended_orf_aa_end_on_transcript   = ( extended_orf_cdna_end   - 1 ) / 3 + 1

							extended_orf_aa_sequence = translate_cdna(extended_orf_cdna_sequence, nt_table)

							hit['extended_orf_aa_sequence'] = extended_orf_aa_sequence
							hit['extended_orf_cdna_sequence'] = extended_orf_cdna_sequence
							hit['extended_orf_cdna_start'] = extended_orf_cdna_start
							hit['extended_orf_cdna_end'] = extended_orf_cdna_end
							hit['extended_orf_aa_start'] = extended_orf_aa_start
							hit['extended_orf_aa_end'] = extended_orf_aa_end
							hit['extended_orf_aa_start'] = extended_orf_aa_start
							hit['extended_orf_aa_end'] = extended_orf_aa_end

							hit['extended_orf_aa_start_on_transcript'] = extended_orf_aa_start_on_transcript
							hit['extended_orf_aa_end_on_transcript']   = extended_orf_aa_end_on_transcript

							if extended_orf_contains_original_orf(hit):
								orf_overlap = overlap_by_orf(hit)
								if orf_overlap < orf_overlap_minimum:
									hit = remove_extended_orf(hit)
							
							else:
								hit = remove_extended_orf(hit)
						else:
							print('Failed to extend orf on {}'.format(hit['header']))

					output_sequences.append(hit)

	if len(output_sequences) > 0:
		output_sequences = sorted(output_sequences, key=lambda d: d['hmm_start'])

		this_out = [orthoid]
		for hit in output_sequences:
			length = len(hit['extended_orf_aa_sequence']) if 'extended_orf_aa_sequence' in hit else len(hit['orf_aa_sequence'])
			this_out.append('{}[{} aa]'.format(hit['header'],length))
		
		summary_out_log.append('\t'.join(this_out))
		
		core_sequences = get_ortholog_group(orthoset_id, orthoid, orthoset_db_con)

		this_aa_out = []
		this_aa_path = os.path.join(aa_out_path,orthoid+'.aa.fa')
		this_aa_out.extend(print_core_sequences(orthoid,core_sequences))
		kicks,output = print_unmerged_sequences(output_sequences, orthoid, 'aa', min_length, species_name, kicks=set())
		this_aa_out.extend(output)

		open(this_aa_path, 'w').writelines(this_aa_out)

		core_sequences_nt = get_ortholog_group_nucleotide(orthoset_id, orthoid, orthoset_db_con)

		this_nt_out = []
		this_nt_path = os.path.join(nt_out_path,orthoid+'.nt.fa')

		this_nt_out.extend(print_core_sequences(orthoid,core_sequences_nt))
		na_kicks,output = print_unmerged_sequences(output_sequences, orthoid, 'nt', min_length, species_name, kicks=kicks)
		this_nt_out.extend(output)

		orthoid_summary_out = []
		for hit in output_sequences:
			if 'orf_aa_sequence' in hit:
				if hit['orf_aa_sequence'] != None:
					header = hit['header']
					length = len(hit['extended_orf_aa_sequence']) if 'extended_orf_aa_sequence' in hit else len(hit['orf_aa_sequence'])
					orthoid_summary_out.append((header,length))

		open(this_nt_path, 'w').writelines(this_nt_out)

	open(summary_out_path,'w').write('\n'.join(summary_out_log))

	if exonerate_verbose:
		print('{} took {:.2f}s. Had {} sequences'.format(orthoid,time()-T_gene_start,len(output_sequences)))
		
def run_internal_filter(this_gene_transcripts, orthoid, min_overlap_internal, score_diff_internal, filter_verbose,input_path,tmp_path,internal_verbose):
    return internal_filter_gene(this_gene_transcripts, orthoid, min_overlap_internal, score_diff_internal, filter_verbose,input_path,tmp_path,internal_verbose)

def internal_filter_gene(this_gene_hits, gene, min_overlap_internal, score_diff_internal, filter_verbose, input_path, tmp_path,internal_verbose):
	if internal_verbose:
		T_internal_start = time()
		print('Checking for internal dupes in {}.'.format(gene))

	descending_hits = sorted(this_gene_hits, key=lambda hit: hit['hmm_score'], reverse=True)
	ascending_hits = sorted(this_gene_hits, key=lambda hit: hit['hmm_score'])
	removed_hits = set()
	already_passed = set()
	filtered_sequences_log = []
	this_gene_passes = this_gene_hits.copy()

	for hit_a in descending_hits:
		for hit_b in ascending_hits:
			if not hit_b['uuid'] in already_passed: #Not removed yet
				if hit_a != hit_b:
					if get_baseheader(hit_a['header']) != get_baseheader(hit_b['header']): #Not the same sequence
						if hit_a['hmm_score'] > hit_b['hmm_score']:
							rangeA = range(hit_a['hmm_start'], hit_a['hmm_end'] + 1)  # Adjusted for range starting at 0
							rangeB = range(hit_b['hmm_start'], hit_b['hmm_end'] + 1)

							overlap = set(rangeA).intersection(set(rangeB))
							amount_of_overlap = len(overlap)
							percentage_of_overlap = amount_of_overlap / len(rangeB)
							if percentage_of_overlap >= min_overlap_internal:
								score_difference = get_difference(hit_a['hmm_score'], hit_b['hmm_score'])
								if score_difference >= score_diff_internal:
									#removed_hits.add(hit_b['uuid'])
									descending_hits.remove(hit_b)
									ascending_hits.remove(hit_b)
									this_gene_passes.remove(hit_b)
									if filter_verbose: 
										filtered_sequences_log.append([hit_b['orthoid'],hit_b['hmmhit'],hit_b['header'],str(hit_b['hmm_score']),str(hit_b['hmm_start']),str(hit_b['hmm_end']),'Internal Overlapped with Lowest Score',hit_a['orthoid'],hit_a['hmmhit'],hit_a['header'],str(hit_a['hmm_score']),str(hit_a['hmm_start']),str(hit_a['hmm_end'])])
								else:
									already_passed.add(hit_b['uuid'])
	
	this_out_data = {'Passes':this_gene_passes, 'Log':filtered_sequences_log, 'gene':gene}

	if internal_verbose:
		print('Found {} dupes in {}. Took {:.2f}s'.format(len(removed_hits),gene,time()-T_internal_start))

	return this_out_data

def run_multi_filter(header,this_hits,filter_verbose,min_overlap_multi,score_diff_multi,tmp_path,multi_verbose):
    return multi_filter_dupes(header,this_hits,filter_verbose,min_overlap_multi,score_diff_multi,tmp_path,multi_verbose)

def multi_filter_dupes(header,this_hits,filter_verbose,min_overlap_multi,score_diff_multi,tmp_path,multi_verbose):
	if multi_verbose:
		T_multi_check = time()
		print('Checking for multi-gene dupes in {}'.format(header))
	
	kick_happend = True
	filtered_sequences_log = []
	this_kicks = []

	this_hits.sort(key=lambda data: (data['hmm_score'],data['orthoid'])) #ascending
	this_hits.reverse() #descending

	while kick_happend:
		kick_happend = False

		master = this_hits[0]
		candidates = this_hits[1:]

		master_env_start = master['env_start']
		master_env_end = master['env_end']

		for candidate in candidates:
			if candidate['orthoid'] == master['orthoid']: #From same gene = Pseudomaster
				# Remove
				this_hits.remove(candidate)
				this_kicks.append(candidate)
				if filter_verbose: 
					filtered_sequences_log.append([candidate['orthoid'],candidate['hmmhit'],candidate['header'],str(candidate['hmm_score']),str(candidate['hmm_start']),str(candidate['hmm_end']),'Pseudomaster',master['orthoid'],master['hmmhit'],master['header'],str(master['hmm_score']),str(master['hmm_start']),str(master['hmm_end'])])
				candidates.remove(candidate)
				# Extend master range
				kick_happend = True
				if candidate['env_start'] < master_env_start:
					master_env_start = candidate['env_start']
				if candidate['env_end'] > master_env_start:
					master_env_end = candidate['env_end']
			else:
				break
		
		miniscule_score = False

		for candidate in candidates:
			rangeA = range(master_env_start, master_env_end + 1)  # Adjusted for range starting at 0
			rangeB = range(candidate['env_start'], candidate['env_end'] + 1)

			overlap = set(rangeA).intersection(set(rangeB))
			amount_of_overlap = len(overlap)
			percentage_of_overlap = amount_of_overlap / len(rangeA)
		
			if percentage_of_overlap >= min_overlap_multi:
				score_difference = get_difference(master['hmm_score'], candidate['hmm_score'])
				
				if score_difference >= score_diff_multi:
					kick_happend = True
					this_hits.remove(candidate)
					this_kicks.append(candidate)
					if filter_verbose: 
						filtered_sequences_log.append([candidate['orthoid'],candidate['hmmhit'],candidate['header'],str(candidate['hmm_score']),str(candidate['env_start']),str(candidate['env_end']),'Multi Overlapped with Lowest Score',master['orthoid'],master['hmmhit'],master['header'],str(master['hmm_score']),str(master['env_start']),str(master['env_end'])])
				else:
					miniscule_score = True
					break
		
		if miniscule_score: #Remove all overlapping candidates if it's score is a miniscule difference of the masters
			for candidate in candidates:
				rangeA = range(master_env_start, master_env_end + 1)  # Adjusted for range starting at 0
				rangeB = range(candidate['env_start'], candidate['env_end'] + 1)
				overlap = set(rangeA).intersection(set(rangeB))
				amount_of_overlap = len(overlap)
				percentage_of_overlap = amount_of_overlap / len(rangeA)
				if percentage_of_overlap >= min_overlap_multi:
					kick_happend = True
					this_hits.remove(candidate)
					this_kicks.append(candidate)
					if filter_verbose: 
						filtered_sequences_log.append([candidate['orthoid'],candidate['hmmhit'],candidate['header'],str(candidate['hmm_score']),str(candidate['hmm_start']),str(candidate['hmm_end']),'Multi Overlapped with Miniscule Score',master['orthoid'],master['hmmhit'],master['header'],str(master['hmm_score']),str(master['hmm_start']),str(master['hmm_end'])])
	
	multi_data = {'Kicks':this_kicks, 'Log':filtered_sequences_log, 'Remaining':this_hits}
	if multi_verbose:
		print('Found {} dupes for {}. Took {:.2f}s.'.format(header,len(this_kicks),time()-T_multi_check))
	return multi_data

def run_reciprocal_search(hmmresults,list_of_wanted_orthoids,taxa_db_con,orthoset_db_path,reference_taxa, score, min_length, aaseqs_in_orthoid, tmp_path,reciprocal_verbose):
	return reciprocal_search(hmmresults,list_of_wanted_orthoids,taxa_db_con,orthoset_db_path,reference_taxa, score, min_length, aaseqs_in_orthoid, tmp_path,reciprocal_verbose)

def reciprocal_search(hmmresults,list_of_wanted_orthoids,taxa_db_con,orthoset_db_path,reference_taxa, score, min_length, aaseqs_in_orthoid, tmp_path,reciprocal_verbose):
	if reciprocal_verbose:
		T_reciprocal_start = time()
		print('Ensuring reciprocal hit for hmmresults in {}'.format(score))
	taxa_db_con = sqlite3.connect(taxa_db_path)
	results = []
	this_fails = []
	for result in hmmresults:
		orthoid = result['orthoid']

		if list_of_wanted_orthoids != []:
			if orthoid not in list_of_wanted_orthoids:
				continue

		result_hmmsearch_id = result['hmmsearch_id']

		blast_results = get_blastresults_for_hmmsearch_ids(taxa_db_con,result_hmmsearch_id)
		x = len(blast_results)
		blast_results = [i for i in blast_results if i['blast_hit'] in aaseqs_in_orthoid]

		this_match = is_reciprocal_match(orthoid, blast_results,taxa_db_path,orthoset_db_path,reference_taxa)

		if this_match == None:
			this_fails.append([result['orthoid'],result['hmmhit'],result['header'],str(result['hmm_score']),str(result['hmm_start']),str(result['hmm_end']),'Reciprocal mismatch'])
		else:
			this_match['orthoid'] = orthoid
			this_match['hmm_score'] = score
			results.append(this_match)

	if reciprocal_verbose:
		print('Checked reciprocal hits for {}. Took {:.2f}s.'.format(score,time()-T_reciprocal_start))

	return {'Results':results,'Kicks':this_fails}

####
# Database Settings
####
db_col_name = 'name'
db_col_target = 'target'
db_col_score = 'score'
db_col_evalue = 'evalue'
db_col_start = 'start'
db_col_end = 'end'
db_col_env_start = 'env_start'
db_col_env_end = 'env_end'
db_col_ali_start = 'ali_start'
db_col_ali_end = 'ali_end'
db_col_hmm_start = 'hmm_start'
db_col_hmm_end = 'hmm_end'
db_col_header = 'header'
db_col_hmmsearch_id = 'hmmsearch_id'
db_col_id = 'id'
db_col_digest = 'digest'
db_col_orthoid = 'ortholog_gene_id'
db_col_taxid = 'taxid'
db_col_query = 'query'
db_col_setid = 'setid'
db_col_sequence = 'sequence'
db_col_aaseq = 'aa_seq'
db_col_seqpair = 'sequence_pair'
db_col_ntseq = 'nt_seq'

orthoset_set_details = "orthograph_set_details"
orthoset_taxa = "orthograph_taxa"
orthoset_seqpairs = "orthograph_sequence_pairs"
orthoset_orthologs = "orthograph_orthologs"
orthoset_aaseqs = "orthograph_aaseqs"
orthoset_ntseqs = "orthograph_ntseqs"

taxa_hmmsearch_table = 'orthograph_hmmsearch'
taxa_blast_table = 'orthograph_blast'
taxa_ests_table = 'orthograph_ests'
taxa_species_info = 'orthograph_species_info'
####
# Misc Settings
####

strict_search_mode = False
orthoid_list_file = None
frameshift_correction = True
extend_orf = True
orf_overlap_minimum = 0.15
clear_output = False

header_seperator = '|'

####

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input',type=str, default='Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa',
		help='Path to directory of Input folder')
	parser.add_argument('-oi','--orthoset_input',type=str, default='Syrphidae/orthosets',
		help='Path to directory of Orthosets folder')
	parser.add_argument('-o','--orthoset',type=str, default='Ortholog_set_Mecopterida_v4',
		help='Orthoset')
	parser.add_argument('-ml','--min_length',type=int, default=30,
		help='Minimum Transcript Length')
	parser.add_argument('-ms','--min_score',type=float, default=40,
		help='Minimum Hit Domain Score')
	parser.add_argument('-sdm','--score_diff_multi',type=float, default=1.05,
		help='Multi-gene Score Difference Adjustment')
	parser.add_argument('-mom','--min_overlap_multi',type=float, default=0.3,
		help='Multi-gene Minimum Overlap Adjustment')
	parser.add_argument('-sdi','--score_diff_internal',type=float, default=1.5,
		help='Internal Score Difference Adjustmen')
	parser.add_argument('-moi','--min_overlap_internal',type=float, default=0.9,
		help='Internal Minimum Overlap Adjustment')
	parser.add_argument('-p', '--processes', type=int, default=1,
        help='Number of threads used to call processes.')
	parser.add_argument('-v', '--verbose', type=int, default=2,
        help='Verbose debug.')
	parser.add_argument('-d', '--debug', type=int, default=0,
        help='Verbose debug.')

	args = parser.parse_args()	

	debug = args.debug != 0

	num_threads = args.processes

	####
	# Filter settings
	####

	min_length = args.min_length
	min_score = args.min_score
	score_diff_multi = args.score_diff_multi
	min_overlap_multi = args.min_overlap_multi
	score_diff_internal = args.score_diff_internal
	min_overlap_internal = args.min_overlap_internal
	verbose = range(0,args.verbose+1)

	####

	input_path = args.input
	species_name = os.path.basename(input_path).split('.')[0]

	print('Doing {}.'.format(species_name))

	if 1 in verbose:
		T_init_db = time()

	if os.path.exists("/run/shm"):
		tmp_path = "/run/shm"
	elif os.path.exists ("/dev/shm"):
		tmp_path = "/dev/shm"
	else:
		tmp_path = os.path.join(input_path,'tmp')
		if not os.path.exists(tmp_path):
			os.mkdir(tmp_path)
	
	if orthoid_list_file:
		list_of_wanted_orthoids = open(orthoid_list_file).read().split('\n')
		wanted_orthoid_only = True
	else:
		list_of_wanted_orthoids = []
		wanted_orthoid_only = False

	orthoset_path = args.orthoset_input
	orthoset = args.orthoset

	orthoset_db_path = os.path.join(orthoset_path,orthoset+'.sqlite')
	orthoset_db_con = sqlite3.connect(orthoset_db_path)

	cache_size = 16000000

	orthoset_db_con.execute('PRAGMA journal_mode = OFF;')
	orthoset_db_con.execute('PRAGMA synchronous = 0;')
	orthoset_db_con.execute(f'PRAGMA cache_size = {cache_size};')
	orthoset_db_con.execute('PRAGMA locking_mode = EXCLUSIVE;')
	orthoset_db_con.execute('PRAGMA temp_store = MEMORY;')

	orthoset_id = get_set_id(orthoset_db_con,orthoset)

	if clear_output:
		clear_output_path(input_path)

	taxa_db_path = os.path.join(input_path,os.path.basename(input_path).split('.')[0]+'.sqlite')
	taxa_db_con = sqlite3.connect(taxa_db_path)

	taxa_db_con.execute('PRAGMA journal_mode = OFF;')
	taxa_db_con.execute('PRAGMA synchronous = 0;')
	taxa_db_con.execute(f'PRAGMA cache_size = {cache_size};')
	taxa_db_con.execute('PRAGMA locking_mode = EXCLUSIVE;')
	taxa_db_con.execute('PRAGMA temp_store = MEMORY;')

	aa_out = 'aa'
	nt_out = 'nt'

	aa_out_path = os.path.join(input_path,aa_out)
	nt_out_path = os.path.join(input_path,nt_out)

	if os.path.exists(aa_out_path):
		shutil.rmtree(aa_out_path)
	os.mkdir(aa_out_path)

	if os.path.exists(nt_out_path):
		shutil.rmtree(nt_out_path)
	os.mkdir(nt_out_path)

	if 2 in verbose:
		T_species_id = time()
		print('Initialized databases. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing species id.'.format(time()-T_global_start,time()-T_init_db))

	species_id = get_species_id(taxa_db_con,taxa_hmmsearch_table)

	if 2 in verbose:
		T_reference_taxa = time()
		print('Got species id. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing reference taxa in set.'.format(time()-T_global_start,time()-T_species_id))

	reference_taxa = get_taxa_in_set(orthoset_id,orthoset_db_con)

	if 2 in verbose:
		T_blast_aaseqs = time()
		print('Got referenca taxa in set. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing blast hit aa seqs'.format(time()-T_global_start,time()-T_reference_taxa))

	aaseqs_in_orthoid = get_orthologs_for_set_hashref(orthoset_id,orthoset_db_con)

	if 2 in verbose:
		T_hmmresults = time()
		print('Got blast hit aa seqs. Elapsed time {:.2f}s. Took {:.2f}s. Grabbing hmmresults.'.format(time()-T_global_start,time()-T_blast_aaseqs))

	gene_based_results,header_based_results,ufr_rows = get_scores_list(min_score,taxa_db_path,orthoset_db_path,min_length,orthoset_id)

	ufr_path = os.path.join(input_path,'unfiltered-hits.csv')

	ufr_out = []
	for row in ufr_rows:
		ufr_out.append(','.join(row)+'\n')

	if debug:
		open(ufr_path,'w').writelines(ufr_out)

	####################################
	if 2 in verbose:
		print('Got hmmresults. Elapsed time {:.2f}s. Took {:.2f}s.'.format(time()-T_global_start,time()-T_hmmresults))

	if 1 in verbose:
		T_filter_start = time()
		if args.verbose != 1:
			print('Retrieved data from DB. Elapsed time {:.2f}s. Took {:.2f}s. Filtering hits.'.format(time()-T_global_start,time()-T_hmmresults))
		else:
			print('Retrieved data from DB. Elapsed time {:.2f}s. Took {:.2f}s. Filtering hits.'.format(time()-T_global_start,time()-T_init_db))


	filter_verbose = True

	if filter_verbose:
		filtered_sequences_log = []	
		filtered_sequences_log_path = os.path.join(input_path,'filtered-hits.csv')

	f_duplicates = {}

	if 2 in verbose:
		T_grab_multi = time()
		print('Grabbing multi-gene headers. Elapsed time {:.2f}s.'.format(time()-T_global_start))

	for orthoid in gene_based_results:
		for hit in gene_based_results[orthoid]:
			if 'revcomp' in hit['header']:
				base_header = get_baseheader(hit['header'])
				raw_length = int(base_header.split('_length_')[1]) / 3
				length = math.floor(raw_length)

				new_env_start = length - int(hit['env_end'])
				new_env_end = length - int(hit['env_start'])
				hit['env_start'] = new_env_start
				hit['env_end'] = new_env_end
				
			if hit['header'] not in f_duplicates:
				f_duplicates[hit['header']] = []

			f_duplicates[hit['header']].append(hit)

	headers = list(f_duplicates.keys())
	for header in headers:
		if len(f_duplicates[header]) > 1:
			unique_genes = list(dict.fromkeys([i['orthoid'] for i in f_duplicates[header]]))
			if len(unique_genes) <= 1: # But the same gene
				f_duplicates.pop(header)
		else: 
			f_duplicates.pop(header)

	for header_left in f_duplicates:
		header_based_results[header_left] = []

	if 2 in verbose:
		T_filter_multi = time()
		print('Got multi-gene headers. Elapsed time {:.2f}s. Took {:.2f}s. Filtering dupes.'.format(time()-T_global_start,time()-T_grab_multi))

	multi_verbose = 4 in verbose

	if num_threads == 1:
		for header in f_duplicates:
			this_hits = f_duplicates[header]
			data = multi_filter_dupes(header,this_hits,filter_verbose,min_overlap_multi,score_diff_multi,tmp_path,multi_verbose)

			filtered_sequences_log.extend(data['Log'])
			header_based_results[data['Remaining'][0]['header']] = data['Remaining']
			
	else:
		arguments = list()
		for header in f_duplicates:
			this_hits = f_duplicates[header]
			arguments.append((header,this_hits,filter_verbose,min_overlap_multi,score_diff_multi,tmp_path,multi_verbose))

		with Pool(num_threads) as pool:
			multi_data = pool.starmap(run_multi_filter, arguments, chunksize=1)

		for data in multi_data:
			filtered_sequences_log.extend(data['Log'])
			header_based_results[data['Remaining'][0]['header']] = data['Remaining']

	if 2 in verbose:
		T_kick_multi = time()
		print('Found multi-gene dupes. Elapsed time {:.2f}s. Took {:.2f}s. Kicking dupes.'.format(time()-T_global_start,time()-T_filter_multi))

	score_based_results = {}

	for header in header_based_results:
		for match in header_based_results[header]:
			if match['hmm_score'] not in score_based_results:
				score_based_results[match['hmm_score']] = []
			score_based_results[match['hmm_score']].append(match)
			
	reciprocal_match_start = time()

	scores = list(score_based_results.keys())
	scores.sort() #Ascending
	scores.reverse() #Descend	

	transcripts_mapped_to = {}
	duplicates = {}

	reciprocal_verbose = 4 in verbose

	if 2 in verbose:
		T_reciprocal_search = time()
		print('Kicked multi-gene dupes. Elapsed time {:.2f}s. Took {:.2f}s. Doing reciprocal hit check.'.format(time()-T_global_start,time()-T_kick_multi))

	if num_threads == 1:
		for score in scores:
			if reciprocal_verbose:
				T_reciprocal_start = time()
				print('Ensuring reciprocal hit for hmmresults in {}'.format(score))
			hmmresults = score_based_results[score]
		
			if len(hmmresults) >= 1:
				for result in hmmresults:
					score_start = time()
					orthoid = result['orthoid']
					result_hmmsearch_id = result['hmmsearch_id']

					if wanted_orthoid_only:
						if orthoid not in list_of_wanted_orthoids:
							continue

					blast_results = get_blastresults_for_hmmsearch_ids(taxa_db_con,result_hmmsearch_id)

					this_match = is_reciprocal_match(orthoid, blast_results,taxa_db_path,orthoset_db_path,reference_taxa,aaseqs_in_orthoid=aaseqs_in_orthoid)

					if this_match == None:
						filtered_sequences_log.append([this_match['orthoid'],this_match['hmmhit'],this_match['header'],str(this_match['hmm_score']),str(this_match['hmm_start']),str(this_match['hmm_end']),'Reciprocal mismatch'])
					else:
						this_match['orthoid'] = orthoid
						this_match['hmm_score'] = score
						this_match['uuid'] = this_match['hmmhit']+'{}{}'.format(this_match['hmm_start'],this_match['hmm_end'])

						if transcript_not_long_enough(this_match, min_length):
							continue
						
						if orthoid not in transcripts_mapped_to:
							transcripts_mapped_to[orthoid] = []

						transcripts_mapped_to[orthoid].append(this_match)
			if reciprocal_verbose:
				print('Checked reciprocal hits for {}. Took {:.2f}'.format(score,time()-T_reciprocal_start))
	else:
		arguments = list()
		for score in scores:
			hmmresults = score_based_results[score]
			arguments.append((hmmresults,list_of_wanted_orthoids,taxa_db_path,orthoset_db_path,reference_taxa, score, min_length, aaseqs_in_orthoid, tmp_path,reciprocal_verbose))
		with Pool(num_threads) as pool:
			reciprocal_data = pool.starmap(run_reciprocal_search, arguments, chunksize=1)

		for data in reciprocal_data:
			filtered_sequences_log.extend(data['Kicks'])
			for this_match in data['Results']:
				this_match['uuid'] = this_match['hmmhit']+'{}{}'.format(this_match['hmm_start'],this_match['hmm_end'])
				orthoid = this_match['orthoid']

				if transcript_not_long_enough(this_match, min_length):
					continue

				if orthoid not in transcripts_mapped_to:
					transcripts_mapped_to[orthoid] = []

				transcripts_mapped_to[orthoid].append(this_match)
		
	if 2 in verbose:
		T_internal_search = time()
		print('Reciprocal check done. Elapsed time {:.2f}s. Took {:.2f}s. Looking for internal dupes.'.format(time()-T_global_start,time()-T_reciprocal_search))

	total_hits = 0
	internal_verbose = 4 in verbose

	if num_threads == 1:
		internal_data = []
		for orthoid in transcripts_mapped_to:
			this_gene_transcripts = transcripts_mapped_to[orthoid]
			internal_data.append(internal_filter_gene(this_gene_transcripts, orthoid, min_overlap_internal, score_diff_internal, filter_verbose,input_path,tmp_path,internal_verbose))
			
	else:
		arguments = list()
		for orthoid in transcripts_mapped_to:
			this_gene_transcripts = transcripts_mapped_to[orthoid]
			arguments.append((this_gene_transcripts, orthoid, min_overlap_internal, score_diff_internal, filter_verbose,input_path,tmp_path,internal_verbose))
		with Pool(num_threads) as pool:
			internal_data = pool.starmap(run_internal_filter, arguments, chunksize=1)

	if 2 in verbose:
		T_internal_kick = time()
		print('Found internal dupes. Elapsed time {:.2f}s. Took {:.2f}s. Kicking internal dupes.'.format(time()-T_global_start,time()-T_internal_search))

	transcripts_mapped_to = {}
	brh_file_out = []
	brh_path = os.path.join(input_path,'best-reciprocal-hits.txt')

	for data in internal_data:
		gene = data['gene']
		if gene not in duplicates.keys():
			duplicates[gene] = set()

		transcripts_mapped_to[gene] = []

		for this_match in data['Passes']:
			if this_match['uuid'] in duplicates[gene]:
				filtered_sequences_log.append([this_match['orthoid'],this_match['hmmhit'],this_match['header'],str(this_match['hmm_score']),str(this_match['hmm_start']),str(this_match['hmm_end']),'Duplicate UUID'])
				continue
			else:
				duplicates[gene].add(this_match['uuid'])
			transcripts_mapped_to[gene].append(this_match)

			brh_file_out.append('\t'.join([this_match['orthoid'],
												this_match['header'], 
												str(this_match['ali_start']),
												str(this_match['ali_end']),
												str(this_match['hmm_score']),
												str(this_match['hmm_evalue']),
												str(this_match['hmm_start']), 
												str(this_match['hmm_end'])])+'\n')

		total_hits += len(data['Passes'])
		filtered_sequences_log.extend(data['Log'])
		
	filtered_sequences_log_out = []
	for line in filtered_sequences_log:
		filtered_sequences_log_out.append(','.join(line)+'\n')

	if debug:
		open(filtered_sequences_log_path,'w').writelines(filtered_sequences_log_out)

	brh_file_out.sort()
	open(brh_path,'w').writelines(brh_file_out)

	if 2 in verbose:
		print('Kicked internal dupes. Elapsed time {:.2f}s. Took {:.2f}s.'.format(time()-T_global_start,time()-T_internal_kick))
		
	if 1 in verbose:
		T_exonerate_genes = time()
		print('Filtering done, found {} hits. Elapsed time {:.2f}s. Took {:.2f}s. Exonerating genes'.format(total_hits,time()-T_global_start,time()-T_filter_start))

	#Clear summary files
	for item in os.listdir(tmp_path):
		item_path = os.path.join(tmp_path,item)
		if '-summary.txt' in item:
			os.remove(item_path)

	exonerate_verbose = 3 in verbose
	summary_file_out = []

	#Disperse into reftaxons
	if num_threads == 1:
		for orthoid in transcripts_mapped_to:
			list_of_hits = transcripts_mapped_to[orthoid]
			exonerate_gene_multi(orthoid,list_of_hits,taxa_db_path,orthoset_db_path,input_path,min_score,orthoset_id,aa_out_path,min_length,species_name,nt_out_path,tmp_path,exonerate_verbose)
			
	else:
		if num_threads > 64:
			this_num_threads = 64
		else:
			this_num_threads = num_threads

		arguments = list()
		for orthoid in transcripts_mapped_to:
			list_of_hits = transcripts_mapped_to[orthoid]
			arguments.append((orthoid,list_of_hits,taxa_db_path,orthoset_db_path,input_path,min_score,orthoset_id,aa_out_path,min_length,species_name,nt_out_path,tmp_path,exonerate_verbose))
		with Pool(this_num_threads) as pool:
			pool.map(run_exonerate, arguments, chunksize=1)
	
	summary_file_out = []
	for item in os.listdir(tmp_path):
		item_path = os.path.join(tmp_path,item)
		if '-summary.txt' in item:
			item_content = open(item_path).read()
			summary_file_out.append(item_content)
			os.remove(item_path)

	summary_file_path = os.path.join(input_path,'summary.txt')
	open(summary_file_path,'w').write('\n'.join(summary_file_out))

	if 1 in verbose:
		print('Done. Final time {:.2f}s. Exonerate took {:.2f}s.'.format(time()-T_global_start,time()-T_exonerate_genes))

	if args.verbose == 0:
		print('Done took {:.2f}s.'.format(time()-T_global_start))