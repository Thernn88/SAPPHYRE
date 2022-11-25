# PhyMMR

## Requirements
	
	Externals
	Mafft 7.489+ https://mafft.cbrc.jp/alignment/software/ sudo apt install mafft
	Exonerate 2.4.0+ sudo apt install Exonerate
	HMMer - 3.3+ https://github.com/EddyRivasLab/hmmer - sudo apt install Hmmer
	Blast - NCBI-Blast 2.2.28+ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ or sudo apt install ncbi-blast+
	SQLite3 - sudo apt install sqlite3 (Will be removed later.)
	sra-toolkit - sudo apt install sra-toolkit
	
	Python - 3.8+ (Recommended 3.9+)
		Python Modules
			wrap-rocks
			tqdm
			numpy - 1.22+ Required
			blosum_distance
			itertools
			bio
			xxhash
			
## Usage

All scripts are located in the phymmr directory. You can call them using

Example commands workflow using Mayer et al. (2021) reference set. Verbose enabled using -v. Multithreading enabled w/ -p <X> threads passed.  

	python3 -m phymmr -p 6 Prepare input/lepidoptera
	python3 -m phymmr -p 16 -v Hmmsearch -o Lepidoptera_orthoDB9_extended PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 16 -v BlastPal -o Lepidoptera_orthoDB9_extended PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 16 -v Reporter -o Lepidoptera_orthoDB9_extended PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 8 -v mafft -o Lepidoptera_orthoDB9_extended PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 16 -v Pal2Nal PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 16 -v FlexCull PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 16 -v OutlierCheck PhyMMR/lepidoptera/*.fa
	
	--------- Breaks below this line as of this moment for this reference set.
	python3 -m phymmr -p 16 -v MergeOverlap PhyMMR/lepidoptera/*.fa

Generic Commands
	
	python3 -m phymmr <args> Prepare <args> input/<DIR>/
	python3 -m phymmr <args> Hmmsearch <args> PhyMMR/<DIR>/*.fa
	python3 -m phymmr <args> BlastPal <args> PhyMMR/<DIR>/*.fa
	python3 -m phymmr <args> Reporter <args> PhyMMR/<DIR>/*.fa
	python3 -m phymmr <args> mafft <args> PhyMMR/<DIR>/*.fa
	python3 -m phymmr <args> Pal2Nal <args> PhyMMR/<DIR>/*.fa
	python3 -m phymmr <args> FlexCull <args> PhyMMR/<DIR>/*.fa
	python3 -m phymmr <args> OutlierCheck <args> PhyMMR/<DIR>/*.fa
	python3 -m phymmr <args> MergeOverlap <args> PhyMMR/<DIR>/*.fa

## Help

Main
		python3 -m phymmr -h
		usage: phymmr [-h] [-v] [-p PROCESSES]
            		  {Prepare,Hmmsearch,BlastPal,Reporter,mafft,Pal2Nal,FlexCull,OutlierCheck,MergeOverlap,MergeGenes,SRADownload}
              ...

	Order: Prepare, Hmmsearch, BlastPal, Reporter, mafft, pal2nal, FlexCull (optional), OutlierCheck, MergeOverlap, MergeGenes

	positional arguments:
  	{Prepare,Hmmsearch,BlastPal,Reporter,mafft,Pal2Nal,FlexCull,OutlierCheck,MergeOverlap,MergeGenes,SRADownload}
	
    Prepare             Loads NT input files (.fa, .fas or .fasta) into a rocksdb database. Unique NT sequences stored
                        with a duplicate count stored for later use. Performs six-fold translation of base NT
                        sequences into AA translation. Unique AA translation stored with duplicate counts stored for
                        later use.
			
    Hmmsearch           Queries protein translations against profile HMMs using the HMMER3 external. Filters HMMER3
                        output using 3 custom filters: MultiGene, InternalMulti & InternalSubpar to remove LQ hits and
                        prevent sequence reuse. Loads hits into RocksDB after filters finish.
			
    BlastPal            Blasts Hmmsearch hits using NCBI-Blast against reference sequences. Loads resulting output
                        into RocksDB.
			
    Reporter            Checks Blast results to ensure a hit is reciprocal. Queries a sequence using exonerate to
                        align it against a target reference and trim it to mapped region. Produces aa and nt output.
    mafft               Aligns AA sequences against existing reference alignment.
    
    Pal2Nal             Mirrors Amino Acid Alignment to Nucleotide data using specified NCBI table. Also performs
                        basic error checking on data.
			
    FlexCull            Adaptive End Trimming algorithm that trims the start and end of candidate reads to remove
                        introns and/or other LQ bases.
			
    OutlierCheck        Calculates a Blosum62 distance matrix which are used to remove outlier sequences above a
                        threshold.
			
    MergeOverlap        Reference-guided De-novo Assembly Algorithm which merges overlapping reads into contiguous
                        segments (Contigs).
			
    MergeGenes          Initial Dataset Construction: Merges AA and NT data files across multiple taxa into a single
                        AA and NT pair per gene.
			
    SRADownload         Download fastq files from www.ncbi.nlm.nih.gov

	options:
 	 -h, --help            show this help message and exit
 	 -v, --verbose         Verbosity level. Repeat for increased verbosity.
 	 -p PROCESSES, --processes PROCESSES
         	               Number of threads used to call processes.

Prepare
	
	python3 -m phymmr Prepare -h
	usage: phymmr Prepare [-h] [-c] [-ml MINIMUM_SEQUENCE_LENGTH] [-sl SEQUENCES_PER_LEVEL] [-k] INPUT

	positional arguments:
	  INPUT                 Path to directory of Input folder

	options:
 	-h, --help            show this help message and exit
  	-c, --clear_database  Overwrite existing rocksdb database.
  	-ml MINIMUM_SEQUENCE_LENGTH, --minimum_sequence_length MINIMUM_SEQUENCE_LENGTH
                        Minimum input sequence length.
  	-sl SEQUENCES_PER_LEVEL, --sequences_per_level SEQUENCES_PER_LEVEL
                        Amount of sequences to store per database entry.
  	-k, --keep_prepared   Writes the prepared input fasta into the output taxa directory.

Hmmsearch
	python3 -m phymmr Hmmsearch -h
	usage: phymmr Hmmsearch [-h] [-oi ORTHOSET_INPUT] [-o ORTHOSET] [-ovw] [-s SCORE] [-e EVALUE]
                        [--excluded-list EXCLUDED_LIST] [--wanted-list WANTED_LIST] [--remake-protfile]
                        [-sdm SCORE_DIFF_MULTI] [-mom MIN_OVERLAP_MULTI] [-momi MINIMUM_OVERLAP_INTERNAL_MULTI]
                        [-sdi SCORE_DIFF_INTERNAL] [-moi MIN_OVERLAP_INTERNAL] [-m MAX_HMM_BATCH_SIZE] [-d DEBUG]
                        [--enable-multi-internal]
                        INPUT [INPUT ...]

	positional arguments:
		INPUT                 Path to input directory.

	options:
	-h, --help            show this help message and exit
	-oi ORTHOSET_INPUT, --orthoset_input ORTHOSET_INPUT
							Path to directory of Orthosets folder
	-o ORTHOSET, --orthoset ORTHOSET
							Orthoset
	-ovw, --overwrite     Remake domtbl files even if previous file exists.
	-s SCORE, --score SCORE
							Score threshold. Defaults to 40
	-e EVALUE, --evalue EVALUE
							Evalue threshold. Defaults to 0
	--excluded-list EXCLUDED_LIST
							File containing names of genes to be excluded
	--wanted-list WANTED_LIST
							File containing list of wanted genes
	--remake-protfile    	Force creation of a new protfile even if one already exists.
	-sdm SCORE_DIFF_MULTI, --score_diff_multi SCORE_DIFF_MULTI
							Multi-gene Score Difference Adjustment
	-mom MIN_OVERLAP_MULTI, --min_overlap_multi MIN_OVERLAP_MULTI
							Multi-gene Minimum Overlap Adjustment
	-momi MINIMUM_OVERLAP_INTERNAL_MULTI, --minimum_overlap_internal_multi MINIMUM_OVERLAP_INTERNAL_MULTI
							Internal Multi Minimum Overlap Adjustment
	-sdi SCORE_DIFF_INTERNAL, --score_diff_internal SCORE_DIFF_INTERNAL
							Internal Score Difference Adjustmen
	-moi MIN_OVERLAP_INTERNAL, --min_overlap_internal MIN_OVERLAP_INTERNAL
							Internal Minimum Overlap Adjustment
	-m MAX_HMM_BATCH_SIZE, --max_hmm_batch_size MAX_HMM_BATCH_SIZE
							Max hits per hmmsearch batch in db. Default: 250 thousand.
	-d DEBUG, --debug DEBUG
							Output debug logs.
	--enable-multi-internal
							Enable Hmmsearch internal multi filter
							
BlastPal
	python3 -m phymmr BlastPal -h
	usage: phymmr BlastPal [-h] [-oi ORTHOSET_INPUT] [-o ORTHOSET] [-bs BLAST_MINIMUM_SCORE] [-be BLAST_MINIMUM_EVALUE]
						[-m MAX_BLAST_BATCH_SIZE] [-ovw]
						INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help            	show this help message and exit
	-oi ORTHOSET_INPUT, --orthoset_input ORTHOSET_INPUT
							Path to directory of Orthosets folder
	-o ORTHOSET, --orthoset ORTHOSET
							Orthoset
	-bs BLAST_MINIMUM_SCORE, --blast_minimum_score BLAST_MINIMUM_SCORE
							Minimum score filter in blast.
	-be BLAST_MINIMUM_EVALUE, --blast_minimum_evalue BLAST_MINIMUM_EVALUE
							Minimum evalue filter in blast.
	-m MAX_BLAST_BATCH_SIZE, --max_blast_batch_size MAX_BLAST_BATCH_SIZE
							Max results per blastpal batch in db. Default: 250 thousand.
	-ovw, --overwrite     Overwrite existing blast results.

Reporter

	python3 -m phymmr Reporter -h
	usage: phymmr Reporter [-h] [-oi ORTHOSET_INPUT] [-o ORTHOSET] [-ml MIN_LENGTH] [-ms MIN_SCORE] [-d DEBUG]
							INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help            show this help message and exit
	-oi ORTHOSET_INPUT, --orthoset_input ORTHOSET_INPUT
							Path to directory of Orthosets folder
	-o ORTHOSET, --orthoset ORTHOSET
							Orthoset
	-ml MIN_LENGTH, --min_length MIN_LENGTH
							Minimum Transcript Length
	-ms MIN_SCORE, --min_score MIN_SCORE
							Minimum Hit Domain Score
	-d DEBUG, --debug DEBUG
							Verbose debug.
							
mafft
	python3 -m phymmr mafft -h
	usage: phymmr mafft [-h] [-oi ORTHOSET_INPUT] [-o ORTHOSET] [-l] INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help            show this help message and exit
	-oi ORTHOSET_INPUT, --orthoset_input ORTHOSET_INPUT
							Path to directory of Orthosets folder
	-o ORTHOSET, --orthoset ORTHOSET
							Orthoset
	-l, --linsi           Enable the use of mafft-linsi.
	
Pal2Nal	
	python3 -m phymmr Pal2Nal -h
	usage: phymmr Pal2Nal [-h] [-t TABLE] INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help            show this help message and exit
	-t TABLE, --table TABLE
							Table ID.
							
FlexCull							
	python3 -m phymmr FlexCull -h
		usage: phymmr FlexCull [-h] [-o OUTPUT] [-aa AMINO_ACID] [-nt NUCLEOTIDE] [-m MATCHES] [-bp BASE_PAIR] [-d]
						INPUT [INPUT ...]
	
	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help            show this help message and exit
	-o OUTPUT, --output OUTPUT
							Output Directory.
	-aa AMINO_ACID, --amino-acid AMINO_ACID
							AA Folder Name.
	-nt NUCLEOTIDE, --nucleotide NUCLEOTIDE
							NT Folder Name.
	-m MATCHES, --matches MATCHES
							Amount of nucleotides that have to match reference.
	-bp BASE_PAIR, --base-pair BASE_PAIR
							Minimum bp after cull.
	-d, --debug           Enable debug. When enabled output log of culls.		

OutlierCheck
		
	python3 -m phymmr OutlierCheck -h
	usage: phymmr OutlierCheck [-h] [-o OUTPUT] [-t THRESHOLD] [--no-references] [-s {cluster,original}] [-d]
                           INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to taxa

	options:
	-h, --help            show this help message and exit
	-o OUTPUT, --output OUTPUT
							Output folder
	-t THRESHOLD, --threshold THRESHOLD
							Greater than reference mean to be counted as an outlier. Default is 1.5x.
	--no-references       Disable output of reference sequences
	-s {cluster,original}, --sort {cluster,original}
							Sort candidate output by cluster and taxa, or preserver original order.
	-d, --debug           Log outliers to csv files
	
MergeOverlap
	python3 -m phymmr MergeOverlap -h
	usage: phymmr MergeOverlap [-h] [-aa AA_INPUT] [-nt NT_INPUT] [-d] [-c COMPARISON] [-io] [-m MAJORITY]
							[-mc MAJORITY_COUNT]
							INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help          	show this help message and exit
	-aa AA_INPUT, --aa_input AA_INPUT
							Path to directory of AA folder
	-nt NT_INPUT, --nt_input NT_INPUT
							Path to directory of NT folder
	-d, --debug         	Enable debug. When enabled displays each component of merged headers.
	-c COMPARISON, --comparison COMPARISON
							Fallback Comparison Taxa. Sequence in which Sequence A and Sequence B is compared to in split
							calculation.
	-io, --ignore_overlap_chunks
							Ignore overlapping chunks and merge all candidates for a reference taxon.
	-m MAJORITY, --majority MAJORITY
							Percentage for majority ruling.
	-mc MAJORITY_COUNT, --majority_count MAJORITY_COUNT
							Percentage for majority ruling.
---------------------


## Development setup

Some dependencies requires a c++ toolset installed, plus various libraries. At the
moment of writing, I do not know which additional libraries are needed (please complete
me later).

	~$ apt install build-essential
	~$ apt install libboost-all-dev  # may not be needed

For regular python:

	~$ virtualenv venv
	~$ source venv/bin/activate
	~$ pip install -r requirements.txt

For pypy:

	~$ virtualenv -p $(command -v pypy) venv
	~$ source venv/bin/activate
	~$ pip install -r requirements.txt
  
