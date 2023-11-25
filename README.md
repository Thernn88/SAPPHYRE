# sapphyre

Sapphyre is an assembly-less solution for processing high-throughput sequencing reads for phylogenetics. Lorem Ipsum,,,

## Requirements
	
	Externals
	Mafft 7.489+ https://mafft.cbrc.jp/alignment/software/ sudo apt install mafft
	Diamond 2.1.6+ - https://github.com/bbuchfink/diamond 
	sra-toolkit - sudo apt install sra-toolkit
	clustalo - sudo apt install clustalo
	
	Python - 3.9+ (Recommended 3.11)
		Python Modules
			beautifulsoup4>=4.12.2
			Bio>=1.6.0
			biopython>=1.81
			msgspec>=0.18.4
			needletail>=0.5.0
			numpy>=1.24.2
			openpyxl>=3.1.2
			pandas>=2.1.3
			parasail>=1.3.4
			phymmr_tools>=0.6.5
			pro2codon>=1.2.4
			pyfastx>=2.0.1
			Requests>=2.31.0
			scipy>=1.11.4
			setuptools>=67.6.1
			tqdm>=4.66.1
			wrap_rocks>=0.3.7
			xxhash>=3.4.1
			
## Usage

Generic command line

	sapphyre <global args> <module> <module args> <FILE>
	sapphyre -p 32 -v align -alm frags test_dataset/*.fa

Example commands workflow using Mayer et al. (2021) reference set. Verbose enabled using -v. Multithreading enabled w/ -p <X> threads passed.  

	python3 -m sapphyre -p 6 Prepare input/lepidoptera
	python3 -m sapphyre -p 16 -v Diamond -o Lepidoptera_orthoDB9_extended sapphyre/lepidoptera/*.fa
	python3 -m sapphyre -p 16 -v Reporter -o Lepidoptera_orthoDB9_extended sapphyre/lepidoptera/*.fa
	python3 -m sapphyre -p 8 -v mafft -o Lepidoptera_orthoDB9_extended sapphyre/lepidoptera/*.fa
	python3 -m sapphyre -p 16 -v Pal2Nal sapphyre/lepidoptera/*.fa
	python3 -m sapphyre -p 16 -v FlexCull sapphyre/lepidoptera/*.fa
	python3 -m sapphyre -p 16 -v OutlierCheck sapphyre/lepidoptera/*.fa
	python3 -m sapphyre -p 16 -v MergeOverlap sapphyre/lepidoptera/*.fa

Generic Commands
	
	python3 -m sapphyre <args> Prepare <args> input/<DIR>/
	python3 -m sapphyre <args> Diamond <args> sapphyre/<DIR>/*.fa
	python3 -m sapphyre <args> Reporter <args> sapphyre/<DIR>/*.fa
	python3 -m sapphyre <args> mafft <args> sapphyre/<DIR>/*.fa
	python3 -m sapphyre <args> Pal2Nal <args> sapphyre/<DIR>/*.fa
	python3 -m sapphyre <args> FlexCull <args> sapphyre/<DIR>/*.fa
	python3 -m sapphyre <args> OutlierCheck <args> sapphyre/<DIR>/*.fa
	python3 -m sapphyre <args> MergeOverlap <args> sapphyre/<DIR>/*.fa

## Help

   **Main**
		
	python3 -m sapphyre -h
	usage: sapphyre [-h] [-v] [-p PROCESSES]
	
	Order: Prepare, Diamond, Reporter, mafft, pal2nal, FlexCull (optional), OutlierCheck, MergeOverlap, MergeGenes

	positional arguments:
  	{Prepare,Diamond,Reporter,mafft,Pal2Nal,FlexCull,OutlierCheck,MergeOverlap,MergeGenes,SRADownload}
	
    Prepare             Loads NT input files (.fa, .fas or .fasta) into a rocksdb database. Unique NT sequences stored
                        with a duplicate count stored for later use. 
			
    Diamond             Translates and queries input NT sequences against protein reference sequences using the Diamond external. Filters 
                        output prevent sequence reuse. Loads hits into RocksDB after filters finish.
			
    Reporter            Queries a sequence using exonerate to align it against a target reference and trim it to mapped region. Produces aa and nt output.
    
    mafft               Aligns AA sequences against existing reference alignment.
    
    Pal2Nal             Mirrors Amino Acid Alignment to Nucleotide data using specified NCBI table. Also performs
                        basic error checking on data.
			
    FlexCull            Adaptive End Trimming algorithm that trims the start and end of candidate reads to remove
                        introns and/or other LQ bases.
			
    OutlierCheck        Calculates a Blosum62 distance matrix used to remove outlier sequences above a
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

**Prepare**
	
	python3 -m sapphyre Prepare -h
	usage: sapphyre Prepare [-h] [-c] [-ml MINIMUM_SEQUENCE_LENGTH] [-sl SEQUENCES_PER_LEVEL] [-k] INPUT

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

**Hmmsearch**
	
	python3 -m sapphyre Hmmsearch -h
	usage: sapphyre Hmmsearch [-h] [-oi ORTHOSET_INPUT] [-o ORTHOSET] [-ovw] [-s SCORE] [-e EVALUE]
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
							
**Reporter**

	python3 -m sapphyre Reporter -h
	usage: sapphyre Reporter [-h] [-oi ORTHOSET_INPUT] [-o ORTHOSET] [-ml MIN_LENGTH] [-ms MIN_SCORE] [-d DEBUG]
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
							
**mafft**
	
	python3 -m sapphyre mafft -h
	usage: sapphyre mafft [-h] [-oi ORTHOSET_INPUT] [-o ORTHOSET] [-l] INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help            show this help message and exit
	-oi ORTHOSET_INPUT, --orthoset_input ORTHOSET_INPUT
							Path to directory of Orthosets folder
	-o ORTHOSET, --orthoset ORTHOSET
							Orthoset
	-l, --linsi           Enable the use of mafft-linsi.
	
**Pal2Nal**
	
	python3 -m sapphyre Pal2Nal -h
	usage: sapphyre Pal2Nal [-h] [-t TABLE] INPUT [INPUT ...]

	positional arguments:
	INPUT                 Path to directory of Input folder

	options:
	-h, --help            show this help message and exit
	-t TABLE, --table TABLE
							Table ID.
							
**FlexCull**
	
	python3 -m sapphyre FlexCull -h
		usage: sapphyre FlexCull [-h] [-o OUTPUT] [-aa AMINO_ACID] [-nt NUCLEOTIDE] [-m MATCHES] [-bp BASE_PAIR] [-d]
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

**OutlierCheck**
		
	python3 -m sapphyre OutlierCheck -h
	usage: sapphyre OutlierCheck [-h] [-o OUTPUT] [-t THRESHOLD] [--no-references] [-s {cluster,original}] [-d]
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
	
**MergeOverlap**
	
	python3 -m sapphyre MergeOverlap -h
	usage: sapphyre MergeOverlap [-h] [-aa AA_INPUT] [-nt NT_INPUT] [-d] [-c COMPARISON] [-io] [-m MAJORITY]
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

Some dependencies require a c++ toolset installed, plus various libraries. At the
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
  
