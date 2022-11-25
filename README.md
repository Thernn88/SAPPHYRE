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

##HELP

#Main
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
                        into RocksDB
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

phymmr  Copyright (C) 2022  PhyMMR Team
License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it under certain conditions.

#Prepare
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
  
