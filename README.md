# PhyMMR

## Requirements
	
	Externals
	Mafft 7.489+ https://mafft.cbrc.jp/alignment/software/ sudo apt install mafft
	Exonerate 2.4.0+ sudo apt install Exonerate
	HMMer - 3.3+ https://github.com/EddyRivasLab/hmmer - sudo apt install Hmmer
	Blast - NCBI-Blast 2.2.28+ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ or sudo apt install ncbi-blast+
	SQLite3 - sudo apt install sqlite3 (Will be removed later.)
	
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
	
	--------- Breaks below this line as of this moment for this reference set.
	
	python3 -m phymmr -p 16 -v Pal2Nal PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 16 -v FlexCull PhyMMR/lepidoptera/*.fa
	python3 -m phymmr -p 16 -v OutlierCheck PhyMMR/lepidoptera/*.fa
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
  
