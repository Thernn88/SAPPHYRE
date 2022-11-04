# PhyMMR

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
python3 -m phymmr -p 6 Prepare input/<DIR>/

	python3 -m phymmr -p 16 Hmmsearch PhyMMR/<DIR>/*.fa
	python3 -m phymmr -p 8 BlastPal PhyMMR/<DIR>/*.fa
	python3 -m phymmr -p 16 Reporter PhyMMR/<DIR>/*.fa
	python3 -m phymmr -p 8 mafft PhyMMR/<DIR>/*.fa
	python3 -m phymmr -p 16 Pal2Nal PhyMMR/<DIR>/*.fa
	python3 -m phymmr -p 16 FlexCull PhyMMR/<DIR>/*.fa
	python3 -m phymmr -p 16 OutlierCheck PhyMMR/<DIR>/*.fa
	python3 -m phymmr -p 64 MergeOverlap PhyMMR/SRR/*.fa

