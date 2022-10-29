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
	Perl - sudo apt install perl (Will be removed later.)
	
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
`python phymmr/<scriptname>.py`.

Data Generation Scripts in order of use

1. PrepareDB.py
2. HmmsearchDB.py
3. BlastPalDB.py
4. ReporterDB.py

Post-Processing Scripts in order of use

1. mafft.py
2. nt_batch.py
3. FlexCull.py (Not Recommended / Optional for Transcriptome input)
4. OutlierCheck.py
5. MergeOverlap.py

TO FINISH
