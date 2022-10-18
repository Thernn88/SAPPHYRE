# PhyMMR

Requirements
	
	Externals
	Mafft 7.489+ https://mafft.cbrc.jp/alignment/software/ sudo apt install mafft
	Exonerate 2.2.0+ sudo apt install Exonerate
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
			
Data Generation Scripts in order of use

	PrepareDB.py

	HmmsearchDB.py

	BlastPalDB.py

	ReporterDB.py


Post-Processing Scripts in order of use

	mafft.py

	nt_batch.py

	FlexCull.py (Not Recommended / Optional for Transcriptome input)

	OutlierCheck.py
	
	MergeOverlap.py

TO FINISH
