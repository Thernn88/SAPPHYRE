# PhyMMR

Requirements
	
	Externals
	Mafft v7.510+ https://mafft.cbrc.jp/alignment/software/ - Currently requires make.
	Exonerate - 2.4.0 https://github.com/nathanweeks/exonerate - sudo apt install exonerate
	HMMer - 3.3+ https://github.com/EddyRivasLab/hmmer - sudo apt install Hmmer
	Cargo - sudo apt install cargo
	SQLite3 - sudo apt install sqlite3 (Will be removed later.)
	Perl - sudo apt install perl (Will be removed later.)
	
	Python - 3.7+ (Recommended 3.9+)
		Python Modules
			wrap-rocks
			tqdm
			numpy - 1.22+ Required
			blosum_distance
			itertools
			
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
