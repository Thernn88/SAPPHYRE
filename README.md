# PhyMMR

Requirements
	
	Externals
	Mafft v7.510+ https://mafft.cbrc.jp/alignment/software/ - Currently requires make.
	Exonerate - 2.4.0 https://github.com/nathanweeks/exonerate - a custom build to utilize modified fastatranslate. Fastatranslate is planed to be replaced with a		rust module w/ thread capabilities.
	HMMer - 3.3+ https://github.com/EddyRivasLab/hmmer - sudo apt install Hmmer
	SQLite3 - sudo apt install sqlite3 (Will be removed later.)
	Perl - sudo apt install perl (Will be removed later.)
	
	Python - 3.8+ (Recommended 3.9+)
		Python Modules
			wrap-rocks
			tqdm
			numpy - 1.22+ Required
			blosum_distance
			itertools
			
	Installing custom exonerate. 
		PreReqs
		sudo apt remove exonerate (if package previously installed)
		sudo apt install autoconf
		sudo apt install libglib2.0-dev
		
		---- All following commands should be run from inside exonerate folder.
		
		cd exonerate
		autoreconf -f -i
		./configure
		
		----Above should be one time only. Below to reinstall. 
		
		cd exonerate(if not already inside the exonerate folder
		make clean
		make
		sudo make install
		restart VM (?)
			
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
