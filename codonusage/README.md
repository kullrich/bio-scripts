Python codonusage
=================

Dependencies
------------
This python script depends on biopython <http://biopython.org/wiki/Download>, namely the module:

* Bio.SeqIO

Usage
-----

Specify your input FASTA file containing CDS sequences and run the script:

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX

This by default extracts the codonusage for each sequence in the FASTA input file and write a table with the number of codons used to OUTPUTPREFIX".codonusage"

If your FASTA file also conatins sequences that length are not divdable by 3 (len % 3 != 0) it will not report these files to the table but will report the sequence ids to std.out.

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX -r

With the 'enc' option one can choose a method to calculate the effective number of codons for each sequence based either on equation (4) of 'Wright (1990). The effective number of codons used in a gene. Gene 87:23-29.' ['eq4Wright'] or equation (3) or equation (5) of 'Sun et al. (2012). An Improved Implementation of Effective Number of Codons (Nc). Mol. Biol. Evol. 30:191-196.' ['eq3Sun','eq4Sun'].

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX -enc eq4Wright

Note: If you specify the option (-six2fourtwo) all sixfold codon groups will be split into one fourfold and one twofold group to account for possible altering tRNA pools. This option will only affect the calculation of the ENC values for 'eq2Sun' and 'eq5Sun'.

Output
------

As Output different tables will be generated.

	Raw Codon counts: OUTPUTPREFIX.codoncnt
	ACTG counts: OUTPUTPREFIX.actgcnt
	First codon position ACTG counts: OUTPUTPREFIX.firstcnt
	Second codon position ACTG counts: OUTPUTPREFIX.secondcnt
	Third codon position ACTG counts: OUTPUTPREFIX.thirdcnt
        Relative Synonymous Codon Usage: OUTPUTPREFIX.rscucnt

If the (-enc) option was used an additional table containing the choosen method or all ENC values (-enc all) will be written to a file.

	Effective Number of Codons: OUTPUTPREFIX.enc

	Wright (1990). The effective number of codons used in a gene. Gene 87:23-29.
	Equation (4) 
	#\widehat{N}_{c} = 2 + GC_{(3)} + (\frac{29}{GC^{2}_{(3)} + (1 - GC^{2}_{(3)})})

	Sun et al. (2012). An Improved Implementation of Effective Number of Codons (Nc). Mol. Biol. Evol. 30:191-196.

	Equation (2)
	#F_{CF} = \sum_{i=1}^{m} (\frac {n_{i}}{n})^{2}

	Equation (3)
	#F_{CF} = \sum_{i=1}^{m} (\frac {n_{i}+1}{n+m})^{2}

	Equation (5)
	

	Relative Synonymous Codon Usage
	#RSCU_{i,j} = \frac{NumberofCodons_{i} \times CodonFrequency_{j}}{\sum_{j=1}^{NumberofCodons_{j}} CodonFrequency_{i}}
