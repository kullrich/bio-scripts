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
