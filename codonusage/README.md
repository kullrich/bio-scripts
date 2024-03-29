Python codonusage
=================

Description
-----------
This python script extracts codonusage from CDS input FASTA file. Output will be raw codon counts (.codoncnt), global ACTG counts (.actgcnt), first (.firstcnt), second (.secondcnt), third (.third) codon position counts and Relative Synonymous Codon Usage (.rscucnt). Optional different methods can be applied to calculate Effective Number of Codons (.enc).

ENC will be calculated based on:
* Wright (1990). The effective number of codons used in a gene. Gene 87:23-29. [Wright 1990](https://doi.org/10.1016/0378-1119(90)90491-9)
* Sun et al. (2012). An Improved Implementation of Effective Number of Codons (Nc). Mol. Biol. Evol. 30:191-196. [Sun 2012](https://doi.org/10.1093/molbev/mss201)


RSCU will be calculated based on:
* Sharp et al. (1986). Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucl. Acids. Res. 14:5125-5143. [Sharp 1986](https://doi.org/10.1093/nar/14.13.5125)

Installation and usage with conda
-----------------------

```
conda create -n codonusage python=3
conda activate codonusage
conda install -c bioconda biopython pandas seaborn ipython
wget https://raw.githubusercontent.com/kullrich/bio-scripts/master/codonusage/codonusage.py
python codonusage.py -h
```

Dependencies
------------
This python script depends on biopython <http://biopython.org/wiki/Download>, namely the module:

* Bio.SeqIO

Usage
-----

	optional arguments:
	-h,	--help	show this help message and exit
	-v,	--verbose	increase output verbosity
	-i	I	specify CDS input file in FASTA format
	-o	O	specify output prefix
	-r		specify if CDS sequences with length modulo 3 unequal to 0 should be removed and reported to std.out
	-enc	{eq4Wright,eq2Sun,eq5Sun,all}	specify equation to calculate ENC. Either equation (4) [eq4Wright] of (Wright. (1990) Gene 87:23-29) or equation (2) [eq2Sun] or equation (5) [eq5Sun] of (Sun et al. (2012) Mol. Biol. Evol. 30:191-196) or [all].
	-six2fourtwo	specify if sixfold codons should be grouped into one fourfold and one twofold group [default: False]. This will only affect calculation of ENC values.

Examples
--------

Specify your input FASTA file containing CDS sequences and run the script:

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX

This by default extracts the codonusage for each sequence in the FASTA input file and write a table with the number of codons used to OUTPUTPREFIX".codonusage"

If your FASTA file also conatins sequences that length are not divdable by 3 (len % 3 != 0) it will not report these files to the table but will report the sequence ids to std.out.

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX -r

With the 'enc' option one can choose a method to calculate the effective number of codons for each sequence based either on equation (4) of '[Wright (1990). The effective number of codons used in a gene. Gene 87:23-29.](https://doi.org/10.1016/0378-1119(90)90491-9)' ['eq4Wright'] or equation (3) or equation (5) of '[Sun et al. (2012). An Improved Implementation of Effective Number of Codons (Nc). Mol. Biol. Evol. 30:191-196.](https://doi.org/10.1093/molbev/mss201)' ['eq3Sun','eq5Sun'].

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX -enc eq4Wright

Note: If you specify the option (-six2fourtwo) all sixfold codon groups will be split into one fourfold and one twofold group to account for possible altering tRNA pools. This option will only affect the calculation of the ENC values for 'eq2Sun' and 'eq5Sun'.

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX -enc eq5Sun -six2fourtwo

To calculate all implemented ENC methods use the following command line

	$ ./codonusage.py -i FASTA -o OUTPUTPREFIX -enc all

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

**Wright (1990). The effective number of codons used in a gene. Gene 87:23-29. Equation (4):**

$$\widehat{N}_{c} = 2 + GC_{(3)} + (\frac{29}{GC^{2}_{(3)} + (1 - GC_{(3)})^{2}})$$

**Sun et al. (2012). An Improved Implementation of Effective Number of Codons (Nc). Mol. Biol. Evol. 30:191-196.**

**Equation (2):**

$$F_{CF} = \sum\limits_{i=1}^{m} (\frac {n_{i}}{n})^{2}$$

**Equation (3):**

$$F_{CF} = \sum\limits_{i=1}^{m} (\frac {n_{i}+1}{n+m})^{2}$$

**Equation (5):**

$$N_{c} = \frac {K_1 \sum\limits_{j}^{K_1} n_{j}}{\sum\limits_{j=1}^{K_1}(n_{j}F_{CF.j})} + \frac {K_2 \sum\limits_{j}^{K_2} n_{j}}{\sum\limits_{j=1}^{K_2}(n_{j}F_{CF.j})} + \frac {K_3 \sum\limits_{j}^{K_3} n_{j}}{\sum\limits_{j=1}^{K_3}(n_{j}F_{CF.j})} + \frac {K_4 \sum\limits_{j}^{K_4} n_{j}}{\sum\limits_{j=1}^{K_4}(n_{j}F_{CF.j})} + (\frac {K_6 \sum\limits_{j}^{K_6} n_{j}}{\sum\limits_{j=1}^{K_6}(n_{j}F_{CF.j})})$$

**Relative Synonymous Codon Usage:**

**Sharp et al. (1986). Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucl. Acids. Res. 14:5125-5143.**

$$RSCU_{ij} = \frac {x_{ij}} {\frac {1}{n_{i}} \sum\limits_{j=1}^{n_{i}} x_{ij}}$$

Plotting ENC distributions
--------------------------



