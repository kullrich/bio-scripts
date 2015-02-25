rbhplus
=======

Description
-----------
Performs reciprocal BLAST+ [-prog blast+] or ghostx [-prog ghostx] searches
between two species by either starting and performing BLAST+/ghostx search or
providing BLAST+/ghostx produced table files. BLAST+/ghostx binary are
supposed to be in the PATH otherwise specify with [pp] option. Mandatory are
FASTA format sequences for two species and their sequence type (see help). For
single species reciprocal BLAST+/ghostx use the same input sequence as query
and database, please "omit the selfblast option in this case!". Optional
parameters can be set, e.g. if you want to add selfbalst (within species)
BLAST+/ghostx comparison [default: only between; can be changed with the
selfblast option], if BLAST+/ghostx should build the necessary BLAST+/ghostx
databases with makeblastdb/ghostx db [default: false] and some BLAST+/ghostx
specific options (see help).

* performs BLAST+ searches for given species
* perform ghostx searches for given species
* performs reciprocal best hit sorting
* extracts selfscore for each query

Dependencies
------------
* python version > 2.6
* rbhplus depends on BLAST+. Please download the current BLAST+ for your system (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
* rbhplus depends on ghostx. Please download the current ghostx for yor system (http://www.bi.cs.titech.ac.jp/ghostx/).

Pre-build binaries
------------------
For the Linux system pre-build binaries exists in the dist folder with all dependent libraries.
The program can be invoked with:
    ./dist/rbhplus/rbhplus

Binaries were produced with PyInstaller (https://github.com/pyinstaller/pyinstaller/wiki).

Usage
-----

    optional arguments:
    -h, --help  show this help message and exit
    -prog   PROG    specify if BLAST+ or ghostx should be used for sequence comparison [default: blast+]
    -pp PP  specify BLAST+/ghostx binary path [default: PATH]
    -selfblast  {True,False}    specify if BLAST+/ghostx search should be run also within each species to extract selfscores. Can also run as standalone option with step option set to [-step self]. This step is necessary if you would like to sort the BLAST+/ghostx result by selfscores.
    -step   {0,1,2,3,4,self}    specify which step should be performed: [1]->BLAST+/ghostx initially builds the necessary BLAST+/ghostx sequence databases. [2]->BLAST+/ghostx search is performed. [3]->BLAST+/ghostx output is parsed. [4]->Reciprocal Best Hits (RBHs) are calculated. Optional: [self]->BLAST+/ghostx initially builds the necessary BLAST+/ghostx sequence databases if it is not already existing and selfscores are obtained for each species. Optional: [0]->BLAST+/ghostx builds the BLAST+/ghostx sequence databases and exists.
    -q  Q   specify query
    -qt {prot,nucl} specify query type [default: prot] or [nucl]
    -db DB  specify database
    -dbt    {prot,nucl} specify database type [default: prot] or [nucl]
    -e  E   specify the evalue cutoff [default: 1E-5]
    -m  M   specify the BLAST+/ghostx matrix [default: BLOSUM62]
    -n  N   specfiy number of threads to use
    -op OP  specify output path [default: .]
    -p  P   specify pident cutoff for BLAST+/ghostx result filtering [default: 30]
    -l  L   specify alignment length cutoff for BLAST+/ghostx result filtering [default: 80]
    -covs   COVS    specify subject coverage for BLAST+ result filtering [default: 0.0]; can only be applied with BLAST+ since ghostx does not provide qlen and slen yet.
    -covq   COVQ    specify query coverage for BLAST+ result filtering [default: 0.0]; can only be applied with BLAST+ since ghostx does not provide qlen and slen yet.
    -pid    {evalue,static,rost1999}    choose pident filtering for BLAST+/ghostx result filtering; either as static [pident+length], evalue or rost1999
    -sort   {pure,evl,pxl,selfscore}    choose reciprocal best hit sorting: either pure [only first hit of blast+/ghostx input will be processed]; evl [hits that fulfill pident cutoff will be sorted by evalue]; pxl [hits that fulfill pident cutoff will be sorted by pident*length] or selfscore [hits that fulfill pident cutoff will be sorted by selfscore provided via the -selfscore OPTION]
    -selfscore  SELFSCORE   SELFSCORE   specify selfscore files for the query and the database species [tab seperated format (ID SCORE)]. You need to have two file one for the query and one for the database. Selfscore files can be produced either by using [step: self] or by including the [selfblast] option, in this case the selfscores will be calculated anyway.
    -pre    {True,False}    specify if pre-existing BLAST+/ghostx output files should be used; e.g. if BLAST+/ghostx output was produced elsewhere [default: False]. Please note that BLAST+ output needs to be in the format "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident gaps score"; ghostx in the current version does not support individual output format so we stick to standard output here. Also there exists a name convention "query input".vs."database input"
    -prep   PREP    specify path with pre-existing BLAST+/ghostx output [default: .]
    -out    OUT specify output file
    -v  {True,False}    verbose output

Pre-build databases
-------------------

Pre-build databases exist for core eukaryotic gene sets:

* Van Bel M, Proost S, Wischnitzki E, Movahedi S, Scheerlinck C, Van de Peer Y, Vandepoele K: Dissecting Plant Genomes with the PLAZA Comparative Genomics Platform. Plant physiology 2012, 158(2):590-600.
* Parra G, Bradnam K, Ning Z, Keane T, Korf I: Assessing the gene space in draft genomes. Nucleic Acids Research 2009, 37(1):289-297.

Examples
--------

To extract reciprocal best hits of two sequence files:

Specify your query FASTA file with additionally given the type of your query sequence data ['prot' or 'nucl'].
Specify your database FASTA file with additionally given the type of your database sequence data ['prot' or 'nucl'].

	$ ./rbhplus.py -prog blast+ -q QFASTA -qt prot -db DBFASTA -dbt prot -sort pure -out out -step 1

or to use ghostx (both programs should be in your PATH, otherwise please provide the folder of located binary files with the [-pp] option):

	$ ./rbhplus.py -prog ghostx -q QFASTA -qt prot -db DBFASTA -dbt prot -sort pure -out out -step 1

Include coverage as filter option:

	$ ./rbhplus.py -q QFASTA -qt prot -db DBFASTA -dbt prot -sort pure -out out -covs 0.5 -covq 0.5 -step 1
	
Filter not by static protein identity and length [-pid static] but by evalue [-pid evalue]:

	$ ./rbhplus.py -q QFASTA -qt prot -db DBFASTA -dbt prot -sort pure -out out -evalue 1e-20 -pid evalue -step 1
	
Filter not by static protein identity and length [-pid static] but by dynamic formula [-pid rost1999] to account for the twilight-zone according to (Rost B: Twilight zone of protein sequence alignments. Protein Engineering 1999, 12(2):85-94.)

To extract reciprocal best hits of two sequence files including reciprocal best hits within each sequence file:

	$ ./rbhplus.py -q QFASTA -qt prot -db DBFASTA -dbt prot -sort pure -selfblast True -out out -step 1

To sort BLAST+/ghostx output according to selfscore [-sort selfscore]:

	$ ./rbhplus.py -q QFASTA -qt prot -db DBFASTA -dbt prot -sort selfscore -selfblast True -out out -step 1
	
To just extract the selfscore for your sequence files either using BLAST+ or ghostx:

	$ ./rbhplus.py -prog blast+ -q QFASTA -qt prot -db DBFASTA -dbt prot -selfblast True -out out -step self
    $ ./rbhplus.py -prog ghostx -q QFASTA -qt prot -db DBFASTA -dbt prot -selfblast True -out out -step self

Speed things up by using more threads:

	$ ./rbhplus.py -prog blast+ -q QFASTA -qt prot -db DBFASTA -dbt prot -selfblast True -out out -step self -n 4

Output
------

As output different tables will be generated:

    RBHs between query and database sequence file: out.qd
    RBHs within query sequence file: out.qq
    RBHs within database sequence file: out.dd

Bugs and Errors
---------------

rbhplus is under active research development at the University of Marburg. Please report any errors or requests to Kristian Ullrich (kristian.ullrich@biologie.uni-marburg.de).
