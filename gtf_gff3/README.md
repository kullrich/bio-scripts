## t2gGTF
extract transcript to gene table for DEG analysis

### Usage

```
usage: t2gGTF <sub-script> [options] [<arguments>...]

extracts transcript to gene table

optional arguments:
  -h, --help  show this help message and exit
  -i I        specify GTF input file
  -o O        specify output file [optional]
  -g          specify if gene names should be appended if they exist
  -b          specify if gene biotype should be appended if they exist
  -p          specify if protein id should be appended if they exist
  -v          specify if gene/transcript/protein version should be appended
  -s          specify if summary should be printed
```

### Example to extract gene, transcript, protein id from [ensembl.org](https://www.ensembl.org/index.html)

```
# download GTF file
wget http://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz
# get transcript to gene table
python t2gGTF.py -i Mus_musculus.GRCm39.105.gtf.gz -o Mus_musculus.GRCm39.105.gtf.t2g.txt -g -b -p -v -s
#55414 gene_id found
#142435 transcript_id found
#142435 protein_id found
#0 duplicate
```

### Output

output: <gene_id> <transcript_id>

output -v: <gene_id>.<v> <transcript_id>.<v>

output -g: <gene_id> <transcript_id> <gene_name>

output -g -v: <gene_id>.<v> <transcript_id>.<v> <gene_name>

output -b: <gene_id> <transcript_id> <gene_biotype>

output -b -v: <gene_id>.<v> <transcript_id>.<v> <gene_biotype>

output -p: <gene_id> <transcript_id> <protein_id>

output -p -v: <gene_id>.<v> <transcript_id>.<v> <protein_id>.<v>

output -g -b -p: <gene_id> <transcript_id> <gene_name> <gene_biotype> <protein_id>

output -g -b -p -v: <gene_id>.<v> <transcript_id>.<v> <gene_name> <gene_biotype> <protein_id>.<v>

```
head Mus_musculus.GRCm39.105.gtf.t2g.txt
#####
<gene_id>.<v>||<transcript_id>.<v>||<gene_name>||<gene_biotype>||<protein_id>.<v>
#####
ENSMUSG00000000001.5	ENSMUST00000000001.5	Gnai3	protein_coding	ENSMUSP00000000001.5
ENSMUSG00000000003.16	ENSMUST00000000003.14	Pbsn	protein_coding	ENSMUSP00000000003.8
ENSMUSG00000000003.16	ENSMUST00000114041.3	Pbsn	protein_coding	ENSMUSP00000109675.3
ENSMUSG00000000028.16	ENSMUST00000000028.14	Cdc45	protein_coding	ENSMUSP00000000028.8
ENSMUSG00000000028.16	ENSMUST00000096990.10	Cdc45	protein_coding	ENSMUSP00000094753.4
ENSMUSG00000000028.16	ENSMUST00000115585.2	Cdc45	protein_coding	ENSMUSP00000111248.2
ENSMUSG00000000028.16	ENSMUST00000231819.2	Cdc45	protein_coding	
ENSMUSG00000000031.17	ENSMUST00000132294.9	H19	lncRNA	
ENSMUSG00000000031.17	ENSMUST00000136359.8	H19	lncRNA	
ENSMUSG00000000031.17	ENSMUST00000140716.2	H19	lncRNA
```
