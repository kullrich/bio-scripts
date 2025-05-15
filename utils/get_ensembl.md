```
mkdir -p cds
mkdir -p cds_longest
```

Get EnsemblVertebrata Data - release-114

```
wget https://ftp.ensembl.org/pub/release-114/species_EnsemblVertebrates.txt
mkdir -p EnsemblVertebrates
cd EnsemblVertebrates
# get EnsemblVertebrates species file
tail -n+2 ../species_EnsemblVertebrates.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *cds.all.fa.gz -pk -e robots=off https://ftp.ensembl.org/pub/release-114/fasta/"$2"/cds/"}' \
> get_EnsemblVertebrates.sh
tail -n+2 ../species_EnsemblVertebrates.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *.114.gtf.gz -pk -e robots=off https://ftp.ensembl.org/pub/release-114/gtf/"$2"/"}' \
> get_EnsemblVertebrates_gtf.sh
sh ./get_EnsemblVertebrates.sh
sh ./get_EnsemblVertebrates_gtf.sh
cd ..
```

Get EnsemblMetazoa Data - release-61

```
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-61/species_EnsemblMetazoa.txt
mkdir -p EnsemblMetazoa
cd EnsemblMetazoa
# get EnsemblMetazoa species file
tail -n+2 ../species_EnsemblMetazoa.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *cds.all.fa.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/metazoa/release-61/fasta/"$2"/cds/"}' \
> get_EnsemblMetazoa.sh
tail -n+2 ../species_EnsemblMetazoa.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *.61.gtf.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/metazoa/release-61/gtf/"$2"/"}' \
> get_EnsemblMetazoa_gtf.sh
sh ./get_EnsemblMetazoa.sh
sh ./get_EnsemblMetazoa_gtf.sh
cd ..
```

Get EnsemblPlants Data - release-61

```
wget http://ftp.ensemblgenomes.org/pub/plants/release-61/species_EnsemblPlants.txt
mkdir -p EnsemblPlants
cd EnsemblPlants
# get EnsemblPlants species file
tail -n+2 ../species_EnsemblPlants.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *cds.all.fa.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/plants/release-61/fasta/"$2"/cds/"}' \
> get_EnsemblPlants.sh
tail -n+2 ../species_EnsemblPlants.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *.61.gtf.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/plants/release-61/gtf/"$2"/"}' \
> get_EnsemblPlants_gtf.sh
sh ./get_EnsemblPlants.sh
sh ./get_EnsemblPlants_gtf.sh
cd ..
```

Get Ensembl Rapid Release Data - 2023-01-26
https://rapid.ensembl.org/info/about/species.html

```
cd EnsemblRapidRelease
tail -n+2 ../EnsemblRapidRelease_Species_2023-01-26.csv \
| awk -F',' '{gsub(" ","_",$1); print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *cds.fa.gz -pk -e robots=off https://ftp.ensembl.org/pub/rapid-release/species/"$1"/"$7"/ensembl/geneset/"}' \
> get_EnsemblRapidRelease.sh
tail -n+2 ../EnsemblRapidRelease_Species_2023-01-26.csv \
| awk -F',' '{gsub(" ","_",$1); print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *cds.fa.gz -pk -e robots=off https://ftp.ensembl.org/pub/rapid-release/species/"$1"/"$7"/refseq/geneset/"}' \
> get_EnsemblRapidRelease_refseq.sh
tail -n+2 ../EnsemblRapidRelease_Species_2023-01-26.csv \
| awk -F',' '{gsub(" ","_",$1); print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *cds.fa.gz -pk -e robots=off https://ftp.ensembl.org/pub/rapid-release/species/"$1"/"$7"/braker/geneset/"}' \
> get_EnsemblRapidRelease_braker.sh
tail -n+2 ../EnsemblRapidRelease_Species_2023-01-26.csv \
| awk -F',' '{gsub(" ","_",$1); print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *genes.gtf.gz -pk -e robots=off https://ftp.ensembl.org/pub/rapid-release/species/"$1"/"$7"/ensembl/geneset/"}' \
> get_EnsemblRapidRelease_gtf.sh
tail -n+2 ../EnsemblRapidRelease_Species_2023-01-26.csv \
| awk -F',' '{gsub(" ","_",$1); print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *genes.gtf.gz -pk -e robots=off https://ftp.ensembl.org/pub/rapid-release/species/"$1"/"$7"/refseq/geneset/"}' \
> get_EnsemblRapidRelease_refseq_gtf.sh
tail -n+2 ../EnsemblRapidRelease_Species_2023-01-26.csv \
| awk -F',' '{gsub(" ","_",$1); print "wget --no-host-directories --cut-dirs=3 -r -l 2 -np -A *genes.gtf.gz -pk -e robots=off https://ftp.ensembl.org/pub/rapid-release/species/"$1"/"$7"/braker/geneset/"}' \
> get_EnsemblRapidRelease_braker_gtf.sh
sh ./get_EnsemblRapidRelease.sh
sh ./get_EnsemblRapidRelease_gtf.sh
sh ./get_EnsemblRapidRelease_refseq.sh
sh ./get_EnsemblRapidRelease_refseq_gtf.sh
sh ./get_EnsemblRapidRelease_braker.sh
sh ./get_EnsemblRapidRelease_braker_gtf.sh
cd ..
```

Get EnsemblFungi Data - release-61

```
wget http://ftp.ensemblgenomes.org/pub/fungi/release-61/species_EnsemblFungi.txt
mkdir -p EnsemblFungi
cd EnsemblFungi
# get EnsemblFungi species file
tail -n+2 ../species_EnsemblFungi.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *cds.all.fa.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/fungi/release-61/fasta/"$2"/cds/"}' \
> get_EnsemblFungi.sh
tail -n+2 ../species_EnsemblFungi.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *.61.gtf.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/fungi/release-61/gtf/"$2"/"}' \
> get_EnsemblFungi_gtf.sh
sh ./get_EnsemblFungi.sh
sh ./get_EnsemblFungi_gtf.sh
tail -n+2 ../species_EnsemblFungi.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=5 -r -l 2 -np -A *cds.all.fa.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/fungi/release-61/fasta/"$14"/"$2"/cds/"}' | sed 's/_core_61_114_1//g' \
> get_EnsemblFungi_core.sh
tail -n+2 ../species_EnsemblFungi.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=5 -r -l 2 -np -A *.61.gtf.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/fungi/release-61/gtf/"$14"/"$2"/"}' | sed 's/_core_61_114_1//g' \
> get_EnsemblFungi_core_gtf.sh
sh ./get_EnsemblFungi_core.sh
sh ./get_EnsemblFungi_core_gtf.sh
cd ..
```

Get EnsemblProtists Data - release-61

```
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/species_EnsemblProtists.txt
mkdir -p EnsemblProtists
cd EnsemblProtists
# get EnsemblProtists species file
tail -n+2 ../species_EnsemblProtists.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *cds.all.fa.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/protists/release-61/fasta/"$2"/cds/"}' \
> get_EnsemblProtists.sh
tail -n+2 ../species_EnsemblProtists.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=4 -r -l 2 -np -A *.61.gtf.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/protists/release-61/gtf/"$2"/"}' \
> get_EnsemblProtists_gtf.sh
sh ./get_EnsemblProtists.sh
sh ./get_EnsemblProtists_gtf.sh
tail -n+2 ../species_EnsemblProtists.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=5 -r -l 2 -np -A *cds.all.fa.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/protists/release-61/fasta/"$14"/"$2"/cds/"}' | sed 's/_core_61_114_1//g' \
> get_EnsemblProtists_core.sh
tail -n+2 ../species_EnsemblProtists.txt \
| awk -F'\t' '{print "wget --no-host-directories --cut-dirs=5 -r -l 2 -np -A *.61.gtf.gz -pk -e robots=off http://ftp.ensemblgenomes.org/pub/protists/release-61/gtf/"$14"/"$2"/"}' | sed 's/_core_61_114_1//g' \
> get_EnsemblProtists_core_gtf.sh
sh ./get_EnsemblProtists_core.sh
sh ./get_EnsemblProtists_core_gtf.sh
cd ..
```

Add taxon ID to each file and process each CDS file to get CDS ID and gene ID mappings and produce longest isoform CDS files

EnsemblVertebrates:

```
library(CRBHits)
library(Biostrings)
library(stringr)
speciesVertebrates <- readLines("species_EnsemblVertebrates.txt")[-1]
speciesVertebrates_taxIDs <- data.frame(
    cbind(
        species=stringr::str_split_fixed(speciesVertebrates, "\t", 6)[,2],
        taxID=stringr::str_split_fixed(speciesVertebrates, "\t", 6)[,4]
    )
)
t2g <- data.frame()
for(i in seq_along(speciesVertebrates_taxIDs$species)){
    sp <- speciesVertebrates_taxIDs$species[i]
    sp_taxID <- speciesVertebrates_taxIDs$taxID[i]
    # load CDS file
    cds_file <- list.files(file.path("EnsemblVertebrates",sp,"cds"),"*cds.all.fa.gz")
    gtf_file <- file.path("EnsemblVertebrates",sp,list.files(file.path("EnsemblVertebrates",sp),"*gtf.gz"))
    sp_cds <- Biostrings::readDNAStringSet(file.path("EnsemblVertebrates",sp,"cds",cds_file))
    # get longest isoform
    sp_cds_longest <- CRBHits::gtf2longest(gtffile=gtf_file, cds=sp_cds, source="ENSEMBL")$cds
    # get t2g
    sp_t2g <- data.frame(
        cbind(
            species=sp,
            taxID=sp_taxID,
            transcript=stringr::word(names(sp_cds)),
            gene=stringr::word(stringr::str_split_fixed(names(sp_cds),"gene:",2)[,2])
        )
    )
    t2g <- rbind(t2g, sp_t2g)
    # write CDS and longest CDS
    Biostrings::writeXStringSet(sp_cds, file=file.path("cds_EnsemblVertebrates",paste0(sp_taxID,".",sp,".cds.fa")))
    Biostrings::writeXStringSet(sp_cds_longest, file=file.path("cds_longest_EnsemblVertebrates",paste0(sp_taxID,".",sp,".cds.fa")))
    cat(i,sp,"\n")
}
write.table(t2g, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="EnsemblVertebrates-release-114.t2g.tsv")

```

EnsemblMetazoa:

```
library(CRBHits)
library(Biostrings)
library(stringr)
speciesMetazoa <- readLines("species_EnsemblMetazoa.txt")[-1]
speciesMetazoa_taxIDs <- data.frame(
    cbind(
        species=stringr::str_split_fixed(speciesMetazoa, "\t", 6)[,2],
        taxID=stringr::str_split_fixed(speciesMetazoa, "\t", 6)[,4]
    )
)
t2g <- data.frame()
for(i in seq_along(speciesMetazoa_taxIDs$species)){
    sp <- speciesMetazoa_taxIDs$species[i]
    sp_taxID <- speciesMetazoa_taxIDs$taxID[i]
    # load CDS file
    cds_file <- list.files(file.path("EnsemblMetazoa",sp,"cds"),"*cds.all.fa.gz")
    gtf_file <- file.path("EnsemblMetazoa",sp,list.files(file.path("EnsemblMetazoa",sp),"*gtf.gz"))
    sp_cds <- Biostrings::readDNAStringSet(file.path("EnsemblMetazoa",sp,"cds",cds_file))
    # get longest isoform
    sp_cds_longest <- CRBHits::gtf2longest(gtffile=gtf_file, cds=sp_cds, source="ENSEMBL")$cds
    # get t2g
    sp_t2g <- data.frame(
        cbind(
            species=sp,
            taxID=sp_taxID,
            transcript=stringr::word(names(sp_cds)),
            gene=stringr::word(stringr::str_split_fixed(names(sp_cds),"gene:",2)[,2])
        )
    )
    t2g <- rbind(t2g, sp_t2g)
    # write CDS and longest CDS
    Biostrings::writeXStringSet(sp_cds, file=file.path("cds_EnsemblMetazoa",paste0(sp_taxID,".",sp,".cds.fa")))
    Biostrings::writeXStringSet(sp_cds_longest, file=file.path("cds_longest_EnsemblMetazoa",paste0(sp_taxID,".",sp,".cds.fa")))
    cat(i,sp,"\n")
}
write.table(t2g, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="EnsemblMetazoa-release-61.t2g.tsv")
```

EnsemblPlants:

```
library(CRBHits)
library(Biostrings)
library(stringr)
speciesPlants <- readLines("species_EnsemblPlants.txt")[-1]
speciesPlants_taxIDs <- data.frame(
    cbind(
        species=stringr::str_split_fixed(speciesPlants, "\t", 6)[,2],
        taxID=stringr::str_split_fixed(speciesPlants, "\t", 6)[,4]
    )
)
t2g <- data.frame()
for(i in seq_along(speciesPlants_taxIDs$species)){
    sp <- speciesPlants_taxIDs$species[i]
    sp_taxID <- speciesPlants_taxIDs$taxID[i]
    # load CDS file
    cds_file <- list.files(file.path("EnsemblPlants",sp,"cds"),"*cds.all.fa.gz")
    gtf_file <- file.path("EnsemblPlants",sp,list.files(file.path("EnsemblPlants",sp),"*gtf.gz"))
    sp_cds <- Biostrings::readDNAStringSet(file.path("EnsemblPlants",sp,"cds",cds_file))
    # get longest isoform
    sp_cds_longest <- CRBHits::gtf2longest(gtffile=gtf_file, cds=sp_cds, source="ENSEMBL")$cds
    # get t2g
    sp_t2g <- data.frame(
        cbind(
            species=sp,
            taxID=sp_taxID,
            transcript=stringr::word(names(sp_cds)),
            gene=stringr::word(stringr::str_split_fixed(names(sp_cds),"gene:",2)[,2])
        )
    )
    t2g <- rbind(t2g, sp_t2g)
    # write CDS and longest CDS
    Biostrings::writeXStringSet(sp_cds, file=file.path("cds_EnsemblPlants",paste0(sp_taxID,".",sp,".cds.fa")))
    Biostrings::writeXStringSet(sp_cds_longest, file=file.path("cds_longest_EnsemblPlants",paste0(sp_taxID,".",sp,".cds.fa")))
    cat(i,sp,"\n")
}
write.table(t2g, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="EnsemblPlants-release-61.t2g.tsv")
```

EnsemblProtists:
```
library(CRBHits)
library(Biostrings)
library(stringr)
speciesProtists <- readLines("species_EnsemblProtists.txt")[-1]
speciesProtists <- speciesProtists[!duplicated(stringr::str_split_fixed(speciesProtists, "\t", 6)[,2])]
speciesProtists_taxIDs <- data.frame(
    cbind(
        species=stringr::str_split_fixed(speciesProtists, "\t", 6)[,2],
        taxID=stringr::str_split_fixed(speciesProtists, "\t", 6)[,4]
    )
)
t2g <- data.frame()
for(i in seq_along(speciesProtists_taxIDs$species)){
    sp <- speciesProtists_taxIDs$species[i]
    sp_taxID <- speciesProtists_taxIDs$taxID[i]
    # load CDS file
    cds_file <- list.files(file.path("EnsemblProtists",sp,"cds"),"*cds.all.fa.gz")
    gtf_file <- file.path("EnsemblProtists",sp,list.files(file.path("EnsemblProtists",sp),"*gtf.gz"))
    sp_cds <- Biostrings::readDNAStringSet(file.path("EnsemblProtists",sp,"cds",cds_file))
    # get longest isoform
    sp_cds_longest <- CRBHits::gtf2longest(gtffile=gtf_file, cds=sp_cds, source="ENSEMBL")$cds
    # get t2g
    sp_t2g <- data.frame(
        cbind(
            species=sp,
            taxID=sp_taxID,
            transcript=stringr::word(names(sp_cds)),
            gene=stringr::word(stringr::str_split_fixed(names(sp_cds),"gene:",2)[,2])
        )
    )
    t2g <- rbind(t2g, sp_t2g)
    # write CDS and longest CDS
    Biostrings::writeXStringSet(sp_cds, file=file.path("cds_EnsemblProtists",paste0(sp_taxID,".",sp,".cds.fa")))
    Biostrings::writeXStringSet(sp_cds_longest, file=file.path("cds_longest_EnsemblProtists",paste0(sp_taxID,".",sp,".cds.fa")))
    cat(i,sp,"\n")
}
write.table(t2g, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="EnsemblProtists-release-61.t2g.tsv")
```

EnsemblFungi:
```
library(CRBHits)
library(Biostrings)
library(stringr)
speciesFungi <- readLines("species_EnsemblFungi.txt")[-1]
speciesFungi_taxIDs <- data.frame(
    cbind(
        species=stringr::str_split_fixed(speciesFungi, "\t", 6)[,2],
        taxID=stringr::str_split_fixed(speciesFungi, "\t", 6)[,4]
    )
)
t2g <- data.frame()
for(i in seq_along(speciesFungi_taxIDs$species)){
    sp <- speciesFungi_taxIDs$species[i]
    sp_taxID <- speciesFungi_taxIDs$taxID[i]
    # load CDS file
    cds_file <- list.files(file.path("EnsemblFungi",sp,"cds"),"*cds.all.fa.gz")
    gtf_file <- file.path("EnsemblFungi",sp,list.files(file.path("EnsemblFungi",sp),"*gtf.gz"))
    sp_cds <- Biostrings::readDNAStringSet(file.path("EnsemblFungi",sp,"cds",cds_file))
    # get longest isoform
    sp_cds_longest <- CRBHits::gtf2longest(gtffile=gtf_file, cds=sp_cds, source="ENSEMBL")$cds
    # get t2g
    sp_t2g <- data.frame(
        cbind(
            species=sp,
            taxID=sp_taxID,
            transcript=stringr::word(names(sp_cds)),
            gene=stringr::word(stringr::str_split_fixed(names(sp_cds),"gene:",2)[,2])
        )
    )
    t2g <- rbind(t2g, sp_t2g)
    # write CDS and longest CDS
    Biostrings::writeXStringSet(sp_cds, file=file.path("cds_EnsemblFungi",paste0(sp_taxID,".",sp,".cds.fa")))
    Biostrings::writeXStringSet(sp_cds_longest, file=file.path("cds_longest_EnsemblFungi",paste0(sp_taxID,".",sp,".cds.fa")))
    cat(i,sp,"\n")
}
write.table(t2g, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="EnsemblFungi-release-58.t2g.tsv")

```

EnsemblRapidRelease:

```

```

Wormbase:

```
library(CRBHits)
library(Biostrings)
library(stringr)
speciesWormbase <- list.files(".","*.gz")
for(i in seq_along(speciesWormbase)){
    sp_cds <- Biostrings::readDNAStringSet(speciesWormbase[i])
    #names(sp_cds) <- gsub("transcript=>","",apply(stringr::str_split_fixed(names(sp_cds),"\t", 3)[,2:3],1,function(x) paste0(x[1]," ",x[2])))
    sp_cds_longest <- CRBHits::isoform2longest(cds=sp_cds, source="WORMBASE")
    Biostrings::writeXStringSet(sp_cds_longest, file=file.path("../cds_longest_WS288",paste0(gsub(".gz","",speciesWormbase[i]))))
    cat(i, "\n")
    cat(length(sp_cds), "\n")
    cat(length(sp_cds_longest), "\n")
}

```
