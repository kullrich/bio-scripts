library(circlize)
library(Biostrings)
library(dplyr)
library(CRBHits)

MASHMAP <- "/opt/mashmap-Linux64-v2.0/mashmap"
REF <- "ref.fasta"
QUERY <- "query.fasta"
THREADS <- 48
OUTPUT <- "ref_query"
WIDTH <- 1200
HEIGHT <- 800
NOSPLIT <- TRUE

if(NOSPLIT){
    system(paste0(MASHMAP," -r ",REF," -q ",QUERY," -t ",THREADS," -o ",OUTPUT,".mashmap.out --noSplit"))
} else{
    system(paste0(MASHMAP," -r ",REF," -q ",QUERY," -t ",THREADS," -o ",OUTPUT,".mashmap.out"))
}

REF.genome <- Biostrings::readDNAStringSet(REF)
QUERY.genome <- Biostrings::readDNAStringSet(QUERY)

REF.genome.sizes <- data.frame(name=paste0("AA1:",stringr::word(names(REF.genome))), start=0, end=width(REF.genome))
QUERY.genome.sizes <- data.frame(name=paste0("AA2:",stringr::word(names(QUERY.genome))), start=0, end=width(QUERY.genome))

genome.df <- rbind(REF.genome.sizes, QUERY.genome.sizes) %>% dplyr::filter(end > 2000000)

bed <- read.table(paste0(OUTPUT,".mashmap.out"),sep=" ")
colnames(bed) <- c("query.chr","query.length","query.start","query.end","strand","ref.chr","ref.length","ref.start","ref.end","identity")
bed$ref.chr <- paste0("AA1:",bed$ref.chr)
bed$query.chr <- paste0("AA2:",bed$query.chr)

bed1 <- bed %>% dplyr::filter(ref.chr %in% genome.df$name) %>% dplyr::filter(query.chr %in% genome.df$name) %>%dplyr::select(ref.chr, ref.start, ref.end, identity)
bed2 <- bed %>% dplyr::filter(ref.chr %in% genome.df$name) %>% dplyr::filter(query.chr %in% genome.df$name) %>%dplyr::select(query.chr, query.start, query.end, identity)
colnames(bed1) <- colnames(bed2) <- c("chr", "start", "end", "identity")

png(paste0(OUTPUT,".mashmap.png"),width=WIDTH,height=HEIGHT)
circos.genomicInitialize(genome.df, labels.cex = 0.5)
circos.genomicLink(bed1, bed2, col = CRBHitsColors(length(levels(as.factor(bed1$chr))), 10)[as.factor(bed1$chr)], border = 0)
dev.off()
