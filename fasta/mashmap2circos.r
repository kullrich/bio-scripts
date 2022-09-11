#!/usr/bin/env Rscript
library(circlize)
library(Biostrings)
library(dplyr)
library(stringr)
library(Redmonder)
require("argparse")

col2transparent <- function(col, alpha.perc=0){
    if(methods::is(col, "character")){
        if(length(col)==1){
            alpha = (100 - alpha.perc) * 255 / 100
            R = col2rgb(col)[1]
            G = col2rgb(col)[2]
            B = col2rgb(col)[3]
            return(rgb(R, G, B, alpha, maxColorValue=255))
        } else {
            return(unlist(lapply(col,
                function(x) {col2transparent(x, alpha.perc)})))
        }
    }
    if(methods::is(col, "palette")){
        return(palette(unlist(lapply(col,
            function(x) {col2transparent(x, alpha.perc)}))))
    }
}

parser <- ArgumentParser()
parser$add_argument("-m", "--mashmap", help="path to mashmap", default="/home/ullrich/opt/ngs/mashmap-Linux64-v2.0/mashmap")
parser$add_argument("-r", "--reference", help="path to reference fasta file")
parser$add_argument("-q", "--query", help="path to query fasta file")
parser$add_argument("-c", "--cpu", type="integer", help = "number of threads (default: 1)", default=1)
parser$add_argument("-p", "--prefix", help="prefix for the outputs", default="mashmap.circos")
parser$add_argument("-pw", "--pwidth", help="png width (default: 1200)", type="integer", default=1200)
parser$add_argument("-ph", "--pheight", help="png height (default: 800)", type="integer", default=800)
parser$add_argument("-n", "--nosplit", help="nosplit", action='store_true')
parser$add_argument("-pi", "--perc_identity", help="threshold for identity (default : 85)", type="integer", default=85)
parser$add_argument("-s", "--segLength", help="mapping segment length (default : 5,000)", type="integer", default=5000)
parser$add_argument("-k", "--kmer", help="kmer size <= 16 (default : 16)", type="integer", default=16)
parser$add_argument("-f", "--filter_mode", help="filter modes in mashmap: 'map', 'one-to-one' or 'none' (default: map)", default="map")
parser$add_argument("-minl", "--minLength", help="minimum scaffold length (default : 5,000)", type="integer", default=5000)
parser$add_argument("-cex", "--cex", help="cex (default : 1.0)", type="double", default=1.0)
parser$add_argument("-alpha", "--alpha", help="alpha (default : 10)", type="integer", default=10)
args <- parser$parse_args()

MASHMAP <- args$mashmap
REF <- args$reference
QUERY <- args$query
THREADS <- args$cpu
OUTPUT <- args$prefix
WIDTH <- args$pwidth
HEIGHT <- args$pheight
NOSPLIT <- args$nosplit
PERCIDENT <- args$perc_identity
SEGLENGTH <- args$segLength
KMER <- args$kmer
FILTERMODE <- args$filter_mode
MINLENGTH <- args$minLength
CEX <- args$cex
ALPHA <- args$alpha

if(NOSPLIT){
    system(paste0(MASHMAP," -r ",REF," -q ",QUERY," -t ",THREADS," -s ",SEGLENGTH," --perc_identity ",PERCIDENT," -k ",KMER," -f ",FILTERMODE," -o ",OUTPUT,".mashmap.out --noSplit"))
} else{
    system(paste0(MASHMAP," -r ",REF," -q ",QUERY," -t ",THREADS," -s ",SEGLENGTH," --perc_identity ",PERCIDENT," -k ",KMER," -f ",FILTERMODE," -o ",OUTPUT,".mashmap.out"))
}

REF.genome <- Biostrings::readDNAStringSet(REF)
QUERY.genome <- Biostrings::readDNAStringSet(QUERY)

REF.genome.sizes <- data.frame(name=paste0("R:",stringr::word(names(REF.genome))), start=0, end=width(REF.genome))
QUERY.genome.sizes <- data.frame(name=paste0("Q:",stringr::word(names(QUERY.genome))), start=0, end=width(QUERY.genome))

genome.df <- rbind(REF.genome.sizes, QUERY.genome.sizes) %>% dplyr::filter(end > MINLENGTH)

bed <- read.table(paste0(OUTPUT,".mashmap.out"),sep=" ")
colnames(bed) <- c("query.chr","query.length","query.start","query.end","strand","ref.chr","ref.length","ref.start","ref.end","identity")
bed$ref.chr <- paste0("R:",bed$ref.chr)
bed$query.chr <- paste0("Q:",bed$query.chr)

bed1 <- bed %>% dplyr::filter(ref.chr %in% genome.df$name) %>% dplyr::filter(query.chr %in% genome.df$name) %>%dplyr::select(ref.chr, ref.start, ref.end, identity)
bed2 <- bed %>% dplyr::filter(ref.chr %in% genome.df$name) %>% dplyr::filter(query.chr %in% genome.df$name) %>%dplyr::select(query.chr, query.start, query.end, identity)
colnames(bed1) <- colnames(bed2) <- c("chr", "start", "end", "identity")

png(paste0(OUTPUT,".mashmap.png"), width=WIDTH, height=HEIGHT)
circos.genomicInitialize(genome.df, labels.cex=CEX)
circos.genomicLink(bed1, bed2, col=col2transparent(Redmonder::redmonder.pal(length(levels(as.factor(bed1$chr))), "qMSOStd"), ALPHA)[as.factor(bed1$chr)], border=0)
dev.off()
