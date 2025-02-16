#!/usr/bin/env Rscript
library(circlize)
library(Biostrings)
library(dplyr)
library(stringr)
#library(Redmonder)
library(intervals)
library(colorspace)
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

CRBHitsColors <- function(n,
    alpha.perc=0
    ){
    palette(c(
        "#CBC106", "#27993C", "#1C6838",
        "#8EBCB5", "#389CA7", "#4D83AB",
        "#CB7B26", "#BF565D", "#9E163C"))
    my.palette <- palette()
    my.n.palette <- colorRampPalette(my.palette)(n)
    my.out.palette <- col2transparent(my.n.palette, alpha.perc)
    return(my.out.palette)
}

file_ext <- function(filename){
    return(rev(strsplit(basename(filename),"\\.")[[1]])[1])
}

merge_intervals <- function(bed, distance){
    df <- bed %>% arrange(query.chr, ref.chr, query.start, ref.start)
    df_split <- df %>%
        group_by(query.chr, ref.chr, strand) %>%
        group_split()
    merged_list <- lapply(df_split, function(sub_df) {
        query_intervals <- Intervals(cbind(sub_df$query.start, sub_df$query.end), closed = c(TRUE, TRUE))
        ref_intervals <- Intervals(cbind(sub_df$ref.start, sub_df$ref.end), closed = c(TRUE, TRUE))
        query_intervals_expanded <- intervals::reduce(Intervals(cbind(sub_df$query.start - distance, sub_df$query.end + distance), closed = c(TRUE, TRUE)))
        ref_intervals_expanded <- intervals::reduce(Intervals(cbind(sub_df$ref.start - distance, sub_df$ref.end + distance), closed = c(TRUE, TRUE)))
        query_overlap <- interval_overlap(query_intervals, query_intervals_expanded)
        ref_overlap <- interval_overlap(ref_intervals, ref_intervals_expanded)
        sub_df$q_overlap <- query_overlap
        sub_df$r_overlap <- ref_overlap
        merged_sub_df <- sub_df %>%
            group_by(q_overlap, r_overlap) %>%
            summarise(
                query.chr = first(query.chr),
                query.length = first(query.length),
                query.start = min(query.start),
                query.end = max(query.end),
                strand = first(strand),
                ref.chr = first(ref.chr),
                ref.length = first(ref.length),
                ref.start = min(ref.start),
                ref.end = max(ref.end),
                V10 = sum(V10),
                V11 = sum(V11),
                V12 = sum(V12),
                identity = mean(identity),
                complexity = mean(complexity),
                .groups = "drop"
            )
        return(merged_sub_df)
        })
    merged_df <- bind_rows(merged_list)
    merged_df <- merged_df %>% select(-q_overlap, -r_overlap)
    return(as.data.frame(merged_df))
}

parser <- ArgumentParser()
parser$add_argument("-m", "--mashmap", help="path to mashmap", default="/Users/ullrich/projects/scientific/books/statisctical_population_genomics_2nd_edition/data/MashMap-3.1.3/build/bin/mashmap")
parser$add_argument("-r", "--reference", help="path to reference fasta file")
parser$add_argument("-q", "--query", help="path to query fasta file")
parser$add_argument("-mo", "--moutput", help="path to pre-existing mashmap output file")
parser$add_argument("-cytoband", "--cytoband", help="path to cytoband file (chr IDs need to match, chr order is kept)")
parser$add_argument("-c", "--cpu", type="integer", help = "number of threads (default: 1)", default=1)
parser$add_argument("-p", "--prefix", help="prefix for the outputs", default="mashmap.circos")
parser$add_argument("-pngw", "--pngwidth", help="png width (default: 1200)", type="integer", default=1200)
parser$add_argument("-pngh", "--pngheight", help="png height (default: 1200)", type="integer", default=1200)
parser$add_argument("-pdfw", "--pdfwidth", help="pdf width (default: 12)", type="integer", default=12)
parser$add_argument("-pdfh", "--pdfheight", help="pdf height (default: 12)", type="integer", default=12)
parser$add_argument("-pr", "--prefix_ref", help="reference chromosome prefix (default: R)", default="R")
parser$add_argument("-pq", "--prefix_query", help="query chromosome prefix (default: Q)", default="Q")
parser$add_argument("-n", "--nosplit", help="nosplit", action='store_true')
parser$add_argument("-mpi", "--mashmap_perc_identity", help="threshold for identity for mashmap command (default: 85)", type="integer", default=85)
parser$add_argument("-pi", "--perc_identity", help="threshold for identity after mashmap (default: 85)", type="integer", default=85)
parser$add_argument("-s", "--segLength", help="mapping segment length (default : 5,000)", type="integer", default=5000)
parser$add_argument("-k", "--kmer", help="kmer size <= 19 (default : 19)", type="integer", default=19)
parser$add_argument("-f", "--filter_mode", help="filter modes in mashmap: 'map', 'one-to-one' or 'none' (default: map)", default="map")
parser$add_argument("-minl", "--minLength", help="minimum scaffold length (default : 5,000)", type="integer", default=5000)
parser$add_argument("-merge", "--merge", help="merge adjacent blocks if distance is smaller than block_gap", action='store_true')
parser$add_argument("-bg", "--block_gap", help="block gap between adjacent blocks as threshold to merge (default : 2,500)", type="integer", default=2500)
parser$add_argument("-cex", "--cex", help="cex (default : 1.0)", type="double", default=1.0)
parser$add_argument("-alpha", "--alpha", help="alpha (default : 10)", type="integer", default=10)
parser$add_argument("-strandcolor", "--strandcolor", help="define how different strand should be colored (lighten, darken, shift, invert) (default : darken)", default="darken")
parser$add_argument("-rev", "--revert", help="revert chromosome order on plot", action='store_true')
args <- parser$parse_args()

MASHMAP <- args$mashmap
REF <- args$reference
QUERY <- args$query
MASHMAPOUT <- args$moutput
CYTOBANDFILE <- args$cytoband
THREADS <- args$cpu
OUTPUT <- args$prefix
PNGWIDTH <- args$pngwidth
PNGHEIGHT <- args$pngheight
PDFWIDTH <- args$pdfwidth
PDFHEIGHT <- args$pdfheight
PREFIXREF <- args$prefix_ref
PREFIXQUERY <- args$prefix_query
NOSPLIT <- args$nosplit
MASHMAPPERCIDENT <- args$mashmap_perc_identity
PERCIDENT <- args$perc_identity
SEGLENGTH <- args$segLength
KMER <- args$kmer
FILTERMODE <- args$filter_mode
MINLENGTH <- args$minLength
MERGE <- args$merge
BLOCKGAP <- args$block_gap
CEX <- args$cex
ALPHA <- args$alpha
STRANDCOLOR <- args$strandcolor
REVERT <- args$revert

if (is.null(MASHMAPOUT)) {
    if (NOSPLIT){
        system(paste0(MASHMAP," -r ",REF," -q ",QUERY," -t ",THREADS," -s ",SEGLENGTH," --perc_identity ",MASHMAPPERCIDENT," -k ",KMER," -f ",FILTERMODE," -o ",OUTPUT,".mashmap.out --noSplit"))
    } else{
        system(paste0(MASHMAP," -r ",REF," -q ",QUERY," -t ",THREADS," -s ",SEGLENGTH," --perc_identity ",MASHMAPPERCIDENT," -k ",KMER," -f ",FILTERMODE," -o ",OUTPUT,".mashmap.out"))
    }
}

if (file_ext(REF) %in% c("fna", "fasta", "fa")) {
    REF.genome <- Biostrings::readDNAStringSet(REF)
    REF.genome.sizes <- data.frame(name=paste0(PREFIXREF,":",stringr::word(names(REF.genome))), start=0, end=width(REF.genome))
} else if (file_ext(REF) %in% c("fai")) {
    REF.genome.fai <- read.table(REF)
    REF.genome.sizes <- data.frame(name=paste0(PREFIXREF,":",REF.genome.fai$V1), start=0, end=REF.genome.fai$V2)
}
if (file_ext(QUERY) %in% c("fna", "fasta", "fa")) {
    QUERY.genome <- Biostrings::readDNAStringSet(QUERY)
    QUERY.genome.sizes <- data.frame(name=paste0(PREFIXQUERY,":",stringr::word(names(QUERY.genome))), start=0, end=width(QUERY.genome))
} else if (file_ext(QUERY) %in% c("fai")) {
    QUERY.genome.fai <- read.table(QUERY)
    QUERY.genome.sizes <- data.frame(name=paste0(PREFIXQUERY,":",QUERY.genome.fai$V1), start=0, end=QUERY.genome.fai$V2)
}

if (REVERT) {
    genome.df <- rbind(QUERY.genome.sizes[dim(QUERY.genome.sizes)[1]:1, ], REF.genome.sizes) %>% dplyr::filter(end > MINLENGTH)
    } else {
    genome.df <- rbind(REF.genome.sizes, QUERY.genome.sizes) %>% dplyr::filter(end > MINLENGTH)
}

if (!is.null(CYTOBANDFILE)) {
    genome.df <- read.table(CYTOBANDFILE)
    colnames(genome.df) <- c("name", "start", "end", "bandname", "stain")
}

if (is.null(MASHMAPOUT)) {
    bed <- read.table(paste0(OUTPUT,".mashmap.out"),sep="\t")
} else {
    bed <- read.table(MASHMAPOUT,sep="\t")
}
colnames(bed) <- c("query.chr","query.length","query.start","query.end","strand","ref.chr","ref.length","ref.start","ref.end","V10","V11","V12","identity","complexity")
bed$ref.chr <- paste0(PREFIXREF,":",bed$ref.chr)
bed$query.chr <- paste0(PREFIXQUERY,":",bed$query.chr)
bed$identity <- as.numeric(gsub("id:f:", "", bed$identity))*100
bed$complexity <- as.numeric(gsub("kc:f:", "", bed$complexity))*100
bed <- bed %>% dplyr::filter(identity >= PERCIDENT)
if (MERGE) {
    bed <- merge_intervals(bed, BLOCKGAP)
}

bed1 <- bed %>% dplyr::filter(ref.chr %in% genome.df$name) %>% dplyr::filter(query.chr %in% genome.df$name) %>%dplyr::select(ref.chr, ref.start, ref.end, identity)
bed2 <- bed %>% dplyr::filter(ref.chr %in% genome.df$name) %>% dplyr::filter(query.chr %in% genome.df$name) %>%dplyr::select(query.chr, query.start, query.end, identity)
colnames(bed1) <- colnames(bed2) <- c("chr", "start", "end", "identity")

chr_colors <- CRBHitsColors(length(levels(as.factor(bed1$chr))), ALPHA)[as.factor(bed1$chr)]
inverted_colors <- desaturate(rev(chr_colors))
chr_colors_hsl <- as(hex2RGB(chr_colors), "HLS")
chr_colors_hsl@coords[, 1] <- (chr_colors_hsl@coords[, 1] + 30) %% 360  # Hue shift
shifted_colors <- hex(HLS(chr_colors_hsl@coords))
light_colors <- lighten(chr_colors, amount = 0.3)
dark_colors <- darken(chr_colors, amount = 0.3)
if (STRANDCOLOR == "lighten") {
    plot_colors <- ifelse(bed$strand == "-", light_colors, chr_colors)
} else if (STRANDCOLOR == "darken") {
    plot_colors <- ifelse(bed$strand == "-", dark_colors, chr_colors)
} else if (STRANDCOLOR == "shift") {
    plot_colors <- ifelse(bed$strand == "-", shifted_colors, chr_colors)
} else if (STRANDCOLOR == "invert") {
    plot_colors <- ifelse(bed$strand == "-", inverted_colors, chr_colors)
} else {
    plot_colors <- chr_colors
}

png(paste0(OUTPUT,".mashmap.png"), width=PNGWIDTH, height=PNGHEIGHT)
if (!is.null(CYTOBANDFILE)) {
    circos.initializeWithIdeogram(genome.df, labels.cex=CEX, sort.chr=FALSE)
} else {
    circos.genomicInitialize(genome.df, labels.cex=CEX)
}
#circos.genomicLink(bed1, bed2, col=col2transparent(Redmonder::redmonder.pal(length(levels(as.factor(bed1$chr))), "qMSOStd"), ALPHA)[as.factor(bed1$chr)], border=0)
#circos.genomicLink(bed1, bed2, col = CRBHits::CRBHitsColors(length(levels(as.factor(bed1$chr))), ALPHA)[as.factor(bed1$chr)], border = CRBHits::CRBHitsColors(length(levels(as.factor(bed1$chr))), 0)[as.factor(bed1$chr)])
circos.genomicLink(bed1, bed2, col = plot_colors, border = NA, inverse = bed$strand=="-")
dev.off()

pdf(paste0(OUTPUT,".mashmap.pdf"), width=PDFWIDTH, height=PDFHEIGHT)
if (!is.null(CYTOBANDFILE)) {
    circos.initializeWithIdeogram(genome.df, labels.cex=CEX, sort.chr=FALSE)
} else {
    circos.genomicInitialize(genome.df, labels.cex=CEX)
}
#circos.genomicLink(bed1, bed2, col=col2transparent(Redmonder::redmonder.pal(length(levels(as.factor(bed1$chr))), "qMSOStd"), ALPHA)[as.factor(bed1$chr)], border=0)
#circos.genomicLink(bed1, bed2, col = CRBHits::CRBHitsColors(length(levels(as.factor(bed1$chr))), ALPHA)[as.factor(bed1$chr)], border = CRBHits::CRBHitsColors(length(levels(as.factor(bed1$chr))), 0)[as.factor(bed1$chr)])
circos.genomicLink(bed1, bed2, col = plot_colors, border = NA, inverse = bed$strand=="-")
dev.off()
