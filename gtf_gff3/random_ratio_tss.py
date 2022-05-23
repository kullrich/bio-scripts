#!/usr/bin/env python3
# SQANTI plugin: random ratio_TSS
# Authors: Kristian K Ullrich (ullrich@evolbio.mpg.de)
# License: MIT


__author__  = "ullrich@evolbio.mpg.de"
__version__ = '0.1'  # Python 3.7


import pdb
import os, re, sys, subprocess, timeit, glob, copy, gzip
import shutil
import distutils.spawn
import itertools
import bisect
import argparse
import math
import pandas
import random
import pybedtools
import numpy as np
import pyranges as pr
from scipy import mean
from collections import defaultdict, Counter, namedtuple
from collections.abc import Iterable
from csv import DictWriter, DictReader
from multiprocessing import Process
from gtfparse import read_gtf


try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Unable to import Biopython! Please make sure Biopython is installed.", file=sys.stderr)
    sys.exit(-1)


def get_Ns(reference, verbose = True):
    if verbose:
        print("parse Ns from %s" % reference)
    ref = SeqIO.parse(reference, "fasta")
    Ndict = {"Chromosome": [], "Start": [], "End": [], "Strand": []}
    for rec in ref:
        if verbose:
            print(rec.id)
        start_pos = 0
        counter = 0
        gap = False
        gap_length = 0
        for nuc in rec.seq:
            if nuc == "N":
                if gap_length == 0:
                    start_pos = counter
                    gap_length = 1
                    gap = True
                else:
                    gap_length += 1
            else:
                if gap:
                    if verbose:
                        print(rec.id + "\t" + str(start_pos) + "\t" + str(start_pos + gap_length))
                    Ndict["Chromosome"].append(rec.id)
                    Ndict["Start"].append(start_pos)
                    Ndict["End"].append(start_pos + gap_length)
                    Ndict["Strand"].append("+")
                    Ndict["Chromosome"].append(rec.id)
                    Ndict["Start"].append(start_pos)
                    Ndict["End"].append(start_pos + gap_length)
                    Ndict["Strand"].append("-")
                    gap_length = 0
                    gap = False
            counter += 1
    Ns = pr.from_dict(Ndict)
    return(Ns)


def get_bam_header2(bam, o_dir):
    out=o_dir + "/chr_order.txt"
    os.system("samtools view -H {b} | grep '^@SQ' | sed 's/@SQ\tSN:\|LN://g'  > {o}".format(b=bam, o=out))
    return(out)


def read_bam_header(bam_header, min_size):
    Bdict = {}
    with open(bam_header) as in_handle:
        for line in in_handle:
            #_, SN, LN = line.strip().split('\t')
            SN, LN = line.strip().split('\t')
            SN = SN.replace("SN:", "")
            LN = int(LN.replace("LN:", ""))
            if LN >= min_size:
                Bdict[SN] = LN
    return(Bdict)


def get_sample_position(sites, random_tss_per_chr):
    Sdict = {"Chromosome": [], "Start": [], "End": [], "Strand": []}
    for i in range(random_tss_per_chr):
        s_tmp = sites.as_df().sample(1)
        s_tmp_chromosome = s_tmp["Chromosome"].values[0]
        s_tmp_start = s_tmp["Start"].values[0]
        s_tmp_end = s_tmp["End"].values[0]
        s_tmp_strand = s_tmp["Strand"].values[0]
        s_tmp_range = range(s_tmp_start, s_tmp_end)
        s_tmp_random = random.sample(s_tmp_range, 1)[0]
        Sdict["Chromosome"].append(s_tmp_chromosome)
        Sdict["Start"].append(s_tmp_random)
        Sdict["End"].append(s_tmp_random + 1)
        Sdict["Strand"].append(s_tmp_strand)
    Sp = pr.from_dict(Sdict)
    return(Sp)


def create_bed(k, gtf_type, pos, o_dir, bam_header_file, seqlen=100):
    inside_tmp = o_dir + "/" + k + "_" + gtf_type + ".inside.bed_tmp"
    outside_tmp = o_dir + "/" + k + "_" + gtf_type + ".outside.bed_tmp"
    with open(inside_tmp, "w") as inside:
        with open(outside_tmp, "w") as outside:
            for i in range(len(pos)):
                CHR, START, END, STRAND = pos.as_df().iloc()[i]
                if STRAND == "+":
                    #out_handle.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n' % (CHR, gtf_type, 'gene', START + 1, END + seqlen - 1, '.', STRAND, '.', 'gene_id "g' + str(i) + '"'))
                    #out_handle.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n' % (CHR, gtf_type, 'transcript', START + 1, END + seqlen - 1, '.', STRAND, '.', 'gene_id "g' + str(i) + '"; transcript_id "t' + str(i) + '"'))
                    #out_handle.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n' % (CHR, gtf_type, 'exon', START + 1, END + seqlen - 1, '.', STRAND, '.', 'gene_id "g' + str(i) + '"; transcript_id "t' + str(i) + '"; exon_id "e' + str(i) + '"'))
                    inside.write('%s\t%i\t%i\t%s\t%i\t%s\n' % (CHR, START + 1, END + seqlen - 1, 't' + str(i), 0, STRAND))
                    outside.write('%s\t%i\t%i\t%s\t%i\t%s\n' % (CHR, START - seqlen, START, 't' + str(i), 0, STRAND))
            else:
                #out_handle.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n' % (CHR, gtf_type, 'gene', START - seqlen, END, '.', STRAND, '.', 'gene_id "g' + str(i) + '"'))
                #out_handle.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n' % (CHR, gtf_type, 'transcript', START - seqlen, END, '.', STRAND, '.', 'gene_id "g' + str(i) + '"; transcript_id "t' + str(i) + '"'))
                #out_handle.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n' % (CHR, gtf_type, 'exon', START - seqlen, END, '.', STRAND, '.', 'gene_id "g' + str(i) + '"; transcript_id "t' + str(i) + '"; exon_id "e' + str(i) + '"'))
                inside.write('%s\t%i\t%i\t%s\t%i\t%s\n' % (CHR, START -seqlen, START, 't' + str(i), 0, STRAND))
                outside.write('%s\t%i\t%i\t%s\t%i\t%s\n' % (CHR, START + 1, END + seqlen -1, 't' + str(i), 0, STRAND))
    i = pybedtools.BedTool(inside_tmp)
    o = pybedtools.BedTool(outside_tmp)
    inside_sorted = o_dir + "/" + k + "_" + gtf_type + ".inside.sorted.bed"
    outside_sorted = o_dir + "/" + k + "_" + gtf_type + ".outside.sorted.bed"
    i.sort(g=bam_header_file, output=inside_sorted)
    o.sort(g=bam_header_file, output=outside_sorted)
    os.system("rm {i} {o}".format(i=inside_tmp , o=outside_tmp))
    return(inside_sorted, outside_sorted)


def get_ratio_TSS(inside_bed, outside_bed, o_dir, input_bam, k, gtf_type, bam_header_file, verbose = True): 
## the idea would be to first calculate the average coverage per sample for in and out beds. Calculate each ratio
## get the maximum the ratios across replicates and return it as a dictionary
    if verbose:
        print("parse gtf_type %s on chromosome %s on BAM file %s" % (gtf_type, k, input_bam))
    out_TSS_tmp = o_dir + "/" + k + "_" + gtf_type + ".ratio_TSS.tmp"    
    out_TSS_file = o_dir + "/" + k + "_" + gtf_type + ".ratio_TSS.tsv"
    in_bed = pybedtools.BedTool(inside_bed)
    out_bed = pybedtools.BedTool(outside_bed)
    in_cov = in_bed.coverage(input_bam, sorted=True, g=bam_header_file)
    out_cov = out_bed.coverage(input_bam, sorted=True, g=bam_header_file)
    inside_df = pandas.DataFrame(columns=['id','inside'])
    for entry in in_cov:
        new_entry = pandas.DataFrame({'id' : [entry.name] , 'inside' : [float(entry[6])]})
        inside_df = pandas.concat([inside_df,new_entry], ignore_index=True)
    outside_df = pandas.DataFrame(columns=['id','outside'])
    for entry in out_cov:
        new_entry = pandas.DataFrame({'id' : [entry.name] , 'outside' : [float(entry[6])]})
        outside_df = pandas.concat([outside_df, new_entry], ignore_index=True)
    merged = pandas.merge(inside_df, outside_df, on="id")
    merged['ratio_TSS'] = (merged['inside']+0.01)/(merged['outside']+0.01)
    merged['ratio_TSS'] = pandas.to_numeric(merged['ratio_TSS'])
    #os.system('rm {i} {o}'.format(i=inside_bed, o=outside_bed))
    merged.to_csv(out_TSS_tmp, sep = "\t", header=False)
    os.system("paste {i} {o} {t} > {m}".format(i=inside_bed, o=outside_bed, t=out_TSS_tmp, m=out_TSS_file))
    os.system('rm {i} {o} {t}'.format(i=inside_bed, o=outside_bed, t=out_TSS_tmp))


def get_random_tss(input_bam, output_dir, reference, reference_gtf, min_size, random_tss_per_chr, mask, gtf_feature, not_parse_reference, write_reference_Ns, verbose):
    # create output directory
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # get mask from BED file
    if mask is not None:
        if verbose:
            print("get mask from %s" % mask)
        mask = pr.read_bed(mask)
        mask.columns = ['Chromosome', 'Start', 'End', 'Strand']
        mask.Strand = mask.Strand.astype("category")
    else:
        mask = pr.PyRanges()
    # get Ns from reference
    if not_parse_reference:
        Ns = get_Ns(reference, verbose)
        if write_reference_Ns:
            Ns.to_csv(path = reference + ".Ns.tsv", sep = "\t")
    else:
        Ns = pr.PyRanges()
    # join mask and Ns
    if verbose:
        print("join mask and Ns")
    mask = mask.join(Ns)
    # get exon chr, start, end position as ranges for possible tss sites
    if verbose:
        print("read GTF %s" % reference_gtf)
    gtf = pr.read_gtf(reference_gtf)
    gtf_tss = gtf.features.tss().merge()
    gtf_tes = gtf.features.tes().merge()
    gtf_introns = gtf.features.introns().merge()
    gtf_subset = gtf[gtf.Feature == gtf_feature]
    gtf_subset_merge = gtf_subset.merge()
    # subtract mask
    gtf_tss = gtf_tss.subtract(mask)
    gtf_tes = gtf_tes.subtract(mask)
    gtf_introns = gtf_introns.subtract(mask)
    gtf_subset_merge = gtf_subset_merge.subtract(mask)
    if verbose:
        print("read BAM header %s" % input_bam)
    bam_header_file = get_bam_header2(input_bam, output_dir)
    bam_header = read_bam_header(bam_header_file, min_size)
    # for each chromosome create bed file
    for k in list(bam_header.keys()):
        if verbose:
            print("create random bed for %s" % k)
        gtf_tss_k = gtf_tss[k]
        gtf_tes_k = gtf_tes[k]
        gtf_introns_k = gtf_introns[k]
        gtf_subset_merge_k = gtf_subset_merge[k]
        gtf_tss_k_pos = get_sample_position(gtf_tss_k, random_tss_per_chr)
        gtf_tes_k_pos = get_sample_position(gtf_tes_k, random_tss_per_chr)
        gtf_introns_k_pos = get_sample_position(gtf_introns_k, random_tss_per_chr)
        gtf_subset_merge_k_pos = get_sample_position(gtf_subset_merge_k, random_tss_per_chr)
        tss_inside_bed, tss_outside_bed = create_bed(k, "random_tss", gtf_tss_k_pos, output_dir, bam_header_file)
        tes_inside_bed, tes_outside_bed = create_bed(k, "random_tes", gtf_tes_k_pos, output_dir, bam_header_file)
        introns_inside_bed, introns_outside_bed = create_bed(k, "random_introns", gtf_introns_k_pos, output_dir, bam_header_file)
        feature_inside_bed, feature_outside_bed = create_bed(k, "random_feature", gtf_subset_merge_k_pos, output_dir, bam_header_file)
        get_ratio_TSS(tss_inside_bed, tss_outside_bed, output_dir, input_bam, k, "random_tss", bam_header_file, verbose)
        get_ratio_TSS(tes_inside_bed, tes_outside_bed, output_dir, input_bam, k, "random_tes", bam_header_file, verbose)
        get_ratio_TSS(introns_inside_bed, introns_outside_bed, output_dir, input_bam, k, "random_introns", bam_header_file, verbose)
        get_ratio_TSS(feature_inside_bed, feature_outside_bed, output_dir, input_bam, k, "random_feature", bam_header_file, verbose)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("random ratio_TSS")
    parser.add_argument("input_bam", help="Input BAM")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("reference", help="Reference FASTA")
    parser.add_argument("reference_gtf", help="Reference GTF")
    parser.add_argument("--min_size", type=int, default=100000, help="Minimum chr/scaffold/contig size in bp (default: 100000)")
    parser.add_argument("--random_tss_per_chr", type=int, default=1000, help="Number of random TSS per chr/scaffold/contig (default: 1000)")
    parser.add_argument("--mask", help="Specify sites to be masked in BED format (CHR, START, END, STRAND)")
    parser.add_argument("--gtf_feature", default="exon", help="Specify GTF feature to be used as possible tss sites (default: exon)")
    parser.add_argument("--not_parse_reference", help="Specify if reference FASTA file parsing for Ns should be turned off (default: True)", action="store_false")
    parser.add_argument("--write_reference_Ns", help="Specify if reference FASTA Ns position should be written to a bed file (default: False)", action="store_true")
    parser.add_argument("--verbose", help="increase output verbosity (default: False)", action="store_true")
    args = parser.parse_args()
    print(args)
    get_random_tss(args.input_bam, args.output_dir, args.reference, args.reference_gtf, args.min_size, args.random_tss_per_chr, args.mask, args.gtf_feature, args.not_parse_reference, args.write_reference_Ns, args.verbose)

