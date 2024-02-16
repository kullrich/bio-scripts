#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian K Ullrich
date: February 2024
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2024 Kristian Ullrich

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


import sys
import os
import argparse
import subprocess
import pysam
import tempfile
import numpy as np
import pandas as pd
from io import BytesIO


def split_by(x, by, getpos):
    return float(str(x).split(by)[getpos])


def as_dist(df):
    keep = np.triu(np.ones(df.shape), k=1).astype('bool').reshape(df.size)
    return df.stack(future_stack=True)[keep]


def subseq(stats_individuals, range_row):
    subseq_dict = {}
    for ind_index, ind_row in stats_individuals.iterrows():
        subseq_dict[ind_row['ind']] = pysam.FastaFile(ind_row['filepath']).fetch(reference=range_row[0], start=int(range_row[1]), end=int(range_row[2]))
    return subseq_dict


def run_hyde(infile, mapfile, outgroup, nind, nsites, ntaxa, threads, pvalue, prefix):
    hyderesults = subprocess.run(['run_hyde_mp.py', '-i', infile, '-m', mapfile, '-o', outgroup, '-n', nind, '-s', nsites, '-t', ntaxa, '-j', threads, '-p', pvalue, '--prefix', prefix], capture_output=True)
    return hyderesults


def run_hyde_individual(infile, mapfile, outgroup, triplets, nind, nsites, ntaxa prefix):
    hyderesults = subprocess.run(['individual_hyde.py', '-i', infile, '-m', mapfile, '-o', outgroup, '-tr', triplets, '-n', nind, '-s', nsites, '-t', ntaxa, '--prefix', prefix], capture_output=True)
    return hyderesults


def fasta2hyde(args, parser):
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease provide file with population assignment')
    if args.m is None:
        parser.print_help()
        sys.exit('\nPlease provide file with population assignment for hyde')
    if args.r is None:
        parser.print_help()
        sys.exit('\nPlease provide file with genomic ranges')
    if args.a and args.e is None:
        parser.print_help()
        sys.exit('\nPlease provide triplet file in case of individual testing')
    outfile = open(args.o, 'w')
    header = ['chr', 'start', 'end',
    'P1', 'Hybrid', 'P2', 'Zscore', 'Pvalue', 'Gamma',
    'AAAA', 'AAAB', 'AABA', 'AABB', 'AABC', 'ABAA', 'ABAB',
    'ABAC', 'ABBA', 'BAAA', 'ABBC', 'CABC', 'BACA', 'BCAA', 'ABCD']
    outfile.write('\t'.join(header)+'\n')
    popassignment = pd.read_csv(args.i, delimiter='\t')
    genomicranges = pd.read_csv(args.r, delimiter='\t', header=None)
    mapfile = args.m
    outgroup = args.g
    nind = str(len(popassignment))
    ntaxa = str(len(set(popassignment['pop'])))
    threads = str(args.t)
    pvalue = str(args.p)
    for range_index, range_row in genomicranges.iterrows():
        range_seq_dict = subseq(popassignment, range_row)
        tmpfile = tempfile.NamedTemporaryFile()
        nsites = str(range_row[2]-range_row[1])
        prefix = str(range_row[0])+'-'+str(range_row[1])+'-'+str(range_row[2])
        with open(tmpfile.name, 'w') as f:
            for k,v in range_seq_dict.items():
                f.write(k+'\t')
                f.write(v+'\n')
        run_hyde(tmpfile.name, mapfile, outgroup, nind, nsites, ntaxa, threads, pvalue, prefix)
        if args.a:
            triplets = args.e
            run_hyde_individual(tmpfile.name, mapfile, outgroup, triplets, nind, nsites, ntaxa, prefix)
        hyde_results = pd.read_csv(prefix+'-out.txt', delimiter='\t')
        for hyde_index, hyde_row in hyde_results.iterrows():
            outfile.write('\t'.join([str(x) for x in [range_row[0], range_row[1], range_row[2]]] + [str(x) for x in hyde_row])+'\n')
        os.remove(prefix+'-out.txt')
        os.remove(prefix+'-out-filtered.txt')
    outfile.close()


def define_parser():
    parser = argparse.ArgumentParser(prog='fasta2hyde', usage='%(prog)s [options] [<arguments>...]',
                                     description='hybridization stats from FASTA sequences')
    parser.add_argument('-i', help='file with population assignment and fasta file path (see example_pop.txt); use full path to avoid file finding')
    parser.add_argument('-m', help='file with population assignment for hyde (see example_map.txt); use full path to avoid file finding')
    parser.add_argument('-r', help='file with genomic ranges to be analysed in bed format (see example_region.txt)')
    parser.add_argument('-g', help='specify outgroup population')
    parser.add_argument('-a', help='specify if individuals should be tested', action='store_true')
    parser.add_argument('-e', help='specify triplet file in case of indvidual testing')
    parser.add_argument('-p', help='specify hyde pvalue [default: 0.05]', default='0.05')
    parser.add_argument('-o', help='output file [default: hyde_stats.txt]', default='hyde_stats.txt')
    parser.add_argument('-t', help='number of threads for hyde [default: 1]', default=1, type=int)
    return parser


def main():
    # parser
    parser = define_parser()
    # get args
    args = parser.parse_args()
    # print args
    sys.stderr.write(str(args))
    # run
    fasta2hyde(args, parser)


if __name__ == '__main__':
    main()
