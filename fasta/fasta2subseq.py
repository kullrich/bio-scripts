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
import argparse
import subprocess
import pysam
import pandas as pd


def subseq(individuals, chr, start, end):
    subseq_dict = {}
    for ind_index, ind_row in individuals.iterrows():
        subseq_dict[ind_row['ind']] = pysam.FastaFile(ind_row['filepath']).fetch(reference=chr, start=int(start), end=int(end))
    return subseq_dict


def fasta2subseq(args, parser):
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease provide file with population assignment')
    if args.chr is None:
        parser.print_help()
        sys.exit('\nPlease provide region:chr')
    if args.chr is None:
        parser.print_help()
        sys.exit('\nPlease provide region:start')
    if args.chr is None:
        parser.print_help()
        sys.exit('\nPlease provide region:end')
    popassignment = pd.read_csv(args.i, delimiter='\t')
    range_seq_dict = subseq(popassignment, args.chr, args.start, args.end)
    if args.o is None:
        for k,v in range_seq_dict.items():
            print('>'+k)
            print(v)
    else:
        with open(args.o, 'w') as outhandle:
            for k,v in range_seq_dict.items():
                    outhandle.write('>'+k+'\n')
                    outhandle.write(v+'\n')


def define_parser():
    parser = argparse.ArgumentParser(prog='fasta2subseq', usage='%(prog)s [options] [<arguments>...]',
                                     description='use samtools indexed FASTA sequences to obtain subsequences for a given region')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-chr', help='region:chr')
    parser.add_argument('-start', help='region:start', type=int)
    parser.add_argument('-end', help='region:end', type=int)
    return parser


def main():
    # parser
    parser = define_parser()
    # get args
    args = parser.parse_args()
    # print args
    #if args.o is None:
    #    sys.stderr.write(str(args))
    #else:
    #    print(args)
    # run
    fasta2subseq(args, parser)


if __name__ == '__main__':
    main()
