#!/usr/bin/python
# -*- coding: UTF-8 -*-

'''
Author: Krisian Ullrich
date: February 2020
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2017 Kristian Ullrich

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

This software relies on the following python packages:
sys
os
gzip
argparse
numpy https://pypi.python.org/pypi/numpy

This software relies on the following external-tools:
iqtree - http://www.iqtree.org/
'''


import sys
import os
import gzip
import argparse
import numpy as np


#main
def main():
    # top-level parser
    parser = argparse.ArgumentParser(prog='phasedvcf2pomo', usage='%(prog)s [options] [<arguments>...]',
                                     description='script to convert phased vcf (only GT field) into iqtree PoMo countfile')
    parser.add_argument('-i', help='specify phased vcf')
    parser.add_argument('-o', help='specify countfile')
    parser.add_argument('-p', help='specify populations as e.g. "-p A:1,2 -p B:3,4"', action='append')
    parser.add_argument('-w', type=int, help='specify window size [default: 50]')
    # get args
    args = parser.parse_args()
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify phased vcf')
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify countfile')
    if args.p is None:
        parser.print_help()
        sys.exit('\nPlease specify populations as e.g. "-p A:1,2 -p B:3,4"')
    print args
    INPUTFILE = args.i
    COUNTFILE = args.o
    POPNAMES = [x.split(':')[0] for x in args.p]
    POPPOS = [[int(y)+8 for y in x.split(':')[1].split(',')] for x in args.p]
    #default window size
    WINDOWSIZE = 50
    if args.w is not None:
        WINDOWSIZE = args.w
    if INPUTFILE.endswith('gz'):
        INPUT = gzip.open(INPUTFILE)
    if INPUTFILE.endswith('vcf'):
        INPUT = open(INPUTFILE)
    curidx = 1
    curpos = 0
    outhandle = open(COUNTFILE+'.'+str(curidx)+'.cf','w')
    outhandle.write('COUNTSFILE\tNPOP\t')
    outhandle.write(str(len(POPNAMES)))
    outhandle.write('\tNSITES\t')
    outhandle.write(str(WINDOWSIZE))
    outhandle.write('\n')
    outhandle.write('CHROM\tPOS')
    for popidx, popname in enumerate(POPNAMES):
        outhandle.write('\t')
        outhandle.write(popname)
    outhandle.write('\n')
    for line in INPUT:
        popdict = {}
        if line[0] != '#':
            curpos += 1
            if curpos > WINDOWSIZE:
                outhandle.close()
                curidx += 1
                outhandle = open(COUNTFILE+'.'+str(curidx)+'.cf','w')
                outhandle.write('COUNTSFILE\tNPOP\t')
                outhandle.write(str(len(POPNAMES)))
                outhandle.write('\tNSITES\t')
                outhandle.write(str(WINDOWSIZE))
                outhandle.write('\n')
                outhandle.write('CHROM\tPOS')
                for popidx, popname in enumerate(POPNAMES):
                    outhandle.write('\t')
                    outhandle.write(popname)
                outhandle.write('\n')
                curpos = 1
            l = line.strip().split('\t')
            CHROM = l[0]
            POS = int(l[1])
            l_ref = l[3]
            l_alt = l[4]
            outhandle.write(CHROM)
            outhandle.write('\t')
            outhandle.write(str(POS))
            for popidx, popname in enumerate(POPNAMES):
                popdict[popname] = {'A':0, 'C':0, 'T':0, 'G':0}
                popGT = [j for i,j in enumerate(l) if i in POPPOS[popidx]]
                popCounts = [[l_ref,l_alt][int(x)] for x in np.array([x.split('|') for x in popGT]).flatten()]
                popdict[popname]['A'] = "".join(popCounts).count('A')
                popdict[popname]['C'] = "".join(popCounts).count('C')
                popdict[popname]['T'] = "".join(popCounts).count('T')
                popdict[popname]['G'] = "".join(popCounts).count('G')
                popout = '%i,%i,%i,%i' % (popdict[popname]['A'],popdict[popname]['C'],popdict[popname]['T'],popdict[popname]['G'])
                outhandle.write('\t')
                outhandle.write(popout)
            outhandle.write('\n')
    outhandle.close()


if __name__ == '__main__':
    main()
