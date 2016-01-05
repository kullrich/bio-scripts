#!/usr/bin/python
# -*- coding: UTF-8 -*-

'''
Author: Krisian Ullrich
date: Dezember 2015
email: kristian.ullrich@biologie.uni-marburg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2015 Kristian Ullrich

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
'''

import sys
import argparse
import numpy
import gzip
from Bio.Data import CodonTable as CT
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqUtils import GC
from itertools import izip

def mainparser(parser):
    return()

def subparser(subparsers):
    # trim static; parser
    parser_static = subparsers.add_parser('static', help='static help')
    parser_static.add_argument('-i', help='input file; single-end mode [in] and paired-end mode [in1,in2]')
    parser_static.add_argument('-o', help='output file; single-end mode [out] and paired-end mode [out1p,out1s,out2p,out2s]')
    parser_static.add_argument('-m', default=50, type=int, help='specify min length [default: 50]')
    parser_static.add_argument('-f', default=0, type=int, help='specify five-prime trim [default: 0]')
    parser_static.add_argument('-t', default=0, type=int, help='specify three-prime trim [default: 0]')
    parser_static.add_argument('-a', default=20.0, type=float, help='specify average quality [default: 20.0]')
    parser_static.add_argument('-p', default=False, type=bool, help='specify if paired-end data [default: False]')
    parser_static.add_argument('-g', default=False, type=bool, help='specify if data is gz [default: False]')
    parser_static.add_argument('-d', default=False, type=bool, help='specify if direction should be added to defline [default: False]')
    parser_static.add_argument('-r', default=False, type=bool, help='specify if output should be reverse complement of read [default: False]')
    parser_static.set_defaults(func=trimstatic)
    # trim dynamic; parser
    parser_dynamic = subparsers.add_parser('dynamic', help='dynamic help')
    parser_dynamic.add_argument('-i', help='input file; single-end mode [in] and paired-end mode [in1,in2]')
    parser_dynamic.add_argument('-o', help='output file; single-end mode [out] and paired-end mode [out1p,out1s,out2p,out2s]')
    parser_dynamic.add_argument('-m', default=50, type=int, help='specify min length [default: 50]')
    parser_dynamic.add_argument('-q', default=20, type=int, help='specify quality cutoff [default: 20]')
    parser_dynamic.add_argument('-a', default=20.0, type=float, help='specify average quality [default: 20.0]')
    parser_dynamic.add_argument('-p', default=False, type=bool, help='specify if paired-end data [default: False]')
    parser_dynamic.add_argument('-g', default=False, type=bool, help='specify if data is gz [default: False]')
    parser_dynamic.add_argument('-d', default=False, type=bool, help='specify if direction should be added to defline [default: False]')
    parser_dynamic.add_argument('-r', default=False, type=bool, help='specify if output should be reverse complement of read [default: False]')
    parser_dynamic.set_defaults(func=trimdynamic)
    # trim mott; parser
    parser_mott = subparsers.add_parser('mott', help='mott help')
    parser_mott.add_argument('-i', help='input file; single-end mode [in] and paired-end mode [in1,in2]')
    parser_mott.add_argument('-o', help='output file; single-end mode [out] and paired-end mode [out1p,out1s,out2p,out2s]')
    parser_mott.add_argument('-m', default=50, type=int, help='specify min length [default: 50]')
    parser_mott.add_argument('-q', default=20, type=int, help='specify quality cutoff [default: 20]')
    parser_mott.add_argument('-a', default=20.0, type=float, help='specify average quality [default: 20.0]')
    parser_mott.add_argument('-p', default=False, type=bool, help='specify if paired-end data [default: False]')
    parser_mott.add_argument('-g', default=False, type=bool, help='specify if data is gz [default: False]')
    parser_mott.add_argument('-d', default=False, type=bool, help='specify if direction should be added to defline [default: False]')
    parser_mott.add_argument('-r', default=False, type=bool, help='specify if output should be reverse complement of read [default: False]')
    parser_mott.set_defaults(func=trimmott)
    # nuc2Nqual; parser
    parser_nuc2Nqual = subparsers.add_parser('nuc2Nqual', help='nuc2Nqual help')
    parser_nuc2Nqual.add_argument('-i', help='input file; single-end mode [in] and paired-end mode [in1,in2]')
    parser_nuc2Nqual.add_argument('-o', help='output file; single-end mode [out] and paired-end mode [out1,out2]')
    parser_nuc2Nqual.add_argument('-q', default=20, type=int, help='specify quality cutoff [default: 20]')
    parser_nuc2Nqual.add_argument('-p', default=False, type=bool, help='specify if paired-end data [default: False]')
    parser_nuc2Nqual.add_argument('-g', default=False, type=bool, help='specify if data is gz [default: False]')
    parser_nuc2Nqual.set_defaults(func=nuc2Nqual)

def trimmott(args, parser):
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    if args.p:
        if len(args.i.split(','))!=2:
            sys.exit('\nPlease specify two comma seperated paired-end input files [in1,in2]')
        if len(args.o.split(','))!=4:
            sys.exit('\nPlease specify four comma seperated paired-end output files [out1p,out1s,out2p,out2s]')
    print(args)
    if not args.p:
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(args.i,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(args.i, "fastq")
        trimmed_fastq = trimmott_se(original1_fastq, args)
        write_trimmed(trimmed_fastq, args)
    if args.p:
        data1 = args.i.split(',')[0]
        data2 = args.i.split(',')[1]
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(data1,'rb'), "fastq")
            original2_fastq = SeqIO.parse(gzip.open(data2,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(data1, "fastq")
            original2_fastq = SeqIO.parse(data2, "fastq")
        trimmed_fastq = trimmott_pe(original1_fastq, original2_fastq, args)
        write_trimmed(trimmed_fastq, args)

def trimmott_se(records, args):
    """Mott trimming reads, trims according to quality cutoff, average quality and minimal read length to keep.
    
    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        record=mottRecord(record, args.q)
        if args.r:
            record=record.reverse_complement()
        if args.d:
            record.name=record.name+'/1'
            record.id=record.id+'/1'
            record.description=record.description+'/1'
        if len(record)>=args.m and numpy.mean(record.letter_annotations['phred_quality'])>=args.a:
            yield record

def trimmott_pe(records1, records2, args):
    """Mott trimming reads, trims according to quality cutoff, average quality and minimal read length to keep.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for rec1, rec2 in izip(records1, records2):
        rec1=mottRecord(rec1,args.q)
        rec2=mottRecord(rec2,args.q)
        if args.r:
            rec1=rec1.reverse_complement()
            rec2=rec2.reverse_complement()
        if args.d:
            rec1.name=rec1.name+'/1'
            rec1.id=rec1.id+'/1'
            rec1.description=rec1.description+'/1'
            rec2.name=rec2.name+'/2'
            rec2.id=rec2.id+'/2'
            rec2.description=rec2.description+'/2'
        y1 = False
        y2 = False
        if len(rec1)>=args.m and numpy.mean(rec1.letter_annotations['phred_quality'])>=args.a:
            y1 = True
        if len(rec2)>=args.m and numpy.mean(rec2.letter_annotations['phred_quality'])>=args.a:
            y2 = True
        if y1 and y2:
            yield rec1, None, rec2, None, 'pe'
        if y1 and not y2:
            yield None, rec1, None, None, 'se1'
        if not y1 and y2:
            yield None, None, None, rec2, 'se2'

def mottRecord(record, cutoff):
    mottScore=0
    index=0
    Scores=[x-cutoff for x in record[::-1].letter_annotations['phred_quality']]
    for i,mS in enumerate(Scores):
        mottScore+=mS
        if mottScore<=0:
            mottScore=0
            index=i+1
    return(record[::-1][index:][::-1])

def nuc2Nqual(args, parser):
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    if args.p:
        if len(args.i.split(','))!=2:
            sys.exit('\nPlease specify two comma seperated paired-end input files [in1,in2]')
        if len(args.o.split(','))!=2:
            sys.exit('\nPlease specify two comma seperated paired-end output files [out1,out2]')
    print(args)
    if not args.p:
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(args.i,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(args.i, "fastq")
        nuc2Nqual_fastq = nuc2Nqual_se(original1_fastq, args)
        write_nuc2Nqual(nuc2Nqual_fastq, args)
    if args.p:
        data1 = args.i.split(',')[0]
        data2 = args.i.split(',')[1]
        rcpe=0
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(data1,'rb'), "fastq")
            original2_fastq = SeqIO.parse(gzip.open(data2,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(data1, "fastq")
            original2_fastq = SeqIO.parse(data2, "fastq")
        nuc2Nqual_fastq = nuc2Nqual_pe(original1_fastq, original2_fastq, args)
        write_nuc2Nqual(nuc2Nqual_fastq, args)

def nuc2Nqual_se(records, args):
    """Exchanges nucleotides given a quality threshold to Ns.
    
    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        changepos = [i for i,x in enumerate(record.letter_annotations['phred_quality']) if x < args.q]
        outread = SeqRecord.SeqRecord(record.seq.tomutable())
        outread.id = record.id
        outread.name = record.name
        outread.description = record.description
        outread.letter_annotations = record.letter_annotations
        for i in changepos:
            outread.seq[i]='N'
        yield outread

def nuc2Nqual_pe(records1, records2, args):
    """Exchanges nucleotides given a quality threshold to Ns.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for rec1, rec2 in izip(records1, records2):
        changepos1 = [i for i,x in enumerate(rec1.letter_annotations['phred_quality']) if x < args.q]
        outread1 = SeqRecord.SeqRecord(rec1.seq.tomutable())
        outread1.id = rec1.id
        outread1.name = rec1.name
        outread1.description = rec1.description
        outread1.letter_annotations = rec1.letter_annotations
        for i in changepos1:
            outread1.seq[i]='N'
        changepos2 = [i for i,x in enumerate(rec2.letter_annotations['phred_quality']) if x < args.q]
        outread2 = SeqRecord.SeqRecord(rec2.seq.tomutable())
        outread2.id = rec2.id
        outread2.name = rec2.name
        outread2.description = rec2.description
        outread2.letter_annotations = rec2.letter_annotations
        for i in changepos2:
            outread2.seq[i]='N'
        yield outread1, outread2

def write_nuc2Nqual(nuc2Nqual_fastq, args):
    if not args.p:
        if args.g:
            count = SeqIO.write(nuc2Nqual_fastq, gzip.open(args.o,'wb'), "fastq")
            print('changed %i reads' % count)
        if not args.g:
            count = SeqIO.write(nuc2Nqual_fastq, args.o, "fastq")
            print('changed %i reads' % count)
    if args.p:
        data1pout = args.o.split(',')[0]
        data2pout = args.o.split(',')[1]
        rcpe=0
        if args.g:
            data1pout_fo = gzip.open(data1pout,'wb')
            data2pout_fo = gzip.open(data2pout,'wb')
            for r in nuc2Nqual_fastq:
                SeqIO.write(r[0], data1pout_fo, "fastq")
                SeqIO.write(r[1], data2pout_fo, "fastq")
                rcpe+=1
        if not args.g:
            data1pout_fo = open(data1pout,'wb')
            data2pout_fo = open(data2pout,'wb')
            for r in nuc2Nqual_fastq:
                SeqIO.write(r[0], data1pout_fo, "fastq")
                SeqIO.write(r[1], data2pout_fo, "fastq")
                rcpe+=1
        data1pout_fo.close()
        data2pout_fo.close()
        print('changed %i paired-end reads' % rcpe)

def trimdynamic(args, parser):
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    if args.p:
        if len(args.i.split(','))!=2:
            sys.exit('\nPlease specify two comma seperated paired-end input files [in1,in2]')
        if len(args.o.split(','))!=4:
            sys.exit('\nPlease specify four comma seperated paired-end output files [out1p,out1s,out2p,out2s]')
    print(args)
    if not args.p:
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(args.i,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(args.i, "fastq")
        trimmed_fastq = trimdynamic_se(original1_fastq, args)
        write_trimmed(trimmed_fastq, args)
    if args.p:
        data1 = args.i.split(',')[0]
        data2 = args.i.split(',')[1]
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(data1,'rb'), "fastq")
            original2_fastq = SeqIO.parse(gzip.open(data2,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(data1, "fastq")
            original2_fastq = SeqIO.parse(data2, "fastq")
        trimmed_fastq = trimdynamic_pe(original1_fastq, original2_fastq, args)
        write_trimmed(trimmed_fastq, args)

def trimdynamic_se(records, args):
    """Dynamic trimming reads, trims according to quality cutoff, average quality and minimal read length to keep.
    
    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        tmp_qual = [0 if x < args.q else 1 for x in record.letter_annotations['phred_quality']]
        tmp_qual.append(0)
        jumps = [i for i,x in enumerate(tmp_qual[:len(tmp_qual)-1]) if [tmp_qual[i],tmp_qual[i+1]]==[1,0]]
        if len(jumps)==0:
            cutpos = 0
        if len(jumps)!=0:
            cutpos = numpy.max(jumps)+1
        record=record[:cutpos]
        if args.r:
            record=record.reverse_complement()
        if args.d:
            record.name=record.name+'/1'
            record.id=record.id+'/1'
            record.description=record.description+'/1'
        if len(record)>=args.m and numpy.mean(record.letter_annotations['phred_quality'])>=args.a:
            yield record

def trimdynamic_pe(records1, records2, args):
    """Dynamic trimming reads, trims according to position, average quality and minimal read length to keep.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for rec1, rec2 in izip(records1, records2):
        tmp_qual1 = [0 if x < args.q else 1 for x in rec1.letter_annotations['phred_quality']]
        tmp_qual1.append(0)
        jumps1 = [i for i,x in enumerate(tmp_qual1[:len(tmp_qual1)-1]) if [tmp_qual1[i],tmp_qual1[i+1]]==[1,0]]
        if len(jumps1)==0:
            cutpos1 = 0
        if len(jumps1)!=0:
            cutpos1 = numpy.max(jumps1)+1
        rec1=rec1[:cutpos1]
        tmp_qual2 = [0 if x < args.q else 1 for x in rec2.letter_annotations['phred_quality']]
        tmp_qual2.append(0)
        jumps2 = [i for i,x in enumerate(tmp_qual2[:len(tmp_qual2)-1]) if [tmp_qual2[i],tmp_qual2[i+1]]==[1,0]]
        if len(jumps2)==0:
            cutpos2 = 0
        if len(jumps2)!=0:
            cutpos2 = numpy.max(jumps2)+1
        rec2=rec2[:cutpos2]
        if args.r:
            rec1=rec1.reverse_complement()
            rec2=rec2.reverse_complement()
        if args.d:
            rec1.name=rec1.name+'/1'
            rec1.id=rec1.id+'/1'
            rec1.description=rec1.description+'/1'
            rec2.name=rec2.name+'/2'
            rec2.id=rec2.id+'/2'
            rec2.description=rec2.description+'/2'
        y1 = False
        y2 = False
        if len(rec1)>=args.m and numpy.mean(rec1.letter_annotations['phred_quality'])>=args.a:
            y1 = True
        if len(rec2)>=args.m and numpy.mean(rec2.letter_annotations['phred_quality'])>=args.a:
            y2 = True
        if y1 and y2:
            yield rec1, None, rec2, None, 'pe'
        if y1 and not y2:
            yield None, rec1, None, None, 'se1'
        if not y1 and y2:
            yield None, None, None, rec2, 'se2'

def trimstatic(args, parser):
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    if args.p:
        if len(args.i.split(','))!=2:
            sys.exit('\nPlease specify two comma seperated paired-end input files [in1,in2]')
        if len(args.o.split(','))!=4:
            sys.exit('\nPlease specify four comma seperated paired-end output files [out1p,out1s,out2p,out2s]')
    print(args)
    if not args.p:
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(args.i,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(args.i, "fastq")
        trimmed_fastq = trimstatic_se(original1_fastq, args)
        write_trimmed(trimmed_fastq, args)
    if args.p:
        data1 = args.i.split(',')[0]
        data2 = args.i.split(',')[1]
        if args.g:
            original1_fastq = SeqIO.parse(gzip.open(data1,'rb'), "fastq")
            original2_fastq = SeqIO.parse(gzip.open(data2,'rb'), "fastq")
        if not args.g:
            original1_fastq = SeqIO.parse(data1, "fastq")
            original2_fastq = SeqIO.parse(data2, "fastq")
        trimmed_fastq = trimstatic_pe(original1_fastq, original2_fastq, args)
        write_trimmed(trimmed_fastq, args)

def trimstatic_se(records, args):
    """Static trimming reads, trims according to position, average quality and minimal read length to keep.
    
    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        record=record[args.f:len(record)-args.t]
        if args.r:
            record=record.reverse_complement()
        if args.d:
            record.name=record.name+'/1'
            record.id=record.id+'/1'
            record.description=record.description+'/1'
        if len(record)>=args.m and numpy.mean(record.letter_annotations['phred_quality'])>=args.a:
            yield record

def trimstatic_pe(records1, records2, args):
    """Static trimming reads, trims according to position, average quality and minimal read length to keep.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for rec1, rec2 in izip(records1, records2):
        rec1=rec1[args.f:len(rec1)-args.t]
        rec2=rec2[args.f:len(rec2)-args.t]
        if args.r:
            rec1=rec1.reverse_complement()
            rec2=rec2.reverse_complement()
        if args.d:
            rec1.name=rec1.name+'/1'
            rec1.id=rec1.id+'/1'
            rec1.description=rec1.description+'/1'
            rec2.name=rec2.name+'/2'
            rec2.id=rec2.id+'/2'
            rec2.description=rec2.description+'/2'
        y1 = False
        y2 = False
        if len(rec1)>=args.m and numpy.mean(rec1.letter_annotations['phred_quality'])>=args.a:
            y1 = True
        if len(rec2)>=args.m and numpy.mean(rec2.letter_annotations['phred_quality'])>=args.a:
            y2 = True
        if y1 and y2:
            yield rec1, None, rec2, None, 'pe'
        if y1 and not y2:
            yield None, rec1, None, None, 'se1'
        if not y1 and y2:
            yield None, None, None, rec2, 'se2'

def write_trimmed(trimmed_fastq, args):
    if not args.p:
        if args.g:
            count = SeqIO.write(trimmed_fastq, gzip.open(args.o,'wb'), "fastq")
            print('retained %i reads' % count)
        if not args.g:
            count = SeqIO.write(trimmed_fastq, args.o, "fastq")
            print('retained %i reads' % count)
    if args.p:
        data1pout = args.o.split(',')[0]
        data1sout = args.o.split(',')[1]
        data2pout = args.o.split(',')[2]
        data2sout = args.o.split(',')[3]
        rcpe=0
        rcse1=0
        rcse2=0
        if args.g:
            data1pout_fo = gzip.open(data1pout,'wb')
            data1sout_fo = gzip.open(data1sout,'wb')
            data2pout_fo = gzip.open(data2pout,'wb')
            data2sout_fo = gzip.open(data2sout,'wb')
            for r in trimmed_fastq:
                if r[0] is not None:
                    SeqIO.write(r[0], data1pout_fo, "fastq")
                if r[1] is not None:
                    SeqIO.write(r[1], data1sout_fo, "fastq")
                if r[2] is not None:
                    SeqIO.write(r[2], data2pout_fo, "fastq")
                if r[3] is not None:
                    SeqIO.write(r[3], data2sout_fo, "fastq")
                if r[4]=='pe':
                    rcpe+=1
                if r[4]=='se1':
                    rcse1+=1
                if r[4]=='se2':
                    rcse2+=1
        if not args.g:
            data1pout_fo = open(data1pout,'wb')
            data1sout_fo = open(data1sout,'wb')
            data2pout_fo = open(data2pout,'wb')
            data2sout_fo = open(data2sout,'wb')
            for r in trimmed_fastq:
                if r[0] is not None:
                    SeqIO.write(r[0], data1pout_fo, "fastq")
                if r[1] is not None:
                    SeqIO.write(r[1], data1sout_fo, "fastq")
                if r[2] is not None:    
                    SeqIO.write(r[2], data2pout_fo, "fastq")
                if r[3] is not None:
                    SeqIO.write(r[3], data2sout_fo, "fastq")
                if r[4]=='pe':
                    rcpe+=1
                if r[4]=='se1':
                    rcse1+=1
                if r[4]=='se2':
                    rcse2+=1
        data1pout_fo.close()
        data1sout_fo.close()
        data2pout_fo.close()
        data2sout_fo.close()
        print('retained %i paired-end reads' % rcpe)
        print('retained %i single-end fwd reads' % rcse1)
        print('retained %i single-end rev reads' % rcse2)

def main():
    # top-level parser
    parser = argparse.ArgumentParser(prog='fastq', usage='%(prog)s <sub-script> [options] [<arguments>...]', description='scripts for FASTQ file handling')
    mainparser(parser)
    subparsers = parser.add_subparsers(title='sub-scripts', description='valid sub-scripts', help='sub-scripts help', dest='cmd')
    # sub-level parser
    subparser(subparsers)
    args = parser.parse_args()
    args.func(args, parser)

if __name__ == '__main__':
    main()
