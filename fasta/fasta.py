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
from Bio.Data import CodonTable as CT
from Bio import SeqIO
from Bio.SeqUtils import GC

transtable_std=CT.CodonTable(forward_table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', '---': 'X', '--A': 'X', '--C': 'X', '--G': 'X',
     '--T': 'X', '-A-': 'X', '-AA': 'X', '-AC': 'X', '-AG': 'X',
     '-AT': 'X', '-C-': 'X', '-CA': 'X', '-CC': 'X', '-CG': 'X',
     '-CT': 'X', '-G-': 'X', '-GA': 'X', '-GC': 'X', '-GG': 'X',
     '-GT': 'X', '-T-': 'X', '-TA': 'X', '-TC': 'X', '-TG': 'X',
     '-TT': 'X', 'A--': 'X', 'A-A': 'X', 'A-C': 'X', 'A-G': 'X',
     'A-T': 'X', 'AA-': 'X', 'AC-': 'X', 'AG-': 'X', 'AT-': 'X',
     'C--': 'X', 'C-A': 'X', 'C-C': 'X', 'C-G': 'X', 'C-T': 'X',
     'CA-': 'X', 'CC-': 'X', 'CG-': 'X', 'CT-': 'X', 'G--': 'X',
     'G-A': 'X', 'G-C': 'X', 'G-G': 'X', 'G-T': 'X', 'GA-': 'X',
     'GC-': 'X', 'GG-': 'X', 'GT-': 'X', 'T--': 'X', 'T-A': 'X',
     'T-C': 'X', 'T-G': 'X', 'T-T': 'X', 'TA-': 'X', 'TC-': 'X',
     'TG-': 'X', 'TT-': 'X', 'NNN': 'X', 'GCN': 'A', 'CGN': 'R',
     'MGR': 'R', 'AAY': 'N', 'GAY': 'D', 'TGY': 'C', 'CAR': 'Q',
     'GAR': 'E', 'GGN': 'G', 'CAY': 'H', 'ATH': 'I', 'YTR': 'L',
     'CTN': 'L', 'AAR': 'K', 'TTY': 'F', 'CCN': 'P', 'TCN': 'S',
     'AGY': 'S', 'ACN': 'T', 'TAY': 'Y', 'GTN': 'V', 'TAR': '*',
     'TRA': '*', },
                   stop_codons = [ 'TAA', 'TAG', 'TGA', ],
                   start_codons = [ 'TTG', 'CTG', 'ATG', ]
)

def mainparser(parser):
    parser.add_argument('-l', help='show sub-script list')

def subparser(subparsers):
    # translate; parser
    parser_translate = subparsers.add_parser('translate', help='translate help')
    parser_translate.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_translate.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_translate.add_argument('-i', help='input file')
    parser_translate.add_argument('-o', help='output file')
    parser_translate.add_argument('-t', default='std', help='transtable [default: std]')
    parser_translate.set_defaults(func=translate)
    # reverse; parser
    parser_reverse = subparsers.add_parser('rev', help='reverse help')
    parser_reverse.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_reverse.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_reverse.add_argument('-i', help='input file')
    parser_reverse.add_argument('-o', help='output file')
    parser_reverse.set_defaults(func=reverse)
    # reverse complement; parser
    parser_reverse_complement = subparsers.add_parser('revcomp', help='reverse complement help')
    parser_reverse_complement.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_reverse_complement.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_reverse_complement.add_argument('-i', help='input file')
    parser_reverse_complement.add_argument('-o', help='output file')
    parser_reverse_complement.set_defaults(func=reverse_complement)
    # complement; parser
    parser_complement = subparsers.add_parser('comp', help='complement help')
    parser_complement.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_complement.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_complement.add_argument('-i', help='input file')
    parser_complement.add_argument('-o', help='output file')
    parser_complement.set_defaults(func=complement)
    # translate nuc26orf; parser
    parser_nuc26orf = subparsers.add_parser('nuc26orf', help='nuc26orf help')
    parser_nuc26orf.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_nuc26orf.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_nuc26orf.add_argument('-i', help='input file')
    parser_nuc26orf.add_argument('-o', help='output file')
    parser_nuc26orf.add_argument('-t', default='std', help='transtable [default: std]')
    parser_nuc26orf.set_defaults(func=nuc26orf)
    # N50 stats; parser
    parser_n50 = subparsers.add_parser('n50', help='n50 stats help')
    parser_n50.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_n50.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_n50.add_argument('-i', help='input file')
    parser_n50.add_argument('-o', help='output file')
    parser_n50.set_defaults(func=n50stats)
    # get length; parser
    parser_getlength = subparsers.add_parser('getlen', help='get length help')
    parser_getlength.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_getlength.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_getlength.add_argument('-i', help='input file')
    parser_getlength.add_argument('-o', help='output file')
    parser_getlength.set_defaults(func=getlen)
    # reduce to min or max length; parser
    parser_redlength = subparsers.add_parser('redlen', help='reduce length help')
    parser_redlength.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_redlength.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_redlength.add_argument('-i', help='input file')
    parser_redlength.add_argument('-o', help='output file')
    parser_redlength.add_argument('-t', default='min', choices=['min','max'], help='specify type [default: min] or [max]')
    parser_redlength.add_argument('-m', default=0, type=int, help='specify length [default: 0]')
    parser_redlength.set_defaults(func=redlen)
    # GC content; parser
    parser_gc = subparsers.add_parser('gc', help='GC content help')
    parser_gc.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_gc.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_gc.add_argument('-i', help='input file')
    parser_gc.add_argument('-o', help='output file')
    parser_gc.add_argument('-t', default='gc', choices=['gc','gc1','gc2','gc3'], help='specify type [default: gc] or [gc1] or [gc2] or [gc3]')
    parser_gc.set_defaults(func=gccontent)
    # reduce names; parser
    parser_redname = subparsers.add_parser('redname', help='reduce IDs help')
    parser_redname.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_redname.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_redname.add_argument('-i', help='input file')
    parser_redname.add_argument('-o', help='output file')
    parser_redname.add_argument('-s', default=' ', help='specify split [default: " "]')
    parser_redname.add_argument('-k', default=0, help='specify keep field [default: 0]')
    parser_redname.set_defaults(func=redname)
    # convert names from strings to numbers, storing name number assignment into a table; parser
    parser_id2num = subparsers.add_parser('id2num', help='id2num help')
    parser_id2num.add_argument('-si', type=bool, default=True, help='specify if input from STDIN')
    parser_id2num.add_argument('-so', type=bool, default=True, help='specify if output to STDOUT')
    parser_id2num.add_argument('-i', help='input file')
    parser_id2num.add_argument('-o', help='output file')
    parser_id2num.set_defaults(func=id2num)
    # convert names from numbers to strings, using a name number assignment table; parser

    # split sequences into short sequences; parser

    # extract longest component trinity; parser

    # extract longest component velvet/oases; parser

def translate(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    if args.t=='std':
        cds2aa_fasta = cds2aa(original_fasta, transtable_std)
    if args.so:
        count = SeqIO.write(cds2aa_fasta, sys.stdout, "fasta")
    if not args.so:
        count = SeqIO.write(cds2aa_fasta, args.o, "fasta")
        print "translated %i sequences" % count

def cds2aa(records, transtable):
    """Translates nucleotide to amino acids assuming that cds is in frame 0.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        aa = SeqIO.SeqRecord(record.seq.translate(transtable),name=record.name,id=record.name,description=record.name)
        yield aa

def reverse(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    rev_fasta = reversefasta(original_fasta)
    if args.so:
        count = SeqIO.write(rev_fasta, sys.stdout, "fasta")
    if not args.so:
        count = SeqIO.write(rev_fasta, args.o, "fasta")
        print "reverted %i sequences" % count

def reversefasta(records):
    """Reverts sequence objects.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        revseq = SeqIO.SeqRecord(record.seq[::-1],name=record.name,id=record.name,description=record.name)
        yield revseq

def reverse_complement(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    revcomp_fasta = reversecomplementfasta(original_fasta)
    if args.so:
        count = SeqIO.write(revcomp_fasta, sys.stdout, "fasta")
    if not args.so:
        count = SeqIO.write(revcomp_fasta, args.o, "fasta")
        print "reverse complemented %i sequences" % count

def reversecomplementfasta(records):
    """Reverse complements sequence objects.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        revcompseq = SeqIO.SeqRecord(record.seq.reverse_complement(),name=record.name,id=record.name,description=record.name)
        yield revcompseq

def complement(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    comp_fasta = complementfasta(original_fasta)
    if args.so:
        count = SeqIO.write(comp_fasta, sys.stdout, "fasta")
    if not args.so:
        count = SeqIO.write(comp_fasta, args.o, "fasta")
        print "reverse complemented %i sequences" % count

def complementfasta(records):
    """Complements sequence objects.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        revcompseq = SeqIO.SeqRecord(record.seq.complement(),name=record.name,id=record.name,description=record.name)
        yield revcompseq

def gccontent(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    gcdict = calc_gc(original_fasta, args)
    gcmean = numpy.mean(gcdict.values())
    gcsd = numpy.std(gcdict.values(),ddof=1)
    if args.so:
        if args.t=='gc':
            print('#GC')
            print('mean GC percent: %s' % (gcmean))
            print('sd GC percent: %s' % (gcsd))
        if args.t=='gc1':
            print('#GC1')
            print('mean GC1 percent: %s' % (gcmean))
            print('sd GC1 percent: %s' % (gcsd))
        if args.t=='gc2':
            print('#GC2')
            print('mean GC2 percent: %s' % (gcmean))
            print('sd GC2 percent: %s' % (gcsd))
        if args.t=='gc3':
            print('#GC3')
            print('mean GC3 percent: %s' % (gcmean))
            print('sd GC3 percent: %s' % (gcsd))
        for k in sorted(gcdict.keys()):
            print('%s\t%s' % (k,gcdict[k]))
    if not args.so:
        with open(args.o,'w') as outhandle:
            if args.t=='gc':
                outhandle.write('#GC\n')
            if args.t=='gc1':
                outhandle.write('#GC1\n')
            if args.t=='gc2':
                outhandle.write('#GC2\n')
            if args.t=='gc3':
                outhandle.write('#GC3\n')
            for k in sorted(gcdict.keys()):
                outhandle.write('%s\t%f\n' % (k,gcdict[k]))
        if args.t=='gc':
            print('#GC')
            print('mean GC percent: %s' % (gcmean))
            print('sd GC percent: %s' % (gcsd))
        if args.t=='gc1':
            print('#GC1')
            print('mean GC1 percent: %s' % (gcmean))
            print('sd GC1 percent: %s' % (gcsd))
        if args.t=='gc2':
            print('#GC2')
            print('mean GC2 percent: %s' % (gcmean))
            print('sd GC2 percent: %s' % (gcsd))
        if args.t=='gc3':
            print('#GC3')
            print('mean GC3 percent: %s' % (gcmean))
            print('sd GC3 percent: %s' % (gcsd))

def calc_gc(records, args):
    gcdict = {}
    for record in records:
        tmp_id = record.id.split()[0]
        if args.t=='gc':
            tmp_gc = GC(record.seq)
        if args.t=='gc1':
            tmp_gc = GC(record.seq[::3])
        if args.t=='gc2':
            tmp_gc = GC(record.seq[1::3])
        if args.t=='gc3':
            tmp_gc = GC(record.seq[2::3])
        if tmp_id in gcdict:
            print('duplicated id: %s; skip GC calculation' % (tmp_id))
        if tmp_id not in gcdict:
            gcdict[tmp_id] = tmp_gc
    return(gcdict)

def redname(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    red_fasta = redids(original_fasta, args)
    if args.so:
        count = SeqIO.write(red_fasta, sys.stdout, "fasta")
    if not args.so:
        count = SeqIO.write(red_fasta, args.o, "fasta")
        print "reduced ids from %i sequences" % count

def redids(records, args):
    """Reduces names of sequence objects.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        redseq = SeqIO.SeqRecord(record.seq,name=record.description.split(args.s)[int(args.k)],id=record.description.split(args.s)[int(args.k)],description=record.description.split(args.s)[int(args.k)])
        yield redseq

def redlen(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    red_fasta = red_len(original_fasta, args)
    if args.so:
        count = SeqIO.write(red_fasta, sys.stdout, "fasta")
    if not args.so:
        count = SeqIO.write(red_fasta, args.o, "fasta")
        print "kept %i sequences" % count

def red_len(records, args):
    """Reduces sequence objects based on min/max length.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        if args.t=='min':
            if len(record)>=args.m:
                yield record
        if args.t=='max':
            if len(record)<=args.m:
                yield record

def getlen(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    lendict = get_len(original_fasta, args)
    lenmean = numpy.mean(lendict.values())
    lensd = numpy.std(lendict.values(),ddof=1)
    if args.so:
        print('#LEN')
        print('mean length: %s' % (lenmean))
        print('sd length: %s' % (lensd))
        for k in sorted(lendict.keys()):
            print('%s\t%s' % (k,lendict[k]))
    if not args.so:
        with open(args.o,'w') as outhandle:
            outhandle.write('#LEN\n')
            for k in sorted(lendict.keys()):
                outhandle.write('%s\t%f\n' % (k,lendict[k]))
            print('#LEN')
            print('mean length: %s' % (lenmean))
            print('sd length: %s' % (lensd))

def get_len(records, args):
    lendict = {}
    for record in records:
        tmp_id = record.id.split()[0]
        tmp_len = len(record)
        if tmp_id in lendict:
            print('duplicated id: %s; skip length calculation' % (tmp_id))
        if tmp_id not in lendict:
            lendict[tmp_id] = tmp_len
    return(lendict)

def n50stats(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    gc_len_summary = calc_n50stats(original_fasta, args)
    if args.so:
        print('#N50stats')
        print('%s\t%s' % ('n',gc_len_summary[0]))
        print('%s\t%s' % ('GCmean',gc_len_summary[1]))
        print('%s\t%s' % ('GCsd',gc_len_summary[2]))
        print('%s\t%s' % ('Sum',gc_len_summary[3]))
        print('%s\t%s' % ('Min',gc_len_summary[4]))
        print('%s\t%s' % ('Max',gc_len_summary[5]))
        print('%s\t%s' % ('Mean',gc_len_summary[6]))
        print('%s\t%s' % ('Median',gc_len_summary[7]))
        print('%s\t%s' % ('N5',gc_len_summary[8]))
        print('%s\t%s' % ('N25',gc_len_summary[9]))
        print('%s\t%s' % ('N50',gc_len_summary[10]))
        print('%s\t%s' % ('N75',gc_len_summary[12]))
        print('%s\t%s' % ('N95',gc_len_summary[12]))
    if not args.so:
        with open(args.o,'w') as outhandle:
            outhandle.write('#N50stats\n')
            outhandle.write('%s\t%s\n' % ('n',gc_len_summary[0]))
            outhandle.write('%s\t%s\n' % ('GCmean',gc_len_summary[1]))
            outhandle.write('%s\t%s\n' % ('GCsd',gc_len_summary[2]))
            outhandle.write('%s\t%s\n' % ('Sum',gc_len_summary[3]))
            outhandle.write('%s\t%s\n' % ('Min',gc_len_summary[4]))
            outhandle.write('%s\t%s\n' % ('Max',gc_len_summary[5]))
            outhandle.write('%s\t%s\n' % ('Mean',gc_len_summary[6]))
            outhandle.write('%s\t%s\n' % ('Median',gc_len_summary[7]))
            outhandle.write('%s\t%s\n' % ('N5',gc_len_summary[8]))
            outhandle.write('%s\t%s\n' % ('N25',gc_len_summary[9]))
            outhandle.write('%s\t%s\n' % ('N50',gc_len_summary[10]))
            outhandle.write('%s\t%s\n' % ('N75',gc_len_summary[12]))
            outhandle.write('%s\t%s\n' % ('N95',gc_len_summary[12]))

def calc_n50stats(records, args):
    gcdict = {}
    lendict = {}
    for record in records:
        tmp_id = record.id.split()[0]
        tmp_gc = GC(record.seq)
        tmp_len = len(record)
        if tmp_id in gcdict:
            print('duplicated id: %s; skip GC and length calculation' % (tmp_id))
        if tmp_id not in gcdict:
            gcdict[tmp_id] = tmp_gc
            if tmp_len in lendict:
                lendict[tmp_len]+=1
            if tmp_len not in lendict:
                lendict[tmp_len]=1
    gcmean = numpy.mean(gcdict.values())
    gcsd = numpy.std(gcdict.values(),ddof=1)
    lenu = numpy.repeat(lendict.keys(),lendict.values())
    lenr = numpy.repeat(lendict.keys(),[x[0]*x[1] for x in lendict.items()])
    minlen = numpy.min(lendict.keys())
    maxlen = numpy.max(lendict.keys())
    meanlen = numpy.mean(lenu)
    medianlen = numpy.median(lenu)
    totallen = len(lenr)
    number_contigs = len(lenu)
    N5 = numpy.percentile(lenr,95)
    N25 = numpy.percentile(lenr,75)
    N50 = numpy.percentile(lenr,50)
    N75 = numpy.percentile(lenr,25)
    N95 = numpy.percentile(lenr,5)
    gc_len_summary = [number_contigs,gcmean,gcsd,totallen,minlen,maxlen,meanlen,medianlen,N5,N25,N50,N75,N95]
    return(gc_len_summary)

def nuc26orf(args, parser):
    if args.si==True and args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')        
    if args.i is not None:
        args.si = False
    if args.si==False and args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    if args.o is not None:
        args.so = False
    if args.so==False and args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    print(args)
    if args.si:
        original_fasta = SeqIO.parse(sys.stdin, "fasta")
    if not args.si:
        original_fasta = SeqIO.parse(args.i, "fasta")
    if args.t=='std':
        nuc26orf_fasta = calc_nuc26orf(original_fasta, transtable_std)
    if args.so:
        count = SeqIO.write(nuc26orf_fasta, sys.stdout, "fasta")
    if not args.so:
        count = SeqIO.write(nuc26orf_fasta, args.o, "fasta")
        print "translated %i sequences" % count

def calc_nuc26orf(records, transtable):
    """Extracting of all 6 open reading frame(ORFs).

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        for strand, nuc in [("+", record.seq), ("-", record.seq.reverse_complement())]:
            for frame in range(3):
                length = 3 * ((len(record)-frame) // 3) #Multiple of three
                pro = SeqIO.SeqRecord(nuc[frame:frame+length].translate(transtable),name=record.name+"_"+str(frame)+str(strand),id=record.name+"_"+str(frame)+str(strand),description=record.name+"_"+str(frame)+str(strand))
                yield pro

def main():
    # top-level parser
    parser = argparse.ArgumentParser(prog='fasta', usage='%(prog)s <sub-script> [options] [<arguments>...]', description='scripts for FASTA file handling')
    mainparser(parser)
    subparsers = parser.add_subparsers(title='sub-scripts', description='valid sub-scripts', help='sub-scripts help', dest='cmd')
    # sub-level parser
    subparser(subparsers)
    args = parser.parse_args()
    args.func(args, parser)

if __name__ == '__main__':
    main()
