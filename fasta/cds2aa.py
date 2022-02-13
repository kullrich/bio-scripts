#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian K Ullrich
date: February 2022
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2022 Kristian Ullrich

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
from Bio import SeqIO
from Bio.Data import CodonTable


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator. Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


transtable = {'std': CodonTable.CodonTable(forward_table={
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
    'TRA': '*'},
    stop_codons=['TAA', 'TAG', 'TGA', ],
    start_codons=['TTG', 'CTG', 'ATG', ])}


def cds2aaRecord(record_iter, args):
    """Translates nucleotide to amino acids assuming that cds is in frame 0.

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            aa = SeqIO.SeqRecord(record.seq.translate(transtable[args.t]), name=record.name, id=record.name, description=record.name)
            yield aa


def cds2aaFasta(args, parser):
    record_iter = None
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    cds2aa_iter = cds2aaRecord(record_iter, args)
    if args.o is None:
        SeqIO.write(cds2aa_iter, sys.stdout, "fasta")
    else:
        count = SeqIO.write(cds2aa_iter, args.o, "fasta")
        print("translated %i sequences" % count)


def define_parser():
    parser = argparse.ArgumentParser(prog='cds2aa', usage='%(prog)s [options] [<arguments>...]',
                                     description='Translates nucleotide to amino acids assuming that CDS is in frame 0')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
    parser.add_argument('-t', help='transtable [default: std]', default='std')
    return parser


def main():
    # parser
    parser = define_parser()
    # get args
    args = parser.parse_args()
    # print args
    if args.o is None:
        sys.stderr.write(str(args))
    else:
        print(args)
    # run
    cds2aaFasta(args, parser)


if __name__ == '__main__':
    main()
