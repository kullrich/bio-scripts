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


def batch_iterator(iterator, batch_size):
    """Converts name of sequence objects to number.

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


def id2intRecord(record_iter, args, iddict):
    """Reduces names of sequence objects.

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    idcounter = 0
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            idcounter += 1
            iddict[record.name] = idcounter
            redseq = SeqIO.SeqRecord(record.seq, name=str(idcounter), id=str(idcounter), description=str(idcounter))
            yield redseq


def id2intFasta(args, parser):
    record_iter = None
    if args.t is None:
        parser.print_help()
        sys.exit('\nPlease specify table output file')
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    iddict = {}
    id2int_iter = id2intRecord(record_iter, args, iddict)
    if args.o is None:
        SeqIO.write(id2int_iter, sys.stdout, "fasta")
    else:
        count = SeqIO.write(id2int_iter, args.o, "fasta")
        print("converted %i sequence IDs" % count)
    with open(args.t, 'w') as outhandle:
        for k in sorted(iddict.keys()):
            outhandle.write('%s\t%s\n' % (k, iddict[k]))


def define_parser():
    parser = argparse.ArgumentParser(prog='fastaID2int', usage='%(prog)s [options] [<arguments>...]',
                                     description='Converts FASTA sequence IDs into number')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
    parser.add_argument('-t', help='table output file')
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
    id2intFasta(args, parser)


if __name__ == '__main__':
    main()
