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


def int2idRecord(record_iter, args, numdict):
    """Converts number of sequence objects to name.

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    idcounter = 0
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            tmp_name = None
            tmp_id = None
            tmp_description = None
            if record.name not in numdict:
                tmp_name = record.name
                tmp_id = tmp_name
                tmp_description = 'not_in_table'
            if record.name in numdict:
                tmp_name = numdict[record.name]
                tmp_id = tmp_name
                tmp_description = ''
            redseq = SeqIO.SeqRecord(record.seq, name=tmp_name, id=tmp_id, description=tmp_description)
            yield redseq


def int2idFasta(args, parser):
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
    numdict = {}
    with open(args.t, 'rU') as inhandle:
        for lines in inhandle:
            line = lines.strip().split('\t')
            numdict[line[1]] = line[0]
    int2id_iter = int2idRecord(record_iter, args, numdict)
    if args.o is None:
        SeqIO.write(int2id_iter, sys.stdout, "fasta")
    else:
        count = SeqIO.write(int2id_iter, args.o, "fasta")
        print("converted %i sequence IDs" % count)


def define_parser():
    parser = argparse.ArgumentParser(prog='fastaINT2id', usage='%(prog)s [options] [<arguments>...]',
                                     description='Converts FASTA sequence number into IDs')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
    parser.add_argument('-t', help='table input file')
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
    int2idFasta(args, parser)


if __name__ == '__main__':
    main()
