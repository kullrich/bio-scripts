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


def complementRecord(record_iter, args):
    """Complements sequence objects.

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            compseq = SeqIO.SeqRecord(record.reverse_complement().seq[::-1], name=record.name, id=record.name, description=record.name)
            yield compseq


def complementFasta(args, parser):
    record_iter = None
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    complement_iter = complementRecord(record_iter, args)
    if args.o is None:
        SeqIO.write(complement_iter, sys.stdout, "fasta")
    else:
        count = SeqIO.write(complement_iter, args.o, "fasta")
        print("complemented %i sequences" % count)


def define_parser():
    parser = argparse.ArgumentParser(prog='fasta2comp', usage='%(prog)s [options] [<arguments>...]',
                                     description='complements FASTA sequences')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
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
    complementFasta(args, parser)


if __name__ == '__main__':
    main()
