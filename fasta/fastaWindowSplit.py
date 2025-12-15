#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian K Ullrich
date: Dece,ber 2025
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2025 Kristian Ullrich

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
import numpy as np
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


def sliding_window_steps_generator(record, w, j):
    '''Sliding Window.

    :param record:
    :param w:
    :param j:
    '''
    s_range = np.arange(0, len(record.seq), j)
    e_range = np.arange(0 + w, len(record.seq), j)
    # Ensure e_range has the same length as s_range
    if len(e_range) < len(s_range):
        e_range = np.concatenate((e_range, np.full(len(s_range) - len(e_range), len(record.seq))))
    s_e_matrix = np.vstack((s_range, e_range))
    return s_e_matrix


def lensplitRecord(record_iter, args):
    """Split sequence object by length.

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            tmp_s_e_matrix = sliding_window_steps_generator(record, args.w, args.j)
            for tmp_s, tmp_e in zip(tmp_s_e_matrix[0], tmp_s_e_matrix[1]):
                tmp_rec = SeqIO.SeqRecord(record.seq[tmp_s:tmp_e], name=record.name, id=record.id, description=record.description + ' ' + str(tmp_s + 1) + '_' + str(tmp_e))
                yield tmp_rec


def lensplitFasta(args, parser):
    record_iter = None
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    lensplit_iter = lensplitRecord(record_iter, args)
    if args.o is None:
        SeqIO.write(lensplit_iter, sys.stdout, "fasta")
    else:
        count = SeqIO.write(lensplit_iter, args.o, "fasta")
        print("splited into %i sequences" % count)


def define_parser():
    parser = argparse.ArgumentParser(prog='fastaLengthSplit', usage='%(prog)s [options] [<arguments>...]',
                                     description='Splits FASTA sequences by length')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
    parser.add_argument('-w', help='window size [default: 1000]', default=1000, type=int)
    parser.add_argument('-j', help='jump size [default: 1000]', default=1000, type=int)
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
    lensplitFasta(args, parser)


if __name__ == '__main__':
    main()
