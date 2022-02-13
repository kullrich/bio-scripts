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
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC


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


def gcRecord(record_iter, args):
    """Length of sequence objects.

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    tmp_gc = None
    gcdict = {}
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            tmp_id = record.id.split()[0]
            if args.gc == 'gc':
                tmp_gc = GC(record.seq)
            if args.gc == 'gc1':
                tmp_gc = GC(record.seq[::3])
            if args.gc == 'gc2':
                tmp_gc = GC(record.seq[1::3])
            if args.gc == 'gc3':
                tmp_gc = GC(record.seq[2::3])
            if tmp_id in gcdict:
                print('duplicated id: %s; skip GC calculation' % tmp_id)
            if tmp_id not in gcdict:
                gcdict[tmp_id] = tmp_gc
    return gcdict


def gcFasta(args, parser):
    record_iter = None
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    gcdict = gcRecord(record_iter, args)
    gcmean = np.mean(list(gcdict.values()))
    gcsd = np.std(list(gcdict.values()), ddof=1)
    if args.o is None:
        if args.gc == 'gc':
            print('#GC')
            print('mean GC percent: %s' % gcmean)
            print('sd GC percent: %s' % gcsd)
        if args.gc == 'gc1':
            print('#GC1')
            print('mean GC1 percent: %s' % gcmean)
            print('sd GC1 percent: %s' % gcsd)
        if args.gc == 'gc2':
            print('#GC2')
            print('mean GC2 percent: %s' % gcmean)
            print('sd GC2 percent: %s' % gcsd)
        if args.gc == 'gc3':
            print('#GC3')
            print('mean GC3 percent: %s' % gcmean)
            print('sd GC3 percent: %s' % gcsd)
        for k in sorted(gcdict.keys()):
            print('%s\t%s' % (k, gcdict[k]))
    else:
        with open(args.o, 'w') as outhandle:
            if args.gc == 'gc':
                outhandle.write('#GC\n')
            if args.gc == 'gc1':
                outhandle.write('#GC1\n')
            if args.gc == 'gc2':
                outhandle.write('#GC2\n')
            if args.gc == 'gc3':
                outhandle.write('#GC3\n')
            for k in sorted(gcdict.keys()):
                outhandle.write('%s\t%f\n' % (k, gcdict[k]))
        if args.gc == 'gc':
            print('#GC')
            print('mean GC percent: %s' % gcmean)
            print('sd GC percent: %s' % gcsd)
        if args.gc == 'gc1':
            print('#GC1')
            print('mean GC1 percent: %s' % gcmean)
            print('sd GC1 percent: %s' % gcsd)
        if args.gc == 'gc2':
            print('#GC2')
            print('mean GC2 percent: %s' % gcmean)
            print('sd GC2 percent: %s' % gcsd)
        if args.gc == 'gc3':
            print('#GC3')
            print('mean GC3 percent: %s' % gcmean)
            print('sd GC3 percent: %s' % gcsd)


def define_parser():
    parser = argparse.ArgumentParser(prog='fasta2GCcontent', usage='%(prog)s [options] [<arguments>...]',
                                     description='GC content of FASTA sequences')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
    parser.add_argument('-gc', default='gc', choices=['gc', 'gc1', 'gc2', 'gc3'], help='specify type [default: gc] or [gc1] or [gc2] or [gc3]')
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
    gcFasta(args, parser)


if __name__ == '__main__':
    main()
