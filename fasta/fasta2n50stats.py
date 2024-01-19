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


def n50statsRecord(record_iter, args):
    """n50 stats from sequence objects.

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    gcdict = {}
    lendict = {}
    n_count = 0
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            tmp_id = record.id.split()[0]
            tmp_gc = GC(record.seq)
            tmp_len = len(record)
            n_count += record.seq.count('N')
            if tmp_id in gcdict:
                print('duplicated id: %s; skip GC and length calculation' % tmp_id)
            if tmp_id not in gcdict:
                gcdict[tmp_id] = tmp_gc
                if tmp_len in lendict:
                    lendict[tmp_len] += 1
                if tmp_len not in lendict:
                    lendict[tmp_len] = 1
    gcmean = np.mean(list(gcdict.values()))
    gcsd = np.std(list(gcdict.values()), ddof=1)
    lenu = np.repeat(list(lendict.keys()), list(lendict.values()))
    lenr = np.repeat(list(lendict.keys()), [x[0] * x[1] for x in lendict.items()])
    minlen = np.min(list(lendict.keys()))
    maxlen = np.max(list(lendict.keys()))
    meanlen = np.mean(lenu)
    medianlen = np.median(lenu)
    totallen = len(lenr)
    number_contigs = len(lenu)
    n5 = np.percentile(lenr, 95)
    n25 = np.percentile(lenr, 75)
    n50 = np.percentile(lenr, 50)
    n75 = np.percentile(lenr, 25)
    n95 = np.percentile(lenr, 5)
    gc_len_summary = [number_contigs, gcmean, gcsd, totallen, minlen, maxlen, meanlen, medianlen, n5, n25, n50, n75, n95, n_count]
    return gc_len_summary


def n50statsFasta(args, parser):
    record_iter = None
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    gc_len_summary = n50statsRecord(record_iter, args)
    if args.o is None:
        print('#N50stats')
        print('%s\t%s' % ('n', gc_len_summary[0]))
        print('%s\t%s' % ('GCmean', gc_len_summary[1]))
        print('%s\t%s' % ('GCsd', gc_len_summary[2]))
        print('%s\t%s' % ('Sum', gc_len_summary[3]))
        print('%s\t%s' % ('Min', gc_len_summary[4]))
        print('%s\t%s' % ('Max', gc_len_summary[5]))
        print('%s\t%s' % ('Mean', gc_len_summary[6]))
        print('%s\t%s' % ('Median', gc_len_summary[7]))
        print('%s\t%s' % ('N5', gc_len_summary[8]))
        print('%s\t%s' % ('N25', gc_len_summary[9]))
        print('%s\t%s' % ('N50', gc_len_summary[10]))
        print('%s\t%s' % ('N75', gc_len_summary[11]))
        print('%s\t%s' % ('N95', gc_len_summary[12]))
        print('%s\t%s' % ('Ns', gc_len_summary[13]))
    else:
        with open(args.o, 'w') as outhandle:
            outhandle.write('#N50stats\n')
            outhandle.write('%s\t%s\n' % ('n', gc_len_summary[0]))
            outhandle.write('%s\t%s\n' % ('GCmean', gc_len_summary[1]))
            outhandle.write('%s\t%s\n' % ('GCsd', gc_len_summary[2]))
            outhandle.write('%s\t%s\n' % ('Sum', gc_len_summary[3]))
            outhandle.write('%s\t%s\n' % ('Min', gc_len_summary[4]))
            outhandle.write('%s\t%s\n' % ('Max', gc_len_summary[5]))
            outhandle.write('%s\t%s\n' % ('Mean', gc_len_summary[6]))
            outhandle.write('%s\t%s\n' % ('Median', gc_len_summary[7]))
            outhandle.write('%s\t%s\n' % ('N5', gc_len_summary[8]))
            outhandle.write('%s\t%s\n' % ('N25', gc_len_summary[9]))
            outhandle.write('%s\t%s\n' % ('N50', gc_len_summary[10]))
            outhandle.write('%s\t%s\n' % ('N75', gc_len_summary[11]))
            outhandle.write('%s\t%s\n' % ('N95', gc_len_summary[12]))
            outhandle.write('%s\t%s\n' % ('Ns', gc_len_summary[13]))


def define_parser():
    parser = argparse.ArgumentParser(prog='fasta2n50stats', usage='%(prog)s [options] [<arguments>...]',
                                     description='n50 stats from FASTA sequences')
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
    n50statsFasta(args, parser)


if __name__ == '__main__':
    main()
