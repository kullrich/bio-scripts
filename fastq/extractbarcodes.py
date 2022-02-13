#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian Ullrich
date: January 2022
email: ullrich@bevolbio.mpg.de
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
import gzip
import re
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from itertools import repeat
import numpy as np


def get_barcodes(i, batch, prog, out):
#def get_barcodes(batch, prog, out):
    filename = "%s_%i.out" % (out, i + 1)
    #filename = "%s_%i.out" % (out, batch[0] + 1)
    with open(filename, 'w') as outhandle:
        batch_total_count = 0
        batch_orig_count = 0
        batch_reverse_count = 0
        batch_complement_count = 0
        batch_reverse_complement_count = 0
        for record in batch:
        #for record in batch[1]:
            batch_total_count += 1
            #orig
            match = prog.search(str(record.seq))
            if match is not None:
                batch_orig_count += 1
                barcode = match.string[match.start():match.end()]
                outhandle.write(record.id+'\to\t'+barcode+'\n')
            #reverse
            match = prog.search(str(record.seq)[::-1])
            if match is not None:
                batch_reverse_count += 1
                barcode = match.string[match.start():match.end()]
                outhandle.write(record.id+'\tr\t'+barcode+'\n')
            #complement
            match = prog.search(str(record.seq.complement()))
            if match is not None:
                batch_complement_count += 1
                barcode = match.string[match.start():match.end()]
                outhandle.write(record.id+'\tc\t'+barcode+'\n')
            #reverse-complement
            match = prog.search(str(record.seq.reverse_complement()))
            if match is not None:
                batch_reverse_complement_count += 1
                barcode = match.string[match.start():match.end()]
                outhandle.write(record.id+'\trc\t'+barcode+'\n')
        return batch_total_count, batch_orig_count, batch_reverse_count, batch_complement_count, batch_reverse_complement_count


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


def main():
    parser = argparse.ArgumentParser(usage='%(prog)s [options]',
                                     description='This script extracts barcodes from FASTQ file given a left and a right border sequence returning only perfect matches. The read sequence will be evaluated in original, reverse, complement and reverse complement orientation. All matches will be returned specifying the hit direction',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument('-fq', help='specify input file in FASTQ format')
    parser.add_argument('-out', help='specify output file')
    parser.add_argument('-l', default="CCGATGTCCACGAAGCTCTCCTACG", help='specify left border [default:CCGATGTCCACGAAGCTCTCCTACG]')
    parser.add_argument('-r', default="CAGTCCAGCGCCAACCAGATAAGTG", help='specify right border [default:CAGTCCAGCGCCAACCAGATAAGTG]')
    parser.add_argument('-bl', type=int, default=25, help='specify length of barcode [default:25]')
    parser.add_argument('-s', type=int, default=1000000, help='specify batch size [default:1000000]')
    #parser.add_argument('-t', type=int, default=1, help='specify number of threads [default:1]')
    args = parser.parse_args()
    if args.fq is None:
        parser.print_help()
        sys.exit('\nPlease specify input FASTQ file')
    if args.out is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    prog = re.compile(args.l+'[A-Z]'+'{'+str(args.bl)+'}'+args.r)
    total_count = 0
    orig_count = 0
    reverse_count = 0
    complement_count = 0
    reverse_complement_count = 0
    if args.fq.split('.')[-1] == 'gz':
        with gzip.open(args.fq, 'rt') as handle:
            #hanlde = gzip.open(args.fq, 'rt')
            record_iter = SeqIO.parse(handle, 'fastq')
            for i, batch in enumerate(batch_iterator(record_iter, args.s)):
                batch_total_count, batch_orig_count, batch_reverse_count, batch_complement_count, batch_reverse_complement_count = get_barcodes(i, batch, prog, args.out)
                total_count += batch_total_count
                orig_count += batch_orig_count
                reverse_count += batch_reverse_count
                complement_count += batch_complement_count
                reverse_complement_count += batch_reverse_complement_count
            #b = zip(enumerate(batch_iterator(record_iter, args.s)), repeat(prog), repeat(args.out))
            #p = Pool(processes=args.t)
            #results = p.starmap_async(get_barcodes, b).get()
            #p.close()
            #handle.close()
            #print(results)
            #total_count, orig_count, reverse_count, complement_count, reverse_complement_count = np.sum(results, 0)
    else:
        with open(args.fq, 'rt') as handle:
            #handle = open(args.fq, 'rt')
            record_iter = SeqIO.parse(handle, 'fastq')
            for i, batch in enumerate(batch_iterator(record_iter, args.s)):
                batch_total_count, batch_orig_count, batch_reverse_count, batch_complement_count, batch_reverse_complement_count = get_barcodes(i, batch, prog, args.out)
                total_count += batch_total_count
                orig_count += batch_orig_count
                reverse_count += batch_reverse_count
                complement_count += batch_complement_count
                reverse_complement_count += batch_reverse_complement_count
            #b = zip(enumerate(batch_iterator(record_iter, args.s)), repeat(prog), repeat(args.out))
            #p = Pool(processes=args.t)
            #results = p.starmap_async(get_barcodes, b).get()
            #p.close()
            #handle.close()
            #print(results)
            #total_count, orig_count, reverse_count, complement_count, reverse_complement_count = np.sum(results, 0)
    print('total parsed reads:\t'+str(total_count))
    print('match orig reads:\t'+str(orig_count))
    print('match reverse reads:\t'+str(reverse_count))
    print('match complement reads:\t'+str(complement_count))
    print('match reverse complement reads:\t'+str(reverse_complement_count))


if __name__ == '__main__':
    main()
