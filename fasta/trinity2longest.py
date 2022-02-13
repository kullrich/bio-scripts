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


def lctrinityRecord(record_iter, args):
    """Extracts longest trinity transcripts component

    :param record_iter:
    :param args:

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    tmp_id = ''
    seq_dict = {}
    for i, batch in enumerate(batch_iterator(record_iter, args.s)):
        for record in batch:
            if args.t == 'gene':
                tmp_id = record.id.split('_i')[0]
            if args.t == 'component':
                tmp_id = record.id.split('_g')[0]
            if tmp_id in seq_dict:
                if len(record) > len(seq_dict[tmp_id]):
                    seq_dict[tmp_id] = record
            if tmp_id not in seq_dict:
                seq_dict[tmp_id] = record
    return seq_dict


def lctrinityFasta(args, parser):
    record_iter = None
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    lctrinity_dict = lctrinityRecord(record_iter, args)
    if args.o is None:
        for k in sorted(lctrinity_dict.keys()):
            SeqIO.write(lctrinity_dict[k], sys.stdout, "fasta")
    else:
        count = SeqIO.write(list(lctrinity_dict.values()), args.o, "fasta")
        print("kept %i sequences" % count)


def define_parser():
    parser = argparse.ArgumentParser(prog='trinity2longest', usage='%(prog)s [options] [<arguments>...]',
                                     description='Extracts longest trinity transcripts component')
    parser.add_argument('-i', help='input file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
    parser.add_argument('-t', default='gene', choices=['gene', 'component'], help='specify type of filter [default: gene] or [component]')
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
    lctrinityFasta(args, parser)


if __name__ == '__main__':
    main()
