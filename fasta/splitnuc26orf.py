#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian Ullrich
date: June 2022
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
import re
from Bio import SeqIO


def splitnuc26orf(records, fout, fmp='\*[^\*]', bmp='[^\*]\*', fbmp='\**[^\*]*\**', min=0):
    front_match_pattern = fmp
    back_match_pattern = bmp
    front_back_match_pattern = fbmp
    front_pattern = re.compile(front_match_pattern)
    back_pattern = re.compile(back_match_pattern)
    front_back_pattern = re.compile(front_back_match_pattern)
    for record in records:
        rec_len = len(record)
        rec_orf = record.id.split('_')[-1] 
        front_iter = front_pattern.finditer('*' + str(record.seq))
        back_iter = back_pattern.finditer(str(record.seq) + '*')
        front_back_iter = front_back_pattern.finditer('*' + str(record.seq) + '*')
        for m, n, l in zip(front_iter, back_iter, front_back_iter):
            out_str = str(l.group(0)).replace('*','')
            out_str = out_str.replace("X", "")
            if len(out_str) >= min:
                if rec_orf=='0+':
                    fout.write('>' + record.id + '_%i_%i' % ((((m.end()-1) * 3) + 0), ((n.end() * 3) + 0) + '\n')
                if rec_orf=='1+':
                    fout.write('>' + record.id + '_%i_%i' % ((((m.end()-1) * 3) + 1), ((n.end() * 3) + 1) + '\n')
                if rec_orf=='2+':
                    fout.write('>' + record.id + '_%i_%i' % ((((m.end()-1) * 3) + 2), ((n.end() * 3) + 2) + '\n')
                if rec_orf=='0-':
                    fout.write('>' + record.id + '_%i_%i' % (rec_len - (((m.end()-1) * 3) + 2), rec_len - ((n.end() * 3) + 2) + '\n')
                if rec_orf=='1-':
                    fout.write('>' + record.id + '_%i_%i' % (rec_len - (((m.end()-1) * 3) + 1), rec_len - ((n.end() * 3) + 1) + '\n')
                if rec_orf=='2-':
                    fout.write('>' + record.id + '_%i_%i' % (rec_len - (((m.end()-1) * 3) + 0), rec_len - ((n.end() * 3) + 0) + '\n')
                fout.write(out_str + '\n')


def main():
    original_fasta = SeqIO.parse(sys.stdin, "fasta")
    splitnuc26orf(original_fasta, sys.stdout)


if __name__ == '__main__':
    main()
