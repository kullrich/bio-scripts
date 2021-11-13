#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Krisian Ullrich
date: November 2021
email: ullrich@evolbio.mpg.de
License: MIT
The MIT License (MIT)
Copyright (c) 2021 Kristian Ullrich
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
import textwrap
import gzip
import re
from collections import Counter


def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in sorted(rep_dict,key=len,reverse=True)]), flags=re.DOTALL)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)


def parse_lines(fin, fou, phased, singletons):
    singletonscount = 0
    totalcount = 0
    equalcount = 0
    for line in fin:
        if line[0] != '#':
            totalcount += 1
            linesplit = line.strip().split('\t')
            CHR = linesplit[0]
            POS = linesplit[1]
            REF = linesplit[3]
            ALT = linesplit[4]
            SNPS = linesplit[9:]
            freq = Counter(multiple_replace(','.join(SNPS),{'|':',','/':','}).split(','))
            if len(freq) == 1:
                equalcount += 1
                continue
            if singletons:
                if freq['0'] == 1 or freq['1'] == 1:
                    singletonscount += 1
                    continue
            if phased:
                changed = multiple_replace(','.join(linesplit[9:]),{
                    '0|0':REF+','+REF,
                    '0|1':REF+','+ALT,
                    '1|0':ALT+','+REF,
                    '1|1':ALT+','+ALT,
                    '0|.':REF+',N',
                    '1|.':ALT+',N',
                    '.|0':'N,'+REF,
                    '.|1':'N,'+ALT,
                    '.|.':'N,N',
                    '0/0':REF+','+REF,
                    '0/1':'.,.',
                    '1/0':'.,.',
                    '1/1':ALT+','+ALT,
                    '0/.':'.,.',
                    '1/.':'.,.',
                    './0':'.,.',
                    './1':'.,.',
                    './.':'N,N'})
            else:
                changed = multiple_replace(','.join(linesplit[9:]),{
                    '0|0':REF,
                    '0|1':'.',
                    '1|0':'.',
                    '1|1':ALT,
                    '0|.':'.',
                    '1|.':'.',
                    '.|0':'.',
                    '.|1':'.',
                    '.|.':'N',
                    '0/0':REF,
                    '0/1':'.',
                    '1/0':'.',
                    '1/1':ALT,
                    '0/.':'.',
                    '1/.':'.',
                    './0':'.',
                    './1':'.',
                    './.':'N'}) 
            fou.write(POS + ',' + changed + '\n')
    print('Parsed ' + str(totalcount) + ' sites.')
    print('Removed ' + str(equalcount) + ' non variable sites.')
    if singletons:
        print('Removed ' + str(singletonscount) + ' singletons.')


def main():
    parser = argparse.ArgumentParser(prog='vcf2hscan', description='Create H-scan input from VCF.', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-vcf', help=textwrap.dedent('''\
specify vcf input file
vcf file should only contain bi-allelic sites and only GT field
bcftools commands to retain only bi-allelic sites and GT field:
(bcftools view -h VCFFILE;
 bcftools query -f
 "%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%QUAL\\t%%FILTER\\t%%INFO\\tGT\\t[%%GT\\t]\\n" VCFFILE)
| cat | bcftools view -m2 -M2 -v snps
 '''))
    parser.add_argument('-out', help='specify output file')
    parser.add_argument('-phased', action='store_true', help='specify if vcf input is phased')
    parser.add_argument('-singletons', action='store_true', help='specify if singletons should be removed from output')
    args = parser.parse_args()
    print(args)
    if args.vcf is None:
        parser.print_help()
        sys.exit('Please specify vcf input file')
    if args.out is None:
        parser.print_help()
        sys.exit('Please specify output file')
    if args.out.endswith('gz'):
        with gzip.open(args.out, 'wt') as fou:
            if args.vcf.endswith('gz'):
                with gzip.open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.phased, args.singletons)
            if not args.vcf.endswith('gz'):
                with open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.phased, args.singletons)
    if not args.out.endswith('gz'):
        with open(args.out, 'wt') as fou:
            if args.vcf.endswith('gz'):
                with gzip.open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.phased, args.singletons)
            if not args.vcf.endswith('gz'):
                with open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.phased, args.singletons)


if __name__ == '__main__':
    main()
