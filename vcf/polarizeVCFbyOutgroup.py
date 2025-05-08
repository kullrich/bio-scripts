#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Krisian Ullrich
date: May 2025
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
import textwrap
import gzip
import re


#def multiple_replace(string, rep_dict):
#    pattern = re.compile("|".join([re.escape(k) for k in sorted(rep_dict,key=len,reverse=True)]), flags=re.DOTALL)
#    return pattern.sub(lambda x: rep_dict[x.group(0)], string)
def swap_genotypes(genotypes_str, rep_dict):
    genotypes = genotypes_str.strip().split('\t')
    swapped = []
    for gt in genotypes:
        if gt in rep_dict:
            swapped.append(rep_dict[gt])
        else:
            print(f"Warning: Genotype '{gt}' not found in replacement dictionary. Keeping as is.")
            swapped.append(gt)
    return '\t'.join(swapped)


def get_chrom_line(fin):
    for line in fin:
        if line.startswith('#') and line.startswith('#CHROM'):
            return line.strip().split('\t')


def parse_lines(fin, fou, ind_names, keep, add):
    switchcount = 0
    removecount = 0
    outmissing = 0
    totalcount = 0
    outgroup_indices = []
    rep_dict = {
        '0/0': '1/1',
        '1/1': '0/0',
        '0/1': '1/0',
        '1/0': '0/1',
        './.': './.'
    }
    for line in fin:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                sample_names = headers[9:]
                missing = [name for name in ind_names if name not in sample_names]
                if missing:
                    raise ValueError(f"Sample name(s) not found in VCF header: {', '.join(missing)}")
                outgroup_indices = [sample_names.index(name) + 9 for name in ind_names]
                if add:
                    fou.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n')
                fou.write(line)
            else:
                fou.write(line)
            continue
        totalcount += 1
        fields = line.strip().split('\t')
        outgroup_gts = [fields[i].split(':')[0] for i in outgroup_indices]
        if any(gt in ['./.', '.|.'] for gt in outgroup_gts):
            outmissing += 1
            continue
        if any(gt in ['0/1', '1/0', '0|1', '1|0'] for gt in outgroup_gts):
            removecount += 1
            if keep:
                fou.write(line)
            continue
        if all(gt in ['0/0', '0|0'] for gt in outgroup_gts):
            if add:
                fields[7] = f'AA={fields[3]};' + fields[7]
            fou.write('\t'.join(fields) + '\n')
        elif all(gt in ['1/1', '1|1'] for gt in outgroup_gts):
            switchcount += 1
            REF = fields[4]
            ALT = fields[3]
            fields[3] = REF
            fields[4] = ALT
            genotypes = '\t'.join(fields[9:])
            swapped_genotypes = swap_genotypes(genotypes, rep_dict)
            if add:
                fields[7] = f'AA={fields[3]};' + fields[7]
            fou.write('\t'.join(fields[:9]) + '\t' + swapped_genotypes + '\n')
        else:
            removecount += 1
            if keep:
                fou.write(line)

    print(f'Parsed {totalcount} sites.')
    print(f'Removed {outmissing} sites due to missing allele info in any outgroup individual.')
    print(f'{"Kept" if keep else "Removed"} {removecount} sites with undefined ancestral state.')
    print(f'Switched REF and ALT allele for {switchcount} sites.')


def main():
    parser = argparse.ArgumentParser(
        prog='polarizeVCFbyOutgroup',
        description='Switch REF and ALT allele of a VCF file if all specified individuals are homozygous ALT.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-vcf', required=True, help='VCF input file (VCF must be bi-allelic and GT-only)')
    parser.add_argument('-out', required=True, help='Output VCF file')
    parser.add_argument('-ind', required=True, nargs='+', help='One or more individual names used as outgroup')
    parser.add_argument('-keep', action='store_true', help='Keep sites with undefined ancestral states')
    parser.add_argument('-add', action='store_true', help='Add ancestral allele info (AA=) to INFO field')
    parser.add_argument('-show_ind', action='store_true', help='Show sample indices and exit')
    args = parser.parse_args()

    if args.out.endswith('.gz'):
        fout = gzip.open(args.out, 'wt')
    else:
        fout = open(args.out, 'wt')

    if args.vcf.endswith('.gz'):
        fin = gzip.open(args.vcf, 'rt')
    else:
        fin = open(args.vcf, 'rt')

    if args.show_ind:
        chrom_line = get_chrom_line(fin)
        for idx, name in enumerate(chrom_line):
            print(f'{idx} : {name}')
        sys.exit(0)

    parse_lines(fin, fout, args.ind, args.keep, args.add)
    fin.close()
    fout.close()


if __name__ == '__main__':
    main()
