#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian K Ullrich
date: Februar 2023
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
import os
import argparse
import tempfile
import subprocess
from Bio import SeqIO
from Bio import AlignIO
from Bio.Data import CodonTable
from Bio.codonalign import CodonSeq
from Bio.codonalign.codonseq import cal_dn_ds


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


transtable = {'std': CodonTable.CodonTable(forward_table={
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', '---': 'X', '--A': 'X', '--C': 'X', '--G': 'X',
    '--T': 'X', '-A-': 'X', '-AA': 'X', '-AC': 'X', '-AG': 'X',
    '-AT': 'X', '-C-': 'X', '-CA': 'X', '-CC': 'X', '-CG': 'X',
    '-CT': 'X', '-G-': 'X', '-GA': 'X', '-GC': 'X', '-GG': 'X',
    '-GT': 'X', '-T-': 'X', '-TA': 'X', '-TC': 'X', '-TG': 'X',
    '-TT': 'X', 'A--': 'X', 'A-A': 'X', 'A-C': 'X', 'A-G': 'X',
    'A-T': 'X', 'AA-': 'X', 'AC-': 'X', 'AG-': 'X', 'AT-': 'X',
    'C--': 'X', 'C-A': 'X', 'C-C': 'X', 'C-G': 'X', 'C-T': 'X',
    'CA-': 'X', 'CC-': 'X', 'CG-': 'X', 'CT-': 'X', 'G--': 'X',
    'G-A': 'X', 'G-C': 'X', 'G-G': 'X', 'G-T': 'X', 'GA-': 'X',
    'GC-': 'X', 'GG-': 'X', 'GT-': 'X', 'T--': 'X', 'T-A': 'X',
    'T-C': 'X', 'T-G': 'X', 'T-T': 'X', 'TA-': 'X', 'TC-': 'X',
    'TG-': 'X', 'TT-': 'X', 'NNN': 'X', 'GCN': 'A', 'CGN': 'R',
    'MGR': 'R', 'AAY': 'N', 'GAY': 'D', 'TGY': 'C', 'CAR': 'Q',
    'GAR': 'E', 'GGN': 'G', 'CAY': 'H', 'ATH': 'I', 'YTR': 'L',
    'CTN': 'L', 'AAR': 'K', 'TTY': 'F', 'CCN': 'P', 'TCN': 'S',
    'AGY': 'S', 'ACN': 'T', 'TAY': 'Y', 'GTN': 'V', 'TAR': '*',
    'TRA': '*'},
    stop_codons=['TAA', 'TAG', 'TGA', ],
    start_codons=['TTG', 'CTG', 'ATG', ])}


def cds2aaRecord(record_iter, batch_size, t_name):
    """Translates nucleotide to amino acids assuming that cds is in frame 0.
    :param record_iter:
    :param args:
    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
        for record in batch:
            aa = SeqIO.SeqRecord(record.seq.translate(transtable[t_name]), name=record.name, id=record.name, description=record.name)
            yield aa


def cds2preprocessRecord(record_iter, batch_size, char2gap_list, remove_gaps=False, char2gap=False):
    """Preprocess nucleotide to substitute char to gap and remove gaps.
    :param record_iter:
    :param args:
    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
        for record in batch:
            if char2gap:
                for cg in char2gap_list:
                    record.seq = record.seq.replace(cg, '-')
            if remove_gaps:
                record.seq = record.seq.replace('-', '')
            yield record


def cds2aaFasta(input, outdir, batch_size, t_name):
    record_iter = None
    if input is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if input is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, 'fasta')
    else:
        record_iter = SeqIO.parse(input, 'fasta')
    cds2aa_iter = cds2aaRecord(record_iter, batch_size, t_name)
    cds2aa_filename = os.path.basename(tempfile.NamedTemporaryFile().name) + '.fa'
    count = SeqIO.write(cds2aa_iter, os.path.join(outdir, cds2aa_filename), 'fasta')
    print('translated %i sequences' % count)
    return cds2aa_filename


def cds2preprocessFasta(input, outdir, batch_size, char2gap_list, remove_gaps=False, char2gap=False):
    record_iter = None
    if input is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if input is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, 'fasta')
    else:
        record_iter = SeqIO.parse(input, 'fasta')
    cds2preprocess_iter = cds2preprocessRecord(record_iter, batch_size, char2gap_list, remove_gaps, char2gap)
    cds2preprocess_filename = os.path.basename(tempfile.NamedTemporaryFile().name) + '.fa'
    count = SeqIO.write(cds2preprocess_iter, os.path.join(outdir, cds2preprocess_filename), 'fasta')
    print('preprocessed %i sequences' % count)
    return cds2preprocess_filename


def align(input, outdir, mafft_options):
    alg_filename = os.path.basename(tempfile.NamedTemporaryFile().name) + '.alg.fa'
    if mafft_options is not None:
        mafft_options.split(',')
        cmd = ['mafft'] + mafft_options + [input]
    else:
        cmd = ['mafft', input]
    alg_fo = open(os.path.join(outdir, alg_filename), 'w')
    subprocess.call(cmd, stderr=None, stdout=alg_fo, shell=False)
    alg_fo.close()
    return alg_filename


def get_codonalg(aa_input, cds_input, outdir):
    codonalg_filename = os.path.basename(tempfile.NamedTemporaryFile().name) + '.codonalg.fa'
    cmd = ['pal2nal.pl', aa_input, cds_input, '-output', 'fasta', '-nogap']
    codonalg_fo = open(os.path.join(outdir, codonalg_filename), 'w')
    subprocess.call(cmd, stderr=None, stdout=codonalg_fo, shell=False)
    codonalg_fo.close()
    return codonalg_filename


def codonalg2tree(input, outdir):
    tmp_filename = os.path.basename(tempfile.NamedTemporaryFile().name)
    phylip_filename =  tmp_filename + '.phylip'
    trimal_cmd = ['trimal', '-in', input, '-out', os.path.join(outdir, phylip_filename), '-phylip']
    subprocess.call(trimal_cmd, stderr=None, stdout=None, shell=False)
    phyml_filename = tmp_filename + '.phylip_phyml_tree.txt'
    phyml_cmd = ['phyml', '-i', os.path.join(outdir, phylip_filename)]
    subprocess.call(phyml_cmd, stderr=None, stdout=None, shell=False)
    return phyml_filename


def twosampletree(codonalg, outdir):
    phyml_filename = os.path.basename(tempfile.NamedTemporaryFile().name) + '.phylip_phyml_tree.txt'
    phyml_fo = open(os.path.join(outdir, phyml_filename), 'w')
    phyml_fo.write('(%s,%s);' % (codonalg[0].name, codonalg[2].name))
    phyml_fo.close()
    return phyml_filename


def run_codeml(input_fasta, input_tree, outdir, model=['M0', 'M1', 'M2'], leaves=False, internals=False, tests='M2,M1', cpu=1):
    if leaves and internals:
        ete3_evol_cmd = ['ete3', 'evol', '-t', input_tree, '--alg', input_fasta, '--models'] + model + ['--leaves', '--internals', '-o', outdir, '--codeml_param', 'cleandata,1', '--cpu', str(cpu)]
    if leaves and not internals:
        ete3_evol_cmd = ['ete3', 'evol', '-t', input_tree, '--alg', input_fasta, '--models'] + model + ['--leaves', '-o', outdir, '--codeml_param', 'cleandata,1', '--cpu', str(cpu)]
    if not leaves and internals:
        ete3_evol_cmd = ['ete3', 'evol', '-t', input_tree, '--alg', input_fasta, '--models'] + model + ['--internals', '-o', outdir, '--codeml_param', 'cleandata,1', '--cpu', str(cpu)]
    if not leaves and not internals:
        ete3_evol_cmd = ['ete3', 'evol', '-t', input_tree, '--alg', input_fasta, '--models'] + model + ['-o', outdir, '--codeml_param', 'cleandata,1', '--cpu', str(cpu)]
    if tests is not None:
        ete3_evol_cmd += ['--tests', tests]
    subprocess.call(ete3_evol_cmd, stderr=None, stdout=None, shell=False)


def define_parser():
    parser = argparse.ArgumentParser(prog='cds2codeml', usage='%(prog)s [options] [<arguments>...]',
                                     description='1. Translates nucleotide to amino acids assuming that CDS is in frame 0.\n'
                                                 '2. Creates amino acid alignments with MAFFT.\n'
                                                 '3. Creates Codon alignments with pal2nal.\n'
                                                 '4. Calculate dNdS with Bio.codonalign.codonseq.cal_dn_ds.\n'
                                                 '5. Create PhyML tree as input for codeml.\n'
                                                 '6. Run codeml.\n\n'
                                                 'Easiest to use if one uses a conda environment as follows:\n\n'
                                                 'conda create -n cds2codeml python=3.8\n'
                                                 'conda activate cds2codeml\n'
                                                 'conda install -c bioconda biopython\n'
                                                 'conda install -c bioconda ete3\n'
                                                 'conda install -c bioconda paml\n'
                                                 'conda install -c bioconda phyml\n'
                                                 'conda install -c bioconda mafft\n'
                                                 'conda install -c bioconda pal2nal\n'
                                                 'conda install -c bioconda trimal\n'
                                                 'conda install -c etetoolkit slr\n',
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', help='input file [mandatory]')
    parser.add_argument('-o', help='output folder [default: results]', default='results')
    parser.add_argument('-m', help='specify models as given here http://etetoolkit.org/documentation/ete-evol/ [default: M0,M1,M2]', default='M0,M1,M2')
    parser.add_argument('-leaves', help='specify to run a branch model for each branch [default: False]', action='store_true')
    parser.add_argument('-internals', help='specify to run a branch model for each internal node [default: False]', action='store_true')
    parser.add_argument('-tests', help='specify tests as given here http://etetoolkit.org/documentation/ete-evol/ [default: M2,M1]', default='M2,M1')
    parser.add_argument('-c', help='specify number of cpu [default: 2]', default=2, type=int)
    parser.add_argument('-s', help='batch size [default: 1000]', default=1000, type=int)
    parser.add_argument('-t', help='transtable [default: std]', default='std')
    parser.add_argument('-cal_dn_ds_method', help='specify method as indicated here https://biopython.org/docs/1.80/api/Bio.codonalign.codonseq.html [default: NG86]', default='NG86')
    parser.add_argument('-mafft_options', help='additional mafft options [default: None]')
    parser.add_argument('-remove_gaps', action='store_true', help='specify if gaps (-) in input CDS file should be removed [default: False]')
    parser.add_argument('-char2gap', action='store_true', help='specify if char should converted into gap char (-) in input CDS file [default: False]')
    parser.add_argument('-char2gap_list', default=['?'], nargs='+', help='specify char that should be converted into gap (-) in input CDS file [default: False]')
    return parser


def main():
    # parser
    parser = define_parser()
    # get args
    args = parser.parse_args()
    # print args
    print(args)
    # mkdir
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    # 0. Preprocess nucleotide
    if args.remove_gaps or args.char2gap:
        cds2preprocess_filename = cds2preprocessFasta(args.i, args.o, args.s, args.char2gap_list, args.remove_gaps, args.char2gap)
    # 1. Translates nucleotide to amino acids assuming that CDS is in frame 0.
        cds2aa_filename = cds2aaFasta(os.path.join(args.o, cds2preprocess_filename), args.o, args.s, args.t)
    else:
        cds2aa_filename = cds2aaFasta(args.i, args.o, args.s, args.t)
    # 2. Creates amino acid alignments with MAFFT
    alg_filename = align(os.path.join(args.o, cds2aa_filename), args.o, args.mafft_options)
    # 3. Creates Codon alignments with pal2nal
    if args.remove_gaps or args.char2gap:
        codonalg_filename = get_codonalg(os.path.join(args.o, alg_filename), os.path.join(args.o, cds2preprocess_filename), args.o)
    else:
        codonalg_filename = get_codonalg(os.path.join(args.o, alg_filename), args.i, args.o)
    # 4. Calculate dNdS with Bio.codonalign.codonseq.cal_dn_ds
    codonalg = AlignIO.read(os.path.join(args.o, codonalg_filename), 'fasta')
    cal_dn_ds_list = []
    for i,i_rec in enumerate(codonalg):
        for j,j_rec in enumerate(codonalg[(i+1):]):
            i_j_dn_ds = cal_dn_ds(CodonSeq(i_rec.seq), CodonSeq(j_rec.seq), method=args.cal_dn_ds_method)
            cal_dn_ds_list += [(i+1, i+2+j, i_rec.name, j_rec.name, i_j_dn_ds)]
    cal_dn_ds_filename = os.path.basename(args.i) + '.cal_dNdS.tsv'
    cal_dn_ds_fo = open(os.path.join(args.o, cal_dn_ds_filename), 'w')
    cal_dn_ds_fo.write('comp1\tcomp2\tname1\tname2\tdN\tdS\tmodel\n')
    for comp in cal_dn_ds_list:
        cal_dn_ds_fo.write(str(comp[0]) + '\t' + str(comp[1]) + '\t' + comp[2] + '\t' + comp[3] + '\t' + str(comp[4][0]) + '\t' + str(comp[4][1]) + '\t' + args.cal_dn_ds_method + '\n')
    cal_dn_ds_fo.close()
    # 5. Create PhyML tree as input for codeml
    if len(codonalg)>2:
        phyml_filename = codonalg2tree(os.path.join(args.o, codonalg_filename), args.o)
    else:
        phyml_filename = twosampletree(codonalg, args.o)
    # 6. Run codeml
    run_codeml(os.path.join(args.o, codonalg_filename), os.path.join(args.o, phyml_filename), args.o, args.m.split(','), args.c)
    print('Done')


if __name__ == '__main__':
    main()
