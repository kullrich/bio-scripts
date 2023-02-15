#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
This script extracts codonusage from CDS input FASTA file and writes a table containing the number of codons used for each sequence in the input FASTA file. In addition to the raw codon counts global ACTG counts, ACTG counts for the first, second and thrid codon position and Relative Synonymous Codon Usage (RSCU) is written to a file. Optional different methods can be applied to calculate Effective Number of Codons (ENC) can be choosen and calculated.

Author: Krisian Ullrich
date: Feb 2023
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2023 Kristian Ullrich

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
import copy
import numpy as np
import pandas as pd
from Bio import SeqIO


class GeneticCode(object):
    def __init__(self, transl_table=1, pseudocount=1, six2fourtwo=False):
        self.transl_table = transl_table
        self.pseudocount = pseudocount
        self.six2fourtwo = six2fourtwo
        self.codons_df = self.get_codons_df()
        self.codon_count_table = self.get_codon_count_table()
        self.codons_missing = 0
        self.actg_count_table = self.get_actg_count_table()
        self.actg_missing = {0:0, 1:0, 2:0, 3:0}
    def get_codons_df(self):
        codons = ['TTT', 'TTC', 'TTA', 'TTG',
              'TCT', 'TCC', 'TCA', 'TCG',
              'TAT', 'TAC', 'TAA', 'TAG',
              'TGT', 'TGC', 'TGA', 'TGG',
              'CTT', 'CTC', 'CTA', 'CTG',
              'CCT', 'CCC', 'CCA', 'CCG',
              'CAT', 'CAC', 'CAA', 'CAG',
              'CGT', 'CGC', 'CGA', 'CGG',
              'ATT', 'ATC', 'ATA', 'ATG',
              'ACT', 'ACC', 'ACA', 'ACG',
              'AAT', 'AAC', 'AAA', 'AAG',
              'AGT', 'AGC', 'AGA', 'AGG',
              'GTT', 'GTC', 'GTA', 'GTG',
              'GCT', 'GCC', 'GCA', 'GCG',
              'GAT', 'GAC', 'GAA', 'GAG',
              'GGT', 'GGC', 'GGA', 'GGG']
        AAs = ['',
           'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
           'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           '',
           '',
           'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
           'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
           '',
           'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           '',
           '',
           '',
           '',
           'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
           'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
           '',
           'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG']
        AAs6242 = ['',
           'FFllSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKss**VVVVAAAADDEEGGGG',
           'FFLLSSSSYY**CCWWttttPPPPHHQQRRRRIIMMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKssssVVVVAAAADDEEGGGG',
           'FFllSSSSYYqqCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           '',
           '',
           'FFllSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKssssVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CC*WLLL?PPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKssggVVVVAAAADDEEGGGG',
           'FFllSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKssssVVVVAAAADDEEGGGG',
           '',
           'FFllSSSSYY*!CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           '',
           '',
           '',
           '',
           'FFllSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKssssVVVVAAAADDEEGGGG',
           'FFllSS*SYY*!CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FF*lSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKsssKVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CCgWLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYY**CC*WLLLaPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           'FFllSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKssrrVVVVAAAADDEEGGGG',
           '',
           'FFllSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKssskVVVVAAAADDEEGGGG']
        NAMEs = ['',
             'Standard',
             'Vertebrate_Mitochondrial',
             'Yeast_Mitochondrial',
             'Mold_Protozoan_and_Coelenterate_Mitochondria_and_Mycoplasma_Spiroplasma',
             'Invertebrate_Mitochondrial',
             'Ciliate_Dasycladacean_and_Hexamita_Nuclear',
             '',
             '',
             'Echinoderm_and_Flatworm_Mitochondrial',
             'Euplotid_Nuclear',
             'Bacterial_Archaeal_and_Plant_Plastid',
             'Alternative_Yeast_Nuclear',
             'Ascidian_Mitochondrial',
             'Alternative_Flatworm_Mitochondrial',
             '',
             'Chlorophycean_Mitochondrial',
             '',
             '',
             '',
             '',
             'Trematode_Mitochondrial',
             'Scenedesmus_obliquus_Mitochondrial',
             'Thraustochytrium Mitochondrial',
             'Rhabdopleuridae Mitochondrial',
             'Candidate_Division_SR1_and_Gracilibacteria',
             'Pachysolen_tannophilus_Nuclear',
             'Karyorelict_Nuclear',
             'Condylostoma_Nuclear',
             'Mesodinium_Nuclear',
             'Peritrich_Nuclear',
             'Blastocrithidia_Nuclear',
             '',
             'Cephalodiscidae_Mitochondrial_UAA_Tyr']
        if self.six2fourtwo:
            codons_df = pd.DataFrame(codons, columns=['codons'])
            for aa_idx, aa in enumerate(AAs6242):
                if aa != '':
                    codons_df[aa_idx] = [x for x in aa]
            codons_df.set_index('codons',inplace=True)
        else:
            codons_df = pd.DataFrame(codons, columns=['codons'])
            for aa_idx, aa in enumerate(AAs):
                if aa != '':
                    codons_df[aa_idx] = [x for x in aa]
            codons_df.set_index('codons',inplace=True)
        return codons_df
    def get_codon_count_table(self):
        count_table_df = pd.DataFrame(list(self.codons_df.index), columns=['codons'])
        count_table_df[self.transl_table] = [[x, self.codons_df[self.transl_table].value_counts()[x], self.pseudocount] for x in self.codons_df[self.transl_table]]
        count_table_df.set_index('codons',inplace=True)
        return count_table_df[self.transl_table].to_dict()
    def get_actg_count_table(self):
        actgtablecount = {0:{
            'A': self.pseudocount,
            'C': self.pseudocount,
            'T': self.pseudocount,
            'G': self.pseudocount}, 1:{
            'A': self.pseudocount,
            'C': self.pseudocount,
            'T': self.pseudocount,
            'G': self.pseudocount}, 2:{
            'A': self.pseudocount,
            'C': self.pseudocount,
            'T': self.pseudocount,
            'G': self.pseudocount}, 3:{
            'A': self.pseudocount,
            'C': self.pseudocount,
            'T': self.pseudocount,
            'G': self.pseudocount}}
        return actgtablecount
    def set_record_id(self, record_id):
        self.record_id = record_id
    def set_record_len(self, record_len):
        self.record_len = record_len
    def set_record_mod3(self, record_mod3):
        self.record_mod3 = record_mod3
    def get_actg_freq(self):
        actg_freq_df = pd.DataFrame(np.array([list(x.values()) for x in self.actg_count_table.values()]), columns=['A', 'C', 'T', 'G'])
        self.actg_freq_df = actg_freq_df.div(actg_freq_df.sum(1), axis='rows')
    def get_codon_freq(self):
        aa_count = pd.DataFrame(self.codon_count_table).transpose().groupby([0])[2].sum().to_dict()
        codon_freq_table = copy.deepcopy(self.codon_count_table)
        for k in codon_freq_table.keys():
            codon_freq_table[k][2] = codon_freq_table[k][2]/aa_count[codon_freq_table[k][0]]
        self.codon_freq_table = codon_freq_table
    def get_eq1sun(self):
        aa_count_df = pd.DataFrame(self.codon_count_table).transpose().groupby([0])[2].sum()
        #only calculate for codon classes with counts
        return pd.DataFrame(self.codon_count_table).transpose()[self.codons_df[self.transl_table].isin(aa_count_df[aa_count_df!=0].index)].groupby([0])[2].apply(lambda x: (((x.sum() * ((x / x.sum())**2)).sum())-1)/(x.sum()-1))






def main():
    parser = argparse.ArgumentParser(usage='%(prog)s [options]',
                                     description='Extracts codonusage from CDS input FASTA file. Output will be raw codon counts (.codoncnt), global ACTG counts (.actgcnt), first (.firstcnt), second (.secondcnt), third (.third) codon position counts and Relative Synonymous Codon Usage (.rscucnt). Optional different methods can be applied to calculate Effective Number of Codons (.enc).')
    parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-i', help='specify CDS input file in FASTA format')
    parser.add_argument('-o', help='specify output prefix')
    parser.add_argument('-gc', help='specify genetic code', default=1)
    parser.add_argument('-pseudocount', help='specify pseudocount', default=1)
    parser.add_argument('-remstop', help='specify if stop codon should be removed prior score calculation', action='store_true')
    parser.add_argument('-mod3', action='store_true',
                        help='specify if CDS sequences with length modulo 3 unequal to 0 should be removed and reported to std.out')
    parser.add_argument('-enc', choices=['eq4Wright', 'eq2Sun', 'eq5Sun', 'all'],
                        help='specify equation to calculate ENC. Either equation (4) [eq4Wright] of (Wright. (1990) Gene 87:23-29) or equation (2) [eq2Sun] or equation (5) [eq5Sun] of (Sun et al. (2012) Mol. Biol. Evol. 30:191-196) or [all].')
    parser.add_argument('-six2fourtwo', action='store_true',
                        help='specify if codons that code for the same amino acid but start with different nucleotides should be treated as distinct codon classes [default: False]. This will only affect calculation of ENC values.')
    args = parser.parse_args()

    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input fasta file')
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output prefix')
    if args.i == args.o:
        sys.exit('\nInput file and output prefix are identical, use "out" as output prefix instead')

    print('\ncommand arguments used:\n')
    print(args)

    infile = args.i
    outfile_codoncount = args.o + '.codoncnt'
    outfile_actgcount = args.o + '.actgcnt'
    outfile_firstcount = args.o + '.firstcnt'
    outfile_secondcount = args.o + '.secondcnt'
    outfile_thirdcount = args.o + '.thirdcnt'
    outfile_rscucount = args.o + '.rscucnt'
    outfile_enc = args.o + '.enc'

    print('read fasta')
    original_fasta = SeqIO.parse(infile, 'fasta')
    print('extract codon counts')
    global_gc = GeneticCode(transl_table=args.gc, pseudocount=0, six2fourtwo=args.six2fourtwo)
    ids_mod3 = []
    with open(outfile_codoncount, 'w') as codonhandle:
        with open(outfile_actgcount, 'w') as actghandle:
            with open(outfile_firstcount, 'w') as firsthandle:
                with open(outfile_secondcount, 'w') as secondhandle:
                    with open(outfile_thirdcount, 'w') as thirdhandle:
                        with open(outfile_rscucount, 'w') as rscuhandle:
                            if args.enc is not None:
                                enchandle = open(outfile_enc, 'w')
                                if args.enc == 'eq4Wright':
                                    enchandle.write('id\tlen\tmo3\tgc3\teq4Wright\n')
                                if args.enc == 'eq2Sun':
                                    enchandle.write('id\tlen\tmo3\tgc3\teq2Sun\n')
                                if args.enc == 'eq5Sun':
                                    enchandle.write('id\tlen\tmo3\tgc3\teq5sun\n')
                                if args.enc == 'all':
                                    enchandle.write('id\tlen\tmo3\tgc3\teq4Wright\teq2Sun\teq5Sun\n')
                            codonhandle.write('id\tlen\tmo3\t' + '\t'.join(sorted(global_codons.keys())) + '\n')
                            actghandle.write('id\tlen\tmo3\t' + '\t'.join(sorted(global_actg.keys())) + '\n')
                            firsthandle.write('id\tlen\tmo3\t' + '\t'.join(sorted(global_first.keys())) + '\n')
                            secondhandle.write('id\tlen\tmo3\t' + '\t'.join(sorted(global_second.keys())) + '\n')
                            thirdhandle.write('id\tlen\tmo3\t' + '\t'.join(sorted(global_third.keys())) + '\n')
                            rscuhandle.write('id\tlen\tmo3\t' + '\t'.join(sorted(global_rscu.keys())) + '\n')
                            c = 0
                            cmod3 = 0
                            for record in original_fasta:
                                c += 1
                                record_gc = GeneticCode(transl_table=args.gc, pseudocount=args.pseudocount, six2fourtwo=args.six2fourtwo)
                                record_gc.set_record_id(record.id)
                                record_gc.set_record_len(len(record))
                                record_gc.set_record_mod3(0)
                                if len(record) % 3 != 0:
                                    record_gc.set_record_mod3(1)
                                    cmod3 += 1
                                    ids_mod3.append(record_gc.record_id)
                                    if args.remstop:
                                        continue
                                for i in range(0, len(record), 3):
                                    codon = str(record[i:i + 3].seq)
                                    for nuc_idx, nuc in enumerate(codon):
                                        if nuc in global_gc.actg_count_table[0]:
                                            global_gc.actg_count_table[0][nuc] += 1
                                            record_gc.actg_count_table[0][nuc] += 1
                                            global_gc.actg_count_table[nuc_idx + 1][nuc] += 1
                                            record_gc.actg_count_table[nuc_idx + 1][nuc] += 1
                                        else:
                                            global_gc.actg_missing[0] += 1
                                            record_gc.actg_missing[0] += 1
                                            global_gc.actg_missing[nuc_idx + 1] += 1
                                            record_gc.actg_missing[nuc_idx + 1] += 1
                                    if len(codon) != 3:
                                        continue
                                    if codon in global_gc.codon_count_table:
                                        global_gc.codon_count_table[codon][2] += 1
                                        record_gc.codon_count_table[codon][2] += 1
                                    if codon not in global_gc.codon_count_table:
                                        global_gc.codons_missing += 1
                                        record_gc.codons_missing += 1
                                if args.enc is not None:
                                    tmp_enc = {}
                                    tmp_gcbypos = gcbypos(tmp_counts, six2fourtwo)
                                    if args.enc == 'eq4Wright':
                                        enchandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + str(
                                            tmp_gcbypos[2]) + '\t' + str(
                                            calc_eq4wright(tmp_gcbypos[2])) + '\n')
                                    if args.enc == 'eq2Sun':
                                        enchandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + str(
                                            tmp_gcbypos[2]) + '\t' + str(
                                            calc_eq2sun(tmp_counts, six2fourtwo)) + '\n')
                                    if args.enc == 'eq5Sun':
                                        enchandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + str(
                                            tmp_gcbypos[2]) + '\t' + str(
                                            calc_eq5sun(tmp_counts, six2fourtwo)) + '\n')
                                    if args.enc == 'all':
                                        enchandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + str(
                                            tmp_gcbypos[2]) + '\t' + str(
                                            calc_eq4wright(tmp_gcbypos[2])) + '\t' + str(
                                            calc_eq2sun(tmp_counts, six2fourtwo)) + '\t' + str(
                                            calc_eq5sun(tmp_counts, six2fourtwo)) + '\n')
                                tmp_rscu = calc_rscu(tmp_counts, codontable)
                                codonhandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + '\t'.join(
                                    [str(tmp_counts[x][2]) for x in sorted(tmp_counts.keys())]) + '\n')
                                actghandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + '\t'.join(
                                    [str(tmp_actgcounts[x][1]) for x in sorted(tmp_actgcounts.keys())]) + '\n')
                                firsthandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + '\t'.join(
                                    [str(tmp_first[x][1]) for x in sorted(tmp_first.keys())]) + '\n')
                                secondhandle.write(
                                    tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + '\t'.join(
                                        [str(tmp_second[x][1]) for x in sorted(tmp_second.keys())]) + '\n')
                                thirdhandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + '\t'.join(
                                    [str(tmp_third[x][1]) for x in sorted(tmp_third.keys())]) + '\n')
                                rscuhandle.write(tmp_id + '\t' + str(tmp_len) + '\t' + str(tmp_mo3) + '\t' + '\t'.join(
                                    [str(tmp_rscu[x][2]) for x in sorted(tmp_rscu.keys())]) + '\n')
                            print('finished codon counting')
                            codonhandle.write('global_count\t' + str(c) + '\t' + str(cmo3) + '\t' + '\t'.join(
                                [str(global_codons[x][2]) for x in sorted(global_codons.keys())]) + '\n')
                            actghandle.write('global_actg\t' + str(c) + '\t' + str(cmo3) + '\t' + '\t'.join(
                                [str(global_actg[x][1]) for x in sorted(global_actg.keys())]) + '\n')
                            firsthandle.write('global_first\t' + str(c) + '\t' + str(cmo3) + '\t' + '\t'.join(
                                [str(global_first[x][1]) for x in sorted(global_first.keys())]) + '\n')
                            secondhandle.write('global_second\t' + str(c) + '\t' + str(cmo3) + '\t' + '\t'.join(
                                [str(global_second[x][1]) for x in sorted(global_second.keys())]) + '\n')
                            thirdhandle.write('global_third\t' + str(c) + '\t' + str(cmo3) + '\t' + '\t'.join(
                                [str(global_third[x][1]) for x in sorted(global_third.keys())]) + '\n')
                            if args.enc is not None:
                                enchandle.close()
    print('finished writing')
    if args.r:
        print('\n'.join(ids_mo3))



# Wright 1990 eqautions
def calc_fa(tmp):
    counts = [x[2] for x in tmp]
    na = sum(counts)
    return [((na * sum([(x / float(na)) ** 2 for x in counts])) - 1) / (na - 1), sum(counts)]


def calc_frc(fa, nrc):
    return 1 / float(nrc) * float(sum(fa))


# Latex function
# \widehat{N}_{c} = 2 + GC_{(3)} + (\frac{29}{GC^{2}_{(3)} + (1 - GC^{2}_{(3)})})
def calc_eq4wright(fgc):
    return 2 + float(fgc) + (29 / ((fgc ** 2) + ((1 - fgc) ** 2)))


# Sun 2012 equations
# FCF according to equation (1) without pseudocounts for eq1Sun


def eq1sun():


# FCF according to equation (2) without pseudocounts for eq2Sun
# F_{CF} = \sum\limits_{i=1}^{m} (\frac {n_{i}}{n})^{2}
def calc_fcf(tmp):
    counts = [x[2] for x in tmp]
    na = sum(counts)
    if na == 0:
        return 'NA', 0
    return sum([(x / float(na)) ** 2 for x in counts]), sum(counts)


# FCF according to equation (3) using pseudocounts for eq5Sun
# F_{CF} = \sum\limits_{i=1}^{m} (\frac {n_{i}+1}{n+m})^{2}
def calc_fcf_(tmp):
    counts = [x[2] for x in tmp]
    pseudocounts = [x + 1 for x in counts]
    na = sum(pseudocounts)
    return sum([(x / float(na)) ** 2 for x in pseudocounts]), sum(pseudocounts)


def calc_eq2sun(codoncounts, six2fourtwo):
    six_fcf = []
    six_nrc = []
    four_fcf = []
    four_nrc = []
    three_fcf = []
    three_nrc = []
    two_fcf = []
    two_nrc = []
    one_fcf = []
    one_nrc = []
    four_aa = []
    three_aa = []
    two_aa = []
    one_aa = []
    ncsix = None
    nc = None
    if not six2fourtwo:
        six_aa = ['S', 'R', 'L']
        four_aa = ['A', 'P', 'T', 'G', 'V']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y']
        one_aa = ['M', 'W']
        for aa in six_aa:
            tmp_six = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
            six = calc_fcf(tmp_six)
            six_fcf.append(six[0])
            six_nrc.append(six[1])
        # setting codons with no counts to average value 1/6
        six_fcf = [float(1) / float(6) if x == 'NA' else x for x in six_fcf]
        ncsix = len(six_fcf) / (sum(six_fcf) / len(six_fcf))
    if six2fourtwo:
        four_aa = ['A', 'P', 'T', 'G', 'V', 'S', 'R', 'L']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y', 'S_', 'R_', 'L_']
        one_aa = ['M', 'W']
    for aa in four_aa:
        tmp_four = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        four = calc_fcf(tmp_four)
        four_fcf.append(four[0])
        four_nrc.append(four[1])
    # setting codons with no counts to average value 1/4
    four_fcf = [float(1) / float(4) if x == 'NA' else x for x in four_fcf]
    ncfour = len(four_fcf) / (sum(four_fcf) / len(four_fcf))
    for aa in three_aa:
        tmp_three = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        three = calc_fcf(tmp_three)
        three_fcf.append(three[0])
        three_nrc.append(three[1])
    # setting codons with no counts to average value 1/3
    three_fcf = [float(1) / float(3) if x == 'NA' else x for x in three_fcf]
    ncthree = len(three_fcf) / (sum(three_fcf) / len(three_fcf))
    for aa in two_aa:
        tmp_two = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        two = calc_fcf(tmp_two)
        two_fcf.append(two[0])
        two_nrc.append(two[1])
    # setting codons with no counts to average value 1/2
    two_fcf = [float(1) / float(2) if x == 'NA' else x for x in two_fcf]
    nctwo = len(two_fcf) / (sum(two_fcf) / len(two_fcf))
    for aa in one_aa:
        tmp_one = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        one = calc_fcf(tmp_one)
        one_fcf.append(one[0])
        one_nrc.append(one[1])
    # setting codons with no counts to average value 1/1
    one_fcf = [float(1) / float(1) if x == 'NA' else x for x in one_fcf]
    ncone = len(one_fcf) / (sum(one_fcf) / len(one_fcf))
    if not six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour + ncsix
    if six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour
    return nc


# N_{c} = \frac {K_1 \sum\limits_{j}^{K_1} n_{j}}{\sum\limits_{j=1}^{K_1}(n_{j}F_{CF.j})}
# + \frac {K_2 \sum\limits_{j}^{K_2} n_{j}}{\sum\limits_{j=1}^{K_2}(n_{j}F_{CF.j})}
# + \frac {K_3 \sum\limits_{j}^{K_3} n_{j}}{\sum\limits_{j=1}^{K_3}(n_{j}F_{CF.j})}
# + \frac {K_4 \sum\limits_{j}^{K_4} n_{j}}{\sum\limits_{j=1}^{K_4}(n_{j}F_{CF.j})}
# + (\frac {K_6 \sum\limits_{j}^{K_6} n_{j}}{\sum\limits_{j=1}^{K_6}(n_{j}F_{CF.j})})
def calc_eq5sun(codoncounts, six2fourtwo):
    six_fcf_nrc = []
    four_fcf_nrc = []
    three_fcf_nrc = []
    two_fcf_nrc = []
    one_fcf_nrc = []
    four_aa = []
    three_aa = []
    two_aa = []
    one_aa = []
    ncsix = None
    nc = None
    if not six2fourtwo:
        six_aa = ['S', 'R', 'L']
        four_aa = ['A', 'P', 'T', 'G', 'V']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y']
        one_aa = ['M', 'W']
        for aa in six_aa:
            tmp_six = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
            six = calc_fcf_(tmp_six)
            six_fcf_nrc.append(six)
        ncsix = len(six_fcf_nrc) / (sum([x[0] * x[1] for x in six_fcf_nrc]) / (sum([x[1] for x in six_fcf_nrc])))
    if six2fourtwo:
        four_aa = ['A', 'P', 'T', 'G', 'V', 'S', 'R', 'L']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y', 'S_', 'R_', 'L_']
        one_aa = ['M', 'W']
    for aa in four_aa:
        tmp_four = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        four = calc_fcf_(tmp_four)
        four_fcf_nrc.append(four)
    ncfour = len(four_fcf_nrc) / (sum([x[0] * x[1] for x in four_fcf_nrc]) / (sum([x[1] for x in four_fcf_nrc])))
    for aa in three_aa:
        tmp_three = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        three = calc_fcf_(tmp_three)
        three_fcf_nrc.append(three)
    ncthree = len(three_fcf_nrc) / (sum([x[0] * x[1] for x in three_fcf_nrc]) / (sum([x[1] for x in three_fcf_nrc])))
    for aa in two_aa:
        tmp_two = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        two = calc_fcf_(tmp_two)
        two_fcf_nrc.append(two)
    nctwo = len(two_fcf_nrc) / (sum([x[0] * x[1] for x in two_fcf_nrc]) / (sum([x[1] for x in two_fcf_nrc])))
    for aa in one_aa:
        tmp_one = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        one = calc_fcf_(tmp_one)
        one_fcf_nrc.append(one)
    ncone = len(one_fcf_nrc) / (sum([x[0] * x[1] for x in one_fcf_nrc]) / (sum([x[1] for x in one_fcf_nrc])))
    if not six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour + ncsix
    if six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour
    return nc


# Relative Synonymous Codon Usage
# RSCU_{ij} = \frac {x_{ij}} {\frac {1}{n_{i}} \sum\limits_{j=1}^{n_{i}} x_{ij}}
def calc_rscu(codoncounts, codontable):
    rscutable = codontable()
    for codon in rscutable.keys():
        aa = rscutable[codon][0]
        tmp_aa = [[x, codoncounts[x][2]] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        counts = [x[1] for x in tmp_aa]
        na = sum(counts)
        if na == 0:
            rscutable[codon][2] = 'NA'
        else:
            tmp_codon_na = [x[1] for x in tmp_aa if x[0] == codon][0]
            tmp_codon_freq = tmp_codon_na / float(na)
            tmp_codon_rscu = tmp_codon_freq * len(counts)
            rscutable[codon][2] = tmp_codon_rscu
    return rscutable


if __name__ == '__main__':
    main()
