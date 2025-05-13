#!/usr/bin/python
# -*- coding: UTF-8 -*-

'''
Author: Krisian Ullrich
date: June 2016
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2016 Kristian Ullrich

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
'''


import sys
import os
import numpy as np
import copy
import argparse


def get_introns(g, tr, tr_strand, feature_array):
    if len(feature_array)==0:
        return []
    feature_starts = [int(x[3]) for x in feature_array]
    feature_type = feature_array[0][2].split('_')[0]
    if tr_strand == "+":
        feature_starts.sort(reverse=False)
    if tr_strand == "-":
        feature_starts.sort(reverse=True)
    feature_array_sorted = []
    for s in feature_starts:
        feature_array_sorted.append([x[:] for x in feature_array if int(x[3]) == s][0])
    if tr_strand=="+":
        intron_array = copy.deepcopy(feature_array_sorted)
        intron_array = intron_array[:len(intron_array)-1][:]
        for i in range(0,len(intron_array)):
            intron_array[i][3]=str(int(feature_array_sorted[i][4])+1)
            intron_array[i][4]=str(int(feature_array_sorted[i+1][3])-1)
            intron_array[i][7]="."
            intron_array[i][2]=feature_type+"Intron_"+intron_array[i][2]
            intron_array[i][8] = intron_array[i][8]+" "+feature_type+"Intron_"+'number "'+str(i+1)+'";'
        return(intron_array)
    if tr_strand=="-":
        intron_array = copy.deepcopy(feature_array_sorted)
        intron_array = intron_array[1:][:]
        for i in range(0,len(intron_array)):
            intron_array[i][3]=str(int(feature_array_sorted[i+1][4])+1)
            intron_array[i][4]=str(int(feature_array_sorted[i][3])-1)
            intron_array[i][7]="."
            intron_array[i][2]=feature_type+"Intron_"+intron_array[i][2]
            intron_array[i][8] = intron_array[i][8]+" "+feature_type+"Intron_"+'number "'+str(i+1)+'";'
        return(intron_array)


def main():
    '''
    Add Intron information to ENSEMBL GTF files.
    '''
    parser = argparse.ArgumentParser(usage='%(prog)s [options]',
                                     description='Given a ENSEMBL GTF file intron information will be added to output.')
    parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-i', help='specify input gtf file')
    parser.add_argument('-o', help='specify output gtf file')
    args = parser.parse_args()
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input gtf file')
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output gtf file')
    print(args)
    gtf={}
    ###fill gtf with genes
    with open(args.i,'rU') as gtf_handle:
        for lines in gtf_handle:
            _seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr = lines.strip().split('\t')
            if _type == 'gene':
                gene_id = _attr.split(';')[0].split()[1].replace('"','')
                gene_biotype = _attr.split('gene_biotype "')[1].split('";')[0]
                if gene_id in gtf:
                    print('duplicated gene_id: %s' % (gene_id))
                if gene_id not in gtf:
                    gtf[gene_id]={}
                    gtf[gene_id]['gene']={}
                    gtf[gene_id]['transcript']={}
                    gtf[gene_id]['gene'][gene_id]=[]
                    _type = _type+'_'+gene_biotype
                    gtf[gene_id]['gene'][gene_id].append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr])
    
    ###fill gtf with transcript
    with open(args.i,'rU') as gtf_handle:
        for lines in gtf_handle:
            _seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr = lines.strip().split('\t')
            if _type == 'transcript':
                gene_id = _attr.split(';')[0].split()[1].replace('"','')
                transcript_id = _attr.split('transcript_id ')[1].split(';')[0].replace('"','')
                transcript_biotype = _attr.split('transcript_biotype "')[1].split('";')[0]
                if transcript_id in gtf[gene_id]['transcript']:
                    print('duplicated transcript_id: %s' % (transcript_id))
                if transcript_id not in gtf[gene_id]['transcript']:
                    gtf[gene_id]['transcript'][transcript_id]={}
                    gtf[gene_id]['transcript'][transcript_id]['transcript']=[]
                    gtf[gene_id]['transcript'][transcript_id]['exon']=[]
                    gtf[gene_id]['transcript'][transcript_id]['exonIntron']=[]
                    gtf[gene_id]['transcript'][transcript_id]['five_prime_utr']=[]
                    gtf[gene_id]['transcript'][transcript_id]['five_prime_utrIntron']=[]
                    gtf[gene_id]['transcript'][transcript_id]['CDS']=[]
                    gtf[gene_id]['transcript'][transcript_id]['CDSIntron']=[]
                    gtf[gene_id]['transcript'][transcript_id]['stop_codon']=[]
                    gtf[gene_id]['transcript'][transcript_id]['stop_codonIntron']=[]
                    gtf[gene_id]['transcript'][transcript_id]['three_prime_utr']=[]
                    gtf[gene_id]['transcript'][transcript_id]['three_prime_utrIntron']=[]
                    _type = _type+'_'+transcript_biotype
                    gtf[gene_id]['transcript'][transcript_id]['transcript'].append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr])
    
    ###fill transcripts
    with open(args.i,'rU') as gtf_handle:
        for lines in gtf_handle:
            _seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr = lines.strip().split('\t')
            if _type == 'exon':
                gene_id = _attr.split(';')[0].split()[1].replace('"','')
                transcript_id = _attr.split('transcript_id ')[1].split(';')[0].replace('"','')
                transcript_biotype = _attr.split('transcript_biotype "')[1].split('";')[0]
                _type = _type+'_'+transcript_biotype
                gtf[gene_id]['transcript'][transcript_id]['exon'].append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr])
            if _type == 'five_prime_utr':
                gene_id = _attr.split(';')[0].split()[1].replace('"','')
                transcript_id = _attr.split('transcript_id ')[1].split(';')[0].replace('"','')
                transcript_biotype = _attr.split('transcript_biotype "')[1].split('";')[0]
                _type = _type+'_'+transcript_biotype
                gtf[gene_id]['transcript'][transcript_id]['five_prime_utr'].append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr])
            if _type == 'CDS':
                gene_id = _attr.split(';')[0].split()[1].replace('"','')
                transcript_id = _attr.split('transcript_id ')[1].split(';')[0].replace('"','')
                transcript_biotype = _attr.split('transcript_biotype "')[1].split('";')[0]
                _type = _type+'_'+transcript_biotype
                gtf[gene_id]['transcript'][transcript_id]['CDS'].append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr])
            if _type == 'stop_codon':
                gene_id = _attr.split(';')[0].split()[1].replace('"','')
                transcript_id = _attr.split('transcript_id ')[1].split(';')[0].replace('"','')
                transcript_biotype = _attr.split('transcript_biotype "')[1].split('";')[0]
                _type = _type+'_'+transcript_biotype
                gtf[gene_id]['transcript'][transcript_id]['stop_codon'].append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr])
            if _type == 'three_prime_utr':
                gene_id = _attr.split(';')[0].split()[1].replace('"','')
                transcript_id = _attr.split('transcript_id ')[1].split(';')[0].replace('"','')
                transcript_biotype = _attr.split('transcript_biotype "')[1].split('";')[0]
                _type = _type+'_'+transcript_biotype
                gtf[gene_id]['transcript'][transcript_id]['three_prime_utr'].append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr])
    
    ###get_introns
    for g in gtf.keys():
        for tr in gtf[g]['transcript'].keys():
            _seqid, _source, _type, _start, _end, _score, _strand, _phase, _attr = gtf[g]['transcript'][tr]['transcript'][0]
            gtf[g]['transcript'][tr]['exonIntron'] = get_introns(g, tr, _strand, gtf[g]['transcript'][tr]['exon'])
            gtf[g]['transcript'][tr]['five_prime_utrIntron'] = get_introns(g, tr, _strand, gtf[g]['transcript'][tr]['five_prime_utr'])
            gtf[g]['transcript'][tr]['CDSIntron'] = get_introns(g, tr, _strand, gtf[g]['transcript'][tr]['CDS'])
            gtf[g]['transcript'][tr]['stop_codonIntron'] = get_introns(g, tr, _strand, gtf[g]['transcript'][tr]['stop_codon'])
            gtf[g]['transcript'][tr]['three_prime_utrIntron'] = get_introns(g, tr, _strand, gtf[g]['transcript'][tr]['three_prime_utr'])

    ###write to output
    with open(args.o,'w') as outhandle:
        for g in gtf.keys():
            outhandle.write('\t'.join(gtf[g]['gene'][g][0]))
            outhandle.write('\n')
            for tr in gtf[g]['transcript'].keys():
                outhandle.write('\t'.join(gtf[g]['transcript'][tr]['transcript'][0]))
                outhandle.write('\n')
                for e in gtf[g]['transcript'][tr]['exon']:
                    outhandle.write('\t'.join(e))
                    outhandle.write('\n')
                for eI in gtf[g]['transcript'][tr]['exonIntron']:
                    outhandle.write('\t'.join(eI))
                    outhandle.write('\n')
                for fpu in gtf[g]['transcript'][tr]['five_prime_utr']:
                    outhandle.write('\t'.join(fpu))
                    outhandle.write('\n')
                for fpuI in gtf[g]['transcript'][tr]['five_prime_utrIntron']:
                    outhandle.write('\t'.join(fpuI))
                    outhandle.write('\n')
                for cds in gtf[g]['transcript'][tr]['CDS']:
                    outhandle.write('\t'.join(cds))
                    outhandle.write('\n')
                for cdsI in gtf[g]['transcript'][tr]['CDSIntron']:
                    outhandle.write('\t'.join(cdsI))
                    outhandle.write('\n')
                for st in gtf[g]['transcript'][tr]['stop_codon']:
                    outhandle.write('\t'.join(st))
                    outhandle.write('\n')
                for stI in gtf[g]['transcript'][tr]['stop_codonIntron']:
                    outhandle.write('\t'.join(stI))
                    outhandle.write('\n')
                for tpu in gtf[g]['transcript'][tr]['three_prime_utr']:
                    outhandle.write('\t'.join(tpu))
                    outhandle.write('\n')
                for tpuI in gtf[g]['transcript'][tr]['three_prime_utrIntron']:
                    outhandle.write('\t'.join(tpuI))
                    outhandle.write('\n')


if __name__ == '__main__':
    main()

