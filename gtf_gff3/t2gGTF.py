#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian Ullrich
date: January 2022
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


import os
import sys
import argparse
import textwrap
import gzip


def parse_gtf(gtf, output, args):
    t2g={}
    t2p={}
    tc=0
    gc=0
    pc=0
    dc=0
    for lines in gtf:
        if len(lines)==0 or lines[0]=='#':
            continue
        line=lines.strip().split('\t')
        if line[2]=='transcript':
            infosplit=line[8].strip().split(';')
            gid=[x for x in infosplit if 'gene_id' in x]
            if len(gid)>0:
                if len(gid)==1:
                    gid=gid[0]
                    gid=gid.replace('gene_id','').replace(' ','').replace('"','')
                else:
                    print('duplicated gene_id field:\t'+lines)
            else:
                print('no gene_id field:\t'+lines)
                continue
            tid=[x for x in infosplit if 'transcript_id' in x]
            if len(tid)>0:
                if len(tid)==1:
                    tid=tid[0]
                    tid=tid.replace('transcript_id','').replace(' ','').replace('"','')
                else:
                    print('duplicated transcript_id field:\t'+lines)
            else:
                print('no transcript_id field:\t'+lines)
                continue
            if args.g:
                gname=[x for x in infosplit if 'gene_name' in x]
                if len(gname)>0:
                    if len(gname)==1:
                        gname=gname[0]
                        gname=gname.replace('gene_name','').replace(' ','').replace('"','')
                    else:
                        print('duplicated gene_name field:\t'+lines)
                else:
                    gname=''
            if args.b:
                gtype=[x for x in infosplit if 'gene_biotype' in x]
                if len(gtype)>0:
                    if len(gtype)==1:
                        gtype=gtype[0]
                        gtype=gtype.replace('gene_biotype','').replace(' ','').replace('"','')
                    else:
                        print('duplicated gene_biotype field:\t'+lines)
                else:
                    gtype=''
            if args.v:
                gv=[x for x in infosplit if 'gene_version' in x and 'havana_gene_version' not in x]
                if len(gv)>0:
                    if len(gv)==1:
                        gv=gv[0]
                        gv=gv.replace('gene_version','').replace(' ','').replace('"','')
                    else:
                        print('duplicated gene_version field:\t'+lines)
                else:
                    print('no gene_version field:\t'+lines)
                    continue
                gid=gid+'.'+gv
                tv=[x for x in infosplit if 'transcript_version' in x and 'havana_transcript_version' not in x]
                if len(tv)>0:
                    if len(tv)==1:
                        tv=tv[0]
                        tv=tv.replace('transcript_version','').replace(' ','').replace('"','')
                    else:
                        print('duplicated transcript_version field:\t'+lines)
                else:
                    print('no transcript_version field:\t'+lines)
                    continue
                tid=tid+'.'+tv
            if gid in t2g:
                if tid in t2g[gid]:
                    dc+=1
                    print('duplicated gid-tid: '+gid+' '+tid)
                    continue
                if tid not in t2g[gid]:
                    tc+=1
                    if args.g and args.b:
                        t2g[gid][tid]=[gid, tid, gname, gtype]
                    if args.g and not args.b:
                        t2g[gid][tid]=[gid, tid, gname]
                    if not args.g and not args.b:
                        t2g[gid][tid]=[gid, tid]
            if gid not in t2g:
                gc+=1
                tc+=1
                t2g[gid]={}
                if args.g and args.b:
                    t2g[gid][tid]=[gid, tid, gname, gtype]
                if args.g and not args.b:
                    t2g[gid][tid]=[gid, tid, gname]
                if not args.g and args.b:
                    t2g[gid][tid]=[gid, tid, gtype]
                if not args.g and not args.b:
                    t2g[gid][tid]=[gid, tid]
        if line[2]=='CDS':
            if args.p:
                infosplit=line[8].strip().split(';')
                tid=[x for x in infosplit if 'transcript_id' in x]
                if len(tid)>0:
                    if len(tid)==1:
                        tid=tid[0]
                        tid=tid.replace('transcript_id','').replace(' ','').replace('"','')
                    else:
                        print('duplicated transcript_id field:\t'+lines)
                else:
                    print('no transcript_id field:\t'+lines)
                    continue
                pid=[x for x in infosplit if 'protein_id' in x]
                if len(pid)>0:
                    if len(pid)==1:
                        pid=pid[0]
                        pid=pid.replace('protein_id','').replace(' ','').replace('"','')
                    else:
                        print('duplicated protein_id field:\t'+lines)
                else:
                    print('no protein_id field:\t'+lines)
                    continue
                if args.v:
                    tv=[x for x in infosplit if 'transcript_version' in x and 'havana_transcript_version' not in x]
                    if len(tv)>0:
                        if len(tv)==1:
                            tv=tv[0]
                            tv=tv.replace('transcript_version','').replace(' ','').replace('"','')
                        else:
                            print('duplicated transcript_version field:\t'+lines)
                    else:
                        print('no transcript_version field:\t'+lines)
                        continue
                    tid=tid+'.'+tv
                    pv=[x for x in infosplit if 'protein_version' in x]
                    if len(pv)>0:
                        if len(pv)==1:
                            pv=pv[0]
                            pv=pv.replace('protein_version','').replace(' ','').replace('"','')
                        else:
                            print('duplicated protein_version field:\t'+lines)
                    else:
                        print('no protein_version field:\t'+lines)
                        continue
                    pid=pid+'.'+pv
                if tid in t2p:
                    continue
                if tid not in t2p:
                    pc+=1
                    t2p[tid]=pid
    if args.p:
        for gidk in sorted(t2g.keys()):
            for tidk in sorted(t2g[gidk].keys()):
                pidk=''
                if tidk in t2p:
                    pidk=t2p[tidk]
                output.write('\t'.join(t2g[gidk][tidk])+'\t'+pidk+'\n')
    else:
        for gidk in sorted(t2g.keys()):
            for tidk in sorted(t2g[gidk].keys()):
                output.write('\t'.join(t2g[gidk][tidk])+'\n')
    if args.s:
        print(str(gc)+' gene_id found')
        print(str(tc)+' transcript_id found')
        print(str(tc)+' protein_id found')
        print(str(dc)+' duplicated')


def main():
    # top-level parser
    parser = argparse.ArgumentParser(prog='t2gGTF', usage='%(prog)s <sub-script> [options] [<arguments>...]',
                                     description='extracts transcript to gene table')
    parser.add_argument('-i', help='specify GTF input file')
    parser.add_argument('-o', help='specify output file [optional]')
    parser.add_argument('-g', help='specify if gene names should be appended if they exist', action='store_true')
    parser.add_argument('-b', help='specify if gene biotype should be appended if they exist', action='store_true')
    parser.add_argument('-p', help='specify if protein id should be appended if they exist', action='store_true')
    parser.add_argument('-v', help='specify if gene/transcript/protein version should be appended', action='store_true')
    parser.add_argument('-s', help='specify if summary should be printed', action='store_true')
    # get args
    args = parser.parse_args()
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify GTF input file')
    if args.i.endswith('gz'):
        gtf = gzip.open(args.i,'rt')
    else:
        gtf = open(args.i,'rt')
    if args.o:
        output = open(args.o, 'w')
    else:
        output = sys.stdout
    parse_gtf(gtf, output, args)
    gtf.close()
    output.close()


if __name__ == '__main__':
    main()
