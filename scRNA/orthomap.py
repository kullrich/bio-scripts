#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian Ullrich
date: October 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import os
import sys
import time
import argparse
import itertools
import pandas as pd
from ete3 import NCBITaxa


def orthomap(args, parser, subparsers):
    """
    orthomap

    :param args:
    :param parser:
    :param subparsers:
    :return:
    """
    if not args.qname:
        subparsers.choices['orthomap'].print_help()
        print('\nError <-qname>: Please specify query species name in orthofinder and taxid')
        sys.exit()
    if not args.qt:
        subparsers.choices['orthomap'].print_help()
        print('\nError <-qt>: Please specify query species taxid')
        sys.exit()
    if not args.sl:
        subparsers.choices['orthomap'].print_help()
        print('\nError <-sl>: Please specify species list as <orthofinder name><tab><species taxid>')
        sys.exit()
    if not args.oc:
        subparsers.choices['orthomap'].print_help()
        print('\nError <-oc>: Please specify orthofinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)')
        sys.exit()
    if not args.og:
        subparsers.choices['orthomap'].print_help()
        print('\nError <-og>: Please specify orthofinder <Orthogroups.tsv> (see Orthogroups directory)')
        sys.exit()
    print(args)
    ncbi = NCBITaxa()
    args.q = None
    query_lineage = ncbi.get_lineage(args.qt)
    query_lineage_names_dict = ncbi.get_taxid_translator(query_lineage)
    query_lineage_names = pd.DataFrame([(x, y, query_lineage_names_dict[y]) for x,y in enumerate(query_lineage)])
    query_lineage_names.columns = ['PSnum', 'PStaxID', 'PSname']
    species_list = pd.read_csv(args.sl, sep='\t', header=None)
    species_list.columns = ['species', 'taxID']
    species_list['lineage']=species_list.apply(lambda x: ncbi.get_lineage(x[1]), axis=1)
    species_list['youngest_common']=[get_youngest_common(query_lineage, x) for x in species_list.lineage]
    species_list['youngest_name']=[list(x.values())[0] for x in [ncbi.get_taxid_translator([x]) for x in list(species_list.youngest_common)]]
    print(args.qname)
    print(args.qt)
    print(species_list)
    oc_og_dict = {}
    with open(args.oc, 'r') as oc_lines:
        oc_species = next(oc_lines).strip().split('\t')
        oc_qidx = [x for x,y in enumerate(oc_species) if y==args.qname]
        if(len(oc_qidx)==0):
            print('\nError <-qname>: query species name not in orthofinder results, please check spelling\ne.g. <head -1 Orthogroups.GeneCounts.tsv>')
            sys.exit()
        for oc_line in oc_lines:
            oc_og = oc_line.strip().split('\t')
            if(int(oc_og[oc_qidx[0]])==0):
                continue
            if(int(oc_og[oc_qidx[0]])>0):
                oc_og_hits = [oc_species[x+1] for x,y in enumerate(oc_og[1::][::-1][1::][::-1]) if int(y)>0]
                # get list of youngest common between query and all other species
                oc_og_hits_youngest_common = list(species_list.youngest_common[[x for x,y in enumerate(species_list.species) if y in oc_og_hits]])
                # evaluate all youngest common nodes to retain the oldest of them and assign as the orthogroup ancestral state (gene age)
                if(len(oc_og_hits_youngest_common)>0):
                    oc_og_oldest_common = get_oldest_common(query_lineage, oc_og_hits_youngest_common)
                    oc_og_dict[oc_og[0]] = oc_og_oldest_common
    with open(args.out, 'w') as outhandle:
        outhandle.write('gene\tPSnum\tPStaxID\tPSname\n')
        with open(args.og, 'r') as og_lines:
            og_species = next(og_lines).strip().split('\t')
            og_qidx = [x for x,y in enumerate(og_species) if y==args.qname]
            if(len(oc_qidx)==0):
                print('\nError <-qname>: query species name not in orthofinder results, please check spelling\ne.g. <head -1 Orthogroups.tsv>')
                sys.exit()
            for og_line in og_lines:
                og_og = og_line.strip().split('\t')
                if(og_og[0] not in oc_og_dict):
                    continue
                else:
                    og_ps = query_lineage_names[query_lineage_names['PStaxID']==oc_og_dict[og_og[0]]].values.tolist()[0]
                    og_ps = '\t'.join([str(x) for x in og_ps])
                    qgenes = [outhandle.write(x.replace(' ','') + '\t' + og_ps + '\n') for x in og_og[og_qidx[0]].split(',')]


def get_youngest_common(ql, tl):
    return [x for x in tl if x in ql][-1]


def get_oldest_common(ql, tl):
    return ql[min([x for x,y in enumerate(ql) if y in tl])]


def get_qtid(args, ncbi):
    if args.q and args.qt:
        taxid2name = ncbi.get_taxid_translator([int(args.qt)])
        name2taxid = ncbi.get_name_translator([taxid2name[int(args.qt)]])
        qname = list(name2taxid.keys())[0]
        qtid = list(taxid2name.keys())[0]
    if not args.q and args.qt:
        taxid2name = ncbi.get_taxid_translator([int(args.qt)])
        name2taxid = ncbi.get_name_translator([taxid2name[int(args.qt)]])
        qname = list(name2taxid.keys())[0]
        qtid = list(taxid2name.keys())[0]
    if args.q and not args.qt:
        name2taxid = ncbi.get_name_translator([args.q])
        taxid2name = ncbi.get_taxid_translator([name2taxid[args.q][0]])
        qname = list(name2taxid.keys())[0]
        qtid = list(taxid2name.keys())[0]
    qlineage = ncbi.get_lineage(qtid)
    qlineagenames = ncbi.get_taxid_translator(qlineage)
    qlineagezip = [(a,b) for a,b, in zip(qlineage,[qlineagenames[x] for x in qlineage])]
    qlineagerev = qlineage[::-1]
    if qlineage[2] == 2:
        qk = 'Bacteria'
    if qlineage[2] == 2157:
        qk = 'Archea'
    if qlineage[2] == 2759:
        qk = 'Eukaryota'
    return [qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk]


def ncbitax(args, parser, subparsers):
    """
    update local ncbi taxonomy database

    :param args:
    :param parser:
    :param subparsers:
    :return:
    """
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()


def qlin(args, parser, subparsers):
    """
    get query lineage based on ncbi taxonomy

    :param args:
    :param parser:
    :param subparsers:
    :return:
    """
    if not args.q and not args.qt:
        subparsers.choices['qlin'].print_help()
        print('\nError <-q> <-qt>: Please specify query species name or taxid')
        sys.exit()
    if args.q and args.qt:
        subparsers.choices['qlin'].print_help()
        print('\nWarning: Since both query species name and taxid are given taxid is used')
        sys.exit()
    print(args)
    ncbi = NCBITaxa()
    qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk = get_qtid(args, ncbi)
    print("query name: %s" % qname)
    print("query taxid: %s" % str(qtid))
    print("query kingdom: %s" % qk)
    print("query lineage names: \n%s" % str([qlineagenames[x] + "(" + str(x) + ")" for x in qlineage]))
    print("query lineage: \n%s" % str(qlineage))


def subparser(subparsers):
    # orthomap; parser
    orthomap_example = '''example:
    
    #
    python orthomap.py orthomap -qname  -qt 10090 -sl -oc -og
    '''
    parser_orthomap = subparsers.add_parser('orthomap', help='extract orthomap (help: <orthomap -h>)',
        description='extract orthomap from orthofinder output for query species',
        epilog=orthomap_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_orthomap.add_argument('-qname', help='query species name in orthofinder (see column names of  <Orthogroups.tsv>)')
    parser_orthomap.add_argument('-qt', help='query species taxid (e.g. use <orthomap qlin -h> to get taxid)')
    parser_orthomap.add_argument('-sl', help='species list as <orthofinder name><tab><species taxid> (only samples in this list will be processed)')
    parser_orthomap.add_argument('-oc', help='specify orthofinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)')
    parser_orthomap.add_argument('-og', help='specify orthofinder <Orthogroups.tsv> (see Orthogroups directory)')
    parser_orthomap.add_argument('-out', help='specify output file <orthomap.tsv> (default: orthomap.tsv)', default='orthomap.tsv')
    parser_orthomap.set_defaults(func=orthomap)
    # ncbitax; parser
    ncbitax_example = '''example:
    
    #update ncbi taxonomy database
    python orthomap.py ncbitax
    '''
    parser_ncbitax = subparsers.add_parser('ncbitax', help='ncbi taxonomy database update (help: <ncbitax -h>)',
        description='update local ncbi taxonomy database',
        epilog=ncbitax_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_ncbitax.set_defaults(func=ncbitax)
    # qlin; parser
    qlin_example = '''example (Mus musculus; 10090):
    
    #get query lineage to be used with orthomap later on using query species taxid
    python orthomap.py qlin -qt 10090
    #using query species name
    python orthomap.py qlin -q "Mus musculus"
    '''
    parser_qlin = subparsers.add_parser('qlin', help='query lineage retrieval (help: <qlin -h>)',
        description='get query lineage based on ncbi taxonomy',
        epilog=qlin_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_qlin.add_argument('-q', help='query species name')
    parser_qlin.add_argument('-qt', help='query species taxid')
    parser_qlin.set_defaults(func=qlin)


def main():
    # top-level parser
    parser = argparse.ArgumentParser(prog='orthomap', usage='%(prog)s <sub-script> [options] [<arguments>...]',
    description='extract orthomap from orthofinder output for query species')
    subparsers = parser.add_subparsers(title='sub-scripts', description='valid sub-scripts', help='sub-scripts help',
    dest='cmd')
    # sub-level parser
    subparser(subparsers)
    # get args
    args = parser.parse_args()
    # call function
    try:
        args.func(args, parser, subparsers)
    except:
        sys.exit(0)


if __name__ == '__main__':
    main()
