#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian K Ullrich
date: February 2024
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2024 Kristian Ullrich

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
import subprocess
import pysam
import tempfile
import numpy as np
import pandas as pd
from io import BytesIO


def split_by(x, by, getpos):
    return float(str(x).split(by)[getpos])


def as_dist(df):
    keep = np.triu(np.ones(df.shape), k=1).astype('bool').reshape(df.size)
    return df.stack(future_stack=True)[keep]


def subseq(stats_individuals, range_row):
    subseq_dict = {}
    for ind_index, ind_row in stats_individuals.iterrows():
        subseq_dict[ind_row['stats_pop']+'-'+ind_row['ind']] = pysam.FastaFile(ind_row['filepath']).fetch(reference=range_row[0], start=int(range_row[1]), end=int(range_row[2]))
    return subseq_dict


def literaldist(args, subseq_file):
    if args.g:
        subseqdist = subprocess.run([args.l, '-b', '-i', '-g', '-z', str(args.z), subseq_file], capture_output=True)
    else:
        subseqdist = subprocess.run([args.l, '-b', '-i', subseq_file], capture_output=True)
    return subseqdist


def get_individuals(stats_type, row, popassignment):
    pop_df = None
    if stats_type == 'dx':
        pop_x = row['x'].split(',')
        pop_xdf = popassignment[popassignment['pop'].isin(pop_x)].copy()
        pop_xdf['stats_pop'] = 'x'
        pop_df = pop_xdf
    if stats_type == 'dxy':
        pop_x = row['x'].split(',')
        pop_xdf = popassignment[popassignment['pop'].isin(pop_x)].copy()
        pop_xdf['stats_pop'] = 'x'
        pop_y = row['y'].split(',')
        pop_ydf = popassignment[popassignment['pop'].isin(pop_y)].copy()
        pop_ydf['stats_pop'] = 'y'
        pop_df = pd.concat([pop_xdf, pop_ydf], ignore_index = True)
    if stats_type == 'dixy':
        pop_i = row['i'].split(',')
        pop_idf = popassignment[popassignment['pop'].isin(pop_i)].copy()
        pop_idf['stats_pop'] = 'i'
        pop_x = row['x'].split(',')
        pop_xdf = popassignment[popassignment['pop'].isin(pop_x)].copy()
        pop_xdf['stats_pop'] = 'x'
        pop_y = row['y'].split(',')
        pop_ydf = popassignment[popassignment['pop'].isin(pop_y)].copy()
        pop_ydf['stats_pop'] = 'y'
        pop_df = pd.concat([pop_idf, pop_xdf, pop_ydf], ignore_index = True)
    if stats_type == 'dxyo':
        pop_x = row['x'].split(',')
        pop_xdf = popassignment[popassignment['pop'].isin(pop_x)].copy()
        pop_xdf['stats_pop'] = 'x'
        pop_y = row['y'].split(',')
        pop_ydf = popassignment[popassignment['pop'].isin(pop_y)].copy()
        pop_ydf['stats_pop'] = 'y'
        pop_o = row['o'].split(',')
        pop_odf = popassignment[popassignment['pop'].isin(pop_o)].copy()
        pop_odf['stats_pop'] = 'o'
        pop_df = pd.concat([pop_xdf, pop_ydf, pop_odf], ignore_index = True)
    if stats_type == 'dixyo':
        pop_i = row['i'].split(',')
        pop_idf = popassignment[popassignment['pop'].isin(pop_i)].copy()
        pop_idf['stats_pop'] = 'i'
        pop_x = row['x'].split(',')
        pop_xdf = popassignment[popassignment['pop'].isin(pop_x)].copy()
        pop_xdf['stats_pop'] = 'x'
        pop_y = row['y'].split(',')
        pop_ydf = popassignment[popassignment['pop'].isin(pop_y)].copy()
        pop_ydf['stats_pop'] = 'y'
        pop_o = row['o'].split(',')
        pop_odf = popassignment[popassignment['pop'].isin(pop_o)].copy()
        pop_odf['stats_pop'] = 'o'
        pop_df = pd.concat([pop_idf, pop_xdf, pop_ydf, pop_odf], ignore_index = True)
    return pop_df


def get_stats_type(row):
    stats_type = None
    if row.isnull()['i'] and row.isnull()['x'] and row.isnull()['y'] and row.isnull()['o']:
        stats_type = None
        sys.exit('\nPlease provide file with population combinations to be analyzed')
    if row.isnull()['i'] and not row.isnull()['x'] and row.isnull()['y'] and row.isnull()['o']:
        stats_type = 'dx'
    if row.isnull()['i'] and not row.isnull()['x'] and not row.isnull()['y'] and row.isnull()['o']:
        stats_type = 'dxy'
    if not row.isnull()['i'] and not row.isnull()['x'] and not row.isnull()['y'] and row.isnull()['o']:
        stats_type = 'dixy'
    if row.isnull()['i'] and not row.isnull()['x'] and not row.isnull()['y'] and not row.isnull()['o']:
        stats_type = 'dxyo'
    if not row.isnull()['i'] and not row.isnull()['x'] and not row.isnull()['y'] and not row.isnull()['o']:
        stats_type = 'dixyo'
    return stats_type


def fasta2dstats(args, parser):
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease provide file with population assignment')
    if args.r is None:
        parser.print_help()
        sys.exit('\nPlease provide file with genomic ranges')
    if args.c is None:
        parser.print_help()
        sys.exit('\nPlease provide file with population combinations to be analyzed')
    if args.l is None:
        parser.print_help()
        sys.exit('\nPlease provide path to literal-dist binary (see https://github.com/kullrich/literal-dists)')
    outfile = open(args.o, 'w')
    header = ['chr', 'start', 'end', 'gaplength',
    'pop_i', 'pop_x', 'pop_y', 'pop_o',
    'len_pop_i', 'len_pop_x', 'len_pop_y', 'len_pop_o',
    'stats_type',
    'dMean_i', 'dSd_i', 'dMin_i', 'dMax_i', 'dSites_i', 'dSum_i',
    'dMean_x', 'dSd_x', 'dMin_x', 'dMax_x', 'dSites_x', 'dSum_x',
    'dMean_y', 'dSd_y', 'dMin_y', 'dMax_y', 'dSites_y', 'dSum_y',
    'dMean_o', 'dSd_o', 'dMin_o', 'dMax_o', 'dSites_o', 'dSum_o',
    'dMean_ix', 'dSd_ix', 'dMin_ix', 'dMax_ix', 'dSites_ix', 'dSum_ix', 'dMean_ix_i_x', 'dF2_ix', 'Gmin_ix',
    'dMean_iy', 'dSd_iy', 'dMin_iy', 'dMax_iy', 'dSites_iy', 'dSum_iy', 'dMean_iy_i_y', 'dF2_iy', 'Gmin_iy',
    'dMean_io', 'dSd_io', 'dMin_io', 'dMax_io', 'dSites_io', 'dSum_io', 'dMean_io_i_o', 'dF2_io', 'Gmin_io',
    'dMean_xy', 'dSd_xy', 'dMin_xy', 'dMax_xy', 'dSites_xy', 'dSum_xy', 'dMean_xy_x_y', 'dF2_xy', 'Gmin_xy',
    'dMean_xo', 'dSd_xo', 'dMin_xo', 'dMax_xo', 'dSites_xo', 'dSum_xo', 'dMean_xo_x_o', 'dF2_xo', 'Gmin_xo',
    'dMean_yo', 'dSd_yo', 'dMin_yo', 'dMax_yo', 'dSites_yo', 'dSum_yo', 'dMean_yo_y_o', 'dF2_yo', 'Gmin_yo',
    'RND_ixy', 'RNDmin_ixy',
    'D3_xy_iy', 'D3_iy_xy',
    'D3_xy_ix', 'D3_ix_xy',
    'D3_ix_iy', 'D3_iy_ix',
    'D3_ixy', 'D3_ixy_abs_min',
    'D3_xy_iy_', 'D3_xy_ix_', 'D3_ix_iy_', 'D3_ixy_',
    'D3_xy_iy__', 'D3_xy_ix__', 'D3_ix_iy__', 'D3_ixy__',
    'D3min_xy_iy', 'D3min_xy_ix', 'D3min_ix_iy', 'D3min_ixy',
    'D3min_xy_iy_', 'D3min_xy_ix_', 'D3min_ix_iy_', 'D3min_ixy_',
    'D3min_xy_iy__', 'D3min_xy_ix__', 'D3min_ix_iy__', 'D3min_ixy__',
    'D3min_xy_min_iy', 'D3min_xy_min_ix', 'D3min_ix_min_iy', 'D3min_ixy_min',
    'D3min_xy_min_iy_', 'D3min_xy_min_ix_', 'D3min_ix_min_iy_', 'D3min_ixy_min_',
    'D3min_xy_min_iy__', 'D3min_xy_min_ix__', 'D3min_ix_min_iy__', 'D3min_ixy_min__',
    'RND_xyo', 'RNDmin_xyo',
    'D3_yo_xo', 'D3_xo_yo',
    'D3_xo_xy', 'D3_xy_xo',
    'D3_yo_xy', 'D3_xy_yo',
    'D3_xyo', 'D3_xyo_abs_min',
    'D3_yo_xo_', 'D3_yo_xy_', 'D3_xo_xy_', 'D3_xyo_',
    'D3_yo_xo__', 'D3_yo_xy__', 'D3_xo_xy__', 'D3_xyo__',
    'D3min_yo_xo', 'D3min_yo_xy', 'D3min_xo_xy', 'D3min_xyo',
    'D3min_yo_xo_', 'D3min_yo_xy_', 'D3min_xo_xy_', 'D3min_xyo_',
    'D3min_yo_xo__', 'D3min_yo_xy__', 'D3min_xo_xy__', 'D3min_xyo__',
    'D3min_yo_min_xo', 'D3min_yo_min_xy', 'D3min_xo_min_xy', 'D3min_xyo_min',
    'D3min_yo_min_xo_', 'D3min_yo_min_xy_', 'D3min_xo_min_xy_', 'D3min_xyo_min_',
    'D3min_yo_min_xo__', 'D3min_yo_min_xy__', 'D3min_xo_min_xy__', 'D3min_xyo_min__']
    outfile.write('\t'.join(header)+'\n')
    popassignment = pd.read_csv(args.i, delimiter='\t')
    genomicranges = pd.read_csv(args.r, delimiter='\t', header=None)
    popcombinations = pd.read_csv(args.c, delimiter='\t', comment='#')
    for pop_index, pop_row in popcombinations.iterrows():
        stats_type = get_stats_type(pop_row)
        stats_individuals = get_individuals(stats_type, pop_row, popassignment)
        pop_i_pos = None
        pop_x_pos = None
        pop_y_pos = None
        pop_o_pos = None
        pop_i_name = None
        pop_x_name = None
        pop_y_name = None
        pop_o_name = None
        len_pop_i = 0
        len_pop_x = 0
        len_pop_y = 0
        len_pop_o = 0
        if stats_type == 'dx':
            pop_x_pos = stats_individuals.index[stats_individuals['stats_pop']=='x']
            len_pop_x = len(pop_x_pos)
            pop_x_name = pop_row['x']
        if stats_type == 'dxy':
            pop_x_pos = stats_individuals.index[stats_individuals['stats_pop']=='x']
            pop_y_pos = stats_individuals.index[stats_individuals['stats_pop']=='y']
            len_pop_x = len(pop_x_pos)
            len_pop_y = len(pop_y_pos)
            pop_x_name = pop_row['x']
            pop_y_name = pop_row['y']
        if stats_type == 'dixy':
            pop_i_pos = stats_individuals.index[stats_individuals['stats_pop']=='i']
            pop_x_pos = stats_individuals.index[stats_individuals['stats_pop']=='x']
            pop_y_pos = stats_individuals.index[stats_individuals['stats_pop']=='y']
            len_pop_i = len(pop_i_pos)
            len_pop_x = len(pop_x_pos)
            len_pop_y = len(pop_y_pos)
            pop_i_name = pop_row['i']
            pop_x_name = pop_row['x']
            pop_y_name = pop_row['y']
        if stats_type == 'dxyo':
            pop_x_pos = stats_individuals.index[stats_individuals['stats_pop']=='x']
            pop_y_pos = stats_individuals.index[stats_individuals['stats_pop']=='y']
            pop_o_pos = stats_individuals.index[stats_individuals['stats_pop']=='o']
            len_pop_x = len(pop_x_pos)
            len_pop_y = len(pop_y_pos)
            len_pop_o = len(pop_o_pos)
            pop_x_name = pop_row['x']
            pop_y_name = pop_row['y']
            pop_o_name = pop_row['o']
        if stats_type == 'dixyo':
            pop_i_pos = stats_individuals.index[stats_individuals['stats_pop']=='i']
            pop_x_pos = stats_individuals.index[stats_individuals['stats_pop']=='x']
            pop_y_pos = stats_individuals.index[stats_individuals['stats_pop']=='y']
            pop_o_pos = stats_individuals.index[stats_individuals['stats_pop']=='o']
            len_pop_i = len(pop_i_pos)
            len_pop_x = len(pop_x_pos)
            len_pop_y = len(pop_y_pos)
            len_pop_o = len(pop_o_pos)
            pop_i_name = pop_row['i']
            pop_x_name = pop_row['x']
            pop_y_name = pop_row['y']
            pop_o_name = pop_row['o']
        for range_index, range_row in genomicranges.iterrows():
            range_seq_dict = subseq(stats_individuals, range_row)
            tmpfile = tempfile.NamedTemporaryFile()
            with open(tmpfile.name, 'w') as f:
                for k,v in range_seq_dict.items():
                    f.write('>'+k+'\n')
                    f.write(v+'\n')
            range_dist = literaldist(args, tmpfile.name)
            range_dist_df = pd.read_csv(BytesIO(range_dist.stdout), delimiter='\t', index_col=0)
            range_dist_df_distances = range_dist_df.map(split_by, by='/', getpos=0)
            range_dist_df_usedsites = range_dist_df.map(split_by, by='/', getpos=1)
            range_dist_df_score = range_dist_df.map(split_by, by='/', getpos=2)
            range_gapsites = np.unique(range_dist_df.map(split_by, by='/', getpos=3))[0]
            range_stats = get_stats(pop_i_pos, pop_x_pos, pop_y_pos, pop_o_pos,
                len_pop_i, len_pop_x, len_pop_y, len_pop_o,
                pop_i_name, pop_x_name, pop_y_name, pop_o_name,
                stats_type, range_dist_df_distances, range_dist_df_usedsites)
            outfile.write('\t'.join([str(x) for x in [range_row[0], range_row[1], range_row[2], range_gapsites,
            pop_i_name, pop_x_name, pop_y_name, pop_o_name,
            len_pop_i, len_pop_x, len_pop_y, len_pop_o,
            stats_type]+range_stats])+'\n')
    outfile.close()


def get_stats(pop_i_pos, pop_x_pos, pop_y_pos, pop_o_pos,
    len_pop_i, len_pop_x, len_pop_y, len_pop_o,
    pop_i_name, pop_x_name, pop_y_name, pop_o_name,
    stats_type, range_dist_df_distances, range_dist_df_usedsites):
    #i
    dMean_i = np.nan
    dSd_i = np.nan
    dMin_i = np.nan
    dMax_i = np.nan
    dSites_i = np.nan
    dSum_i = np.nan
    #x
    dMean_x = np.nan
    dSd_x = np.nan
    dMin_x = np.nan
    dMax_x = np.nan
    dSites_x = np.nan
    dSum_x = np.nan
    #y
    dMean_y = np.nan
    dSd_y = np.nan
    dMin_y = np.nan
    dMax_y = np.nan
    dSites_y = np.nan
    dSum_y = np.nan
    #o
    dMean_o = np.nan
    dSd_o = np.nan
    dMin_o = np.nan
    dMax_o = np.nan
    dSites_o = np.nan
    dSum_o = np.nan
    #ix
    dMean_ix = np.nan
    dSd_ix = np.nan
    dMin_ix = np.nan
    dMax_ix = np.nan
    dSites_ix = np.nan
    dSum_ix = np.nan
    dMean_ix_i_x = np.nan # dMean_ix - (dMean_i/2) - (dMean_x/2)
    dF2_ix = np.nan # dMean_ix - ((dMean_i + dMean_x)/2)
    Gmin_ix = np.nan # dMin_ix / dMean_ix
    #iy
    dMean_iy = np.nan
    dSd_iy = np.nan
    dMin_iy = np.nan
    dMax_iy = np.nan
    dSites_iy = np.nan
    dSum_iy = np.nan
    dMean_iy_i_y = np.nan # dMean_iy - (dMean_i/2) - (dMean_y/2)
    dF2_iy = np.nan # dMean_iy - ((dMean_i + dMean_y)/2)
    Gmin_iy = np.nan # dMin_iy / dMean_iy
    #io
    dMean_io = np.nan
    dSd_io = np.nan
    dMin_io = np.nan
    dMax_io = np.nan
    dSites_io = np.nan
    dSum_io = np.nan
    dMean_io_i_o = np.nan # dMean_io - (dMean_i/2) - (dMean_o/2)
    dF2_io = np.nan # dMean_io - ((dMean_i + dMean_o)/2)
    Gmin_io = np.nan # dMin_io / dMean_io
    #xy
    dMean_xy = np.nan
    dSd_xy = np.nan
    dMin_xy = np.nan
    dMax_xy = np.nan
    dSites_xy = np.nan
    dSum_xy = np.nan
    dMean_xy_x_y = np.nan # dMean_xy - (dMean_x/2) - (dMean_y/2)
    dF2_xy = np.nan # dMean_xy - ((dMean_x + dMean_y)/2)
    Gmin_xy = np.nan # dMin_xy / dMean_xy
    #xo
    dMean_xo = np.nan
    dSd_xo = np.nan
    dMin_xo = np.nan
    dMax_xo = np.nan
    dSites_xo = np.nan
    dSum_xo = np.nan
    dMean_xo_x_o = np.nan # dMean_xo - (dMean_x/2) - (dMean_o/2)
    dF2_xo = np.nan # dMean_xo - ((dMean_x + dMean_o)/2)
    Gmin_xo = np.nan # dMin_xo / dMean_xo
    #yo
    dMean_yo = np.nan
    dSd_yo = np.nan
    dMin_yo = np.nan
    dMax_yo = np.nan
    dSites_yo = np.nan
    dSum_yo = np.nan
    dMean_yo_y_o = np.nan # dMean_yo - (dMean_y/2) - (dMean_o/2)
    dF2_yo = np.nan # dMean_yo - ((dMean_y + dMean_o)/2)
    Gmin_yo = np.nan # dMin_yo / dMean_yo
    #ixy
    RND_ixy = np.nan # dMean_xy - ((dMean_ix + dMean_iy)/2)
    RNDmin_ixy = np.nan # dMin_xy - ((dMean_ix + dMean_iy)/2)
    D3_xy_iy = np.nan # (dMean_xy - dMean_iy)/(dMean_xy + dMean_iy)
    D3_iy_xy = np.nan # (dMean_iy - dMean_xy)/(dMean_iy + dMean_xy)
    D3_xy_ix = np.nan # (dMean_xy - dMean_ix)/(dMean_xy + dMean_ix)
    D3_ix_xy = np.nan # (dMean_ix - dMean_xy)/(dMean_ix + dMean_xy)
    D3_ix_iy = np.nan # (dMean_ix - dMean_iy)/(dMean_ix + dMean_iy)
    D3_iy_ix = np.nan # (dMean_iy - dMean_ix)/(dMean_iy + dMean_ix)
    D3_ixy = np.nan
    D3_ixy_abs_min = np.nan
    D3_xy_iy_ = np.nan # (dMean_xy - dMean_iy)/(dMin_xy + dMin_iy)
    D3_xy_ix_ = np.nan # (dMean_xy - dMean_ix)/(dMin_xy + dMin_ix)
    D3_ix_iy_ = np.nan # (dMean_ix - dMean_iy)/(dMin_ix + dMin_iy)
    D3_ixy_ = np.nan
    D3_xy_iy__ = np.nan # (dMean_xy - dMean_iy)/(dMax_xy + dMax_iy)
    D3_xy_ix__ = np.nan # (dMean_xy - dMean_ix)/(dMax_xy + dMax_ix)
    D3_ix_iy__ = np.nan # (dMean_ix - dMean_iy)/(dMax_ix + dMax_iy)
    D3_ixy__ = np.nan
    D3min_xy_iy = np.nan # (dMin_xy - dMean_iy)/(dMean_xy + dMean_iy)
    D3min_xy_ix = np.nan # (dMin_xy - dMean_ix)/(dMean_xy + dMean_ix)
    D3min_ix_iy = np.nan # (dMin_ix - dMean_iy)/(dMean_ix + dMean_iy)
    D3min_ixy = np.nan
    D3min_xy_iy_ = np.nan # (dMin_xy - dMean_iy)/(dMin_xy + dMin_iy)
    D3min_xy_ix_ = np.nan # (dMin_xy - dMean_ix)/(dMin_xy + dMin_ix)
    D3min_ix_iy_ = np.nan # (dMin_ix - dMean_iy)/(dMin_ix + dMin_iy)
    D3min_ixy_ = np.nan
    D3min_xy_iy__ = np.nan # (dMin_xy - dMean_iy)/(dMax_xy + dMax_iy)
    D3min_xy_ix__ = np.nan # (dMin_xy - dMean_ix)/(dMax_xy + dMax_ix)
    D3min_ix_iy__ = np.nan # (dMin_ix - dMean_iy)/(dMax_ix + dMax_iy)
    D3min_ixy__ = np.nan
    D3min_xy_min_iy = np.nan # (dMin_xy - dMin_iy)/(dMean_xy + dMean_iy)
    D3min_xy_min_ix = np.nan # (dMin_xy - dMin_ix)/(dMean_xy + dMean_ix)
    D3min_ix_min_iy = np.nan # (dMin_ix - dMin_iy)/(dMean_ix + dMean_iy)
    D3min_ixy_min = np.nan
    D3min_xy_min_iy_ = np.nan # (dMin_xy - dMin_iy)/(dMin_xy + dMin_iy)
    D3min_xy_min_ix_ = np.nan # (dMin_xy - dMin_ix)/(dMin_xy + dMin_ix)
    D3min_ix_min_iy_ = np.nan # (dMin_ix - dMin_iy)/(dMin_ix + dMin_iy)
    D3min_ixy_min_ = np.nan
    D3min_xy_min_iy__ = np.nan # (dMin_xy - dMin_iy)/(dMax_xy + dMax_iy)
    D3min_xy_min_ix__ = np.nan # (dMin_xy - dMin_ix)/(dMax_xy + dMax_ix)
    D3min_ix_min_iy__ = np.nan # (dMin_ix - dMin_iy)/(dMax_ix + dMax_iy)
    D3min_ixy_min__ = np.nan
    #xyo
    RND_xyo = np.nan # dMean_xy / ((dMean_xo + dMean_yo)/2)
    RNDmin_xyo = np.nan # dMin_xy / ((dMean_xo + dMean_yo)/2)
    D3_yo_xo = np.nan # (dMean_yo - dMean_xo)/(dMean_yo + dMean_xo)
    D3_xo_yo = np.nan # (dMean_xo - dMean_yo)/(dMean_xo + dMean_yo)
    D3_yo_xy = np.nan # (dMean_yo - dMean_xy)/(dMean_yo + dMean_xy)
    D3_xy_yo = np.nan # (dMean_xy - dMean_yo)/(dMean_xy + dMean_yo)
    D3_xo_xy = np.nan # (dMean_xo - dMean_xy)/(dMean_xo + dMean_xy)
    D3_xy_xo = np.nan # (dMean_xy - dMean_xo)/(dMean_xy + dMean_xo)
    D3_xyo = np.nan
    D3_xyo_abs_min = np.nan
    D3_yo_xo_ = np.nan # (dMean_yo - dMean_xo)/(dMin_yo + dMin_xo)
    D3_yo_xy_ = np.nan # (dMean_yo - dMean_xy)/(dMin_yo + dMin_xy)
    D3_xo_xy_ = np.nan # (dMean_xo - dMean_xy)/(dMin_xo + dMin_xy)
    D3_xyo_ = np.nan
    D3_yo_xo__ = np.nan # (dMean_yo - dMean_xo)/(dMax_yo + dMax_xo)
    D3_yo_xy__ = np.nan # (dMean_yo - dMean_xy)/(dMax_yo + dMax_xy)
    D3_xo_xy__ = np.nan # (dMean_xo - dMean_xy)/(dMax_xo + dMax_xy)
    D3_xyo__ = np.nan
    D3min_yo_xo = np.nan # (dMin_yo - dMean_xo)/(dMean_yo + dMean_xo)
    D3min_yo_xy = np.nan # (dMin_yo - dMean_xy)/(dMean_yo + dMean_xy)
    D3min_xo_xy = np.nan # (dMin_xo - dMean_xy)/(dMean_xo + dMean_xy)
    D3min_xyo = np.nan
    D3min_yo_xo_ = np.nan # (dMin_yo - dMean_xo)/(dMin_yo + dMin_xo)
    D3min_yo_xy_ = np.nan # (dMin_yo - dMean_xy)/(dMin_yo + dMin_xy)
    D3min_xo_xy_ = np.nan # (dMin_xo - dMean_xy)/(dMin_xo + dMin_xy)
    D3min_xyo_ = np.nan
    D3min_yo_xo__ = np.nan # (dMin_yo - dMean_xo)/(dMax_yo + dMax_xo)
    D3min_yo_xy__ = np.nan # (dMin_yo - dMean_xy)/(dMax_yo + dMax_xy)
    D3min_xo_xy__ = np.nan # (dMin_xo - dMean_xy)/(dMax_xo + dMax_xy)
    D3min_xyo__ = np.nan
    D3min_yo_min_xo = np.nan # (dMin_yo - dMin_xo)/(dMean_yo + dMean_xo)
    D3min_yo_min_xy = np.nan # (dMin_yo - dMin_xy)/(dMean_yo + dMean_xy)
    D3min_xo_min_xy = np.nan # (dMin_xo - dMin_xy)/(dMean_xo + dMean_xy)
    D3min_xyo_min = np.nan
    D3min_yo_min_xo_ = np.nan # (dMin_yo - dMin_xo)/(dMin_yo + dMin_xo)
    D3min_yo_min_xy_ = np.nan # (dMin_yo - dMin_xy)/(dMin_yo + dMin_xy)
    D3min_xo_min_xy_ = np.nan # (dMin_xo - dMin_xy)/(dMin_xo + dMin_xy)
    D3min_xyo_min_ = np.nan
    D3min_yo_min_xo__ = np.nan # (dMin_yo - dMin_xo)/(dMax_yo + dMax_xo)
    D3min_yo_min_xy__ = np.nan # (dMin_yo - dMin_xy)/(dMax_yo + dMax_xy)
    D3min_xo_min_xy__ = np.nan # (dMin_xo - dMin_xy)/(dMax_xo + dMax_xy)
    D3min_xyo_min__ = np.nan
    #ixyo
    if stats_type == 'dx':
        gaplength = range_gapsites
        if len_pop_x > 1:
            dMean_x = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSd_x = np.std(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMin_x = np.min(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMax_x = np.max(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSites_x = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSum_x = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
        else:
            dMean_x = 0.0
            dSd_x = 0.0
            dMin_x = 0.0
            dMax_x = 0.0
            dSites_x = 0.0
            dSum_x = 0.0
    if stats_type == 'dxy':
        if len_pop_x > 1:
            dMean_x = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSd_x = np.std(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMin_x = np.min(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMax_x = np.max(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSites_x = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSum_x = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
        else:
            dMean_x = 0.0
            dSd_x = 0.0
            dMin_x = 0.0
            dMax_x = 0.0
            dSites_x = 0.0
            dSum_x = 0.0
        if len_pop_y > 1:
            dMean_y = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSd_y = np.std(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMin_y = np.min(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMax_y = np.max(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSites_y = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSum_y = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
        else:
            dMean_y = 0.0
            dSd_y = 0.0
            dMin_y = 0.0
            dMax_y = 0.0
            dSites_y = 0.0
            dSum_y = 0.0
        #xy
        dMean_xy = np.mean(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSd_xy = np.std(np.array(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)]))
        dMin_xy = np.min(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dMax_xy = np.max(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSites_xy = np.mean(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSum_xy = np.sum(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)].sum())
        dMean_xy_x_y = dMean_xy - (dMean_x/2) - (dMean_y/2)
        dF2_xy = dMean_xy - ((dMean_x + dMean_y)/2)
        Gmin_xy = dMin_xy / dMean_xy
    if stats_type == 'dixy':
        if len_pop_i > 1:
            dMean_i = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dSd_i = np.std(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dMin_i = np.min(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dMax_i = np.max(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dSites_i = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dSum_i = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_i_pos)]))
        else:
            dMean_i = 0.0
            dSd_i = 0.0
            dMin_i = 0.0
            dMax_i = 0.0
            dSites_i = 0.0
            dSum_i = 0.0
        if len_pop_x > 1:
            dMean_x = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSd_x = np.std(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMin_x = np.min(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMax_x = np.max(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSites_x = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSum_x = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
        else:
            dMean_x = 0.0
            dSd_x = 0.0
            dMin_x = 0.0
            dMax_x = 0.0
            dSites_x = 0.0
            dSum_x = 0.0
        if len_pop_y > 1:
            dMean_y = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSd_y = np.std(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMin_y = np.min(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMax_y = np.max(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSites_y = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSum_y = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_y_pos)]))
        else:
            dMean_y = 0.0
            dSd_y = 0.0
            dMin_y = 0.0
            dMax_y = 0.0
            dSites_y = 0.0
            dSum_y = 0.0
        #ix
        dMean_ix = np.mean(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)])
        dSd_ix = np.std(np.array(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)]))
        dMin_ix = np.min(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)])
        dMax_ix = np.max(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)])
        dSites_ix = np.mean(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_x_pos)])
        dSum_ix = np.sum(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_x_pos)].sum())
        dMean_ix_i_x = dMean_ix - (dMean_i/2) - (dMean_x/2)
        dF2_ix = dMean_ix - ((dMean_i + dMean_x)/2)
        Gmin_ix = dMin_ix / dMean_ix
        #iy
        dMean_iy = np.mean(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)])
        dSd_iy = np.std(np.array(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)]))
        dMin_iy = np.min(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)])
        dMax_iy = np.max(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)])
        dSites_iy = np.mean(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_y_pos)])
        dSum_iy = np.sum(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_y_pos)].sum())
        dMean_iy_i_y = dMean_iy - (dMean_i/2) - (dMean_y/2)
        dF2_iy = dMean_iy - ((dMean_i + dMean_y)/2)
        Gmin_iy = dMin_iy / dMean_iy
        #xy
        dMean_xy = np.mean(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSd_xy = np.std(np.array(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)]))
        dMin_xy = np.min(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dMax_xy = np.max(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSites_xy = np.mean(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSum_xy = np.sum(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)].sum())
        dMean_xy_x_y = dMean_xy - (dMean_x/2) - (dMean_y/2)
        dF2_xy = dMean_xy - ((dMean_x + dMean_y)/2)
        Gmin_xy = dMin_xy / dMean_xy
        #ixy
        RND_ixy = dMean_xy - ((dMean_ix + dMean_iy)/2)
        RNDmin_ixy = dMin_xy - ((dMean_ix + dMean_iy)/2)
        D3_xy_iy = (dMean_xy - dMean_iy)/(dMean_xy + dMean_iy)
        D3_iy_xy = (dMean_iy - dMean_xy)/(dMean_iy + dMean_xy)
        D3_xy_ix = (dMean_xy - dMean_ix)/(dMean_xy + dMean_ix)
        D3_ix_xy = (dMean_ix - dMean_xy)/(dMean_ix + dMean_xy)
        D3_ix_iy = (dMean_ix - dMean_iy)/(dMean_ix + dMean_iy)
        D3_iy_ix = (dMean_iy - dMean_ix)/(dMean_iy + dMean_ix)
        D3_ixy = [D3_xy_iy, D3_xy_ix, D3_ix_iy][np.argmin([abs(x) for x in [D3_xy_iy, D3_xy_ix, D3_ix_iy]])]
        D3_ixy_abs_min = ['((I,X),Y)', '((I,Y),X)', '((X,Y),I)'][np.argmin([abs(x) for x in [D3_xy_iy, D3_xy_ix, D3_ix_iy]])]
        D3_xy_iy_ = (dMean_xy - dMean_iy)/(dMin_xy + dMin_iy)
        D3_xy_ix_ = (dMean_xy - dMean_ix)/(dMin_xy + dMin_ix)
        D3_ix_iy_ = (dMean_ix - dMean_iy)/(dMin_ix + dMin_iy)
        D3_ixy_ = [D3_xy_iy_, D3_xy_ix_, D3_ix_iy_][np.argmin([abs(x) for x in [D3_xy_iy_, D3_xy_ix_, D3_ix_iy_]])]
        D3_xy_iy__ = (dMean_xy - dMean_iy)/(dMax_xy + dMax_iy)
        D3_xy_ix__ = (dMean_xy - dMean_ix)/(dMax_xy + dMax_ix)
        D3_ix_iy__ = (dMean_ix - dMean_iy)/(dMax_ix + dMax_iy)
        D3_ixy__ = [D3_xy_iy__, D3_xy_ix__, D3_ix_iy__][np.argmin([abs(x) for x in [D3_xy_iy__, D3_xy_ix__, D3_ix_iy__]])]
        D3min_xy_iy = (dMin_xy - dMean_iy)/(dMean_xy + dMean_iy)
        D3min_xy_ix = (dMin_xy - dMean_ix)/(dMean_xy + dMean_ix)
        D3min_ix_iy = (dMin_ix - dMean_iy)/(dMean_ix + dMean_iy)
        D3min_ixy = [D3min_xy_iy, D3min_xy_ix, D3min_ix_iy][np.argmin([abs(x) for x in [D3min_xy_iy, D3min_xy_ix, D3min_ix_iy]])]
        D3min_xy_iy_ = (dMin_xy - dMean_iy)/(dMin_xy + dMin_iy)
        D3min_xy_ix_ = (dMin_xy - dMean_ix)/(dMin_xy + dMin_ix)
        D3min_ix_iy_ = (dMin_ix - dMean_iy)/(dMin_ix + dMin_iy)
        D3min_ixy_ = [D3min_xy_iy_, D3min_xy_ix_, D3min_ix_iy_][np.argmin([abs(x) for x in [D3min_xy_iy_, D3min_xy_ix_, D3min_ix_iy_]])]
        D3min_xy_iy__ = (dMin_xy - dMean_iy)/(dMax_xy + dMax_iy)
        D3min_xy_ix__ = (dMin_xy - dMean_ix)/(dMax_xy + dMax_ix)
        D3min_ix_iy__ = (dMin_ix - dMean_iy)/(dMax_ix + dMax_iy)
        D3min_ixy__ = [D3min_xy_iy__, D3min_xy_ix__, D3min_ix_iy__][np.argmin([abs(x) for x in [D3min_xy_iy__, D3min_xy_ix__, D3min_ix_iy__]])]
        D3min_xy_min_iy = (dMin_xy - dMin_iy)/(dMean_xy + dMean_iy)
        D3min_xy_min_ix = (dMin_xy - dMin_ix)/(dMean_xy + dMean_ix)
        D3min_ix_min_iy = (dMin_ix - dMin_iy)/(dMean_ix + dMean_iy)
        D3min_ixy_min = [D3min_xy_min_iy, D3min_xy_min_ix, D3min_ix_min_iy][np.argmin([abs(x) for x in [D3min_xy_min_iy, D3min_xy_min_ix, D3min_ix_min_iy]])]
        D3min_xy_min_iy_ = (dMin_xy - dMin_iy)/(dMin_xy + dMin_iy)
        D3min_xy_min_ix_ = (dMin_xy - dMin_ix)/(dMin_xy + dMin_ix)
        D3min_ix_min_iy_ = (dMin_ix - dMin_iy)/(dMin_ix + dMin_iy)
        D3min_ixy_min_ = [D3min_xy_min_iy_, D3min_xy_min_ix_, D3min_ix_min_iy_][np.argmin([abs(x) for x in [D3min_xy_min_iy_, D3min_xy_min_ix_, D3min_ix_min_iy_]])]
        D3min_xy_min_iy__ = (dMin_xy - dMin_iy)/(dMax_xy + dMax_iy)
        D3min_xy_min_ix__ = (dMin_xy - dMin_ix)/(dMax_xy + dMax_ix)
        D3min_ix_min_iy__ = (dMin_ix - dMin_iy)/(dMax_ix + dMax_iy)
        D3min_ixy_min__ = [D3min_xy_min_iy__, D3min_xy_min_ix__, D3min_ix_min_iy__][np.argmin([abs(x) for x in [D3min_xy_min_iy__, D3min_xy_min_ix__, D3min_ix_min_iy__]])]
    if stats_type == 'dxyo':
        if len_pop_x > 1:
            dMean_x = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSd_x = np.std(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMin_x = np.min(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMax_x = np.max(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSites_x = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSum_x = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
        else:
            dMean_x = 0.0
            dSd_x = 0.0
            dMin_x = 0.0
            dMax_x = 0.0
            dSites_x = 0.0
            dSum_x = 0.0
        if len_pop_y > 1:
            dMean_y = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSd_y = np.std(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMin_y = np.min(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMax_y = np.max(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSites_y = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSum_y = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_y_pos)]))
        else:
            dMean_y = 0.0
            dSd_y = 0.0
            dMin_y = 0.0
            dMax_y = 0.0
            dSites_y = 0.0
            dSum_y = 0.0
        if len_pop_o > 1:
            dMean_o = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dSd_o = np.std(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dMin_o = np.min(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dMax_o = np.max(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dSites_o = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dSum_o = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_o_pos), list(pop_o_pos)]))
        else:
            dMean_o = 0.0
            dSd_o = 0.0
            dMin_o = 0.0
            dMax_o = 0.0
            dSites_o = 0.0
            dSum_o = 0.0
        #xy
        dMean_xy = np.mean(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSd_xy = np.std(np.array(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)]))
        dMin_xy = np.min(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dMax_xy = np.max(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSites_xy = np.mean(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSum_xy = np.sum(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)].sum())
        dMean_xy_x_y = dMean_xy - (dMean_x/2) - (dMean_y/2)
        dF2_xy = dMean_xy - ((dMean_x + dMean_y)/2)
        Gmin_xy = dMin_xy / dMean_xy
        #xo
        dMean_xo = np.mean(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)])
        dSd_xo = np.std(np.array(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)]))
        dMin_xo = np.min(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)])
        dMax_xo = np.max(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)])
        dSites_xo = np.mean(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_o_pos)])
        dSum_xo = np.sum(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_o_pos)].sum())
        dMean_xo_x_o = dMean_xo - (dMean_x/2) - (dMean_o/2)
        dF2_xo = dMean_xo - ((dMean_x + dMean_o)/2)
        Gmin_xo = dMin_xo / dMean_xo
        #yo
        dMean_yo = np.mean(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)])
        dSd_yo = np.std(np.array(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)]))
        dMin_yo = np.min(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)])
        dMax_yo = np.max(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)])
        dSites_yo = np.mean(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_o_pos)])
        dSum_yo = np.sum(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_o_pos)].sum())
        dMean_yo_y_o = dMean_yo - (dMean_y/2) - (dMean_o/2)
        dF2_yo = dMean_yo - ((dMean_y + dMean_o)/2)
        Gmin_yo = dMin_yo / dMean_yo
        #xyo
        RND_xyo = dMean_xy / ((dMean_xo + dMean_yo)/2)
        RNDmin_xyo = dMin_xy / ((dMean_xo + dMean_yo)/2)
        D3_yo_xo = (dMean_yo - dMean_xo)/(dMean_yo + dMean_xo)
        D3_xo_yo = (dMean_xo - dMean_yo)/(dMean_xo + dMean_yo)
        D3_yo_xy = (dMean_yo - dMean_xy)/(dMean_yo + dMean_xy)
        D3_xy_yo = (dMean_xy - dMean_yo)/(dMean_xy + dMean_yo)
        D3_xo_xy = (dMean_xo - dMean_xy)/(dMean_xo + dMean_xy)
        D3_xy_xo = (dMean_xy - dMean_xo)/(dMean_xy + dMean_xo)
        D3_xyo = [D3_yo_xo, D3_yo_xy, D3_xo_xy][np.argmin([abs(x) for x in [D3_yo_xo, D3_yo_xy, D3_xo_xy]])]
        D3_xyo_abs_min = ['((X,Y),O)', '((O,X),Y)', '((Y,O),X)'][np.argmin([abs(x) for x in [D3_yo_xo, D3_yo_xy, D3_xo_xy]])]
        D3_yo_xo_ = (dMean_yo - dMean_xo)/(dMin_yo + dMin_xo)
        D3_yo_xy_ = (dMean_yo - dMean_xy)/(dMin_yo + dMin_xy)
        D3_xo_xy_ = (dMean_xo - dMean_xy)/(dMin_xo + dMin_xy)
        D3_xyo_ = [D3_yo_xo_, D3_yo_xy_, D3_xo_xy_][np.argmin([abs(x) for x in [D3_yo_xo_, D3_yo_xy_, D3_xo_xy_]])]
        D3_yo_xo__ = (dMean_yo - dMean_xo)/(dMax_yo + dMax_xo)
        D3_yo_xy__ = (dMean_yo - dMean_xy)/(dMax_yo + dMax_xy)
        D3_xo_xy__ = (dMean_xo - dMean_xy)/(dMax_xo + dMax_xy)
        D3_xyo__ = [D3_yo_xo__, D3_yo_xy__, D3_xo_xy__][np.argmin([abs(x) for x in [D3_yo_xo__, D3_yo_xy__, D3_xo_xy__]])]
        D3min_yo_xo = (dMin_yo - dMean_xo)/(dMean_yo + dMean_xo)
        D3min_yo_xy = (dMin_yo - dMean_xy)/(dMean_yo + dMean_xy)
        D3min_xo_xy = (dMin_xo - dMean_xy)/(dMean_xo + dMean_xy)
        D3min_xyo = [D3min_yo_xo, D3min_yo_xy, D3min_xo_xy][np.argmin([abs(x) for x in [D3min_yo_xo, D3min_yo_xy, D3min_xo_xy]])]
        D3min_yo_xo_ = (dMin_yo - dMean_xo)/(dMin_yo + dMin_xo)
        D3min_yo_xy_ = (dMin_yo - dMean_xy)/(dMin_yo + dMin_xy)
        D3min_xo_xy_ = (dMin_xo - dMean_xy)/(dMin_xo + dMin_xy)
        D3min_xyo_ = [D3min_yo_xo_, D3min_yo_xy_, D3min_xo_xy_][np.argmin([abs(x) for x in [D3min_yo_xo_, D3min_yo_xy_, D3min_xo_xy_]])]
        D3min_yo_xo__ = (dMin_yo - dMean_xo)/(dMax_yo + dMax_xo)
        D3min_yo_xy__ = (dMin_yo - dMean_xy)/(dMax_yo + dMax_xy)
        D3min_xo_xy__ = (dMin_xo - dMean_xy)/(dMax_xo + dMax_xy)
        D3min_xyo__ = [D3min_yo_xo__, D3min_yo_xy__, D3min_xo_xy__][np.argmin([abs(x) for x in [D3min_yo_xo__, D3min_yo_xy__, D3min_xo_xy__]])]
        D3min_yo_min_xo = (dMin_yo - dMin_xo)/(dMean_yo + dMean_xo)
        D3min_yo_min_xy = (dMin_yo - dMin_xy)/(dMean_yo + dMean_xy)
        D3min_xo_min_xy = (dMin_xo - dMin_xy)/(dMean_xo + dMean_xy)
        D3min_xyo_min = [D3min_yo_min_xo, D3min_yo_min_xy, D3min_xo_min_xy][np.argmin([abs(x) for x in [D3min_yo_min_xo, D3min_yo_min_xy, D3min_xo_min_xy]])]
        D3min_yo_min_xo_ = (dMin_yo - dMin_xo)/(dMin_yo + dMin_xo)
        D3min_yo_min_xy_ = (dMin_yo - dMin_xy)/(dMin_yo + dMin_xy)
        D3min_xo_min_xy_ = (dMin_xo - dMin_xy)/(dMin_xo + dMin_xy)
        D3min_xyo_min_ = [D3min_yo_min_xo_, D3min_yo_min_xy_, D3min_xo_min_xy_][np.argmin([abs(x) for x in [D3min_yo_min_xo_, D3min_yo_min_xy_, D3min_xo_min_xy_]])]
        D3min_yo_min_xo__ = (dMin_yo - dMin_xo)/(dMax_yo + dMax_xo)
        D3min_yo_min_xy__ = (dMin_yo - dMin_xy)/(dMax_yo + dMax_xy)
        D3min_xo_min_xy__ = (dMin_xo - dMin_xy)/(dMax_xo + dMax_xy)
        D3min_xyo_min__ = [D3min_yo_min_xo__, D3min_yo_min_xy__, D3min_xo_min_xy__][np.argmin([abs(x) for x in [D3min_yo_min_xo__, D3min_yo_min_xy__, D3min_xo_min_xy__]])]
    if stats_type == 'dixyo':
        if len_pop_i > 1:
            dMean_i = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dSd_i = np.std(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dMin_i = np.min(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dMax_i = np.max(as_dist(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dSites_i = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_i_pos)]))
            dSum_i = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_i_pos)]))
        else:
            dMean_i = 0.0
            dSd_i = 0.0
            dMin_i = 0.0
            dMax_i = 0.0
            dSites_i = 0.0
            dSum_i = 0.0
        if len_pop_x > 1:
            dMean_x = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSd_x = np.std(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMin_x = np.min(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dMax_x = np.max(as_dist(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSites_x = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
            dSum_x = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_x_pos)]))
        else:
            dMean_x = 0.0
            dSd_x = 0.0
            dMin_x = 0.0
            dMax_x = 0.0
            dSites_x = 0.0
            dSum_x = 0.0
        if len_pop_y > 1:
            dMean_y = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSd_y = np.std(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMin_y = np.min(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dMax_y = np.max(as_dist(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSites_y = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_y_pos)]))
            dSum_y = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_y_pos)]))
        else:
            dMean_y = 0.0
            dSd_y = 0.0
            dMin_y = 0.0
            dMax_y = 0.0
            dSites_y = 0.0
            dSum_y = 0.0
        if len_pop_o > 1:
            dMean_o = np.mean(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dSd_o = np.std(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dMin_o = np.min(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dMax_o = np.max(as_dist(range_dist_df_distances.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dSites_o = np.mean(as_dist(range_dist_df_usedsites.iloc[list(pop_o_pos), list(pop_o_pos)]))
            dSum_o = np.sum(as_dist(range_dist_df_usedsites.iloc[list(pop_o_pos), list(pop_o_pos)]))
        else:
            dMean_o = 0.0
            dSd_o = 0.0
            dMin_o = 0.0
            dMax_o = 0.0
            dSites_o = 0.0
            dSum_o = 0.0
        #ix
        dMean_ix = np.mean(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)])
        dSd_ix = np.std(np.array(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)]))
        dMin_ix = np.min(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)])
        dMax_ix = np.max(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_x_pos)])
        dSites_ix = np.mean(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_x_pos)])
        dSum_ix = np.sum(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_x_pos)].sum())
        dMean_ix_i_x = dMean_ix - (dMean_i/2) - (dMean_x/2)
        dF2_ix = dMean_ix - ((dMean_i + dMean_x)/2)
        Gmin_ix = dMin_ix / dMean_ix
        #iy
        dMean_iy = np.mean(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)])
        dSd_iy = np.std(np.array(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)]))
        dMin_iy = np.min(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)])
        dMax_iy = np.max(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_y_pos)])
        dSites_iy = np.mean(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_y_pos)])
        dSum_iy = np.sum(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_y_pos)].sum())
        dMean_iy_i_y = dMean_iy - (dMean_i/2) - (dMean_y/2)
        dF2_iy = dMean_iy - ((dMean_i + dMean_y)/2)
        Gmin_iy = dMin_iy / dMean_iy
        #io
        dMean_io = np.mean(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_o_pos)])
        dSd_io = np.std(np.array(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_o_pos)]))
        dMin_io = np.min(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_o_pos)])
        dMax_io = np.max(range_dist_df_distances.iloc[list(pop_i_pos), list(pop_o_pos)])
        dSites_io = np.mean(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_o_pos)])
        dSum_io = np.sum(range_dist_df_usedsites.iloc[list(pop_i_pos), list(pop_o_pos)].sum())
        dMean_io_i_o = dMean_io - (dMean_i/2) - (dMean_o/2)
        dF2_io = dMean_io - ((dMean_i + dMean_o)/2)
        Gmin_io = dMin_io / dMean_io
        #xy
        dMean_xy = np.mean(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSd_xy = np.std(np.array(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)]))
        dMin_xy = np.min(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dMax_xy = np.max(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSites_xy = np.mean(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)])
        dSum_xy = np.sum(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_y_pos)].sum())
        dMean_xy_x_y = dMean_xy - (dMean_x/2) - (dMean_y/2)
        dF2_xy = dMean_xy - ((dMean_x + dMean_y)/2)
        Gmin_xy = dMin_xy / dMean_xy
        #xo
        dMean_xo = np.mean(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)])
        dSd_xo = np.std(np.array(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)]))
        dMin_xo = np.min(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)])
        dMax_xo = np.max(range_dist_df_distances.iloc[list(pop_x_pos), list(pop_o_pos)])
        dSites_xo = np.mean(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_o_pos)])
        dSum_xo = np.sum(range_dist_df_usedsites.iloc[list(pop_x_pos), list(pop_o_pos)].sum())
        dMean_xo_x_o = dMean_xo - (dMean_x/2) - (dMean_o/2)
        dF2_xo = dMean_xo - ((dMean_x + dMean_o)/2)
        Gmin_xo = dMin_xo / dMean_xo
        #yo
        dMean_yo = np.mean(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)])
        dSd_yo = np.std(np.array(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)]))
        dMin_yo = np.min(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)])
        dMax_yo = np.max(range_dist_df_distances.iloc[list(pop_y_pos), list(pop_o_pos)])
        dSites_yo = np.mean(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_o_pos)])
        dSum_yo = np.sum(range_dist_df_usedsites.iloc[list(pop_y_pos), list(pop_o_pos)].sum())
        dMean_yo_y_o = dMean_yo - (dMean_y/2) - (dMean_o/2)
        dF2_yo = dMean_yo - ((dMean_y + dMean_o)/2)
        Gmin_yo = dMin_yo / dMean_yo
        #ixy
        RND_ixy = dMean_xy - ((dMean_ix + dMean_iy)/2)
        RNDmin_ixy = dMin_xy - ((dMean_ix + dMean_iy)/2)
        D3_xy_iy = (dMean_xy - dMean_iy)/(dMean_xy + dMean_iy)
        D3_iy_xy = (dMean_iy - dMean_xy)/(dMean_iy + dMean_xy)
        D3_xy_ix = (dMean_xy - dMean_ix)/(dMean_xy + dMean_ix)
        D3_ix_xy = (dMean_ix - dMean_xy)/(dMean_ix + dMean_xy)
        D3_ix_iy = (dMean_ix - dMean_iy)/(dMean_ix + dMean_iy)
        D3_iy_ix = (dMean_iy - dMean_ix)/(dMean_iy + dMean_ix)
        D3_ixy = [D3_xy_iy, D3_xy_ix, D3_ix_iy][np.argmin([abs(x) for x in [D3_xy_iy, D3_xy_ix, D3_ix_iy]])]
        D3_ixy_abs_min = ['((I,X),Y)', '((I,Y),X)', '((X,Y),I)'][np.argmin([abs(x) for x in [D3_xy_iy, D3_xy_ix, D3_ix_iy]])]
        D3_xy_iy_ = (dMean_xy - dMean_iy)/(dMin_xy + dMin_iy)
        D3_xy_ix_ = (dMean_xy - dMean_ix)/(dMin_xy + dMin_ix)
        D3_ix_iy_ = (dMean_ix - dMean_iy)/(dMin_ix + dMin_iy)
        D3_ixy_ = [D3_xy_iy_, D3_xy_ix_, D3_ix_iy_][np.argmin([abs(x) for x in [D3_xy_iy_, D3_xy_ix_, D3_ix_iy_]])]
        D3_xy_iy__ = (dMean_xy - dMean_iy)/(dMax_xy + dMax_iy)
        D3_xy_ix__ = (dMean_xy - dMean_ix)/(dMax_xy + dMax_ix)
        D3_ix_iy__ = (dMean_ix - dMean_iy)/(dMax_ix + dMax_iy)
        D3_ixy__ = [D3_xy_iy__, D3_xy_ix__, D3_ix_iy__][np.argmin([abs(x) for x in [D3_xy_iy__, D3_xy_ix__, D3_ix_iy__]])]
        D3min_xy_iy = (dMin_xy - dMean_iy)/(dMean_xy + dMean_iy)
        D3min_xy_ix = (dMin_xy - dMean_ix)/(dMean_xy + dMean_ix)
        D3min_ix_iy = (dMin_ix - dMean_iy)/(dMean_ix + dMean_iy)
        D3min_ixy = [D3min_xy_iy, D3min_xy_ix, D3min_ix_iy][np.argmin([abs(x) for x in [D3min_xy_iy, D3min_xy_ix, D3min_ix_iy]])]
        D3min_xy_iy_ = (dMin_xy - dMean_iy)/(dMin_xy + dMin_iy)
        D3min_xy_ix_ = (dMin_xy - dMean_ix)/(dMin_xy + dMin_ix)
        D3min_ix_iy_ = (dMin_ix - dMean_iy)/(dMin_ix + dMin_iy)
        D3min_ixy_ = [D3min_xy_iy_, D3min_xy_ix_, D3min_ix_iy_][np.argmin([abs(x) for x in [D3min_xy_iy_, D3min_xy_ix_, D3min_ix_iy_]])]
        D3min_xy_iy__ = (dMin_xy - dMean_iy)/(dMax_xy + dMax_iy)
        D3min_xy_ix__ = (dMin_xy - dMean_ix)/(dMax_xy + dMax_ix)
        D3min_ix_iy__ = (dMin_ix - dMean_iy)/(dMax_ix + dMax_iy)
        D3min_ixy__ = [D3min_xy_iy__, D3min_xy_ix__, D3min_ix_iy__][np.argmin([abs(x) for x in [D3min_xy_iy__, D3min_xy_ix__, D3min_ix_iy__]])]
        D3min_xy_min_iy = (dMin_xy - dMin_iy)/(dMean_xy + dMean_iy)
        D3min_xy_min_ix = (dMin_xy - dMin_ix)/(dMean_xy + dMean_ix)
        D3min_ix_min_iy = (dMin_ix - dMin_iy)/(dMean_ix + dMean_iy)
        D3min_ixy_min = [D3min_xy_min_iy, D3min_xy_min_ix, D3min_ix_min_iy][np.argmin([abs(x) for x in [D3min_xy_min_iy, D3min_xy_min_ix, D3min_ix_min_iy]])]
        D3min_xy_min_iy_ = (dMin_xy - dMin_iy)/(dMin_xy + dMin_iy)
        D3min_xy_min_ix_ = (dMin_xy - dMin_ix)/(dMin_xy + dMin_ix)
        D3min_ix_min_iy_ = (dMin_ix - dMin_iy)/(dMin_ix + dMin_iy)
        D3min_ixy_min_ = [D3min_xy_min_iy_, D3min_xy_min_ix_, D3min_ix_min_iy_][np.argmin([abs(x) for x in [D3min_xy_min_iy_, D3min_xy_min_ix_, D3min_ix_min_iy_]])]
        D3min_xy_min_iy__ = (dMin_xy - dMin_iy)/(dMax_xy + dMax_iy)
        D3min_xy_min_ix__ = (dMin_xy - dMin_ix)/(dMax_xy + dMax_ix)
        D3min_ix_min_iy__ = (dMin_ix - dMin_iy)/(dMax_ix + dMax_iy)
        D3min_ixy_min__ = [D3min_xy_min_iy__, D3min_xy_min_ix__, D3min_ix_min_iy__][np.argmin([abs(x) for x in [D3min_xy_min_iy__, D3min_xy_min_ix__, D3min_ix_min_iy__]])]
        #xyo
        RND_xyo = dMean_xy / ((dMean_xo + dMean_yo)/2)
        RNDmin_xyo = dMin_xy / ((dMean_xo + dMean_yo)/2)
        D3_yo_xo = (dMean_yo - dMean_xo)/(dMean_yo + dMean_xo)
        D3_xo_yo = (dMean_xo - dMean_yo)/(dMean_xo + dMean_yo)
        D3_yo_xy = (dMean_yo - dMean_xy)/(dMean_yo + dMean_xy)
        D3_xy_yo = (dMean_xy - dMean_yo)/(dMean_xy + dMean_yo)
        D3_xo_xy = (dMean_xo - dMean_xy)/(dMean_xo + dMean_xy)
        D3_xy_xo = (dMean_xy - dMean_xo)/(dMean_xy + dMean_xo)
        D3_xyo = [D3_yo_xo, D3_yo_xy, D3_xo_xy][np.argmin([abs(x) for x in [D3_yo_xo, D3_yo_xy, D3_xo_xy]])]
        D3_xyo_abs_min = ['((X,Y),O)', '((O,X),Y)', '((Y,O),X)'][np.argmin([abs(x) for x in [D3_yo_xo, D3_yo_xy, D3_xo_xy]])]
        D3_yo_xo_ = (dMean_yo - dMean_xo)/(dMin_yo + dMin_xo)
        D3_yo_xy_ = (dMean_yo - dMean_xy)/(dMin_yo + dMin_xy)
        D3_xo_xy_ = (dMean_xo - dMean_xy)/(dMin_xo + dMin_xy)
        D3_xyo_ = [D3_yo_xo_, D3_yo_xy_, D3_xo_xy_][np.argmin([abs(x) for x in [D3_yo_xo_, D3_yo_xy_, D3_xo_xy_]])]
        D3_yo_xo__ = (dMean_yo - dMean_xo)/(dMax_yo + dMax_xo)
        D3_yo_xy__ = (dMean_yo - dMean_xy)/(dMax_yo + dMax_xy)
        D3_xo_xy__ = (dMean_xo - dMean_xy)/(dMax_xo + dMax_xy)
        D3_xyo__ = [D3_yo_xo__, D3_yo_xy__, D3_xo_xy__][np.argmin([abs(x) for x in [D3_yo_xo__, D3_yo_xy__, D3_xo_xy__]])]
        D3min_yo_xo = (dMin_yo - dMean_xo)/(dMean_yo + dMean_xo)
        D3min_yo_xy = (dMin_yo - dMean_xy)/(dMean_yo + dMean_xy)
        D3min_xo_xy = (dMin_xo - dMean_xy)/(dMean_xo + dMean_xy)
        D3min_xyo = [D3min_yo_xo, D3min_yo_xy, D3min_xo_xy][np.argmin([abs(x) for x in [D3min_yo_xo, D3min_yo_xy, D3min_xo_xy]])]
        D3min_yo_xo_ = (dMin_yo - dMean_xo)/(dMin_yo + dMin_xo)
        D3min_yo_xy_ = (dMin_yo - dMean_xy)/(dMin_yo + dMin_xy)
        D3min_xo_xy_ = (dMin_xo - dMean_xy)/(dMin_xo + dMin_xy)
        D3min_xyo_ = [D3min_yo_xo_, D3min_yo_xy_, D3min_xo_xy_][np.argmin([abs(x) for x in [D3min_yo_xo_, D3min_yo_xy_, D3min_xo_xy_]])]
        D3min_yo_xo__ = (dMin_yo - dMean_xo)/(dMax_yo + dMax_xo)
        D3min_yo_xy__ = (dMin_yo - dMean_xy)/(dMax_yo + dMax_xy)
        D3min_xo_xy__ = (dMin_xo - dMean_xy)/(dMax_xo + dMax_xy)
        D3min_xyo__ = [D3min_yo_xo__, D3min_yo_xy__, D3min_xo_xy__][np.argmin([abs(x) for x in [D3min_yo_xo__, D3min_yo_xy__, D3min_xo_xy__]])]
        D3min_yo_min_xo = (dMin_yo - dMin_xo)/(dMean_yo + dMean_xo)
        D3min_yo_min_xy = (dMin_yo - dMin_xy)/(dMean_yo + dMean_xy)
        D3min_xo_min_xy = (dMin_xo - dMin_xy)/(dMean_xo + dMean_xy)
        D3min_xyo_min = [D3min_yo_min_xo, D3min_yo_min_xy, D3min_xo_min_xy][np.argmin([abs(x) for x in [D3min_yo_min_xo, D3min_yo_min_xy, D3min_xo_min_xy]])]
        D3min_yo_min_xo_ = (dMin_yo - dMin_xo)/(dMin_yo + dMin_xo)
        D3min_yo_min_xy_ = (dMin_yo - dMin_xy)/(dMin_yo + dMin_xy)
        D3min_xo_min_xy_ = (dMin_xo - dMin_xy)/(dMin_xo + dMin_xy)
        D3min_xyo_min_ = [D3min_yo_min_xo_, D3min_yo_min_xy_, D3min_xo_min_xy_][np.argmin([abs(x) for x in [D3min_yo_min_xo_, D3min_yo_min_xy_, D3min_xo_min_xy_]])]
        D3min_yo_min_xo__ = (dMin_yo - dMin_xo)/(dMax_yo + dMax_xo)
        D3min_yo_min_xy__ = (dMin_yo - dMin_xy)/(dMax_yo + dMax_xy)
        D3min_xo_min_xy__ = (dMin_xo - dMin_xy)/(dMax_xo + dMax_xy)
        D3min_xyo_min__ = [D3min_yo_min_xo__, D3min_yo_min_xy__, D3min_xo_min_xy__][np.argmin([abs(x) for x in [D3min_yo_min_xo__, D3min_yo_min_xy__, D3min_xo_min_xy__]])]
    results = [dMean_i, dSd_i, dMin_i, dMax_i, dSites_i, dSum_i,
        dMean_x, dSd_x, dMin_x, dMax_x, dSites_x, dSum_x,
        dMean_y, dSd_y, dMin_y, dMax_y, dSites_y, dSum_y,
        dMean_o, dSd_o, dMin_o, dMax_o, dSites_o, dSum_o,
        dMean_ix, dSd_ix, dMin_ix, dMax_ix, dSites_ix, dSum_ix, dMean_ix_i_x, dF2_ix, Gmin_ix,
        dMean_iy, dSd_iy, dMin_iy, dMax_iy, dSites_iy, dSum_iy, dMean_iy_i_y, dF2_iy, Gmin_iy,
        dMean_io, dSd_io, dMin_io, dMax_io, dSites_io, dSum_io, dMean_io_i_o, dF2_io, Gmin_io,
        dMean_xy, dSd_xy, dMin_xy, dMax_xy, dSites_xy, dSum_xy, dMean_xy_x_y, dF2_xy, Gmin_xy,
        dMean_xo, dSd_xo, dMin_xo, dMax_xo, dSites_xo, dSum_xo, dMean_xo_x_o, dF2_xo, Gmin_xo,
        dMean_yo, dSd_yo, dMin_yo, dMax_yo, dSites_yo, dSum_yo, dMean_yo_y_o, dF2_yo, Gmin_yo,
        RND_ixy, RNDmin_ixy,
        D3_xy_iy, D3_iy_xy,
        D3_xy_ix, D3_ix_xy,
        D3_ix_iy, D3_iy_ix,
        D3_ixy, D3_ixy_abs_min,
        D3_xy_iy_, D3_xy_ix_, D3_ix_iy_, D3_ixy_,
        D3_xy_iy__, D3_xy_ix__, D3_ix_iy__, D3_ixy__,
        D3min_xy_iy, D3min_xy_ix, D3min_ix_iy, D3min_ixy,
        D3min_xy_iy_, D3min_xy_ix_, D3min_ix_iy_, D3min_ixy_,
        D3min_xy_iy__, D3min_xy_ix__, D3min_ix_iy__, D3min_ixy__,
        D3min_xy_min_iy, D3min_xy_min_ix, D3min_ix_min_iy, D3min_ixy_min,
        D3min_xy_min_iy_, D3min_xy_min_ix_, D3min_ix_min_iy_, D3min_ixy_min_,
        D3min_xy_min_iy__, D3min_xy_min_ix__, D3min_ix_min_iy__, D3min_ixy_min__,
        RND_xyo, RNDmin_xyo,
        D3_yo_xo, D3_xo_yo,
        D3_yo_xy, D3_xy_yo,
        D3_xo_xy, D3_xy_xo,
        D3_xyo, D3_xyo_abs_min,
        D3_yo_xo_, D3_yo_xy_, D3_xo_xy_, D3_xyo_,
        D3_yo_xo__, D3_yo_xy__, D3_xo_xy__, D3_xyo__,
        D3min_yo_xo, D3min_yo_xy, D3min_xo_xy, D3min_xyo,
        D3min_yo_xo_, D3min_yo_xy_, D3min_xo_xy_, D3min_xyo_,
        D3min_yo_xo__, D3min_yo_xy__, D3min_xo_xy__, D3min_xyo__,
        D3min_yo_min_xo, D3min_yo_min_xy, D3min_xo_min_xy, D3min_xyo_min,
        D3min_yo_min_xo_, D3min_yo_min_xy_, D3min_xo_min_xy_, D3min_xyo_min_,
        D3min_yo_min_xo__, D3min_yo_min_xy__, D3min_xo_min_xy__, D3min_xyo_min__]
    return results


def define_parser():
    parser = argparse.ArgumentParser(prog='fasta2dstats', usage='%(prog)s [options] [<arguments>...]',
                                     description='dxy/dxyi/dxyo/dxyio stats from FASTA sequences')
    parser.add_argument('-i', help='file with population assignment and fasta file path (see example_pop.txt); use full path to avoid file finding')
    parser.add_argument('-r', help='file with genomic ranges to be analysed in bed format (see example_region.txt)')
    parser.add_argument('-c', help='file with population combinations to be used for stats (see example_combn.txt)')
    parser.add_argument('-l', help='path to literal-dist binary')
    parser.add_argument('-o', help='output file [default: dxy_stats.txt]', default='dxy_stats.txt')
    parser.add_argument('-g', help='remove sites due to gap frequency threshold', action='store_true')
    parser.add_argument('-z', help='gap frequency threshold', type=float, default=0.0)
    return parser


def main():
    # parser
    parser = define_parser()
    # get args
    args = parser.parse_args()
    # print args
    sys.stderr.write(str(args))
    # run
    fasta2dstats(args, parser)


if __name__ == '__main__':
    main()
