#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Author: Krisian Ullrich
date: April 2025
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


import os
import gzip
import sys
import argparse
import numpy as np
import scipy.sparse as sp
import pickle
from collections import defaultdict


def define_parser():
    parser = argparse.ArgumentParser(prog='orthodb2sample', usage='%(prog)s [options] [<arguments>...]',
                                     description='extract sample specific Orthogroups.GeneCount and Orthogroups from OrthoDB')
    parser.add_argument('-species_id', help='specify species id from OrthoDB species e.g. 3702_0')
    parser.add_argument('-species', help='input file: odb12v1_species.tab.gz')
    parser.add_argument('-og2genes', help='input file: odb12v1_OG2genes.tab.gz')
    parser.add_argument('-og_pairs', help='input file: odb12v1_OG_pairs.tab.gz')
    parser.add_argument('-genes', help='input file: odb12v1_genes.tab.gz')
    parser.add_argument('-sl', help='output file: species_list.tsv', default='species_list.tsv')
    parser.add_argument('-oc', help='output file: Orthogroups.GeneCount.tsv', default='Orthogroups.GeneCount.tsv')
    parser.add_argument('-og', help='output file: Orthogroups.tsv', default='Orthogroups.tsv')
    parser.add_argument('-use_filtered_ogs', action='store_true', help='use only the filtered set of orthogroups (most ancestral + isolated)')
    parser.add_argument('-use_original_protein_id', action='store_true', help='not use Ortho DB unique gene id but instead protein original sequence id')
    return parser


def get_samples(args, parser):
    samples = []
    with gzip.open(args.species, 'rt') as f:
        for line in f:
            tax_id, species_id, *_ = line.strip().split('\t')
            samples.append(species_id)
    return samples


def get_sample_ogs(args, parser):
    sample_ogs = set()
    with gzip.open(args.og2genes, 'rt') as f:
        for line in f:
            og_id, gene_id = line.strip().split('\t')
            species_id = gene_id.split(":")[0]
            if species_id == args.species_id:
                sample_ogs.add(og_id)
    return sample_ogs


def get_filtered_sample_ogs(args, parser, sample_ogs):
    descendants = set()
    ancestors = set()
    with gzip.open(args.og_pairs, 'rt') as f:
        for line in f:
            desc, anc = line.strip().split('\t')
            descendants.add(desc)
            ancestors.add(anc)
    most_ancestral = ancestors - descendants
    involved_ogs = descendants.union(ancestors)
    completely_isolated = sample_ogs - involved_ogs
    filtered_sample_ogs = (sample_ogs & most_ancestral).union(completely_isolated)
    return filtered_sample_ogs


def get_sample_oc_og(args, parser, sample_to_index, sample_ogs):
    sample = args.species_id  # The sample of interest
    sample_idx = sample_to_index[sample]
    # Prepare data structures
    og_counts = defaultdict(lambda: [0] * len(sample_to_index))
    og_genes = defaultdict(lambda: [''] * len(sample_to_index))
    with gzip.open(args.og2genes, 'rt') as f:
        for line in f:
            og_id, gene_id = line.strip().split('\t')
            species_id = gene_id.split(":")[0]
            if og_id not in sample_ogs:
                continue
            if species_id not in sample_to_index:
                continue
            idx = sample_to_index[species_id]
            og_counts[og_id][idx] += 1
            # Add gene_id to the correct place in the gene list, with comma separation
            if og_genes[og_id][idx] == '':
                og_genes[og_id][idx] = gene_id
            else:
                og_genes[og_id][idx] += f",{gene_id}"
    return og_counts, og_genes


def filter_empty_columns(data_dict, samples, is_counts=True):
    num_cols = len(samples)
    keep_indices = set(range(num_cols))
    for idx in range(num_cols):
        all_empty = True
        for values in data_dict.values():
            val = values[idx]
            if is_counts and val != 0:
                all_empty = False
                break
            if not is_counts and val != '':
                all_empty = False
                break
        if all_empty:
            keep_indices.discard(idx)
    keep_indices = sorted(keep_indices)
    filtered_samples = [samples[i] for i in keep_indices]
    filtered_data = {
        og: [values[i] for i in keep_indices]
        for og, values in data_dict.items()
    }
    return filtered_data, filtered_samples


def create_species_list(args, parser, samples):
    with open(args.sl, 'w') as fo:
        with gzip.open(args.species, 'rt') as f:
            for line in f:
                tax_id, species_id, *_ = line.strip().split('\t')
                if species_id in samples:
                    fo.write(species_id + '\t' + tax_id + '\n')


def get_gene_to_protein_mapping(args, parser):
    mapping = {}
    with gzip.open(args.genes, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            gene_id, species_id, protein_id, *_ = parts
            if species_id == args.species_id:
                mapping[gene_id] = protein_id
    return mapping


def write_tsv(data_dict, samples, output_path, is_counts=True, gene_to_protein=None, target_species_idx=None):
    with open(output_path, 'w') as f:
        header = ['Orthogroup'] + samples + ['Total'] if is_counts else ['Orthogroup'] + samples
        f.write('\t'.join(header) + '\n')
        for og, values in sorted(data_dict.items()):
            if is_counts:
                total = str(sum(values))
                row = [og] + [str(val) for val in values] + [total]
            else:
                new_values = []
                for idx, val in enumerate(values):
                    if gene_to_protein and idx == target_species_idx and val:
                        ids = val.split(',')
                        mapped = [gene_to_protein.get(gid, gid) for gid in ids]
                        new_values.append(','.join(mapped))
                    else:
                        new_values.append(val)
                row = [og] + new_values
            f.write('\t'.join(row) + '\n')


def main():
    # parser
    parser = define_parser()
    # get args
    args = parser.parse_args()
    # print args
    if args.species_id is None:
        parser.print_help()
        sys.exit('Please specify species id from OrthoDB species e.g. 3702_0')
    if args.species is None:
        parser.print_help()
        sys.exit('Please specify odb12v1_species.tab.gz file path')
    if args.og2genes is None:
        parser.print_help()
        sys.exit('Please specify odb12v1_OG2genes.tab.gz file path')
    if args.og_pairs is None:
        parser.print_help()
        sys.exit('Please specify odb12v1_OG_pairs.tab.gz file path')
    if args.genes is None:
        parser.print_help()
        sys.exit('Please specify odb12v1_genes.tab.gz file path')
    gene_to_protein = None
    if args.use_original_protein_id:
        gene_to_protein = get_gene_to_protein_mapping(args, parser)
    samples = get_samples(args, parser)
    sample_to_index = {name: idx for idx, name in enumerate(samples)}
    sample_ogs = get_sample_ogs(args, parser)
    filtered_sample_ogs = get_filtered_sample_ogs(args, parser, sample_ogs)
    og_set_to_use = filtered_sample_ogs if args.use_filtered_ogs else sample_ogs
    sample_oc, sample_og = get_sample_oc_og(args, parser, sample_to_index, og_set_to_use)
    sample_oc_filtered, samples_filtered_counts = filter_empty_columns(sample_oc, samples, is_counts=True)
    sample_og_filtered, samples_filtered_genes = filter_empty_columns(sample_og, samples, is_counts=False)
    assert samples_filtered_counts == samples_filtered_genes, "Mismatch in filtered sample columns"
    #Write species list
    create_species_list(args, parser, samples_filtered_counts)
    # Write zipped gene counts
    sample_to_index_filtered = {name: idx for idx, name in enumerate(samples_filtered_genes)}
    target_species_idx = sample_to_index_filtered.get(args.species_id)
    #write_tsv(sample_oc, samples, args.oc, is_counts=True)
    write_tsv(sample_oc_filtered, samples_filtered_counts, args.oc, is_counts=True)
    # Write zipped gene identifiers
    #write_tsv(sample_og, samples, args.og, is_counts=False)
    write_tsv(sample_og_filtered, samples_filtered_genes, args.og, is_counts=False, gene_to_protein=gene_to_protein, target_species_idx=target_species_idx)
    print(f"Files written: {args.oc}, {args.og}")


if __name__ == '__main__':
    main()
