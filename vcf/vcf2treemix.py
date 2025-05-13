#!/usr/bin/python


import sys
import gzip
from collections import OrderedDict


vcf_file = sys.argv[1]
pop_file = sys.argv[2]
treemix_file = sys.argv[3]


def get_pops(pop_file):
    """
    Returns a dictionary with pop identifier as key and taxa as a list of
    strings. In the pop file, each populations should be in one line, starting
    withe pop name, a colon and the corresponding taxa separated by whitespace.
    E.g.:
    pop1: taxon1 taxon2 taxon3
    """
    pops = OrderedDict()
    with open(pop_file) as fh:

        for line in fh:
            fields = line.strip().split(':')
            pops[fields[0]] = fields[1].split()
    return pops


pop_obj = get_pops(pop_file)


with gzip.open(vcf_file,'rt') as inhandle:
    with open(treemix_file,'w') as outhandle:
        outhandle.write("{}\n".format(" ".join([x for x in list(pop_obj.keys())])))
        for line in inhandle:
            if(line[0]=="#"):
                if line.strip().split('\t')[0]=="#CHROM":
                    taxa_pos = line.strip().split('\t')
                else:
                    continue
            if(line[0]!="#"):
                fields = line.strip().split('\t')
                temp_pop = OrderedDict((x, [0,0]) for x in list(pop_obj.keys()))
                for pop, taxa in list(pop_obj.items()):
                    for taxon in taxa:
                        gt = fields[taxa_pos.index(taxon)]
                        if gt==".|.":
                            continue
                        temp_pop[pop][0] += gt.count("0")
                        temp_pop[pop][1] += gt.count("1")
                outhandle.write("{}\n".format(" ".join([str(x[0]) +  "," + str(x[1]) for x in list(temp_pop.values())])))
