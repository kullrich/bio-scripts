#!/usr/bin/python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
import sys
import argparse
import math

def main():
  parser = argparse.ArgumentParser(usage='%(prog)s [options]',description='Extracts codonusage from CDS input FASTA file')
  parser.add_argument("-v", "--verbose", help="increase output verbosity",action="store_true")
  parser.add_argument('-i', help='specify CDS input file in FASTA format')
  parser.add_argument('-o', help='specify output prefix')
  parser.add_argument('-r', action='store_true', help='specify if CDS sequences with length modulo 3 unequal to 0 should be removed and reported to std.out')
  args = parser.parse_args()

  if args.i is None:
    parser.print_help()
    sys.exit('\nPlease specify input fasta file')
  if args.o is None:
    parser.print_help()
    sys.exit('\nPlease specify output prefix')
  if args.i == args.o:
    sys.exit('\nInput file and output prefix are identical, use "out" as output prefix instead')

  print '\ncommand arguments used:\n'
  print args

  infile = args.i
  outfile = args.o

  print("read fasta")
  original_fasta = SeqIO.parse(infile, "fasta")
  print("extract codon counts")
  global_codons = codontable()
  ids_mo3 = []
  with open(outfile,"w") as handle:
    handle.write("id\tlen\tmo3\t" + "\t".join(sorted(global_codons.keys())) + "\n")
    c = 0
    cmo3 = 0
    for record in original_fasta:
      c += 1
      tmp_counts = codontable()
      tmp_id = record.id
      tmp_len = len(record)
      tmp_mo3 = 0
      if len(record)%3 != 0:
        tmp_mo3 = 1
        cmo3 += 1
        ids_mo3.append(tmp_id)
        if args.r:
          continue
      for i in range(0, len(record), 3):
        codon = str(record[i:i+3].seq)
        if len(codon) != 3:
          continue
        if codon in global_codons:
          global_codons[codon][1]+=1
          tmp_counts[codon][1]+=1
        if codon not in global_codons:
          global_codons['XXX'][1]+=1
          tmp_counts['XXX'][1]+=1
      handle.write(tmp_id + "\t" + str(tmp_len) + "\t" + str(tmp_mo3) + "\t" + "\t".join([str(tmp_counts[x][1]) for x in sorted(tmp_counts.keys())]) + "\n")
    handle.write("global_count\t" + str(c) + "\t" + str(cmo3) + "\t" + "\t".join([str(global_codons[x][1]) for x in sorted(global_codons.keys())]) + "\n")
  if args.r:
    print '\n'.join(ids_mo3)

def inversetable():
  inversetablecount = {
    'A': {'GCT':0,'GCC':0,'GCA':0,'GCG':0},
    'R': {'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
    'N': {'AAT':0,'AAC':0},
    'D': {'GAT':0,'GAC':0},
    'C': {'TGT':0,'TGC':0},
    'Q': {'CAA':0,'CAG':0},
    'E': {'GAA':0,'GAG':0},
    'G': {'GGT':0,'GGC':0,'GGA':0,'GGG':0},
    'H': {'CAT':0,'CAC':0},
    'I': {'ATT':0,'ATC':0,'ATA':0},
    'M': {'ATG':0},
    'L': {'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
    'K': {'AAA':0,'AAG':0},
    'F': {'TTT':0,'TTC':0},
    'P': {'CCT':0,'CCC':0,'CCA':0,'CCG':0},
    'S': {'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
    'T': {'ACT':0,'ACC':0,'ACA':0,'ACG':0},
    'W': {'TGG':0},
    'Y': {'TAT':0,'TAC':0},
    'V': {'GTT':0,'GTC':0,'GTA':0,'GTG':0},
    '*': {'TAA':0,'TGA':0,'TAG':0}
  }
  return inversetablecount

def codontable():
  codontablecount = {
    'GCT':['A',0],
    'GCC':['A',0],
    'GCA':['A',0],
    'GCG':['A',0],
    'CGT':['R',0],
    'CGC':['R',0],
    'CGG':['R',0],
    'CGA':['R',0],
    'AGA':['R',0],
    'AGG':['R',0],
    'AAT':['N',0],
    'AAC':['N',0],
    'GAT':['D',0],
    'GAC':['D',0],
    'TGT':['C',0],
    'TGC':['C',0],
    'CAA':['Q',0],
    'CAG':['Q',0],
    'GAA':['E',0],
    'GAG':['E',0],
    'GGT':['G',0],
    'GGC':['G',0],
    'GGA':['G',0],
    'GGG':['G',0],
    'CAT':['H',0],
    'CAC':['H',0],
    'ATT':['I',0],
    'ATC':['I',0],
    'ATA':['I',0],
    'ATG':['M',0],
    'TTA':['L',0],
    'TTG':['L',0],
    'CTT':['L',0],
    'CTC':['L',0],
    'CTA':['L',0],
    'CTG':['L',0],
    'AAA':['K',0],
    'AAG':['K',0],
    'TTT':['F',0],
    'TTC':['F',0],
    'CCT':['P',0],
    'CCC':['P',0],
    'CCA':['P',0],
    'CCG':['P',0],
    'TCT':['S',0],
    'TCC':['S',0],
    'TCA':['S',0],
    'TCG':['S',0],
    'AGT':['S',0],
    'AGC':['S',0],
    'ACT':['T',0],
    'ACC':['T',0],
    'ACA':['T',0],
    'ACG':['T',0],
    'TGG':['W',0],
    'TAT':['Y',0],
    'TAC':['Y',0],
    'GTT':['V',0],
    'GTC':['V',0],
    'GTA':['V',0],
    'GTG':['V',0],
    'TAA':['*',0],
    'TGA':['*',0],
    'TAG':['*',0],
    'XXX':['_missing',0]
  }
  return codontablecount

if __name__ == '__main__':
    main()
