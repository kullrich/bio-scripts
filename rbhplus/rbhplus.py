#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Reciprocal BLAST+ to find homologs.
Copyright 2015 Kristian Ullrich
author Kristian Ullrich
author email kristian.ullrich@biologie.uni-marburg.de
LICENCSE MIT - http://www.opensource.org/licenses/mit-license.php
"""

import os
import sys
import argparse
import math

RBHPLUS_version = '0.1'

###FOR TESTING###
#prog='ghostx'
#oprog='ghostx'
#oblastpath=''
#oquery_species='VanBel_2928.fa'
#odatabase_species='CEGMA_248_ARATH.fa'
#query_species='VanBel_2928.fa'
#database_species='CEGMA_248_ARATH.fa'
#oquery_species_type='prot'
#odatabase_species_type='prot'
#oselfblast='True'
#outpath='.'
#if prog == 'blast+':
#  outfile_qd = outpath+'/'+query_species+'.vs.'+database_species+'.blast.out'
#  outfile_dq = outpath+'/'+database_species+'.vs.'+query_species+'.blast.out'
#  outfile_qq = outpath+'/'+query_species+'.vs.'+query_species+'.blast.out'
#  outfile_dd = outpath+'/'+database_species+'.vs.'+database_species+'.blast.out'
#  selfscore_q = outpath+'/'+query_species+'.blast.selfscore'
#  selfscore_d = outpath+'/'+database_species+'.blast.selfscore'
#if prog == 'ghostx':
#  outfile_qd = outpath+'/'+query_species+'.vs.'+database_species+'.ghostx.out'
#  outfile_dq = outpath+'/'+database_species+'.vs.'+query_species+'.ghostx.out'
#  outfile_qq = outpath+'/'+query_species+'.vs.'+query_species+'.ghostx.out'
#  outfile_dd = outpath+'/'+database_species+'.vs.'+database_species+'.ghostx.out'
#  selfscore_q = outpath+'/'+query_species+'.ghostx.selfscore'
#  selfscore_d = outpath+'/'+database_species+'.ghostx.selfscore'
#ocomparison_qd, ocomparison_dq, ocomparison_qq, ocomparison_dd = get_comparison_by_type(oquery_species_type, odatabase_species_type, oselfblast)
#oselfscore=''
#ooutfile='test'
#ooption_pident=float(30)
#ooption_length=float(80)
#ooption_covs=float(0.5)
#ooption_covq=float(0.5)
#opident_method='static'
#osort_method='evl'
#onum_targets=5
#ostep='1'
#ooption_evalue=float(1e-5)
#omatrix='BLOSUM62'
#onum_threads=int(1)
#ooutpath='.'
#oprecheck='False'
#oprepath='.'
#ooutfile_qd=outfile_qd
#ooutfile_dq=outfile_dq
#ooutfile_qq=outfile_qq
#ooutfile_dd=outfile_dd
#oselfscore_q=selfscore_q
#oselfscore_d=selfscore_d
###

def get_species_type_by_comparison(comparison_qd):
  if comparison_qd == 'blastp':
    query_species_type = 'prot'
    database_species_type = 'prot'
    return [query_species_type, database_species_type]
  if comparison_qd == 'blastx':
    query_species_type = 'nucl'
    database_sepcies_type = 'prot'
    return [query_species_type, database_species_type]
  if comparison_qd == 'tblastn':
    query_species_type = 'prot'
    database_sepcies_type = 'nucl'
    return [query_species_type, database_species_type]
  if comparison_qd == 'tblastx':
    query_species_type = 'nucl'
    database_sepcies_type = 'nucl'
    return [query_species_type, database_species_type]

def get_comparison_by_type(qt, dbt, selfblast):
  if qt == 'prot' and dbt == 'prot' and selfblast=='False':
    comparison_qd = 'blastp'
    comparison_dq = 'blastp'
    comparison_qq = None
    comparison_dd = None
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
  if qt == 'prot' and dbt == 'prot' and selfblast=='True':
    comparison_qd = 'blastp'
    comparison_dq = 'blastp'
    comparison_qq = 'blastp'
    comparison_dd = 'blastp'
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
  if qt == 'prot' and dbt == 'nucl' and selfblast=='False':
    comparison_qd = 'tblastn'
    comparison_dq = 'blastx'
    comparison_qq = None
    comparison_dd = None
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
  if qt == 'prot' and dbt == 'nucl' and selfblast=='True':
    comparison_qd = 'tblastn'
    comparison_dq = 'blastx'
    comaprison_qq = 'blastp'
    comparison_dd = 'tblastx'
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
  if qt == 'nucl' and dbt == 'prot' and selfblast=='False':
    comparison_qd = 'blastx'
    comparison_dq = 'tblastn'
    comparison_qq = None
    comparison_dd = None
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
  if qt == 'nucl' and dbt == 'prot' and selfblast=='True':
    comparison_qd = 'blastx'
    comaprison_dq = 'tblastn'
    comparison_qq = 'tblastx'
    comparison_dd = 'blastp'
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
  if qt == 'nucl' and dbt == 'nucl' and selfblast=='False':
    comparison_qd = 'tblastx'
    compariosn_dq = 'tblastx'
    comparison_qq = None
    comparison_dd = None
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
  if qt == 'nucl' and dbt == 'nucl' and selfblast=='True':
    comparison_qd = 'tblastx'
    comaprison_dq = 'tblastx'
    comparison_qq = 'tblastx'
    comparison_dd = 'tblastx'
    return [comparison_qd, comparison_dq, comparison_qq, comparison_dd]
    
def commandoptions():
  parser = argparse.ArgumentParser(prog='rbhplus', usage='%(prog)s [options] -prog blast+ -q AFASTTA -qt prot -d BFASTA -dt prot -step 1 -out AFASTA.BFASTA.out', description='Performs conditional reciprocal BLAST+ [-prog blast+] or ghostx [-prog ghostx] or [-prog diamond] searches between two species by either starting and performing BLAST+/ghostx/diamond search or providing BLAST+/ghostx/diamond produced table files. BLAST+/ghostx/diamond binary are supposed to be in the PATH otherwise specify with [pp] option. Mandatory are FASTA format sequences for two species and their sequence type (see help). For single species reciprocal BLAST+/ghostx/diamond use the same input sequence as query and database, please "omit the selfblast option in this case!". Optional parameters can be set, e.g. if you want to add selfblast (within species) BLAST+/ghostx/diamond comparison [default: only between; can be changed with the selfblast option], if BLAST+/ghostx/dimaond should build the necessary BLAST+/ghostx/diamond databases with makeblastdb/ghostx db/diamond makedb [default: false] and some BLAST+/ghostx/diamond specific options (see help).')
  parser.add_argument('-prog', default='blast+', help='specify if BLAST+ or ghostx should be used for sequence comparison [default: blast+]')
  parser.add_argument('-pp', default='', help='specify BLAST+/ghostx binary path [default: PATH]')
  parser.add_argument('-selfblast', default='False', choices=['True','False'], help='specify if BLAST+/ghostx search should be run also within each species to extract selfscores. Can also run as standalone option with step option set to [-step self]. This step is necessary if you would like to sort the BLAST+/ghostx result by selfscores.')
  parser.add_argument('-step', default='1', choices=['0','1','2','3','4','self'], help='specify which step should be performed: [1]->BLAST+/ghostx/diamond initially builds the necessary BLAST+/ghostx/diamond sequence databases. [2]->BLAST+/ghostx/diamond search is performed. [3]->BLAST+/ghostx/diamond output is parsed. [4]->Reciprocal Best Hits (RBHs) are calculated. Optional: [self]->BLAST+/ghostx/diamond initially builds the necessary BLAST+/ghostx/diamond sequence databases if it is not already existing and selfscores are obtained for each species. Optional: [0]->BLAST+/ghostx/diamond builds the BLAST+/ghostx sequence databases and exists.')
  parser.add_argument('-q', help='specify query')
  parser.add_argument('-qt', default='prot', choices=['prot','nucl'], help='specify query type [default: prot] or [nucl]')
  parser.add_argument('-db', help='specify database')
  parser.add_argument('-dbt', default='prot', choices=['prot','nucl'], help='specify database type [default: prot] or [nucl]')
  parser.add_argument('-e', default=float(1E-5), type=float, help='specify the evalue cutoff [default: 1E-5]')
  parser.add_argument('-m', default='BLOSUM62', help='specify the BLAST+/ghostx/diamond matrix [default: BLOSUM62]')
  parser.add_argument('-n',default=int(1), help='specfiy number of threads to use')
  parser.add_argument('-op', default='.', help='specify output path [default: .]')
  parser.add_argument('-p', default=float(30), type=float, help='specify pident cutoff for BLAST+/ghostx/diamond result filtering [default: 30]')
  parser.add_argument('-l', default=float(80), type=float, help='specify alignment length cutoff for BLAST+/ghostx/diamond result filtering [default: 80]')
  parser.add_argument('-covs', default=float(0.0), type=float, help='specify subject coverage [alignment length / subject length] for BLAST+ result filtering [default: 0.0]; can only be applied with BLAST+ and diamond [see help file how to compile diamond to also provide qlen and slen]; ghostx does not provide qlen and slen yet.')
  parser.add_argument('-covq', default=float(0.0), type=float, help='specify query coverage [alignment length / query length] for BLAST+ result filtering [default: 0.0]; can only be applied with BLAST+ and diamond [see help file how to compile diamond to also provied qlen and slen]; ghostx does not provide qlen and slen yet.')
  parser.add_argument('-nt', default=int(250), type=int, help='specify number of max. targets to be considered by BLAST+/ghostx/diamond [default: 250].')
  parser.add_argument('-pid', default='static', choices=['evalue', 'static', 'rost1999'], help='choose pident filtering for BLAST+/ghostx/diamond result filtering; either as static [pident+length], evalue or rost1999')
  parser.add_argument('-sort', default='pure', choices=['pure', 'evl', 'pxl', 'selfscore'], help='choose reciprocal best hit sorting: either [default: pure] (number of max. targets is restricted to 50 and only first hit of those which fulfill the pid criteria will be processed); [evl] (hits that fulfill pid criteria will be sorted by evalue, any number of max. targets possible); [pxl] (hits that fulfill pid criteria will be sorted by pident*length, any number of max. targets possible) or [selfscore] (hits that fulfill pid criteria will be sorted by normalized selfscore provided via the -selfscore OPTION or directly calculated via the -selfblast option, any number of max. targets possible).')
  parser.add_argument('-selfscore', nargs=2, help='specify selfscore files for the query and the database species [tab seperated format (ID\tSCORE)]. You need to have two file one for the query and one for the database. Selfscore files can be produced either by using [step: self] or by including the [selfblast] option, in this case the selfscores will be calculated anyway.')
  parser.add_argument('-pre', default='False', choices=['True','False'], help='specify if pre-existing BLAST+/ghostx/diamond output files should be used; e.g. if BLAST+/ghostx/diamond output was produced elsewhere [default: False]. Please note that BLAST+ output needs to be in the format "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident gaps score"; diamond output needs to be in the format ""; ghostx in the current version does not support individual output format so we stick to standard output here. Also there exists a name convention "query input".vs."database input".')
  parser.add_argument('-prep', default='.', help='specify path with pre-existing BLAST+/ghostx/diamond output [default: .]')
  parser.add_argument('-out', default='out', help='specify output file')
  parser.add_argument('-v', default='False', choices=['True','False'], help='verbose output')
  args = parser.parse_args()
  if args.q is None and args.db is None:
    parser.print_help()
    sys.exit('\n#####\nExit program: Please specify sequence files.')
  if args.q is not None and args.db is None and args.selfblast=='False':
    parser.print_help()
    sys.exit('\n#####\nExit program: If you want to perform a selfblast within the query species please use the query species also as the database species.\nOtherwise please specify the sequence file for the database species.')
  if args.q is None and args.db is not None and args.selfblast=='False':
    parser.print_help()
    sys.exit('\n#####\nExit program: If you want to perform a selfblast within the database species please use the database species also as the query species.\nOtherwise please specify the sequence file for the query species.')
  if args.q is not None and args.db is not None and args.selfblast=='False':
    query_species = args.q
    database_species = args.db
  if args.q is not None and args.db is not None and args.selfblast=='True' and args.q == args.db:
    parser.print_help()
    sys.exit('\n#####\nExit program: If you want to perform a selfblast within one species please use the query species also as the database species and "omit the selfblast option in this case!".\nOtherwise please specify the sequence file for the query species.')
  if args.q is not None and args.db is not None and args.selfblast=='True':
    print 'within option used -> including selfblast\n'
    query_species = args.q
    database_species = args.db
  prog = args.prog
  blastpath = args.pp
  comparison_qd, comparison_dq, comparison_qq, comparison_dd = get_comparison_by_type(args.qt, args.dbt, args.selfblast)
  query_species_type = args.qt
  database_species_type = args.dbt
  selfblast = args.selfblast
  if args.out is None and args.step != '0' and args.step != 'self':
    parser.print_help()
    sys.exit('\nPlease specify output file')
  outfile = args.out
  option_pident = args.p
  option_length = args.l
  option_covs = args.covs
  option_covq = args.covq
  pident_method = args.pid
  sort_method = args.sort
  num_targets = args.nt
  if sort_method == 'pure':
    num_targets = 50
  step = args.step
  option_evalue = float(args.e)
  matrix = args.m
  num_threads = int(args.n)
  outpath = args.op
  precheck = args.pre
  prepath = args.prep
  option_verbose = args.v
  if prog == 'blast+':
    outfile_qd = outpath+'/'+query_species+'.vs.'+database_species+'.blast.out'
    outfile_dq = outpath+'/'+database_species+'.vs.'+query_species+'.blast.out'
    outfile_qq = outpath+'/'+query_species+'.vs.'+query_species+'.blast.out'
    outfile_dd = outpath+'/'+database_species+'.vs.'+database_species+'.blast.out'
    selfscore_q = outpath+'/'+query_species+'.blast.selfscore'
    selfscore_d = outpath+'/'+database_species+'.blast.selfscore'
  if prog == 'ghostx':
    outfile_qd = outpath+'/'+query_species+'.vs.'+database_species+'.ghostx.out'
    outfile_dq = outpath+'/'+database_species+'.vs.'+query_species+'.ghostx.out'
    outfile_qq = outpath+'/'+query_species+'.vs.'+query_species+'.ghostx.out'
    outfile_dd = outpath+'/'+database_species+'.vs.'+database_species+'.ghostx.out'
    selfscore_q = outpath+'/'+query_species+'.ghostx.selfscore'
    selfscore_d = outpath+'/'+database_species+'.ghostx.selfscore'
  if prog == 'diamond':
    outfile_qd = outpath+'/'+query_species+'.vs.'+database_species+'.diamond.out'
    outfile_dq = outpath+'/'+database_species+'.vs.'+query_species+'.diamond.out'
    outfile_qq = outpath+'/'+query_species+'.vs.'+query_species+'.diamond.out'
    outfile_dd = outpath+'/'+database_species+'.vs.'+database_species+'.diamond.out'
    selfscore_q = outpath+'/'+query_species+'.diamond.selfscore'
    selfscore_d = outpath+'/'+database_species+'.diamond.selfscore'
  if args.selfscore is None:
    selfscore = ''
  if args.selfscore is not None:
    selfscore = args.selfscore
    selfscore_q = args.selfscore[0]
    selfscore_d = args.selfscore[1]
  print '\ncommand arguments used:\n'
  print args
  return[prog, blastpath, query_species, database_species, query_species_type, database_species_type, selfblast, comparison_qd, comparison_dq, comparison_qq, comparison_dd, selfscore, outfile, option_pident, option_length, option_covs, option_covq, pident_method, sort_method, num_targets, step, option_evalue, matrix, num_threads, outpath, precheck, prepath, outfile_qd, outfile_dq, outfile_qq, outfile_dd, selfscore_q, selfscore_d, option_verbose]

def main():
  print 'start rbhplus\n'
  print 'parsing command options:\n'
  oprog, oblastpath, oquery_species, odatabase_species, oquery_species_type, odatabase_species_type, oselfblast, ocomparison_qd, ocomparison_dq, ocomparison_qq, ocomparison_dd, oselfscore, ooutfile, ooption_pident, ooption_length, ooption_covs, ooption_covq, opident_method, osort_method, onum_targets, ostep, ooption_evalue, omatrix, onum_threads, ooutpath, oprecheck, oprepath, ooutfile_qd, ooutfile_dq, ooutfile_qq, ooutfile_dd, oselfscore_q, oselfscore_d, ooption_verbose = commandoptions()
  abs_oquery_species=os.path.abspath(oquery_species)
  abs_odatabase_species=os.path.abspath(odatabase_species)
  abs_ooutfile_qd=os.path.abspath(ooutfile_qd)
  abs_ooutfile_dq=os.path.abspath(ooutfile_dq)
  abs_ooutfile_qq=os.path.abspath(ooutfile_qq)
  abs_ooutfile_dd=os.path.abspath(ooutfile_dd)
  abs_oselfscore_q=os.path.abspath(oselfscore_q)
  abs_oselfscore_d=os.path.abspath(oselfscore_d)
  abs_ooutfile=os.path.abspath(os.path.abspath(ooutpath)+'/'+ooutfile)
  abs_oprepath=os.path.abspath(oprepath)
  print 'step: %s' % (str(ostep))
  print 'query_species: %s' % (oquery_species)
  print 'database_species: %s'  % (odatabase_species)
  print 'comparison_qd: %s'  % (ocomparison_qd)
  print 'comparison_dq: %s'  % (ocomparison_dq)
  if oselfblast=='True':
    print 'comparison_qq: %s'  % (ocomparison_qq)
    print 'comparison_dd: %s'  % (ocomparison_dd)
  print 'evalue: %s'  % (str(ooption_evalue))
  print 'matrix: %s'  % (omatrix)
  print 'num_threads: %i'  % (onum_threads)
  print 'outpath: %s'  % (ooutpath)
  print 'prepath: %s'  % (oprepath)
  print 'precheck: %s'  % (str(oprecheck))
  print 'outfile_qd: %s'  % (ooutfile_qd)
  print 'outfile_dq: %s'  % (ooutfile_dq)
  if oselfblast=='True':
    print 'outfile_qq: %s'  % (ooutfile_qq)
    print 'outfile_dd: %s'  % (ooutfile_dd)
    print 'selfscore_q: %s' % (oselfscore_q)
    print 'selfscore_d: %s' % (oselfscore_d)
  if oselfblast=='False':
    if ostep=='self' and oquery_species==odatabase_species:
      print 'outfile_qq: %s'  % (ooutfile_qq)
      print 'selfscore_q: %s' % (oselfscore_q)
    if ostep=='self' and oquery_species!=odatabase_species:
      print 'outfile_qq: %s'  % (ooutfile_qq)
      print 'outfile_dd: %s'  % (ooutfile_dd)
      print 'selfscore_q: %s' % (oselfscore_q)
      print 'selfscore_d: %s' % (oselfscore_d)
  if oselfblast=='False' and oquery_species==odatabase_species and osort_method == 'selfscore':
    print 'selfscore_q: %s' % (oselfscore_q)
  print 'pident: %s' % (str(ooption_pident))
  print 'length: %s' % (str(ooption_length))
  print 'covs: %s' % (str(ooption_covs))
  print 'covq: %s' % (str(ooption_covq))
  print 'pident filter: %s' % (opident_method)
  print 'sort: %s' % (osort_method)
  print 'blastpath: %s'  % (oblastpath)
  print 'verbose: %s\n' % (ooption_verbose)
  if ostep=='0':
    print '\n#####\nSTEP0: building %s sequence databases:\n' % (oprog)
    if oquery_species==odatabase_species and oprog=='blast+':
      build_blast_db(oblastpath,abs_oquery_species,oquery_species_type,abs_oprepath,oprecheck)
    if oquery_species==odatabase_species and oprog=='ghostx':
      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,abs_oprepath,oprecheck)
#    if oquery_species==odatabase_species and oprog=='diamond':
#      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,abs_oprepath,oprecheck)
    if oquery_species!=odatabase_species and oprog=='blast+':
      build_blast_db(oblastpath,abs_oquery_species,oquery_species_type,abs_oprepath,oprecheck)
      build_blast_db(oblastpath,abs_odatabase_species,odatabase_species_type,abs_oprepath,oprecheck)
#    if oquery_species!=odatabase_species and oprog=='ghostx':
#      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,abs_oprepath,oprecheck)
#      build_ghostx_db(oblastpath,abs_odatabase_species,odatabase_species_type,abs_oprepath,oprecheck)
    if oquery_species!=odatabase_species and oprog=='diamond':
      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,abs_oprepath,oprecheck)
      build_ghostx_db(oblastpath,abs_odatabase_species,odatabase_species_type,abs_oprepath,oprecheck)
    ostep='exit'
  if ostep=='1':
    print '\n#####\nSTEP1: building %s sequence databases:\n' % (oprog)
    if oquery_species==odatabase_species and oprog=='blast+':
      build_blast_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
    if oquery_species==odatabase_species and oprog=='ghostx':
      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
#    if oquery_species==odatabase_species and oprog=='diamond':
#      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
    if oquery_species!=odatabase_species and oprog=='blast+':
      build_blast_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
      build_blast_db(oblastpath,abs_odatabase_species,odatabase_species_type,oprepath,oprecheck)
    if oquery_species!=odatabase_species and oprog=='ghostx':
      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
      build_ghostx_db(oblastpath,abs_odatabase_species,odatabase_species_type,oprepath,oprecheck)
    if oquery_species!=odatabase_species and oprog=='diamond':
      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
      build_ghostx_db(oblastpath,abs_odatabase_species,odatabase_species_type,oprepath,oprecheck)
    ostep='2'
  if ostep=='2':
    print '\n#####\nSTEP2: %s searches\n' % (oprog)
    if oselfblast=='False':
      if oquery_species==odatabase_species and oprog=='blast+':
        print '\n###\nBLAST+ q vs d'
        _qd = blast(oblastpath, ocomparison_qd, abs_oquery_species, abs_odatabase_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
        if _qd.check() and oprecheck=='True':
          print 'BLAST+ q vs d - already exists'
        else:
          print 'start BLAST+ q vs d'
          _qd.process()
          print 'finished BLAST+ q vs d'
      if oquery_species==odatabase_species and oprog=='ghostx':
        print '\n###\nghostx q vs d'
        _qd = ghostx(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
        if _qd.check() and oprecheck=='True':
          print 'ghostx q vs d - already exists'
        else:
          print 'start ghostx q vs d'
          _qd.process()
          print 'finished ghostx q vs d'
#      if oquery_species==odatabase_species and oprog=='diamond':
#        print '\n###\ndiamond q vs d'
#        _qd = diamond(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
#        if _qd.check() and oprecheck=='True':
#          print 'diamond q vs d - already exists'
#        else:
#          print 'start diamond q vs d'
#          _qd.process()
#          print 'finished diamond q vs d'
      if oquery_species!=odatabase_species and oprog=='blast+':
        print '\n###\nBLAST+ q vs d'
        _qd = blast(oblastpath, ocomparison_qd, abs_oquery_species, abs_odatabase_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
        if _qd.check() and oprecheck=='True':
          print 'BLAST+ q vs d - already exists'
        else:
          print 'start BLAST+ q vs d'
          _qd.process()
          print 'finished BLAST+ q vs d'
        print '\n###\nBLAST+ d vs q'
        _dq = blast(oblastpath, ocomparison_dq, abs_odatabase_species, abs_oquery_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_dq, abs_oprepath, oprecheck, ooutfile_dq)
        if _dq.check() and oprecheck=='True':
          print 'BLAST+ d vs q - already exists'
        else:
          print 'start BLAST+ d vs q'
          _dq.process()
          print 'finished BLAST+ d vs q'
      if oquery_species!=odatabase_species and oprog=='ghostx':
        print '\n###\nghostx q vs d'
        _qd = ghostx(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
        if _qd.check() and oprecheck=='True':
          print 'ghostx q vs d - already exists'
        else:
          print 'start ghostx q vs d'
          _qd.process()
          print 'finished ghostx q vs d'
        print '\n###\nghostx d vs q'
        _dq = ghostx(oblastpath, oprog, abs_odatabase_species, abs_oquery_species, odatabase_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_dq, abs_oprepath, oprecheck, ooutfile_dq)
        if _dq.check() and oprecheck=='True':
          print 'ghostx d vs q - already exists'
        else:
          print 'start ghostx d vs q'
          _dq.process()
          print 'finished ghostx d vs q'
#      if oquery_species!=odatabase_species and oprog=='diamond':
#        print '\n###\nghostx q vs d'
#        _qd = diamond(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
#        if _qd.check() and oprecheck=='True':
#          print 'diamond q vs d - already exists'
#        else:
#          print 'start diamond q vs d'
#          _qd.process()
#          print 'finished diamond q vs d'
#        print '\n###\ndiamond d vs q'
#        _dq = diamond(oblastpath, oprog, abs_odatabase_species, abs_oquery_species, odatabase_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_dq, abs_oprepath, oprecheck, ooutfile_dq)
#        if _dq.check() and oprecheck=='True':
#          print 'diamond d vs q - already exists'
#        else:
#          print 'start diamond d vs q'
#          _dq.process()
#          print 'finished diamond d vs q'
      ostep='3'
    if oselfblast=='True':
      if oprog=='blast+':
        print '\n###\nBLAST+ q vs d'
        _qd = blast(oblastpath, ocomparison_qd, abs_oquery_species, abs_odatabase_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
        if _qd.check() and oprecheck=='True':
          print 'BLAST+ q vs d - already exists'
        else:
          print 'start BLAST+ q vs d'
          _qd.process()
          print 'finished BLAST+ q vs d'
        print '\n###\nBLAST+ d vs q'
        _dq = blast(oblastpath, ocomparison_dq, abs_odatabase_species, abs_oquery_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_dq, abs_oprepath, oprecheck, ooutfile_dq)
        if _dq.check() and oprecheck=='True':
          print 'BLAST+ d vs q - already exists'
        else:
          print 'start BLAST+ d vs q'
          _dq.process()
          print 'finished BLAST+ d vs q'
        print '\n###\nBLAST+ q vs q'
        _qq = blast(oblastpath, ocomparison_qq, abs_oquery_species, abs_oquery_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_qq, abs_oprepath, oprecheck, ooutfile_qq)
        if _qq.check() and oprecheck=='True':
          print 'BLAST+ q vs q - already exists'
        else:
          print 'start BLAST+ q vs q'
          _qq.process()
          print 'finished BLAST+ q vs q'
        print '\n###\nBLAST+ d vs d'
        _dd = blast(oblastpath, ocomparison_dd, abs_odatabase_species, abs_odatabase_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_dd, abs_oprepath, oprecheck, ooutfile_dd)
        if _dd.check() and oprecheck=='True':
          print 'BLAST+ d vs d - already exists'
        else:
          print 'start BLAST+ d vs d'
          _dd.process()
          print 'finished BLAST+ d vs d'
      if oprog=='ghostx':
        print '\n###\nghostx q vs d'
        _qd = ghostx(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
        if _qd.check() and oprecheck=='True':
          print 'ghostx q vs d - already exists'
        else:
          print 'start ghostx q vs d'
          _qd.process()
          print 'finished ghostx q vs d'
        print '\n###\nghostx d vs q'
        _dq = ghostx(oblastpath, oprog, abs_odatabase_species, abs_oquery_species, odatabase_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_dq, abs_oprepath, oprecheck, ooutfile_dq)
        if _dq.check() and oprecheck=='True':
          print 'ghostx d vs q - already exists'
        else:
          print 'start ghostx d vs q'
          _dq.process()
          print 'finished ghostx d vs q'
        print '\n###\nghostx q vs q'
        _qq = ghostx(oblastpath, oprog, abs_oquery_species, abs_oquery_species, oquery_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_qq, abs_oprepath, oprecheck, ooutfile_qq)
        if _qq.check() and oprecheck=='True':
          print 'ghostx q vs q - already exists'
        else:
          print 'start ghostx q vs q'
          _qq.process()
          print 'finished ghostx q vs q'
        print '\n###\nghostx d vs d'
        _dd = ghostx(oblastpath, oprog, abs_odatabase_species, abs_odatabase_species, odatabase_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_dd, abs_oprepath, oprecheck, ooutfile_dd)
        if _dd.check() and oprecheck=='True':
          print 'ghostx d vs d - already exists'
        else:
          print 'start ghostx d vs d'
          _dd.process()
          print 'finished ghostx d vs d'
#      if oprog=='diamond':
#        print '\n###\nghostx q vs d'
#        _qd = diamond(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
#        if _qd.check() and oprecheck=='True':
#          print 'diamond q vs d - already exists'
#        else:
#          print 'start diamond q vs d'
#          _qd.process()
#          print 'finished diamond q vs d'
#        print '\n###\ndiamond d vs q'
#        _dq = diamond(oblastpath, oprog, abs_odatabase_species, abs_oquery_species, odatabase_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_dq, abs_oprepath, oprecheck, ooutfile_dq)
#        if _dq.check() and oprecheck=='True':
#          print 'diamond d vs q - already exists'
#        else:
#          print 'start diamond d vs q'
#          _dq.process()
#          print 'finished diamond d vs q'
#        print '\n###\ndiamond q vs q'
#        _qq = diamond(oblastpath, oprog, abs_oquery_species, abs_oquery_species, oquery_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_qq, abs_oprepath, oprecheck, ooutfile_qq)
#        if _qq.check() and oprecheck=='True':
#          print 'diamond q vs q - already exists'
#        else:
#          print 'start diamond q vs q'
#          _qq.process()
#          print 'finished diamond q vs q'
#        print '\n###\ndiamond d vs d'
#        _dd = diamond(oblastpath, oprog, abs_odatabase_species, abs_odatabase_species, odatabase_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_dd, abs_oprepath, oprecheck, ooutfile_dd)
#        if _dd.check() and oprecheck=='True':
#          print 'diamond d vs d - already exists'
#        else:
#          print 'start diamond d vs d'
#          _dd.process()
#          print 'finished diamond d vs d'
      ostep='3'
  if ostep=='3':
    print '\n#####\nSTEP3: parse %s output' % (oprog)
    if osort_method=='selfscore':
      if oselfscore!='':
        print("start reading selfscore table query species")
        selfscore_q_dict = selfscore_dict(abs_oselfscore_q, abs_oselfscore_q, 'score')
        selfscore_q_dict.parse()
        selfscore_q_dict_dict = selfscore_q_dict.dict
        print("start reading selfscore table database species")
        selfscore_d_dict = selfscore_dict(abs_oselfscore_d, abs_oselfscore_d, 'score')
        selfscore_d_dict.parse()
        selfscore_d_dict_dict = selfscore_d_dict.dict
      if oselfscore=='':
        if oselfblast=='False':
          if oquery_species==odatabase_species:
            print 'query species (%s) equals database species (%s)\n' % (oquery_species, odatabase_species)
            print 'extract %s selfscores' % (oprog)
            if oprog=='blast+':
              selfscore_q_dict = selfscore_dict(abs_ooutfile_qd, abs_oselfscore_q, 'blast+')
              selfscore_q_dict.parse()
              selfscore_q_dict.write()
              if oquery_species!=odatabase_species:
                parser.print_help()
                sys.exit('\n#####\nExit program: selfscores can not be extracted, if you want to do this please use the "-selfblast" option or directly "-step self".')
            if oprog=='ghostx':
              selfscore_q_dict = selfscore_dict(abs_ooutfile_qd, abs_oselfscore_q, 'ghostx')
              selfscore_q_dict.parse()
              selfscore_q_dict.write()
              if oquery_species!=odatabase_species:
                parser.print_help()
                sys.exit('\n#####\nExit program: selfscores can not be extracted, if you want to do this please use the "-selfblast" option or directly "-step self".')
#            if oprog=='diamond':
#              selfscore_q_dict = selfscore_dict(abs_ooutfile_qd, abs_oselfscore_q, 'diamond')
#              selfscore_q_dict.parse()
#              selfscore_q_dict.write()
#              if oquery_species!=odatabase_species:
#                parser.print_help()
#                sys.exit('\n#####\nExit program: selfscores can not be extracted, if you want to do this please use the "-selfblast" option or directly "-step self".')
            selfscore_q_dict_dict = selfscore_q_dict.dict
        if oselfblast=='True':
          print 'since selfblast was running the selfscores can now be extracted\n'
          if oprog=='blast+':
            print 'extract %s selfscores q' % (oprog)
            selfscore_q_dict = selfscore_dict(abs_ooutfile_qq, abs_oselfscore_q,'blast+')
            selfscore_q_dict.parse()
            selfscore_q_dict.write()
            print 'extract %s selfscores d' % (oprog)
            selfscore_d_dict = selfscore_dict(abs_ooutfile_dd, abs_oselfscore_d,'blast+')
            selfscore_d_dict.parse()
            selfscore_d_dict.write()
          if oprog=='ghostx':
            print 'extract %s selfscores q' % (oprog)
            selfscore_q_dict = selfscore_dict(abs_ooutfile_qq, abs_oselfscore_q,'ghostx')
            selfscore_q_dict.parse()
            selfscore_q_dict.write()
            print 'extract %s selfscores d' % (oprog)
            selfscore_d_dict = selfscore_dict(abs_ooutfile_dd, abs_oselfscore_d,'ghostx')
            selfscore_d_dict.parse()
            selfscore_d_dict.write()
#          if oprog=='diamond':
#            print 'extract %s selfscores q' % (oprog)
#            selfscore_q_dict = selfscore_dict(abs_ooutfile_qq, abs_oselfscore_q,'diamond')
#            selfscore_q_dict.parse()
#            selfscore_q_dict.write()
#            print 'extract %s selfscores d' % (oprog)
#            selfscore_d_dict = selfscore_dict(abs_ooutfile_dd, abs_oselfscore_d,'diamond')
#            selfscore_d_dict.parse()
#            selfscore_d_dict.write()
          selfscore_q_dict_dict = selfscore_q_dict.dict
          selfscore_d_dict_dict = selfscore_d_dict.dict
    if osort_method!='selfscore':
      selfscore_q_dict_dict = {}
      selfscore_d_dict_dict = {}
    #blast_dict_qd
    blast_dict_qd=blast_dict(abs_ooutfile_qd, osort_method, opident_method, ooption_pident, ooption_length, ooption_evalue, ooption_covs, ooption_covq, ocomparison_qd, oprog, selfscore_q_dict_dict, ooption_verbose)
    blast_dict_qd.parse()
    blast_dict_qd_dict=blast_dict_qd.dict
    if oquery_species==odatabase_species:
      ostep='4'
    if oquery_species!=odatabase_species:
      #blast_dict_dq
      blast_dict_dq=blast_dict(abs_ooutfile_dq, osort_method, opident_method, ooption_pident, ooption_length, ooption_evalue, ooption_covs, ooption_covq, ocomparison_dq, oprog, selfscore_d_dict_dict, ooption_verbose)
      blast_dict_dq.parse()
      blast_dict_dq_dict=blast_dict_dq.dict
      if oselfblast=='True':
        #blast_dict_qq
        blast_dict_qq=blast_dict(abs_ooutfile_qq, osort_method, opident_method, ooption_pident, ooption_length, ooption_evalue, ooption_covs, ooption_covq, ocomparison_qq, oprog, selfscore_q_dict_dict, ooption_verbose)
        blast_dict_qq.parse()
        blast_dict_qq_dict=blast_dict_qq.dict
        #blast_dict_dd
        blast_dict_dd=blast_dict(abs_ooutfile_dd, osort_method, opident_method, ooption_pident, ooption_length, ooption_evalue, ooption_covs, ooption_covq, ocomparison_dd, oprog, selfscore_d_dict_dict, ooption_verbose)
        blast_dict_dd.parse()
        blast_dict_dd_dict=blast_dict_dd.dict
      ostep='4'
  if ostep=='4':
    print '\n#####\nSTEP4: calculate RBHs\n'
    reciprocal_besthit_pairs_qq=[]
    reciprocal_besthit_pairs_dd=[]
    reciprocal_besthit_pairs_qd=[]
    if oselfblast=='False':
      if oquery_species==odatabase_species:
        print 'calculate query RBHs due to given options (p: %s, l: %s, evalue: %s, pid: %s, sort: %s, selfscore: %s)\n' % (ooption_pident,ooption_length,ooption_evalue,opident_method,osort_method,oselfscore)
        for aquery in blast_dict_qd_dict:
          blast_query_a = aquery
          blast_query_a_besthit = blast_dict_qd_dict[blast_query_a][0]
          blast_query_a_pident = blast_dict_qd_dict[blast_query_a][1]
          blast_query_a_length = blast_dict_qd_dict[blast_query_a][2]
          blast_query_a_evalue = blast_dict_qd_dict[blast_query_a][3]
          blast_query_a_bitscore = blast_dict_qd_dict[blast_query_a][4]
          blast_query_a_score = blast_dict_qd_dict[blast_query_a][5]
          blast_query_a_normscore = blast_dict_qd_dict[blast_query_a][6]
          blast_query_a_qlen = blast_dict_qd_dict[blast_query_a][7]
          blast_query_a_slen = blast_dict_qd_dict[blast_query_a][8]
          blast_query_a_nident = blast_dict_qd_dict[blast_query_a][9]
          blast_query_a_covs = blast_dict_qd_dict[blast_query_a][10]
          blast_query_a_covS = blast_dict_qd_dict[blast_query_a][11]
          blast_query_a_covns = blast_dict_qd_dict[blast_query_a][12]
          blast_query_a_covq = blast_dict_qd_dict[blast_query_a][13]
          blast_query_a_covQ = blast_dict_qd_dict[blast_query_a][14]
          blast_query_a_covnq = blast_dict_qd_dict[blast_query_a][15]
          reciprocal_besthit_pairs_qq.append([blast_query_a,blast_query_a_besthit,blast_query_a_pident,blast_query_a_length,blast_query_a_evalue,blast_query_a_bitscore,blast_query_a_score,blast_query_a_normscore,blast_query_a_qlen,blast_query_a_slen,blast_query_a_nident,blast_query_a_covs,blast_query_a_covS,blast_query_a_covns,blast_query_a_covq,blast_query_a_covQ,blast_query_a_covnq])
        with open(abs_ooutfile+'.qq','w') as handle:
          handle.write("blast_query_a\tblast_query_a_besthit\tblast_query_a_pident\tblast_query_a_length\tblast_query_a_evalue\tblast_query_a_bitscore\tblast_query_a_score\tblast_query_a_normscore\tblast_query_a_qlen\tblast_query_a_slen\tblast_query_a_nident\tblast_query_a_covs\tblast_query_a_covS\tblast_query_a_covns\tblast_query_a_covq\tblast_query_a_covQ\tblast_query_a_covnq\n")
          for record in reciprocal_besthit_pairs_qq:
            handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record[0],record[1],record[2],record[3],record[4],record[5],record[6],record[7],record[8],record[9],record[10],record[11],record[12],record[13],record[14],record[15],record[16])) 
      if oquery_species!=odatabase_species:
        print 'calculate query vs. database RBHs due to given options (p: %s, l: %s, evalue: %s, pid: %s, sort: %s, selfscore: %s)\n' % (ooption_pident,ooption_length,ooption_evalue,opident_method,osort_method,oselfscore)
        for aquery in blast_dict_qd_dict:
          blast_query_a = aquery
          blast_query_a_besthit = blast_dict_qd_dict[blast_query_a][0]
          blast_query_a_pident = blast_dict_qd_dict[blast_query_a][1]
          blast_query_a_length = blast_dict_qd_dict[blast_query_a][2]
          blast_query_a_evalue = blast_dict_qd_dict[blast_query_a][3]
          blast_query_a_bitscore = blast_dict_qd_dict[blast_query_a][4]
          blast_query_a_score = blast_dict_qd_dict[blast_query_a][5]
          blast_query_a_normscore = blast_dict_qd_dict[blast_query_a][6]
          blast_query_a_qlen = blast_dict_qd_dict[blast_query_a][7]
          blast_query_a_slen = blast_dict_qd_dict[blast_query_a][8]
          blast_query_a_nident = blast_dict_qd_dict[blast_query_a][9]
          blast_query_a_covs = blast_dict_qd_dict[blast_query_a][10]
          blast_query_a_covS = blast_dict_qd_dict[blast_query_a][11]
          blast_query_a_covns = blast_dict_qd_dict[blast_query_a][12]
          blast_query_a_covq = blast_dict_qd_dict[blast_query_a][13]
          blast_query_a_covQ = blast_dict_qd_dict[blast_query_a][14]
          blast_query_a_covnq = blast_dict_qd_dict[blast_query_a][15]
          blast_query_b = blast_query_a_besthit
          if ooption_verbose=='True':
            if blast_query_b not in blast_dict_dq_dict:
              print "%s not in dictionary" % (blast_query_b)
          if blast_query_b in blast_dict_dq_dict:
            blast_query_b_besthit = blast_dict_dq_dict[blast_query_b][0]
            blast_query_b_pident = blast_dict_dq_dict[blast_query_b][1]
            blast_query_b_length = blast_dict_dq_dict[blast_query_b][2]
            blast_query_b_evalue = blast_dict_dq_dict[blast_query_b][3]
            blast_query_b_bitscore = blast_dict_dq_dict[blast_query_b][4]
            blast_query_b_score = blast_dict_dq_dict[blast_query_b][5]
            blast_query_b_normscore = blast_dict_dq_dict[blast_query_b][6]
            blast_query_b_qlen = blast_dict_dq_dict[blast_query_b][7]
            blast_query_b_slen = blast_dict_dq_dict[blast_query_b][8]
            blast_query_b_nident = blast_dict_dq_dict[blast_query_b][9]
            blast_query_b_covs = blast_dict_dq_dict[blast_query_b][10]
            blast_query_b_covS = blast_dict_dq_dict[blast_query_b][11]
            blast_query_b_covns = blast_dict_dq_dict[blast_query_b][12]
            blast_query_b_covq = blast_dict_dq_dict[blast_query_b][13]
            blast_query_b_covQ = blast_dict_dq_dict[blast_query_b][14]
            blast_query_b_covnq = blast_dict_dq_dict[blast_query_b][15]
            if blast_query_b_besthit == blast_query_a:
              reciprocal_besthit_pairs_qd.append([blast_query_a,blast_query_a_besthit,blast_query_a_pident,blast_query_b_pident,blast_query_a_length,blast_query_b_length,blast_query_a_evalue,blast_query_b_evalue,blast_query_a_bitscore,blast_query_b_bitscore,blast_query_a_score,blast_query_b_score,blast_query_a_normscore,blast_query_b_normscore,blast_query_a_qlen,blast_query_b_qlen,blast_query_a_slen,blast_query_b_slen,blast_query_a_nident,blast_query_b_nident,blast_query_a_covs,blast_query_b_covs,blast_query_a_covS,blast_query_b_covS,blast_query_a_covns,blast_query_b_covns,blast_query_a_covq,blast_query_b_covq,blast_query_a_covQ,blast_query_b_covQ,blast_query_a_covnq,blast_query_b_covnq])
        with open(abs_ooutfile+'.qd','w') as handle:
          handle.write("blast_query_a\tblast_query_a_besthit\tblast_query_a_pident\tblast_query_b_pident\tblast_query_a_length\tblast_query_b_length\tblast_query_a_evalue\tblast_query_b_evalue\tblast_query_a_bitscore\tblast_query_b_bitscore\tblast_query_a_score\tblast_query_b_score\tblast_query_a_normscore\tblast_query_b_normscore\tblast_query_a_qlen\tblast_query_b_qlen\tblast_query_a_slen\tblast_query_b_slen\tblast_query_a_covs\tblast_query_b_covs\tblast_query_a_covS\tblast_query_b_covS\tblast_query_a_covns\tblast_query_b_covns\tblast_query_a_covq\tblast_query_b_covq\tblast_query_a_covQ\tblast_query_b_covQ\tblast_query_a_covnq\tblast_query_b_covnq\n")
          for record in reciprocal_besthit_pairs_qd:
            handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record[0],record[1],record[2],record[3],record[4],record[5],record[6],record[7],record[8],record[9],record[10],record[11],record[12],record[13],record[14],record[15],record[16],record[17],record[18],record[19],record[20],record[21],record[22],record[23],record[24],record[25],record[26],record[27],record[28],record[29],record[30],record[31]))
    if oselfblast=='True':
      print 'calculate query vs. query RBHs due to given options (p: %s, l: %s, evalue: %s, pid: %s, sort: %s, selfscore: %s)\n' % (ooption_pident,ooption_length,ooption_evalue,opident_method,osort_method,oselfscore)
      for aquery in blast_dict_qq_dict:
        blast_query_a = aquery
        blast_query_a_besthit = blast_dict_qq_dict[blast_query_a][0]
        blast_query_a_pident = blast_dict_qq_dict[blast_query_a][1]
        blast_query_a_length = blast_dict_qq_dict[blast_query_a][2]
        blast_query_a_evalue = blast_dict_qq_dict[blast_query_a][3]
        blast_query_a_bitscore = blast_dict_qq_dict[blast_query_a][4]
        blast_query_a_score = blast_dict_qq_dict[blast_query_a][5]
        blast_query_a_normscore = blast_dict_qq_dict[blast_query_a][6]
        blast_query_a_qlen = blast_dict_qq_dict[blast_query_a][7]
        blast_query_a_slen = blast_dict_qq_dict[blast_query_a][8]
        blast_query_a_nident = blast_dict_qq_dict[blast_query_a][9]
        blast_query_a_covs = blast_dict_qq_dict[blast_query_a][10]
        blast_query_a_covS = blast_dict_qq_dict[blast_query_a][11]
        blast_query_a_covns = blast_dict_qq_dict[blast_query_a][12]
        blast_query_a_covq = blast_dict_qq_dict[blast_query_a][13]
        blast_query_a_covQ = blast_dict_qq_dict[blast_query_a][14]
        blast_query_a_covnq = blast_dict_qq_dict[blast_query_a][15]
        reciprocal_besthit_pairs_qq.append([blast_query_a,blast_query_a_besthit,blast_query_a_pident,blast_query_a_length,blast_query_a_evalue,blast_query_a_bitscore,blast_query_a_score,blast_query_a_normscore,blast_query_a_qlen,blast_query_a_slen,blast_query_a_nident,blast_query_a_covs,blast_query_a_covS,blast_query_a_covns,blast_query_a_covq,blast_query_a_covQ,blast_query_a_covnq])
      with open(abs_ooutfile+'.qq','w') as handle:
        handle.write("blast_query_a\tblast_query_a_besthit\tblast_query_a_pident\tblast_query_a_length\tblast_query_a_evalue\tblast_query_a_bitscore\tblast_query_a_score\tblast_query_a_normscore\tblast_query_a_qlen\tblast_query_a_slen\tblast_query_a_nident\tblast_query_a_covs\tblast_query_a_covS\tblast_query_a_covns\tblast_query_a_covq\tblast_query_a_covQ\tblast_query_a_covnq\n")
        for record in reciprocal_besthit_pairs_qq:
          handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record[0],record[1],record[2],record[3],record[4],record[5],record[6],record[7],record[8],record[9],record[10],record[11],record[12],record[13],record[14],record[15],record[16]))
      print 'calculate database vs. database RBHs due to given options (p: %s, l: %s, evalue: %s, pid: %s, sort: %s, selfscore: %s)\n' % (ooption_pident,ooption_length,ooption_evalue,opident_method,osort_method,oselfscore)
      for bquery in blast_dict_dd_dict:
        blast_query_b = bquery
        blast_query_b_besthit = blast_dict_dd_dict[blast_query_b][0]
        blast_query_b_pident = blast_dict_dd_dict[blast_query_b][1]
        blast_query_b_length = blast_dict_dd_dict[blast_query_b][2]
        blast_query_b_evalue = blast_dict_dd_dict[blast_query_b][3]
        blast_query_b_bitscore = blast_dict_dd_dict[blast_query_b][4]
        blast_query_b_score = blast_dict_dd_dict[blast_query_b][5]
        blast_query_b_normscore = blast_dict_dd_dict[blast_query_b][6]
        blast_query_b_qlen = blast_dict_dd_dict[blast_query_b][7]
        blast_query_b_slen = blast_dict_dd_dict[blast_query_b][8]
        blast_query_b_nident = blast_dict_dd_dict[blast_query_b][9]
        blast_query_b_covs = blast_dict_dd_dict[blast_query_b][10]
        blast_query_b_covS = blast_dict_dd_dict[blast_query_b][11]
        blast_query_b_covns = blast_dict_dd_dict[blast_query_b][12]
        blast_query_b_covq = blast_dict_dd_dict[blast_query_b][13]
        blast_query_b_covQ = blast_dict_dd_dict[blast_query_b][14]
        blast_query_b_covnq = blast_dict_dd_dict[blast_query_b][15]
        reciprocal_besthit_pairs_dd.append([blast_query_b,blast_query_b_besthit,blast_query_b_pident,blast_query_b_length,blast_query_b_evalue,blast_query_b_bitscore,blast_query_b_score,blast_query_b_normscore,blast_query_b_qlen,blast_query_b_slen,blast_query_b_nident,blast_query_b_covs,blast_query_b_covS,blast_query_b_covns,blast_query_b_covq,blast_query_b_covQ,blast_query_b_covnq])
      with open(abs_ooutfile+'.dd','w') as handle:
        handle.write("blast_query_b\tblast_query_b_besthit\tblast_query_b_pident\tblast_query_b_length\tblast_query_b_evalue\tblast_query_b_bitscore\tblast_query_b_score\tblast_query_b_normscore\tblast_query_b_qlen\tblast_query_b_slen\tblast_query_b_nident\tblast_query_b_covs\tblast_query_b_covS\tblast_query_b_covns\tblast_query_b_covq\tblast_query_b_covQ\tblast_query_b_covnq\n")
        for record in reciprocal_besthit_pairs_dd:
          handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record[0],record[1],record[2],record[3],record[4],record[5],record[6],record[7],record[8],record[9],record[10],record[11],record[12],record[13],record[14],record[15],record[16]))
      print 'calculate query vs. database RBHs due to given options (p: %s, l: %s, evalue: %s, pid: %s, sort: %s, selfscore: %s)\n' % (ooption_pident,ooption_length,ooption_evalue,opident_method,osort_method,oselfscore)
      for aquery in blast_dict_qd_dict:
        blast_query_a = aquery
        blast_query_a_besthit = blast_dict_qd_dict[blast_query_a][0]
        blast_query_a_pident = blast_dict_qd_dict[blast_query_a][1]
        blast_query_a_length = blast_dict_qd_dict[blast_query_a][2]
        blast_query_a_evalue = blast_dict_qd_dict[blast_query_a][3]
        blast_query_a_bitscore = blast_dict_qd_dict[blast_query_a][4]
        blast_query_a_score = blast_dict_qd_dict[blast_query_a][5]
        blast_query_a_normscore = blast_dict_qd_dict[blast_query_a][6]
        blast_query_a_qlen = blast_dict_qd_dict[blast_query_a][7]
        blast_query_a_slen = blast_dict_qd_dict[blast_query_a][8]
        blast_query_a_nident = blast_dict_qd_dict[blast_query_a][9]
        blast_query_a_covs = blast_dict_qd_dict[blast_query_a][10]
        blast_query_a_covS = blast_dict_qd_dict[blast_query_a][11]
        blast_query_a_covns = blast_dict_qd_dict[blast_query_a][12]
        blast_query_a_covq = blast_dict_qd_dict[blast_query_a][13]
        blast_query_a_covQ = blast_dict_qd_dict[blast_query_a][14]
        blast_query_a_covnq = blast_dict_qd_dict[blast_query_a][15]
        blast_query_b = blast_query_a_besthit
        if ooption_verbose=='True':
          if blast_query_b not in blast_dict_dq_dict:
            print "%s not in dictionary" % (blast_query_b)
        if blast_query_b in blast_dict_dq_dict:
          blast_query_b_besthit = blast_dict_dq_dict[blast_query_b][0]
          blast_query_b_pident = blast_dict_dq_dict[blast_query_b][1]
          blast_query_b_length = blast_dict_dq_dict[blast_query_b][2]
          blast_query_b_evalue = blast_dict_dq_dict[blast_query_b][3]
          blast_query_b_bitscore = blast_dict_dq_dict[blast_query_b][4]
          blast_query_b_score = blast_dict_dq_dict[blast_query_b][5]
          blast_query_b_normscore = blast_dict_dq_dict[blast_query_b][6]
          blast_query_b_qlen = blast_dict_dq_dict[blast_query_b][7]
          blast_query_b_slen = blast_dict_dq_dict[blast_query_b][8]
          blast_query_b_nident = blast_dict_dq_dict[blast_query_b][9]
          blast_query_b_covs = blast_dict_dq_dict[blast_query_b][10]
          blast_query_b_covS = blast_dict_dq_dict[blast_query_b][11]
          blast_query_b_covns = blast_dict_dq_dict[blast_query_b][12]
          blast_query_b_covq = blast_dict_dq_dict[blast_query_b][13]
          blast_query_b_covQ = blast_dict_dq_dict[blast_query_b][14]
          blast_query_b_covnq = blast_dict_dq_dict[blast_query_b][15]
          if blast_query_b_besthit == blast_query_a:
            reciprocal_besthit_pairs_qd.append([blast_query_a,blast_query_a_besthit,blast_query_a_pident,blast_query_b_pident,blast_query_a_length,blast_query_b_length,blast_query_a_evalue,blast_query_b_evalue,blast_query_a_bitscore,blast_query_b_bitscore,blast_query_a_score,blast_query_b_score,blast_query_a_normscore,blast_query_b_normscore,blast_query_a_qlen,blast_query_b_qlen,blast_query_a_slen,blast_query_b_slen,blast_query_a_nident,blast_query_b_nident,blast_query_a_covs,blast_query_b_covs,blast_query_a_covS,blast_query_b_covS,blast_query_a_covns,blast_query_b_covns,blast_query_a_covq,blast_query_b_covq,blast_query_a_covQ,blast_query_b_covQ,blast_query_a_covnq,blast_query_b_covnq])
      with open(abs_ooutfile+'.qd','w') as handle:
        handle.write("blast_query_a\tblast_query_a_besthit\tblast_query_a_pident\tblast_query_b_pident\tblast_query_a_length\tblast_query_b_length\tblast_query_a_evalue\tblast_query_b_evalue\tblast_query_a_bitscore\tblast_query_b_bitscore\tblast_query_a_score\tblast_query_b_score\tblast_query_a_normscore\tblast_query_b_normscore\tblast_query_a_qlen\tblast_query_b_qlen\tblast_query_a_slen\tblast_query_b_slen\tblast_query_a_covs\tblast_query_b_covs\tblast_query_a_covS\tblast_query_b_covS\tblast_query_a_covns\tblast_query_b_covns\tblast_query_a_covq\tblast_query_b_covq\tblast_query_a_covQ\tblast_query_b_covQ\tblast_query_a_covnq\tblast_query_b_covnq\n")
        for record in reciprocal_besthit_pairs_qd:
          handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record[0],record[1],record[2],record[3],record[4],record[5],record[6],record[7],record[8],record[9],record[10],record[11],record[12],record[13],record[14],record[15],record[16],record[17],record[18],record[19],record[20],record[21],record[22],record[23],record[24],record[25],record[26],record[27],record[28],record[29],record[30],record[31]))
  if ostep=='self':
    print '\n#####\nSTEPself: selfblast BLAST+/ghostx search with selfscore extraction\n'
    print '\n#####\nSTEP1: building %s sequence databases:\n' % (oprog)
    if oquery_species==odatabase_species and oprog=='blast+':
      build_blast_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
    if oquery_species==odatabase_species and oprog=='ghostx':
      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
#    if oquery_species==odatabase_species and oprog=='diamond':
#      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
    if oquery_species!=odatabase_species and oprog=='blast+':
      build_blast_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
      build_blast_db(oblastpath,abs_odatabase_species,odatabase_species_type,oprepath,oprecheck)
    if oquery_species!=odatabase_species and oprog=='ghostx':
      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
      build_ghostx_db(oblastpath,abs_odatabase_species,odatabase_species_type,oprepath,oprecheck)
#    if oquery_species!=odatabase_species and oprog=='diamond':
#      build_ghostx_db(oblastpath,abs_oquery_species,oquery_species_type,oprepath,oprecheck)
#      build_ghostx_db(oblastpath,abs_odatabase_species,odatabase_species_type,oprepath,oprecheck)
    print '\n#####\nSTEP2: %s searches\n' % (oprog)
    if oquery_species==odatabase_species and oprog=='blast+':
      print '\n###\nBLAST+ q vs d'
      _qd = blast(oblastpath, ocomparison_qd, abs_oquery_species, abs_odatabase_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
      if _qd.check() and oprecheck=='True':
        print 'BLAST+ q vs d - already exists'
      else:
        print 'start BLAST+ q vs d'
        _qd.process()
        print 'finished BLAST+ q vs d'
    if oquery_species==odatabase_species and oprog=='ghostx':
      print '\n###\nghostx q vs d'
      _qd = ghostx(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
      if _qd.check() and oprecheck=='True':
        print 'ghostx q vs d - already exists'
      else:
        print 'start ghostx q vs d'
        _qd.process()
        print 'finished ghostx q vs d'
#    if oquery_species==odatabase_species and oprog=='diamond':
#      print '\n###\ndiamond q vs d'
#      _qd = diamond(oblastpath, oprog, abs_oquery_species, abs_odatabase_species, oquery_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_qd, abs_oprepath, oprecheck, ooutfile_qd)
#      if _qd.check() and oprecheck=='True':
#        print 'diamond q vs d - already exists'
#      else:
#        print 'start diamond q vs d'
#        _qd.process()
#        print 'finished diamond q vs d'
    if oquery_species!=odatabase_species and oprog=='blast+':
      print '\n###\nBLAST+ q vs q'
      _qq = blast(oblastpath, ocomparison_qq, abs_oquery_species, abs_oquery_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_qq, abs_oprepath, oprecheck, ooutfile_qq)
      if _qq.check() and oprecheck=='True':
        print 'BLAST+ q vs q - already exists'
      else:
        print 'start BLAST+ q vs q'
        _qq.process()
        print 'finished BLAST+ q vs q'
      print '\n###\nBLAST+ d vs d'
      _dd = blast(oblastpath, ocomparison_dd, abs_odatabase_species, abs_odatabase_species, ooption_evalue, omatrix, onum_targets, onum_threads, abs_ooutfile_dd, abs_oprepath, oprecheck, ooutfile_dd)
      if _dd.check() and oprecheck=='True':
        print 'BLAST+ d vs d - already exists'
      else:
        print 'start BLAST+ d vs d'
        _dd.process()
        print 'finished BLAST+ d vs d'
    if oquery_species!=odatabase_species and oprog=='ghostx':
      print '\n###\nghostx q vs q'
      _qq = ghostx(oblastpath, oprog, abs_oquery_species, abs_oquery_species, oquery_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_qq, abs_oprepath, oprecheck, ooutfile_qq)
      if _qq.check() and oprecheck=='True':
        print 'ghostx q vs q - already exists'
      else:
        print 'start ghostx q vs q'
        _qq.process()
        print 'finished ghostx q vs q'
      print '\n###\nghostx d vs d'
      _dd = ghostx(oblastpath, oprog, abs_odatabase_species, abs_odatabase_species, odatabase_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_dd, abs_oprepath, oprecheck, ooutfile_dd)
      if _dd.check() and oprecheck=='True':
        print 'ghostx d vs d - already exists'
      else:
        print 'start ghostx d vs d'
        _dd.process()
        print 'finished ghostx d vs d'
#    if oquery_species!=odatabase_species and oprog=='diamond':
#      print '\n###\ndiamond q vs q'
#      _qq = diamond(oblastpath, oprog, abs_oquery_species, abs_oquery_species, oquery_species_type, oquery_species_type, onum_targets, onum_threads, abs_ooutfile_qq, abs_oprepath, oprecheck, ooutfile_qq)
#      if _qq.check() and oprecheck=='True':
#        print 'diamond q vs q - already exists'
#      else:
#        print 'start diamond q vs q'
#        _qq.process()
#        print 'finished diamond q vs q'
#      print '\n###\ndiamond d vs d'
#      _dd = diamond(oblastpath, oprog, abs_odatabase_species, abs_odatabase_species, odatabase_species_type, odatabase_species_type, onum_targets, onum_threads, abs_ooutfile_dd, abs_oprepath, oprecheck, ooutfile_dd)
#      if _dd.check() and oprecheck=='True':
#        print 'diamond d vs d - already exists'
#      else:
#        print 'start diamond d vs d'
#        _dd.process()
#        print 'finished diamond d vs d'
    if oquery_species==odatabase_species:
      print 'query species (%s) equals database species (%s)\n' % (oquery_species, odatabase_species)
      print 'extract %s selfscores' % (oprog)
      if oprog=='blast+':
        selfscore_q_dict = selfscore_dict(abs_ooutfile_qd, abs_oselfscore_q, 'blast+')
        selfscore_q_dict.parse()
        selfscore_q_dict.write()
      if oprog=='ghostx':
        selfscore_q_dict = selfscore_dict(abs_ooutfile_qd, abs_oselfscore_q, 'ghostx')
        selfscore_q_dict.parse()
        selfscore_q_dict.write()
#      if oprog=='diamond':
#        selfscore_q_dict = selfscore_dict(abs_ooutfile_qd, abs_oselfscore_q, 'diamond')
#        selfscore_q_dict.parse()
#        selfscore_q_dict.write()
    if oquery_species!=odatabase_species:
      if oprog=='blast+':
        print 'extract %s selfscores q' % (oprog)
        selfscore_q_dict = selfscore_dict(abs_ooutfile_qq, abs_oselfscore_q,'blast+')
        selfscore_q_dict.parse()
        selfscore_q_dict.write()
        print 'extract %s selfscores d' % (oprog)
        selfscore_d_dict = selfscore_dict(abs_ooutfile_dd, abs_oselfscore_d,'blast+')
        selfscore_d_dict.parse()
        selfscore_d_dict.write()
      if oprog=='ghostx':
        print 'extract %s selfscores q' % (oprog)
        selfscore_q_dict = selfscore_dict(abs_ooutfile_qq, abs_oselfscore_q,'ghostx')
        selfscore_q_dict.parse()
        selfscore_q_dict.write()
        print 'extract selfscores d'
        selfscore_d_dict = selfscore_dict(abs_ooutfile_dd, abs_oselfscore_d,'ghostx')
        selfscore_d_dict.parse()
        selfscore_d_dict.write()
#      if oprog=='diamond':
#        print 'extract %s selfscores q' % (oprog)
#        selfscore_q_dict = selfscore_dict(abs_ooutfile_qq, abs_oselfscore_q,'diamond')
#        selfscore_q_dict.parse()
#        selfscore_q_dict.write()
#        print 'extract selfscores d'
#        selfscore_d_dict = selfscore_dict(abs_ooutfile_dd, abs_oselfscore_d,'diamond')
#        selfscore_d_dict.parse()
#        selfscore_d_dict.write()
    ostep='finished'
  if ostep=='finished':
    print '\n#####\nFinished'

def build_blast_db(path, species, species_type, prepath, precheck):
  if precheck=='False':
    print 'overwrite existing BLAST+ sequeunce files since "-pre" option is not set\n'
    print 'building %s sequence database\n' % (species)
    os.system(path+'makeblastdb -in '+species+' -dbtype '+species_type)
  if precheck=='True':
    if species_type == 'prot':
      if os.path.isfile(species+'.phr') and os.path.isfile(species+'.pin') and os.path.isfile(species+'.psq'):
        print 'skip building %s sequence database - already exists\n' % (species)                    
      else:
        print 'building %s sequence database\n' % (species)
        os.system(path+'makeblastdb -in '+species+' -dbtype '+species_type)
    if species_type == 'nucl':
      if os.path.isfile(species+'.nhr') and os.path.isfile(species+'.nin') and os.path.isfile(species+'.nsq'):
        print 'skip building %s sequence database - already exists\n' % (species)
      else:
        print 'building %s sequence database - already exists\n' % (species)
        os.system(path+'makeblastdb -in '+species+' -dbtype '+species_type)

def get_species_type_ghostx(species_type):
  if species_type=='prot':
    return 'p'
  if species_type=='nucl':
    return 'd'

#def get_species_type_diamond(species_type):
#  if species_type=='prot':
#    return 'p'
#  if species_type=='nucl':
#    return 'd'

def build_ghostx_db(path, species, species_type, prepath, precheck):
  if precheck=='False':
    print 'overwrite existing ghostx sequeunce files since "-pre" option is not set\n'
    print 'building %s sequence database\n' % (species)
    os.system(path+'ghostx db -i '+species+' -t '+get_species_type_ghostx(species_type)+' -o '+species)
  if precheck=='True':
    if species_type == 'prot':
      if os.path.isfile(species+'.inf'):
        print 'skip building %s sequence database - already exists\n' % (species)
      else:
        print 'building %s sequence database\n' % (species)
        os.system(path+'ghostx db -i '+species+' -t '+get_species_type_ghostx(species_type)+' -o '+species)
    if species_type == 'nucl':
      if os.path.isfile(species+'.inf'):
        print 'skip building %s sequence database\n' % (species)
      else:
        print 'building %s sequence database - already exists\n' % (species)
        os.system(path+'ghostx db -i '+species+' -t '+get_species_type_ghostx(species_type)+' -o '+species)

#def build_diamond_db(path, species, species_type, prepath, precheck):
#  if precheck=='False':
#    print 'overwrite existing diamond sequeunce files since "-pre" option is not set\n'
#    print 'building %s sequence database\n' % (species)
#    os.system(path+'diamond makedb -i '+species+' -t '+get_species_type_diamond(species_type)+' -o '+species)
#  if precheck=='True':
#    if species_type == 'prot':
#      if os.path.isfile(species+'.inf'):
#        print 'skip building %s sequence database - already exists\n' % (species)
#      else:
#        print 'building %s sequence database\n' % (species)
#        os.system(path+'diamond db -i '+species+' -t '+get_species_type_ghostx(species_type)+' -o '+species)
#    if species_type == 'nucl':
#      if os.path.isfile(species+'.inf'):
#        print 'skip building %s sequence database\n' % (species)
#      else:
#        print 'building %s sequence database - already exists\n' % (species)
#        os.system(path+'diamond makedb -i '+species+' -t '+get_species_type_ghostx(species_type)+' -o '+species)

def get_pident_by_length(x):
  if x<=11:
    return float(100)
  if x<=450:
    return 480*float(x)**(-0.32*(1+math.exp(-float(x)/float(1000))))
  if x>450:
    return 19.5

#TODO: each blast_dict needs two selfscore files [query and database selfscore]
class blast_dict(object):
  def __init__(self, infile, sort_method, pident_method, option_pident, option_length, option_evalue, option_covs, option_covq, comparison, form, selfscore, verbose):
    self.infile = infile
    self.sort_method = sort_method
    self.pident_method = pident_method
    self.option_pident = float(option_pident)
    self.option_length = float(option_length)
    self.option_evalue = float(option_evalue)
    self.option_covs = float(option_covs)
    self.option_covq = float(option_covq)
    self.comparison = comparison
    self.form = form
    self.selfscore = selfscore
    self.verbose = verbose
    self.dict = {}
    
  def parse(self):
    with open(self.infile,'rU') as inhandle:
        for line in inhandle:
          parts = line.strip().split('\t')
          if self.form=='blast+':
            qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,nident,gaps,score = parts
            score = float(score)
            evalue = float(evalue)
            bitscore = float(bitscore)
            pident = float(pident)
            length = float(length)
            if self.comparison=='blastp' or self.comparison=='blastn':
              qlen = float(qlen)
              slen = float(slen)
            if self.comparison=='blastx':
              qlen = float(qlen)/3
              slen = float(slen)
            if self.comparison=='tblastn':
              qlen = float(qlen)
              slen = float(slen)/3
            if self.comparison=='tblastx':
              qlen = float(qlen)/3
              slen = float(slen)/3
            nident = float(nident)
            gaps = float(gaps)
            score = float(score)
            covs = length/slen
            covS = qlen/slen
            covns = nident/slen
            covq = length/qlen
            covQ = slen/qlen
            covnq = nident/qlen
          if self.form=='ghostx':
            qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = parts
            evalue = float(evalue)
            bitscore = float(bitscore)
            pident = float(pident)
            length = float(length)
            score = bitscore
            #if self.comparison=='blastp' or self.comparison=='blastn':
              #qlen = float(qlen)
              #slen = float(slen)
            #if self.comparison=='blastx':
              #qlen = float(qlen)/3
              #slen = float(slen)
            #if self.comparison=='tblastn':
              #qlen = float(qlen)
              #slen = float(slen)/3
            #if self.comparison=='tblastx':
              #qlen = float(qlen)/3
              #slen = float(slen)/3
            qlen = float(1)
            slen = float(1)
            nident = float(-999)
            covs = float(999)
            covS = float(999)
            covns = float(999)
            covq = float(999)
            covQ = float(999)
            covnq = float(999)
#          if self.form=='diamond':
#            qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = parts
#            evalue = float(evalue)
#            bitscore = float(bitscore)
#            pident = float(pident)
#            length = float(length)
#            score = bitscore
            #if self.comparison=='blastp' or self.comparison=='blastn':
              #qlen = float(qlen)
              #slen = float(slen)
            #if self.comparison=='blastx':
              #qlen = float(qlen)/3
              #slen = float(slen)
            #if self.comparison=='tblastn':
              #qlen = float(qlen)
              #slen = float(slen)/3
            #if self.comparison=='tblastx':
              #qlen = float(qlen)/3
              #slen = float(slen)/3
#            qlen = float(1)
#            slen = float(1)
#            nident = float(-999)
#            covs = float(999)
#            covS = float(999)
#            covns = float(999)
#            covq = float(999)
#            covQ = float(999)
#            covnq = float(999)
          if qseqid==sseqid:
            continue
          if self.sort_method=="pure":
            if self.pident_method=="evalue":
              if evalue<self.option_evalue and qseqid in self.dict:
                continue
              if evalue<self.option_evalue and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
            if self.pident_method=="static":
              if pident>=self.option_pident and length>=self.option_length and qseqid in self.dict:
                continue
              if pident>=self.option_pident and length>=self.option_length and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
            if self.pident_method=="rost1999":
              if pident>=get_pident_by_length(length) and qseqid in self.dict:
                continue
              if pident>=get_pident_by_length(length) and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
          if self.sort_method=="pxl":
              if self.pident_method=="evalue":
                if evalue<self.option_evalue and qseqid in self.dict and self.dict[qseqid][1]*self.dict[qseqid][2]<pident*length and covs>=self.option_covs and covq>=self.option_covq:
                  if self.verbose=='True':
                    print qseqid
                    print "%s %f is larger than %s %f" % (sseqid,pident*length,self.dict[qseqid][0],self.dict[qseqid][1]*self.dict[qseqid][2])
                  self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
                if evalue<self.option_evalue and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                  self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if self.pident_method=="static":
                if pident>=self.option_pident and length>=self.option_length and qseqid in self.dict and self.dict[qseqid][1]*self.dict[qseqid][2]<pident*length and covs>=self.option_covs and covq>=self.option_covq:
                  if self.verbose=='True':
                    print qseqid
                    print "%s %f is larger than %s %f" % (sseqid,pident*length,self.dict[qseqid][0],self.dict[qseqid][1]*self.dict[qseqid][2])
                  self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
                if pident>=self.option_pident and length>=self.option_length and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                  self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if self.pident_method=="rost1999":
                if pident>=get_pident_by_length(length) and qseqid in self.dict and self.dict[qseqid][1]*self.dict[qseqid][2]<pident*length and covs>=self.option_covs and covq>=self.option_covq:
                  if self.verbose=='True':
                    print qseqid
                    print "%s %f is larger than %s %f" % (sseqid,pident*length,self.dict[qseqid][0],self.dict[qseqid][1]*self.dict[qseqid][2])
                  self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
                if pident>=get_pident_by_length(length) and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                  self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
          if self.sort_method=="evl":
            if self.pident_method=="evalue":
              if evalue<self.option_evalue and qseqid in self.dict and self.dict[qseqid][3]>evalue and covs>=self.option_covs and covq>=self.option_covq:
                if self.verbose=='True':
                  print qseqid
                  print "%s %f is smaller than %s %f" % (sseqid,evalue,self.dict[qseqid][0],self.dict[qseqid][3])
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if evalue<self.option_evalue and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
            if self.pident_method=="static":
              if pident>=self.option_pident and length>=self.option_length and qseqid in self.dict and self.dict[qseqid][3]>evalue and covs>=self.option_covs and covq>=self.option_covq:
                if self.verbose=='True':
                  print qseqid
                  print "%s %f is smaller than %s %f" % (sseqid,evalue,self.dict[qseqid][0],self.dict[qseqid][3])
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if pident>=self.option_pident and length>=self.option_length and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
            if self.pident_method=="rost1999":
              if pident>=get_pident_by_length(length) and qseqid in self.dict and self.dict[qseqid][3]>evalue and covs>=self.option_covs and covq>=self.option_covq:
                if self.verbose=='True':
                  print qseqid
                  print "%s %f is smaller than %s %f" % (sseqid,evalue,self.dict[qseqid][0],self.dict[qseqid][3])
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if pident>=get_pident_by_length(length) and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,score,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
          if self.sort_method=="selfscore":
            normscore=score/self.selfscore[qseqid]
            if self.pident_method=="evalue":
              if evalue<self.option_evalue and qseqid in self.dict and self.dict[qseqid][6]<normscore and covs>=self.option_covs and covq>=self.option_covq:
                if self.verbose=='True':
                  print qseqid
                  print "%s %f is larger than %s %f" % (sseqid,normscore,self.dict[qseqid][0],self.dict[qseqid][6])
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,normscore,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if evalue<self.option_evalue and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,normscore,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
            if self.pident_method=="static":
              if pident>=self.option_pident and length>=self.option_length and qseqid in self.dict and self.dict[qseqid][6]<normscore and covs>=self.option_covs and covq>=self.option_covq:
                if self.verbose=='True':
                  print qseqid
                  print "%s %f is larger than %s %f" % (sseqid,normscore,self.dict[qseqid][0],self.dict[qseqid][6])
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,normscore,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if pident>=self.option_pident and length>=self.option_length and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,normscore,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
            if self.pident_method=="rost1999":
              if pident>=get_pident_by_length(length) and qseqid in self.dict and self.dict[qseqid][6]<normscore and covs>=self.option_covs and covq>=self.option_covq:
                if self.verbose=='True':
                  print qseqid
                  print "%s %f is larger than %s %f" % (sseqid,normscore,self.dict[qseqid][0],self.dict[qseqid][6])
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,normscore,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]
              if pident>=get_pident_by_length(length) and qseqid not in self.dict and covs>=self.option_covs and covq>=self.option_covq:
                self.dict[qseqid]=[sseqid,pident,length,evalue,bitscore,score,normscore,qlen,slen,nident,covs,covS,covns,covq,covQ,covnq]

class selfscore_dict(object):
  def __init__(self, infile, outfile, form):
    self.infile = infile
    self.outfile = outfile
    self.form = form
    self.dict = {}
        
  def parameters(self):
    print(self.infile, self.outfile)

  def parse(self):
    if self.form == 'blast+':
      with open(self.infile,'rU') as inhandle:
        for line in inhandle:
          parts = line.strip().split('\t')
          qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,nident,gaps,score = parts
          score = float(score)
          if qseqid not in self.dict and qseqid == sseqid:
            self.dict[qseqid] = float(score)
    if self.form == 'ghostx':
      with open(self.infile,'rU') as inhandle:
        for line in inhandle:
          parts = line.strip().split('\t')
          qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = parts
          score = float(bitscore)
          if qseqid not in self.dict and qseqid == sseqid:
            self.dict[qseqid] = float(bitscore)
#    if self.form == 'diamond':
#      with open(self.infile,'rU') as inhandle:
#        for line in inhandle:
#          parts = line.strip().split('\t')
#          qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = parts
#          score = float(bitscore)
#          if qseqid not in self.dict and qseqid == sseqid:
#            self.dict[qseqid] = float(bitscore)
    if self.form == 'score':
      with open(self.infile,'rU') as inhandle:
        inhandle.next()
        for lines in inhandle:
          line = lines.strip().split('\t')
          self.dict[line[0]]=float(line[1])

  def write(self):
    with open(self.outfile,'w') as outhandle:
      outhandle.write('query\tselfscore\n')
      for record in self.dict:
        outhandle.write('%s\t%s\n' % (record,str(self.dict[record])))
   
class blast(object):
  def __init__(self, path, prog, query, db, evalue, matrix, num_targets, num_threads, abs_out, prepath, precheck, out):
    self.path = path
    self.prog = prog
    self.query = query
    self.db = db
    self.evalue = str(evalue)
    self.matrix = matrix
    self.num_targets = str(num_targets)
    self.num_threads = str(num_threads)
    self.outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident gaps score"
    self.abs_out = abs_out
    self.out = out
    self.prepath = prepath
    self.precheck = precheck

  def parameters(self):
    print(self.path, self.prog, self.query, self.db, self.evalue, self.matrix, self.num_targets, self.num_threads, self.outfmt, self.abs_out, self.prepath, self.precheck, self.out)

  def process(self):
    print self.path+self.prog+' -query '+self.query+' -db '+self.db+' -evalue '+self.evalue+' -matrix '+self.matrix+' -max_target_seqs '+self.num_targets+' -num_threads '+self.num_threads+' -outfmt "'+self.outfmt+'" -out '+self.out
    os.system(self.path+self.prog+' -query '+self.query+' -db '+self.db+' -evalue '+self.evalue+' -matrix '+self.matrix+' -max_target_seqs '+self.num_targets+' -num_threads '+self.num_threads+' -outfmt "'+self.outfmt+'" -out '+self.abs_out)
        
  def check(self):
    if self.precheck:
      if os.path.isfile(self.prepath+'/'+self.out):
        return True

class ghostx(object):
  def __init__(self, path, prog, query, db, query_type, db_type, num_targets, num_threads, abs_out, prepath, precheck, out):
    self.path = path
    self.prog = prog
    self.query = query
    self.db = db
    self.query_type = query_type
    self.db_type = db_type
    self.num_targets = str(num_targets)
    self.num_threads = str(num_threads)
    self.abs_out = abs_out
    self.out = out
    self.prepath = prepath
    self.precheck = precheck

  def parameters(self):
    print(self.path, self.prog, self.query, self.db, self.query_type, self.db_type, self.num_targets, self.num_threads, self.abs_out, self.prepath, self.precheck, self.out)

  def process(self):
    print self.path+self.prog+' aln -i '+self.query+' -d '+self.db+' -b '+self.num_targets+' -a '+self.num_threads+' -o '+self.abs_out+' -q '+get_species_type_ghostx(self.query_type)+' -t '+get_species_type_ghostx(self.db_type)
    os.system(self.path+self.prog+' aln -i '+self.query+' -d '+self.db+' -b '+self.num_targets+' -a '+self.num_threads+' -o '+self.abs_out+' -q '+get_species_type_ghostx(self.query_type)+' -t '+get_species_type_ghostx(self.db_type))

  def check(self):
    if self.precheck:
      if os.path.isfile(self.prepath+'/'+self.out):
        return True

#class diamond(object):
#  def __init__(self, path, prog, query, db, query_type, db_type, num_targets, num_threads, abs_out, prepath, precheck, out):
#    self.path = path
#    self.prog = prog
#    self.query = query
#    self.db = db
#    self.query_type = query_type
#    self.db_type = db_type
#    self.num_targets = str(num_targets)
#    self.num_threads = str(num_threads)
#    self.abs_out = abs_out
#    self.out = out
#    self.prepath = prepath
#    self.precheck = precheck
#
#  def parameters(self):
#    print(self.path, self.prog, self.query, self.db, self.query_type, self.db_type, self.num_targets, self.num_threads, self.abs_out, self.prepath, self.precheck, self.out)
#
#  def process(self):
#    print self.path+self.prog+' aln -i '+self.query+' -d '+self.db+' -b '+self.num_targets+' -a '+self.num_threads+' -o '+self.abs_out+' -q '+get_species_type_ghostx(self.query_type)+' -t '+get_species_type_ghostx(self.db_type)
#    os.system(self.path+self.prog+' aln -i '+self.query+' -d '+self.db+' -b '+self.num_targets+' -a '+self.num_threads+' -o '+self.abs_out+' -q '+get_species_type_ghostx(self.query_type)+' -t '+get_species_type_ghostx(self.db_type))
#
#  def check(self):
#    if self.precheck:
#      if os.path.isfile(self.prepath+'/'+self.out):
#        return True

if __name__ == '__main__':
    main()
