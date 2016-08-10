#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import string
import argparse
import sys
import traceback

from os import path
from assignment.noesy import ResonanceList
from fasta import read_fasta
import library
from amino_acids import NO_resonances, NO_bb_resonances
from PDB.Polypeptide import one_to_three

parser = argparse.ArgumentParser(prog=path.basename(__file__), description="to check the completeness of CS of each residue")
parser.add_argument("-prot", help="chemical shift file");
parser.add_argument("-select", help="all, bb(backbone) or sc(sidechain)",default='all');
library.add_standard_args( parser )

args = parser.parse_args()
resonance_list=ResonanceList.read_from_stream( open(args.prot,'r') )

seq=resonance_list.sequence()
for i,aa in enumerate(seq):
	resid=i+1
	count=0
	if args.select=='all':
		try:
			for resonance in resonance_list.by_residue(resid):
				if resonance.atom().elem() not in 'HCN': continue
				elif resonance.name()=='C': continue
				else: count+=1
			print "all:   %5d %5s   %5.3f"%(resid,aa,float(count)/NO_resonances[one_to_three(aa)])
		except KeyError:
			print "all:   %5d %5s   %5.3f"%(resid,aa,0.0)
	elif args.select=='bb':
		try:
			for resonance in resonance_list.by_residue(resid):
				if resonance.atom().elem() not in 'HCN': continue
				elif resonance.name() in ['N','H','HA','HA3','HA2','HA','CA','CB']: count+=1
			print "backbone:  %5d %5s   %5.3f"%(resid,aa,float(count)/NO_bb_resonances[one_to_three(aa)])
		except KeyError:
			print "backbone:  %5d %5s   %5.3f"%(resid,aa,0.0)

	elif args.select=='sc':
		if aa=='G':
			print "sidechain:  %5d %5s   %5.3f"%(resid,aa,1.0)
		else:
			try:
				for resonance in resonance_list.by_residue(resid):
					if resonance.atom().elem() not in 'HCN': continue
					elif resonance.name() not in ['C','N','H','HA','HA3','HA2','HA','CA','CB']: count+=1
				print "sidechain:  %5d %5s   %5.3f"%(resid,aa,float(count)/(NO_resonances[one_to_three(aa)]-NO_bb_resonances[one_to_three(aa)]))
			except KeyError:
				print "sidechain:  %5d %5s   %5.3f"%(resid,aa,0.0)
