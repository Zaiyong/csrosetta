#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os.path import basename
import os
import string
import argparse
import sys
#sys.path.append('/home/zak/Downloads/numpy-1.6.2')
import library
import StringIO
from assignment import noesy
#from noesy import CrossPeakInfo
#from noesy import CrossPeakList
#from noesy import Resonance
#from noesy import ResonanceList
from assignment.noesy import SpinSystem, get_strips,spinsystem_evaluation
import fasta
import copy

parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-resid",nargs='*', help='only show crosspeaks assigned to these residues',type=int, default=None);
parser.add_argument("-min_vc", default=None, type=float )
parser.add_argument("-min_overall_vc", default=None, type=float )
parser.add_argument("-min_seq", default=None, type=int )
library.add_standard_args( parser )
args = parser.parse_args()

if args.prot:
	resonance_list=noesy.ResonanceList.read_from_stream( open(args.prot,'r') )

crosspeaks=noesy.read_peak_files(args.peaks)

filtered=noesy.CrossPeakList()

for p in crosspeaks:
	use_peak=True
	if args.resid:
		use_peak=False
		for a in p:
		  match=True
			for r in args.resid:
			  if not a.has_resid(r):
					match=False
					break
			if match:
				use_peak=True
				break
	if use_peak:
		has_stats=False
		if args.min_vc:
			has_stats=True
			op=copy.copy(p)
			p.clear_assignments()
			best_decoy_comp=1000.0
			best_seq_sep=1000
			min_seq_sep=10000
			max_vc=0
			sym_max_vc=0
			cs_max_vc=0
			for a in op:
				vc = a.volume_contribution()
				if vc < args.min_vc:
					continue
				if max_vc < vc:
					max_vc = vc
					sym_max_vc = a.weight(1)
					cs_max_vc = a.weight(0)
				p.add_assignment(a)
				seq_sep=abs( a.atom( p.spin_info(1).atom_col ).resid() - a.atom( p.spin_info(2).atom_col ).resid() )
				if a.weight(6) < best_decoy_comp:
					best_decoy_comp=a.weight(6)
					best_seq_sep=seq_sep
					best_vc=a.volume_contribution()
					best_a=a
				if seq_sep < min_seq_sep:
					min_seq_sep=seq_sep
					min_vc=a.volume_contribution()
					min_a=a
	ambig=p.nassign()
	if use_peak and has_stats:
		if best_seq_sep < args.min_seq:
			use_peak=False
		if best_decoy_comp<-0.5:
			use_peak=False
	if use_peak and ambig>0:
		filtered.add_crosspeak(p)
		if has_stats: print "STATS: %5d %5.3f %5.2f   %5.2f %5.2f %5.2f   %5d %5.2f   %5d %5.2f  ----- %d %5s %s %s"%(ambig,best_decoy_comp,p.cum_peak_vol,max_vc,sym_max_vc, cs_max_vc, best_seq_sep,best_vc, min_seq_sep, min_vc, p.id(), p.info().experiment_id(), best_a.atom_str(), min_a.atom_str())
	if use_peak:
		min_seq_sep = 1000
		for a in p:
			seq_sep=abs( a.atom( p.spin_info(1).atom_col ).resid() - a.atom( p.spin_info(2).atom_col ).resid() )
			if seq_sep < min_seq_sep and a.volume_contribution() > 0.01:
				min_seq_sep=seq_sep
		try:
			print "CP_INFO: %5d %8.3f %8.3f %4d %20s %8.3f"%( min_seq_sep, p.cum_peak_vol,p.prob,p.class_id,p.class_str, p.min_viol)
		except AttributeError:
			pass


filtered.write_to_stream( sys.stdout )
