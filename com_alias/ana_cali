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
parser.add_argument("-min_vc", default=0.5, type=float )
parser.add_argument("-min_seq", default=None, type=int )
library.add_standard_args( parser )
args = parser.parse_args()

if args.prot:
	resonance_list=noesy.ResonanceList.read_from_stream( open(args.prot,'r') )

crosspeaks=noesy.read_peak_files(args.peaks)

filtered=noesy.CrossPeakList()
cp_map={}

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

	#now we have only min_vc assignments retained.
	if use_peak and p.nassign()>0:
		filtered.append( p )
		labelled_col = p.spin_info( 1 ).atom_col
		unlab_col = p.spin_info( 2 ).atom_col
		cp_map.setdefault(p.front().atom(labelled_col),{}).setdefault(p.front().atom(unlab_col),[]).append( p )

is_done={}
for p in filtered:
	labelled_col = p.spin_info( 1 ).atom_col
	unlab_col = p.spin_info( 2 ).atom_col
	indir=p.front().atom( labelled_col )
	dir=p.front().atom( unlab_col )
#	if p.front().atom( unlab_col ) in is_done: continue
#	if p.front().atom( labelled_col ) in is_done: continue
#	if dir.resid() < indir.resid(): continue
	res_dir = resonance_list.by_atom( dir )
	res_indir = resonance_list.by_atom( indir )

	try:
		sym_peaks = cp_map[ dir ][ indir  ]
		if len(sym_peaks)>1: continue
		sp = sym_peaks[ 0 ]
		if sp.info()!=p.info(): continue
		str=" %2s %2s  "%( resonance_list.aa_from_atom( dir ), resonance_list.aa_from_atom( indir ) )
		str+=" ".join([ "%s"%sp.front().atom( i ) for i in [1, 3, 2] ])
		print "%50s -- %8.3f %8.3f %8.3f %5d %5d"%(str,pow(1e-7*sp.volume(),-1.0/6),pow(1e-7*p.volume(),-1.0/6),pow(sp.volume()/p.volume(),-1.0/6),res_dir.intensity(), res_indir.intensity() )

#		is_done[p.front().atom( unlab_col ) ]=True
#		is_done[p.front().atom( labelled_col ) ]=True


	except KeyError:
		pass


