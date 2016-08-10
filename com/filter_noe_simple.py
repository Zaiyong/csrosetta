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
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-min_vc", default=None, type=float )
parser.add_argument("-cheats", default=None, help='file with peak-id, assigned residue pair' );

library.add_standard_args( parser )
args = parser.parse_args()

crosspeaks=noesy.read_peak_files(args.peaks)
filtered=noesy.CrossPeakList()

for p in crosspeaks:
	if args.min_vc:
		op=copy.copy(p)
		p.clear_assignments()
		for a in op:
			vc = a.volume_contribution()
			if vc < args.min_vc:
				continue
			p.add_assignment(a)
	filtered.add_crosspeak(p)

def to_string( pair ):
	s="%3d - "%pair[0]
	if len(pair)>1:
		s+='%-3d'%pair[1]
	else:
		s+='%-3s'%' '
	return s

if args.cheats:
	cheat_by_peak={}
	cheat_by_pair={}
	for line in open(args.cheats,'r'):
		tags=line.split()
		if not len(tags)==3: continue

		pair=tuple(sorted(list(set([int(tags[1]),int(tags[2])]))))
		if 0==pair[0]:
			continue
		id = int(tags[0])
		cheat_by_peak[id]=pair
		cheat_by_pair[pair]=id


	found_pairs={} #map pair to VC max
	mismatched_cheat_pairs=set()
	mismatched_found_pairs={}
	for p in filtered:
		found = None
		try:
			for a in p:
				my_pair=a.resid_pair()
				vc=a.volume_contribution()
				maxvc=max(found_pairs.setdefault(my_pair,vc),vc)
				found_pairs[my_pair]=maxvc
				if my_pair == cheat_by_peak[p.id()]:
					found = my_pair
					break
			if not found:
				s='mismatch for peak-id %5d'%p.id()+' cheat: '+to_string(cheat_by_peak[p.id()])+'  found: '
				indent=len(s)
				ct=1
				for a in p:
					if ct>=2:
						s+='\n'+' '*indent
					s+=to_string( a.resid_pair() )
					ct+=1
				print s

				mismatched_cheat_pairs.add(cheat_by_peak[p.id()])
				vc=a.volume_contribution()
				maxvc=max(mismatched_found_pairs.setdefault(my_pair,vc),vc)
				mismatched_found_pairs[my_pair]=maxvc

		except KeyError:
			print 'peak %5d is not assigned in cheats'%p.id()
			for a in p:
				my_pair=a.resid_pair()
				vc=a.volume_contribution()
				maxvc=max(mismatched_found_pairs.setdefault(my_pair,vc),vc)
				mismatched_found_pairs[my_pair]=maxvc

	for p in sorted(mismatched_found_pairs.keys()):
		if not p in cheat_by_pair:
			vc = mismatched_found_pairs[p]
			print 'not in cheat: vc %5.3f %s'%(vc,to_string(p))
		else:
			print 'mismatched but later found in cheat: %s'%to_string(p)

	for p in mismatched_cheat_pairs:
		if not p in found_pairs:
			print 'only in cheat: %s'%to_string(p)
		else:
			print 'mismatched but found somewhere else in autoNOE: %s'%to_string(p)

filtered.write_split_files()

