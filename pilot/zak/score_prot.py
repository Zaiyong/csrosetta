#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#standard modules
import argparse

#csrosetta modules
from assignment import noesy
from assignment import PeakCollection
from assignment.scoring import ScoreFunction
from assignment import AssignmentCollection
from assignment import random_items
from chemical import AtomTree, Atom
from assignment.rules import ScoreDistanceMatcher
from assignment.scoring.methods import FragDistanceScore, ConformationDistanceScore
import library
from basic.options import ApplicationArgumentParser
from basic import Tracer
import math
import time
import BmrbAtomNames
import fasta

start=time.clock()
parser = ApplicationArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-prot", help="chemical shift file",default=None);
#parser.add_argument("-ref_prot", help="ref chemical shift file",default=None);
parser.add_argument("-fasta", help="fasta file",default=None);
parser.add_argument("-output", help="how to call the output files", default="default" )
library.add_standard_args( parser )
args=parser.parse_args()
library.init( args )

tr=Tracer('main')



# def assign_peak_consistently( peak, state, resonances, scorefxn ):
# 	#known_fm=KnownFreqMatcher( resonances )
# 	fm=SelfUpdateFreqMatcher( resonances, state )
# 	known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
# 	known_dist.max_sequence_separation=2
# #	print '%20s'%'assign peak: ', peak
# 	fm.default_answer=False
# 	for match in peak.matches( molecule, random_samples=1, frequency_matcher=fm, distance_matcher=known_dist ):
# #		print '%20s'%'FIX:   ', match
# 		return match
# 	fm.default_answer=True
# 	match = select_scored_match( peak.matches( molecule, frequency_matcher=fm, distance_matcher=known_dist ), state, scorefxn )
# #	if match:
# #		print '%20s'%'LUCKY: ', match
# 	return match

# def dump_assigned_freqs( fd, state, resonances, molecule, scorefxn ):
# 	import numpy
# 	import math
# 	for resonance in resonances.itervalues():
# 		try:
# 			atom = molecule.atom( Atom( resonance.name(), resonance.resid() ) )
# 		except KeyError:
# 			try:
# 				atom = molecule.atom( Atom( resonance.name(), resonance.resid() ).methyl_atom() )
# 			except KeyError:
# 				continue
# 		assignments = state.by_atom.get( atom )
# 		if assignments:
# 			x = [ match.freq for match in assignments ]
# 			freq = sum( x )/len(assignments)
# 			ref=resonance.freq() #should use resonances methods to account for folding. no folding in GmR137
# 			if resonance.atom().elem()=='H':
# 				error=0.03
# 			else:
# 				error=0.3
# 			min_dist = min( abs(w-ref) for w in x )/error
# 			matched=min_dist<2
# 			ids = [ match.peak.id for match in assignments ]
# 			score_str = scorefxn.local_score_str( state, atom )
# 			fd.write('%15s %5s %8.3f %8.3f %8.3f %20s\n'%(atom.long_str(), matched, ref, freq, min_dist, score_str, ) )
# #str(ids) ) )
# 		else:
# 			fd.write('%15s %8s\n'%(atom.long_str(), 'nan') )

# def assign_more_peaks( state, molecule, resonances, scorefxn ):
# 	from assignment import random_items
# 	trm=Tracer('mover',tr)
# 	ct=0
# 	total=0
# 	for peak in random_items( state.unassigned_peaks, 50 ):
# 		total+=1
# 		match = assign_peak_consistently( peak, state, resonances, scorefxn )
# #		trm.Debug('peak: ',peak )
# #		trm.Debug('match :',match)
# 		if match:
# 			ct+=1
# 			state.add( match )

# 	trm.Info('more_peaks: matched %5d (%5d) additional peaks...'%(ct,total))


# def filter_resonances( res_list ):
# 	bb = noesy.ResonanceList()
# 	bb.set_sequence( res_list.sequence() )
# 	for resonance in res_list.itervalues():
# 		n=resonance.name()
# 		if n=='H' or \
# 					n=='HA' or\
# 					n=='C' or \
# 					n=='CA' or \
# 					n=='N' or \
# 					n=='HA2' or \
# 					n=='HA3':
# 			bb.add_resonance( resonance )
# 	return bb


# class FreqMatcher:
# 	def __init__(self, resonances):
# 		self.resonances=resonances
# 		self.default_answer=False

# 	def __call__(self, atom, freq, tol, folder):
# 		try:
# 			return self.evaluate(atom,freq,tol,folder)
# 		except KeyError:
# 			return self.default_answer

# 	def evaluate(self, atom, freq, tol, folder):
# 		reso=self.resonances.by_chemical_atom( atom )
# 		return reso.match( freq, tol, folder )

def methyl_atom_pool_detect(aa,proton_pool):
	if aa=='A' and proton_pool=='QB':
		return True
	elif aa=='M' and proton_pool=='QE':
		return True
	elif aa=='T' and proton_pool=='QG2':
		return True
	elif aa=='V' and proton_pool in ['QG1','QG2']:
		return True
	elif aa=='L' and proton_pool in ['QD1','QD2']:
		return True
	elif aa=='I' and proton_pool in ['QD1','QG2']:
		return True
	return False

def unpack_unmethyl_atom_pool(resonances):
	seq=resonances.sequence()
	res_to_add=[]
	for res in resonances.iter_by_atom_elem('H'):
		aa=seq[res.resid()-1]
		if res.atom().name[0]=='Q':
			if methyl_atom_pool_detect(aa,res.atom().name):
				continue
			else:
				nmr_names=BmrbAtomNames.get_nmr_names(aa,res.atom().name)
				if nmr_names:
					res.atom().set_name(nmr_names[0])
					for name in nmr_names[1:]:
						res_to_add.append(noesy.Resonance(atom=Atom(name,res.resid()),freq=res.freq(),error=res.error()))
	for res in res_to_add:
		#print res.atom().resid()
		resonances.add_resonance(res)

def initial_assign(peak,molecule,fm,known_dist):
	for match in random_items( peak.matches( molecule, frequency_matcher=fm, distance_matcher=known_dist ), 1 ):
		if match:
			return match
	return None

#ref_resonances = noesy.ResonanceList.read_from_stream( open(args.ref_prot,'r') )
resonances = noesy.ResonanceList.read_from_stream( open(args.prot,'r') )
unpack_unmethyl_atom_pool(resonances)
#resonances = filter_resonances( ref_resonances )
sequence=fasta.read_fasta(args.fasta)
resonances.set_sequence(sequence)
peaks = PeakCollection.from_peak_files( args.peaks, ignore=True )
molecule=AtomTree.from_sequence( resonances.sequence() )
state=AssignmentCollection( peaks, molecule )
scorefxn=ScoreFunction(bmrb=1,consistency=1,symmetry=1)

import random
from assignment import ConstantFreqMatcher
peak_order = [ peak for peak in peaks ]
#random.shuffle( peak_order )
fm=ConstantFreqMatcher( resonances )
known_dist=ScoreDistanceMatcher( ConformationDistanceScore(), abs(math.log(0.3)), 0 )
#known_dist.max_sequence_separation=9
count=0
for peak in peak_order:
	count+=1
#	tr.Debug('main:match %s'%peak)
	match=initial_assign(peak,molecule,fm,None)
	if match:
#		tr.Debug('main:match %s'%match)
		state.add( match )
	if count%100==0:
		tr.Info('%d peaks have been assigned'%count)
#print 'the program takes %d seconds'%(time.clock()-start)
scorefxn.print_scores( state )
#state.write_to_stream( open(args.output+'_.state','w') )
#dump_assigned_freqs( open(args.output+'_.prot','w'), state, ref_resonances, molecule, scorefxn )

