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

import library
from basic.options import ApplicationArgumentParser
from basic import Tracer

parser = ApplicationArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-state", help='recover state from this file', default=None )
parser.add_argument("-sample", action='store_true', default=False )
parser.add_argument("-print_residues", action='store_true', default=False )
parser.add_argument("-output", help="how to call the output files", default="default" )
library.add_standard_args( parser )
args=parser.parse_args()
library.init( args )

tr=Tracer('main')

def filter_resonances( res_list ):
	bb = noesy.ResonanceList()
	bb.set_sequence( res_list.sequence() )
	for resonance in res_list.itervalues():
		n=resonance.name()
		if n=='H' or \
					n=='HA' or\
					n=='C' or \
					n=='CA' or \
					n=='N' or \
					n=='HA2' or \
					n=='HA3':
			bb.add_resonance( resonance )
	return bb

class KnownFreqMatcher(object):
	def __init__(self,resonances):
		self.resonances=resonances
		self.default_answer=False

	def __call__(self, atom, freq, tol, folder):
		try:
			reso=self.resonances.by_chemical_atom( atom )
			return reso.match( freq, tol, folder )
		except KeyError:
			return self.default_answer # if nothing is known we allow it

class SelfUpdateFreqMatcher:
	def __init__(self, resonances, state):
		self.resonances=resonances
		self.default_answer=True
		self.state=state

	def __call__(self, atom, freq, tol, folder):
		try:
			return self.evaluate(atom,freq,tol,folder)
		except KeyError:
			return self.default_answer

	def evaluate(self, atom, freq, tol, folder):
		try:
			reso=self.resonances.by_chemical_atom( atom )
			return reso.match( freq, tol, folder )
		except KeyError:
			assignments = self.state.by_atom[ atom ]
			import numpy as np
			import math
			if assignments and len(assignments)>0:
				avg_freq=np.mean( [match.freq for match in assignments] )
				return abs(avg_freq-freq)<=tol/2
			else:
				raise

def freq_rmsd( state, resonances ):
	import numpy
	import math
	rmsd=0
	known=0
	for resonance in resonances.itervalues():
#		atom = Atom( resonance.name(), resonance.resid() )
		try:
			atom = molecule.atom( Atom( resonance.name(), resonance.resid() ) )
		except KeyError:
			try:
				atom = molecule.atom( Atom( resonance.name(), resonance.resid() ).methyl_atom() )
			except KeyError:
				continue
		assignments = state.by_atom.get( atom )
		if assignments:
			x = [ match.freq for match in assignments ]
			freq = sum( x )/len(assignments)
			delta=freq-resonance.freq() #should use resonances methods to account for folding. no folding in GmR137
			known+=min( w-resonance.freq() for w in x )<resonance.error()
#			if atom.name=='CA':
#				print '%10s %8.3f'%(atom, abs(delta) )
			rmsd+=min(1,delta*delta)
		else:
#			print '%10s %8s'%(atom, ' nan')
			rmsd+=1
	return math.sqrt(rmsd/len(resonances)),known

def dump_assigned_freqs( fd, state, resonances, molecule, scorefxn ):
	import numpy
	import math
	for resonance in resonances.itervalues():
		try:
			atom = molecule.atom( Atom( resonance.name(), resonance.resid() ) )
		except KeyError:
			try:
				atom = molecule.atom( Atom( resonance.name(), resonance.resid() ).methyl_atom() )
			except KeyError:
				continue
		assignments = state.by_atom.get( atom )
		if assignments:
			x = [ match.freq for match in assignments ]
			freq = sum( x )/len(assignments)
			ref=resonance.freq() #should use resonances methods to account for folding. no folding in GmR137
			min_dist = min( abs(w-ref) for w in x )/resonance.error()
			matched=min_dist<2
			ids = [ match.peak.id for match in assignments ]
			score_str = scorefxn.local_score_str( state, atom )
			fd.write('%15s %5s %8.3f %8.3f %8.3f %20s %20s\n'%(atom.long_str(), matched, ref, freq, min_dist, score_str,
																															'['+','.join(['%8.3f'%w for w in x])+']' ) )
#str(ids) ) )
		else:
			fd.write('%15s %8s\n'%(atom.long_str(), 'nan') )

#### assignment methods
def assign_partial( peak_collection, molecule, freq_matcher ):
	import random
	all_possible_assignments=AssignmentCollection(peak_collection,molecule)
	for peak in peak_collection.iterpeaks():
		for mask in peak.rule.spinsystem_match_masks():
			for match in peak.matches( molecule, None, frequency_matcher=freq_matcher, match_mask=mask ):
				print match
				all_possible_assignments.add(match)
	return all_possible_assignments


def select_scored_match( match_generator, state, scorefxn ):
	import math
	import random
	import bisect
	from assignment.scoring.methods import BmrbShiftScore
	totals=[]
	scored_matches=[]
	old_score=scorefxn( state )
	T=1
	running_total=0
	ct=0
#	bmrb_score=BmrbShiftScore()
	for match in match_generator:
		if match in scored_matches: continue
#		score=bmrb_score.score_assignment( match )
		state.add( match )
		w=math.exp( (scorefxn(state)-old_score)/T )
		#w=math.exp( score/T )
		state.remove( match )
		if w<0.8: continue
		running_total+=w
		ct+=1
#		print '%5d %5.3f %5.3f %8.3f %s'%(ct, w, running_total, scorefxn(state), match)
		scored_matches.append( match )
		totals.append( running_total )

	if len(scored_matches)==0: return None
#	scored_matches=sorted( scored_matches, key=lambda x: x[0] ) #sorting reduces run-time in the following list
	p=random.random()*running_total
	select_i = bisect.bisect_right( totals, p )
#	print 'selected: %d'%(select_i+1)
	return scored_matches[ select_i ]

#	for w, match in scored_matches:
#		ct+=1
#		total+=w/Z
#		print '%5d %5.3f %5.3f %5.3f %5s'%(ct, w/Z, total, p, match)
#		if p<total: return match
#	assert False

def assign_peak_consistently( peak, state, resonances, scorefxn ):
	#known_fm=KnownFreqMatcher( resonances )
	fm=SelfUpdateFreqMatcher( resonances, state )
	known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
	known_dist.max_sequence_separation=2
#	print '%20s'%'assign peak: ', peak
	fm.default_answer=False
	tr.Trace('main:assign_peak_consistently')
	for match in peak.matches( molecule, random_samples=1, frequency_matcher=fm, distance_matcher=known_dist ):
#		print '%20s'%'FIX:   ', match
		return match
	fm.default_answer=True
	tr.Trace('main:assign_peak_consistently_2nd')
	match = select_scored_match( peak.matches( molecule, frequency_matcher=fm, distance_matcher=known_dist ), state, scorefxn )
#	if match:
#		print '%20s'%'LUCKY: ', match
	return match


def reassign_atoms( state, molecule, resonances, scorefxn ):
	fm=SelfUpdateFreqMatcher( resonances, state )
	known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
	known_dist.max_sequence_separation=2
	import copy

	removed_peaks=set()
	matches_to_remove=set()

	for ct in range(0,10):
		heavy_atom=random.choice( [atom for atom in state.by_atom.iterkeys() if atom.element!='H' and not resonances.has_chemical_atom( atom )] )
		matches_to_remove.update( match.peak_match for match in state.by_atom[ heavy_atom ] )
		for proton in molecule.proton_partners( heavy_atom ):
			matches_to_remove.update( match.peak_match for match in state.by_atom[ proton ] )

	for match in sorted( matches_to_remove ):
#		print 'match to remove: ', match
		state.remove( match )
		removed_peaks.add( match.peak )

	ct_new=0
	ct_tot=0
	ct_remove=0
	for peak in sorted( removed_peaks ):
		match = assign_peak_consistently( peak, state, resonances, scorefxn )
		ct_tot+=1
		if not match:
#			print 'nothing'
			ct_remove+=1
			continue
		state.add( match )
		if match not in matches_to_remove:
			ct_new+=1
#			print 'new match: ', match
#		else:
#			print 'old match: ', match
	trm=Tracer('mover',tr)
	trm.Info('reassign_atoms: reconsidered %5d, %5.2f%% new %5.2f%% removed'%(ct_tot, float(ct_new)/ct_tot*100, float(ct_remove)/ct_tot*100 ))
#	scorefxn.print_scores( state )

def assign_more_peaks( state, molecule, resonances, scorefxn ):
	from assignment import random_items
	trm=Tracer('mover',tr)
	ct=0
	total=0
	for peak in random_items( state.unassigned_peaks, 50 ):
		total+=1
		match = assign_peak_consistently( peak, state, resonances, scorefxn )
#		trm.Debug('peak: ',peak )
#		trm.Debug('match :',match)
		if match:
			ct+=1
			state.add( match )

	trm.Info('more_peaks: matched %5d (%5d) additional peaks...'%(ct,total))

def initial_assign( molecule, resonances, peaks ):
	scorefxn=ScoreFunction(bmrb=1)
	state=AssignmentCollection( peaks, molecule )
	import random
	peak_order = [ peak for peak in peaks ]
	random.shuffle( peak_order )
	for peak in peak_order:
		match = assign_peak_consistently( peak, state, resonances, scorefxn )
		if match:
			state.add( match )
	return state

def extract_residue_matches( state, molecule, resid ):
	matches=set()
	for atom in molecule.iter_residue( resid ):
		if atom.element!='H':
			for match in state.by_atom[ atom ]:
				matches.add( match.peak_match )

	return sorted(matches, key=lambda match: match.peak_match)



#### main program
## setup
from assignment.scoring.methods import FragDistanceScore
from assignment.rules import ScoreDistanceMatcher
from assignment.strips import collect_strips_from_peak_collection
import math
from copy import copy
import random
import sys
import time

ref_resonances = noesy.ResonanceList.read_from_stream( open(args.prot,'r') )
resonances = filter_resonances( ref_resonances )
peaks = PeakCollection.from_peak_files( args.peaks, ignore=True )
molecule=AtomTree.from_sequence( resonances.sequence() )

#known_fm=KnownFreqMatcher( filter_resonances( resonances ) )
#known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
#known_dist.max_sequence_separation=2
#dump_assigned_freqs( open('prot_test.assigned','w'), state, resonances, molecule, scorefxn )
#state.write_to_stream( open('state_test.assigned','w') )

if args.state:
	state = AssignmentCollection.from_saved_state( peaks, molecule, open( args.state, 'r' ) )
#	scorefxn=ScoreFunction(bmrb=1,consistency=1,frag=1)
#	scorefxn( state )
#	state.write_to_stream( open('state_read.assigned','w') )
#	dump_assigned_freqs( open('prot_read.assigned','w'), state, ref_resonances, molecule, scorefxn )

else:
	start = time.clock()
	state = initial_assign( molecule, resonances, peaks )
	elapsed = (time.clock() - start)
	print 'assignment took %5.1f seconds.'%elapsed

	print 'dump...'
	scorefxn=ScoreFunction(bmrb=1,consistency=1,frag=1,symmetry=1)
	scorefxn( state )
	dump_assigned_freqs( open(args.output+'_init.prot','w'), state, ref_resonances, molecule, scorefxn )
	state.write_to_stream( open(args.output+'_init.state','w') )


if args.print_residues:
	sequence=molecule.sequence
	for resid in range(molecule.first_residue,molecule.last_residue+1):
		matches=extract_residue_matches( state, molecule, resid )
		print 'Residue: %5d %1s'%(resid, sequence[resid-1])
		for match in matches:
			print '%-80s'%match
		print

scorefxn=ScoreFunction(bmrb=1,consistency=1,frag=1,symmetry=1)
scorefxn.print_scores( state )
rmsd,known=freq_rmsd(state,ref_resonances)
print '%8.3f %5d'%(rmsd, known)

#perturb_state( state, molecule, resonances, scorefxn )

#rmsd,known=freq_rmsd(state,ref_resonances)
#print '%8.3f %5d'%(rmsd, known)
#dump_assigned_freqs( open('prot_perturb.assigned','w'), state, ref_resonances, molecule, scorefxn )
#state.write_to_stream( open('state_perturb.assigned','w') )

if not args.sample: exit(0)

last_accepted_score=scorefxn(state)
last_accepted_state=copy(state)

temp=1
print 'score        rmsd'
Nacc_window=50
from collections import deque
acceptance_stats=deque([1.0]*Nacc_window)
acceptance_rate=1
trscore=Tracer('score',tr)
for i in xrange( 0, 10000):
	reassign_atoms( state, molecule, resonances, scorefxn )
	assign_more_peaks( state, molecule, resonances, scorefxn )

	score=scorefxn(state)
	deltaV=last_accepted_score-score
	if random.random() < min( 1, math.exp( - deltaV / temp ) ):
		last_accepted_score=score
		last_accepted_state=state
		acceptance_stats.append(1.0)
		acceptance_rate+=1.0/Nacc_window
	else:
		acceptance_stats.append(0.0)
	state=copy(last_accepted_state)
	acceptance_rate-=acceptance_stats.popleft()/Nacc_window
	rmsd,known=freq_rmsd(state,ref_resonances)
	trscore.Info('ROUND:  %5d %8.3f %8.3f %8.3f %5d %4.0f%%'%(i, deltaV, scorefxn(state), rmsd, known , acceptance_rate*100 ))
	trscore.Info(scorefxn.score_str( state ))
	dump_freq=10
	if i%dump_freq==0:
		state.write_to_stream( open(args.output+'_%03d.state'%(i/dump_freq),'w') )
		dump_assigned_freqs( open(args.output+'_%03d.prot'%(i/dump_freq),'w'), state, ref_resonances, molecule, scorefxn )
state.write_to_stream( open(args.output+'_final.state','w') )
dump_assigned_freqs( open(args.output+'_final.prot','w'), state, ref_resonances, molecule, scorefxn )


# score...
# SCORE:      total      frag       bmrb consistency
# WEIGHT:               1.000      1.000      1.000
# SCORE:   3196.866   914.327   1780.689    501.850
# dump...
#   0.546   458


# score...
# SCORE:      total      frag       bmrb consistency
# WEIGHT:               1.000      1.000      1.000
# SCORE:   3078.939   887.961   1689.736    501.243
# dump...
#    0.547   469
# assignment took 1333.5 seconds.

# via strips
# score...
# SCORE:      total      frag       bmrb consistency
# WEIGHT:               1.000      1.000      1.000
# SCORE:   1968.223   566.496   1001.696    400.031
# dump...
#    0.636   381
# assignment took 1188.6 seconds.
