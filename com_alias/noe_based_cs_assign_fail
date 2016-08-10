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


parser = ApplicationArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-prot", help="chemical shift file",default=None);

library.add_standard_args( parser )
args=parser.parse_args()
library.init( args )
resonances = noesy.ResonanceList.read_from_stream( open(args.prot,'r') )
peaks = PeakCollection.from_peak_files( args.peaks, ignore=True )

def filter_resonances( res_list ):
	bb = noesy.ResonanceList()
	for resonance in res_list.itervalues():
		n=resonance.name()
		if n=='H' or \
					n=='HA' or\
					n=='C' or \
					n=='CA' or \
					n=='N':
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

def assign_randomly( peak_collection, molecule, freq_matcher, dist_matcher ):
	import random
	state=AssignmentCollection(peak_collection,molecule)
	for peak_list in peak_collection.experiments.itervalues():
		for peak in peak_list:
#			print peak
			if random.random() > 0.0: #assign 50% of peaks
				for match in peak.matches( molecule, None, frequency_matcher=freq_matcher, distance_matcher=dist_matcher ):
#					print match
					state.add(match)
	return state

def assign_partial( peak_collection, molecule, freq_matcher ):
	import random
	all_possible_assignments=AssignmentCollection(peak_collection,molecule)
	for peak in peak_collection.iterpeaks():
		for mask in peak.rule.spinsystem_match_masks():
			for match in peak.matches( molecule, None, frequency_matcher=freq_matcher, match_mask=mask ):
#				print match
				all_possible_assignments.add(match)
	return all_possible_assignments

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

class DiscoveredFreqMatcher:
	def __init__( self, state ):
		self.harvest_freq_data( state )

	def __call__(self, atom, freq, tol, folder):
		try:
			known=self.freqs[atom]
			best_match=min( (freq-w)**2/tol**2 for w in known )
#			print atom, best_match, freq, known
			return best_match < 2
		except KeyError:
			if atom.element == 'H': return False
			else: return True


	def harvest_freq_data( self, state ):
		import numpy
		import math
		self.freqs={}
		for atom, assignments in state.by_atom.iteritems():
			if assignments and len(assignments)>0:
				self.freqs[atom]=[ match.freq for match in assignments ]

def perturb_state( state, npeaks, freq_matcher, dist_matcher ):
	for peak in random_items( state.peak_collection.iterpeaks(), npeaks ):
		assignments=state.by_peak.get( peak )
		if assignments:
			for match in assignments:
#				print 'remove: ', match
				state.remove( match )
		for match in peak.matches( molecule, 1, frequency_matcher=freq_matcher, distance_matcher=dist_matcher ):
#			print 'add: ', match
			state.add(match)

def perturb_state_from_pool( pool, state, npeaks, freq_matcher, dist_matcher ):
#	print '\npertub...'
	ct=npeaks
	while ct>0:
		for partial_matches in random_items( pool.by_peak.itervalues(), npeaks ):
			peak=partial_matches[0].peak
			assignments=state.by_peak.get( peak )
			if assignments:
				for match in assignments:
#					print 'remove: ', match
					state.remove( match )
			for partial_match in random_items( partial_matches ):
				for match in peak.matches( molecule,
																	 1, #random_samples
																	 frequency_matcher=freq_matcher,
																	 distance_matcher=dist_matcher,
																	 partial_match=partial_match ):
#					print 'add: ', match
					if state.add( match ):
						ct-=1

def fill_state( pool, freq_matcher, dist_matcher ):
	state=AssignmentCollection( pool.peak_collection, pool.molecule )
	for partial_matches in pool.by_peak.itervalues():
		for partial_match in random_items( partial_matches ):
			for match in partial_match.peak.matches( pool.molecule,
																						 1,
																						 frequency_matcher=freq_matcher,
																						 distance_matcher=dist_matcher,
																						 partial_match=partial_match.peak_match ):
				state.add( match )
	return state

def assign_unassigned_peaks( state, freq_matcher, dist_matcher ):
	for peak in state.peak_collection.iterpeaks():
		if not peak in state.by_peak:
#			print peak
			for match in peak.matches( molecule, None, frequency_matcher=freq_matcher, distance_matcher=dist_matcher ):
	#			print match
				state.add(match)
	return state


def dump_assigned_freqs( fd, state, resonances, molecule ):
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
			min_dist = abs(min( w-ref for w in x ))/resonance.error()
			matched=min_dist<2
			fd.write('%15s %5s %8.3f %8.3f %5.3f %s\n'%(atom.long_str(), matched, ref, freq, min_dist, str(x) ) )
		else:
			fd.write('%15s %8s\n'%(atom.long_str(), 'nan') )


class MetaFreqMatcher(list):
#	def _init__(list):
#		self.default_answer=False

	def __call__(self, atom, freq, tol, folder ):
		for fm in self:
			try:
				return fm.evaluate(atom,freq,tol,folder)
			except KeyError:
				pass
		return self.default_answer

#### assignment methods
def assign_strips( strips, molecule, resonances ):
	import random
	known_fm=KnownFreqMatcher( resonances )
	assignments=AssignmentCollection(strips,molecule)
	for strip in strips:
		for match in strip.matches( molecule, frequency_matcher=known_fm ):
			print match
			assignments.add(match)
	assignments.commit()
	return assignments


def assign_strips_consistently( strips, molecule, resonances, scorefxn ):
	assignments=AssignmentCollection( strips, molecule )
	known_fm=KnownFreqMatcher( resonances )
	fm=SelfUpdateFreqMatcher( resonances, assignments )
	for strip in strips:
		has_match=False
		for match in strip.matches( molecule, 1, frequency_matcher=known_fm ):
			assignments.add( match )
			has_match=True
		if has_match: continue
		match = select_scored_match( strip.matches( molecule, 10, frequency_matcher=fm ), assignments, scorefxn )
		if match:
#			print match
			assignments.add( match )
	return assignments

def assign_peaks_of_strips( strip_assignments, peaks, molecule, resonances, scorefxn ):
	assignments=AssignmentCollection( peaks, molecule )
	fm1=SelfUpdateFreqMatcher( resonances, strip_assignments )
	fm2=SelfUpdateFreqMatcher( resonances, assignments )
	fm=MetaFreqMatcher([fm1,fm2])
	fm.default_answer=False
	known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
	known_dist.max_sequence_separation=2
	for strip in strip_assignments:
		print 'STRIP: ', strip
		for peak in strip.strip:
			has_match=False
			fm.default_answer=False
			for match in peak.matches( molecule, 1, frequency_matcher=fm, distance_matcher=known_dist, partial_match=strip._peak_match ):
				print 'FIXED: ', match
				assignments.add( match )
				has_match=True
			if has_match: continue
			fm.default_answer=True
			match = select_scored_match( peak.matches( molecule, 10, frequency_matcher=fm, distance_matcher=known_dist, partial_match=strip._peak_match ), assignments, scorefxn )
			if match:
				print 'LUCKY: ', match
				assignments.add( match )

	#		print match
#			assignments._add( match )
	return assignments


def initial_assign_strips( molecule, resonances, peaks ):
	from assignment.strips import collect_strips_from_peak_collection
	print 'cluster frequencies to obtain strips...'
	strips=collect_strips_from_peak_collection( peaks )

	print 'assign partially with known data...'
	known_fm.default_answer=False
	scorefxn=ScoreFunction(bmrb=1)
	matched_strips=assign_strips_consistently( strips, molecule, filter_resonances( resonances ), scorefxn )
	scorefxn.print_scores( matched_strips )

	rmsd,known=freq_rmsd(matched_strips,resonances)
	print '%8.3f %5d'%(rmsd, known)
	dump_assigned_freqs( open('prot_strip.assigned','w'), matched_strips, resonances, molecule, scorefxn )
	matched_strips.write_to_stream( open('state_strip.assigned','w') )

	print 'expand to full peak assignments...'
	scorefxn=ScoreFunction(bmrb=1,consistency=1,frag=1)
	state=assign_peaks_of_strips( matched_strips, peaks, molecule, filter_resonances( resonances ), scorefxn )
	return state




molecule=AtomTree.from_sequence( resonances.sequence() )
known_fm=KnownFreqMatcher( filter_resonances( resonances ) )
from assignment.scoring.methods import FragDistanceScore
from assignment.rules import ScoreDistanceMatcher
import math
from copy import copy
import random
import sys
known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
known_dist.max_sequence_separation=1
scorefxn=ScoreFunction(bmrb=1,frag=1,consistency=1)

print 'assign partially with known data...'
known_fm.default_answer=False
partial_matched_pool=assign_partial( peaks, molecule, known_fm )

print 'create initial state by randomly assigning the rest...'
known_fm.default_answer=True
state=fill_state( partial_matched_pool, known_fm, known_dist )

discovered_fm=DiscoveredFreqMatcher( state )
assign_unassigned_peaks( state, discovered_fm,  known_dist )

state.write_to_stream( open('state_init.assigned','w') )
dump_assigned_freqs( open('prot_init.assigned','w'), state, resonances, molecule )


last_accepted_score=scorefxn(state)
last_accepted_state=copy(state)

temp=1
print 'score        rmsd'
Nacc_window=50
from collections import deque
acceptance_stats=deque([1.0]*Nacc_window)
acceptance_rate=1
for i in xrange( 0, 10000):
	known_fm.default_answer=True
	perturb_state_from_pool( partial_matched_pool, state, 1, known_fm, known_dist )
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
	rmsd,known=freq_rmsd(state,resonances)
	print '%5d %8.3f %8.3f %8.3f %5d %4.0f%%'%(i, deltaV, scorefxn(state), rmsd, known , acceptance_rate*100 )
	dump_freq=100
	if i%dump_freq==0:
		state.write_to_stream( open('state_%03d.assigned'%(i/dump_freq),'w') )
		dump_assigned_freqs( open('prot_%03d.assigned'%(i/dump_freq),'w'), state, resonances, molecule )
state.write_to_stream( open('state_final.assigned','w') )
dump_assigned_freqs( open('prot_final.assigned','w'), state, resonances, molecule )

