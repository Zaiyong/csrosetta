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
parser.add_argument("-states", nargs='*', help='recover state from this file', default=None )
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

class ResonanceCollector(object):
	def __init__( self, resonances, molecule ):
		self.reference_resonances={}
		self.molecule=molecule
		self.resonances={}
		self._init_atoms(resonances)

	def _init_atoms(self, resonances):
		for resonance in resonances.itervalues():
			n=resonance.name()
			if n=='H' or \
					n=='HA' or\
					n=='C' or \
					n=='CA' or \
					n=='N' or \
					n=='HA2' or \
					n=='HA3':
				continue
			try:
				atom = molecule.atom( Atom( resonance.name(), resonance.resid() ) )
			except KeyError:
				try:
					atom = molecule.atom( Atom( resonance.name(), resonance.resid() ).methyl_atom() )
				except KeyError:
					continue
			self.resonances[atom]=[]
			self.reference_resonances[atom]=resonance

	def append( self, state ):
		import numpy
		import math
		for atom, freqs in self.resonances.iteritems():
			assignments = state.by_atom.get( atom )
			if assignments:
				x = [ match.freq for match in assignments ]
				freq = sum( x )/len(assignments)
				freqs.append( freq )

	def dump_freqs( self, fd ):
		import numpy
		import math
		for atom, freqs in self.resonances.iteritems():
			if len(freqs)==0:
				fd.write('%15s %8s\n'%(atom.long_str(), 'nan') )
				continue
			ref_reso=self.reference_resonances[atom] #should use resonances methods to account for folding. no folding in GmR137
			ref=ref_reso.freq()
			mean_freq = sum( freqs )/len( freqs )
			if len(freqs)>1:
				std_freq = math.sqrt(sum( [ (x-mean_freq)*(x-mean_freq) for x in freqs ] )/(len(freqs)-1))
			else:
				std_freq = 0
			min_dist = min( abs(w-ref) for w in freqs )/ref_reso.error()
			avg_dist = abs( mean_freq - ref )/ref_reso.error()
			good_count = sum([ 1 for x in freqs if (abs( x - ref )/ref_reso.error()) < 1.5 ])
			count = len(freqs)
			matched=min_dist<2
			fd.write('%15s %5s %8.3f %8.3f %8.3f %8.3f %8.3f %5d %5d\n'%(atom.long_str(), matched, ref, mean_freq, std_freq, min_dist, avg_dist, good_count, count) )


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

states=[]
freqs=ResonanceCollector( ref_resonances, molecule )

for saved_state in args.states:
	state = AssignmentCollection.from_saved_state( peaks, molecule, open( saved_state, 'r' ) )
	states.append( state )
	freqs.append( state )

freqs.dump_freqs( open( args.output, 'w' ) )

