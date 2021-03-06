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
from assignment.scoring.methods import FragDistanceScore
import library
from basic.options import ApplicationArgumentParser
from basic import Tracer
import math


parser = ApplicationArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-fasta", help="fasta file",default=None);
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


class SelfUpdateFreqMatcher:
	def __init__(self, resonances, state):
		self.resonances=resonances
		self.default_answer=False
		self.state=state

	def __call__(self, atom, freq, tol, folder):
		try:
			return self.evaluate(atom,freq,tol,folder)
		except KeyError:
			return self.default_answer

	def evaluate(self, atom, freq, tol, folder):
		reso=self.resonances.by_chemical_atom( atom )
		return reso.match( freq, tol, folder )

def initial_assign(peak,molecule,fm,known_dist):
	for match in random_items( peak.matches( molecule, frequency_matcher=fm, distance_matcher=known_dist ), 1 ):
		if match:
			return match
	return None

resonances = noesy.ResonanceList.read_from_stream( open(args.prot,'r') )
#resonances = filter_resonances( ref_resonances )
peaks = PeakCollection.from_peak_files( args.peaks, ignore=True )
molecule=AtomTree.from_sequence( resonances.sequence() )
state=AssignmentCollection( peaks, molecule )
scorefxn=ScoreFunction(bmrb=1,consistency=1,symmetry=1)

import random
peak_order = [ peak for peak in peaks ]
#random.shuffle( peak_order )
fm=SelfUpdateFreqMatcher( resonances, state )
#known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
#known_dist.max_sequence_separation=2
for peak in peak_order:
	tr.Debug('main:match %s'%peak)
	match=initial_assign(peak,molecule,fm,None)
	if match:
		tr.Debug('main:match %s'%match)
		state.add( match )
scorefxn.print_scores( state )
state.write_to_stream( open(args.output+'_final.state','w') )
