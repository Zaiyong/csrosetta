#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments
import assignment.scoring

import string
import argparse
import sys
import library
import fasta
import traceback
from os import path
from chemical import AtomTree
from assignment import PeakList, PeakCollection
import assignment
#from assignment.score import methods
from basic.options import ApplicationArgumentParser
from assignment import scoring
from assignment import AssignmentCollection

from assignment.scoring import methods
from basic import Tracer
#############################
#if len(argv) <=1:
#    Help()

#file = argv[1]
parser =ApplicationArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*',help="files with peak-lists in xeasy format")
parser.add_argument("-fasta", help='sequence information', default=None );
parser.add_argument("-s", help='structural information', default=None );
#parser.add_argument("-ignore", help='ignore assignments in peak file', action='store_true', dest='ignore', default=True )
parser.add_argument("-noignore", help='ignore assignments in peak file', action='store_false', dest='ignore', default=True )

library.add_standard_args( parser )
args = parser.parse_args()
library.init( args )

def assign_randomly(peak_collection,molecule):
	state=AssignmentCollection(peak_collection,molecule)
	for peak_list in peak_collection:
		for peak in peak_list:
			for match in peak.matches( molecule, 1 ):
				state.add(match)
	return state

def visit_frequencies(assignments):
	for atom in assignments.molecule:
		try:
			freqs=[w.freq for w in assignments.by_atom[atom]]
		except KeyError:
			pass

def remove_some(assignments,some_to_remove):
	for match in some_to_remove:
#		print 'remove ', match
		assignments.remove(match)

def symmetry_score(assignments):
	tr=Tracer('scoring.symmetry')
	from assignment import PeakMatch
	total_score=0
	total_ct=0
	total_sym=0
	for assignment in assignments:
		pair=assignment.distance_atoms()
		tr.Debug('')
		tr.Debug('for: ',assignment)
		ct=1
		ct_sym=1
		for experiment in assignments.peak_collection:
			rule=experiment[0].rule
			for expected in rule.expected_symmetric_peaks(molecule, pair[0]):
				if assignment.peak_match==expected: continue
				tr.Debug('expect: '+' '*30,expected)
				ct+=1
				try:
					match=assignments[expected]
					if match:
						ct_sym+=1
						tr.Debug('found ',', '.join(['%s'%m for m in match]))
				except KeyError as exc:
					pass
		total_ct+=ct
		total_sym+=ct_sym
		total_score+=float(ct_sym)/ct
	tr.Info('total peaks: %6d with symmetry: %6d'%(total_ct,total_sym))
	return total_score


def dump_score(assignments):
	offset=0
	total_score=0
	keys=sorted(assignments.scores.assignments.iterkeys(),key=lambda ass: ass.peak)
	for i,key in enumerate(keys):
		try:
			sym_group=assignments.scores.assignments[key]
			a=sym_group['symmetry']
			total_score+=a.score()
			print '%d'%(i+offset),sym_group['symmetry']
		except KeyError:
			offset-=1
		#	print 'no entry for ',ass
	print 'total= ',total_score
try:
	library.hello( __file__ )
	if args.fasta:
		sequence=fasta.read_fasta( args.fasta )

	molecule=AtomTree.from_sequence( sequence )

	peak_collection=PeakCollection()
	for file in args.peaks:
		name=path.splitext(path.basename(file))[0]
		print 'read peak-list %s from %s...'%(name,file)
		peak_collection.add_experiment( PeakList.read_from_stream( name, open(file,'r'), args.ignore ) )
	print 'done reading... '

#	from assignment.scoring.methods.FragDistanceScore import FragDistanceScore
#	frag_dist=FragDistanceScore()
#	frag_dist.fill_expected_peaks( peak_collection, molecule )
#	print peak_collection.expected

#	for assignment in peak_collection.expected:
#		print assignment
#	frag_dist_score2=score.methods.FragDistanceScore()
#	frag_data=FragDistanceDistribtionLibrary.get_from_fragdist_table(
#	print 'start timing.. '
	from timeit import Timer
	from assignment.scoring.methods.ScoreMethod import ScoreCache
	from assignment import random_items
#	t = Timer(lambda : assign_randomly(peak_collection,molecule))
#	print t.timeit(10)/10
#	assignments=assign_randomly(peak_collection,molecule)
	assignments=AssignmentCollection.from_hard_assignments(peak_collection,molecule)
# 	few_assignments=AssignmentCollection(peak_collection,molecule)
# 	small=2000
# 	large=1300
 	remove=100
# 	for i,x in enumerate(assignments):
# 		if i<small and i>large:
# 			few_assignments.add(x)
# 	assignments=few_assignments

	from assignment.scoring import ScoreFunction
	scorefxn=ScoreFunction(bmrb=1,frag=0,consistency=1,symmetry=1,expected=0)
	print 'scorefxn: ',scorefxn
#	t=Timer(lambda: scorefxn(assignments))
#	scorefxn.print_scores(assignments)
#	print t.timeit(1)
#	print t.timeit(1)

#	symmetry=symmetry_score(assignments)
#	print 'old symmetry: ',symmetry

#	assignments.write_to_stream(sys.stdout)
#	frag=frag_score(assignments)
	some_to_remove=list([ x for i,x in enumerate( assignments ) if i<remove])
	scorefxn.print_scores(assignments)
	scorefxn.print_scores(assignments)
#	t2=Timer(lambda: remove_some(assignments,some_to_remove))
	import copy
	def copy_and_score(ass):
		new_ass=copy.copy(ass)
		remove_some(new_ass,some_to_remove)
		scorefxn.print_scores(new_ass)
#		print scorefxn(new_ass)

	tc=Timer(lambda: copy_and_score(assignments))
	#print tc.timeit(10)/10
	print 'make copy'
	old_assignments=copy.copy(assignments)
	copy_and_score(assignments)
	scorefxn.print_scores(old_assignments)

	print 'remove some...'
	remove_some(assignments,some_to_remove)
	print 'score new assignments...'
	scorefxn.print_scores(assignments)
	assignments.scores=ScoreCache()
	scorefxn.print_scores(assignments)
	print 'score copy of original assignents...'
	scorefxn.print_scores(old_assignments)
	print 'now add assignments again'
	for match in some_to_remove:
#		print 'remove ', match
		assignments.add(match)
	print 'score new assignments...'
	scorefxn.print_scores(assignments)
	assignments.scores=ScoreCache()
	scorefxn.print_scores(assignments)

#	dump_score(assignments)
	print 're-read original list'
	assignments=AssignmentCollection.from_hard_assignments(peak_collection,molecule)
	few_assignments=AssignmentCollection(peak_collection,molecule)
	for i,x in enumerate(assignments):
		if i<small and i>large:
			few_assignments.add(x)
	assignments=few_assignments
	scorefxn.print_scores(assignments)
	print 'removed same 100'
	remove_some(assignments,some_to_remove)
#	dump_score(assignments)
	scorefxn.print_scores(assignments)
#	dump_score(assignments)
	exit(1)

#	consistency=consistency_score(assignments)
#	symmetry=symmetry_score(assignments)
	print 'total scores: '
#	print '%-30s: %8.3f'%('frag_score',frag)
	print '%-30s: %8.3f'%('bmrb_score',bmrb)
#	print '%-30s: %8.3f'%('consistency_score',consistency)
#	print '%-30s: %8.3f'%('symmetry_score',symmetry)
	exit(1)
	cum_score=0
#	symmetry_score(assignments)
	t=Timer(lambda : 	symmetry_score(assignments) )
	print t.timeit(1)
	print 'visit freqs: ',Timer(lambda : visit_frequencies(assignments)).timeit(100)
	exit(1)
	for atom in molecule:
		result='no assignment'
		try:
			freqs=[w.freq for w in atomic_assignment[atom]]
			cons_score=consistency_score( freqs )
			cum_score+=cons_score
			result='['+', '.join(['%8.3f'%w for w in freqs])+']'
		except KeyError:
			pass
		print atom,'%8.3f'%cons_score,result
	print 'total_score: ', cum_score
except Exception:
	library.print_exception( args.traceback )
