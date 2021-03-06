#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os import path
import copy
import unittest
import assignment
from ScoreFunction import ScoreFunction
from methods.FragDistanceScore import parser as frag_dist_parser
from assignment import PeakCollection, PeakList, AssignmentCollection
from chemical import AtomTree

import libraries
import fasta

import basic
data_path = basic.get_unittest_data()
basic.fix_unittest_args()

from basic.Tracer import init_tracer_from_cmdline

frag_dist_parser.set_testing_args(dist=data_path+'dist_table.dat')
init_tracer_from_cmdline(['assignment.collection:error','assignment.peaks:warning','scoring.frag:warning'])
#consistency_parser.set_testing_args()
#bmrb_parser.set_testing_args()
#symmetry_parser.set_testing_args()

def remove_some(assignments,some_to_remove):
	for match in some_to_remove:
		try:
			assignments.remove(match)
		except ValueError:
			pass

class ScoreCachingTestCase(unittest.TestCase):

	def __init__(self,name):
		super(ScoreCachingTestCase,self).__init__(name)

	@classmethod
	def setUpClass(cls):
	#setup code that really should only run once
		print 'initialize ScoreCachingTestCase...'
		cls.sequence=fasta.read_fasta( data_path+'gmr137.fasta' )
		cls.molecule=AtomTree.from_sequence( cls.sequence )
		cls.peak_collection=PeakCollection()
		peaks=['aliC.peaks','aroC.peaks','n.peaks']
		for p in peaks:
			file=data_path+'assigned/'+p
			name=path.splitext(path.basename(file))[0]
#			print 'read peak-list %s from %s...'%(name,file)
			cls.peak_collection.add_experiment( PeakList.read_from_stream( name, open(file,'r'), False ) )
#		print 'done reading... '
		cls.assignments=AssignmentCollection.from_hard_assignments( cls.peak_collection, cls.molecule )
		cls.scorefxn=ScoreFunction(bmrb=1, frag=1, symmetry=1, consistency=1)
		cls.some_to_remove=list([ x for i,x in enumerate( cls.assignments ) if i<100])
		cls.some_to_remove2=list([ x for i,x in enumerate( cls.assignments ) if i<200])

	def setUp(self):
		#restart score-caching
		self.assignments.clear_scores()

	def test_copy_assignments(self):
		scorefxn = self.scorefxn
		orig_assignments = copy.copy( self.assignments )
		other_assignments = copy.copy( self.assignments )
		score_orig = scorefxn( orig_assignments )
		score_other = scorefxn( other_assignments )
		self.assertEqual( score_orig, score_other )
		#now remove some from other, and see if score of orig is stable
		remove_some( other_assignments, self.some_to_remove )
		score_other = scorefxn( other_assignments )
		self.assertEqual( score_orig, scorefxn( orig_assignments ) ) #should just reload scores from cache
		self.assertAlmostEqual( score_orig, scorefxn( self.assignments ) ) #now we score the original assignments for comparison
		self.assertNotEqual( score_orig, score_other ) #score of other_assignment should be different
		other_assignments.clear_scores()
		self.assertEqual( score_other, scorefxn( other_assignments ) ) #re-calculate from scratch, should be the same
		copy_scored=copy.copy( other_assignments )
		remove_some( copy_scored, self.some_to_remove2 )
		score_double_removed=scorefxn( copy_scored )
		copy_scored.clear_scores()
		self.assertEqual( score_double_removed, scorefxn( copy_scored ) )

		#remove and add continuously
		from assignment import random_items
		state=copy_scored
		for peak in random_items( state.peak_collection.iterpeaks(), 300 ):
			assignments=state.by_peak.get( peak )
			if assignments:
				for match in assignments:
					state.remove( match )
			try:
				for match in peak.matches( state.molecule, 1 ): #, distance_matcher=dist_matcher ):
					state.add(match)
			except KeyError:
				pass #some of the aromatics have CE, CD pre-assignments which mess us up here #print 'keyerror for ', peak.hard_assignments
	  score_after_confusion = scorefxn( state )
		state.clear_scores()
		self.assertAlmostEqual( score_after_confusion, scorefxn( state ) )

	def test_scores_add_remove(self):
		from assignment.rules import ScoreDistanceMatcher
		from assignment.scoring.methods import FragDistanceScore
		from assignment import random_items
		import math
		known_dist=ScoreDistanceMatcher( FragDistanceScore(), abs(math.log(0.3)), 0 )
		known_dist.max_sequence_separation=2
		state=copy.copy( self.assignments )
		new_matches=[]
		has_diag=0
		for peak in random_items( state.peak_collection, 50 ):
			npeak=copy.copy( peak )
			npeak.hard_assignments=[]
			for match in npeak.matches( state.molecule, 1, distance_matcher=known_dist ):
				if match.peak_match[1]==match.peak_match[2]:
					new_matches.append( match )
					has_diag+=1
			if has_diag > 10: break
		for peak in random_items( state.peak_collection, 30 ):
			npeak=copy.copy( peak )
			npeak.hard_assignments=[]
			for match in npeak.matches( state.molecule, 1, distance_matcher=known_dist ):
				if match.peak_match[1]!=match.peak_match[2]:
					new_matches.append( match )

		def test_individual_scores_add_remove(state,matches,score_type):
			scorefxn = ScoreFunction().set_weight(score_type,1.0)
			state.clear_scores()
			state=AssignmentCollection(state.peak_collection, state.molecule)
			orig_score = scorefxn( state )
			for match in matches:
				state.add( match )
			new_score = scorefxn( state )
			score_cache=state.scores
			state.clear_scores()
			new_score2 = scorefxn( state )
			state.scores=score_cache
			for match in matches:
				try:
					state.remove( match )
					state.commit()
				except:
					import library
					library.augment_exception_msg('ScoreMethod: %s cannot handle removal of match %s'%(score_type, match ) )

			final_score = scorefxn( state )
			state.clear_scores()
			final_score2 = scorefxn( state )
			try:
				self.assertNotEqual( orig_score, new_score )
				self.assertEqual( new_score, new_score2 )
				self.assertEqual( final_score, final_score2 )
				self.assertEqual( orig_score, final_score )
			except:
				import library
				match_str=''
				for match in matches:
					match_str+='%s\n'%match
				library.augment_exception_msg('failure in ScoreMethod: %s\n%s'%(score_type,match_str))

		from methods.ScoreMethodManager import score_methods
		for type in score_methods.iterkeys():
			if type=='expected': continue #these usually have 0 for the selected matches
			if type=='symmetry': continue #the error tested here was coming from diagonal peaks which have no particular meaning for these anyways
			test_individual_scores_add_remove(state,new_matches,type)

	def test_cache_is_equal( self ):
		scorefxn = self.scorefxn
		score_first = scorefxn( self.assignments )
		score_second = scorefxn( self.assignments )
		self.assertEqual( score_first, score_second )

	def test_cache_can_remove(self):
		scorefxn = self.scorefxn
		assignments = copy.copy( self.assignments )
		#score to fill cache
		score = scorefxn( assignments )
		remove_some( assignments, self.some_to_remove )
		#score should come mostly from cached items
		score_cached = scorefxn( assignments )
		assignments.clear_scores()
		#now this is a clean score from scratch, should be the same
		score_direct = scorefxn( assignments )
		self.assertEqual( score_cached, score_direct )
		#can we also add the same assignments and get back to the original score?
		for match in self.some_to_remove:
			assignments.add(match)
		score_cached = scorefxn( assignments )
		assignments.clear_scores()
		#now this is a clean score from scratch, should be the same
		score_direct = scorefxn( assignments )
		self.assertEqual( score_cached, score_direct )
		self.assertAlmostEqual( score, score_direct )

	def test_cache_saves_time(self):
		scorefxn=self.scorefxn
		assignments=self.assignments
		from timeit import Timer
		t=Timer(lambda: scorefxn(assignments))
		t1=t.timeit(1)
		t2=t.timeit(1)
		self.assertLess(t2,t1*0.1)

		def test_individual_score(assignments,score_type,fails):
			assignments.clear_scores()
			t_sym=Timer(lambda: ScoreFunction().set_weight(score_type,1.0)(assignments))
			t1=t_sym.timeit(1)
			t2=t_sym.timeit(1)
			try:
				self.assertLess(t2,t1*0.2)
#				print 'score: %s, speedup %8.5f'%( score_type, t1/t2)
			except AssertionError as exc:
				fails.append(score_type)

		def test_individual_scores_copy_cache(assignments,score_type,fails):
			scorefxn=ScoreFunction().set_weight(score_type,1.0)
			scorefxn( assignments )
			copy_of = copy.copy( assignments )
			t_sym=Timer(lambda: ScoreFunction().set_weight(score_type,1.0)(copy_of))
			t1=t_sym.timeit(1)
			t2=t_sym.timeit(1)
			try:
#				print 'score: %s, copy_slow_down %8.5f'%( score_type, t1/t2)
				self.assertGreater(t2,t1*0.2)
			except AssertionError as exc:
				fails.append(score_type)

		from methods.ScoreMethodManager import score_methods
		fails=[]
		ignore='expected'
		for type in score_methods.iterkeys():
			test_individual_score(assignments,type,fails)
		test_ok=True
		if len(fails):
			import sys
			sys.stderr.write('NO SPEEDUP FOR: '+' '.join(fails)+' ')
			for st in fails:
				if not st in ignore:
					test_ok=False
		self.assertTrue(test_ok)

		cache_fails=[]
		for type in score_methods.iterkeys():
			if type in fails: continue
			test_individual_scores_copy_cache(assignments,type,cache_fails)
		test_ok=True
		if len(cache_fails):
			import sys
			sys.stderr.write(' <=>  CACHE FAILED TO COPY FOR: '+' '.join(cache_fails)+' ')
			for st in cache_fails:
				if not st in ignore:
					test_ok=False
		self.assertTrue(test_ok)


