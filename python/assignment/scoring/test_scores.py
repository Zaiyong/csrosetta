#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os import path
import unittest
import assignment
from ScoreFunction import ScoreFunction
import libraries
import fasta
from assignment import PeakCollection, PeakList, AssignmentCollection
from chemical import AtomTree


#########
##  general unit-test setups
#get data path
import basic
data_path = basic.get_unittest_data()
basic.fix_unittest_args()
#
#fix the FragDistance Score
from methods.FragDistanceScore import parser as frag_dist_parser
frag_dist_parser.set_testing_args(dist=data_path+'dist_table.dat')
#
#fix Tracers
from basic.Tracer import init_tracer_from_cmdline
init_tracer_from_cmdline(['assignment.collection:error','assignment.peaks:warning','scoring.methods:warning'])

######
## The Test
class ScoreFunctionTestCase(unittest.TestCase):

	def __init__(self,name):
		super(ScoreFunctionTestCase,self).__init__(name)

	@classmethod
	def setUpClass(cls):
	#setup code that really should only run once
		print 'initialize ScoreFunctionTestCase...'
		cls.sequence=fasta.read_fasta( data_path+'gmr137.fasta' )
		cls.molecule=AtomTree.from_sequence( cls.sequence )
		cls.peak_collection=PeakCollection()
		peaks=['aroC.peaks','n.peaks']
		for p in peaks:
			file=data_path+'assigned/'+p
			name=path.splitext(path.basename(file))[0]
			cls.peak_collection.add_experiment( PeakList.read_from_stream( name, open(file,'r'), False ) )
		cls.assignments=AssignmentCollection.from_hard_assignments( cls.peak_collection, cls.molecule )
		cls.some_to_remove=list([ x for i,x in enumerate( cls.assignments ) if i<100])

	def setUp(self):
		pass

  def test_score_weights(self):
		ass=self.assignments
		bmrb1=ScoreFunction(bmrb=1)(ass)
		sym1=ScoreFunction(symmetry=1)(ass)
		self.assertEqual(bmrb1*2.3,ScoreFunction(bmrb=2.3)(ass))
		self.assertEqual(bmrb1*0.1,ScoreFunction(bmrb=0.1)(ass))
		self.assertEqual(sym1*0.1,ScoreFunction(symmetry=0.1)(ass))
		self.assertEqual(sym1*0.1+bmrb1*3,ScoreFunction(symmetry=0.1,bmrb=3)(ass))
		self.assertEqual(sym1*0.1+-3*bmrb1,ScoreFunction(symmetry=0.1,bmrb=-3)(ass))

	def test_local_weights(self):
		ass=self.assignments
		score_bmrb1=ScoreFunction(bmrb=1)
		score_sym1=ScoreFunction(symmetry=1)
		score_cons1=ScoreFunction(consistency=1)
		#apply scores to fill local score-info
		score_bmrb1(ass)
		score_sym1(ass)
		score_cons1(ass)

		cache=ass.scores.atomic
		key=cache.keys()[0]
		self.assertEqual(0,score_sym1.evaluate_local_scores(cache[key]))
		self.assertNotEqual(0,score_bmrb1.evaluate_local_scores(cache[key]))
		cache=ass.scores.assignments
		#find a non-zero symmetry entry
		for key in cache.iterkeys():
			if cache[key]['symmetry'].score!=0:
				break
		self.assertNotEqual(0,score_sym1.evaluate_local_scores(cache[key]))
		self.assertEqual(0,score_bmrb1.evaluate_local_scores(cache[key]))

		cache=ass.scores.atomic
		key=cache.keys()[0]
		bmrb1=score_bmrb1.evaluate_local_scores(cache[key])
		cons1=score_cons1.evaluate_local_scores(cache[key])
		self.assertEqual(bmrb1*0.2+cons1*3.2,ScoreFunction(bmrb=0.2,consistency=3.2).evaluate_local_scores(cache[key]))
		self.assertEqual(-bmrb1+cons1*3.2,ScoreFunction(bmrb=-1,consistency=3.2).evaluate_local_scores(cache[key]))

	def test_score_option_change(self):
		ass=self.assignments
		score_bmrb1=ScoreFunction(bmrb=1)
		#what happens if we change the score options
		orig_score=score_bmrb1(ass)

		#prepare change of options
		from methods.BmrbShiftScore import BmrbShiftScoreOptions
		x0=score_bmrb1.methods['bmrb'].options.x0
		new_options=BmrbShiftScoreOptions(x0/2)

		#change options in place
		score_bmrb1.methods['bmrb'].options=new_options
		self.assertNotEqual(orig_score,score_bmrb1(ass))
		score_with_x0half=score_bmrb1(ass)

		#try illegal change of options in place
		def do_bad():
			score_bmrb1.methods['bmrb'].options.x0=new_options.x0
		self.assertRaises(TypeError, do_bad )

		#now change it by differently
		from methods.BmrbShiftScore import BmrbShiftScore
		new_method=BmrbShiftScore(BmrbShiftScoreOptions(x0*2))
		score_bmrb1.append(new_method,1.0)
		self.assertNotEqual(orig_score,score_bmrb1(ass))
		self.assertNotEqual(score_with_x0half,score_bmrb1(ass))

	def test_score_init(self):
		self.assertRaises(KeyError, ScoreFunction,gogleplex=1)
		self.assertIsNotNone(ScoreFunction(gogleplex=0))
		self.assertIsNotNone(ScoreFunction(bmrb=1,symmetry=1))


