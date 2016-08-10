#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#score method based on fragment distance distribution

import math
import library
from basic import options
from basic import Tracer
from assignment.scoring import score_definitions
from assignment.scoring import libraries
from ScoreMethod import ScoreMethod,ScoreOptions,ScoreItem
from assignment import PeakList, Peak

SCORE_TYPE='frag'
parser = options.ModuleArgumentParser("%s-score"%SCORE_TYPE,
																			prefix=SCORE_TYPE,
																			description='options for scoring against distance of structures',
																			add_help=False )
parser.add_argument("-dist",help="distance table generated from ", default=None )
parser.add_argument("-norm", help='expect x-percent of frag_dist below noesy_cut_off [default 0.3]', type=float, default=0.3 )

tr=Tracer('scoring.%s'%SCORE_TYPE)

class FragDistOptions(ScoreOptions):
	def __init__(self,norm=0):
		if not norm:
			args=parser.parse_args()
			norm=args.norm
	  self.prob_norm=1-norm

#######################
##
## class FragDistanceScore
##
## reads a dist-table generated with FragsToAtomDist.default.linuxgccrelease
## NOESY assignments to atom-pairs with distances small according to fragments
## get favorable score
## we evaluate the integral of distance distribution above the noesy cutoff
## if almost all distances are above the noesy cutoff it is unlikely that this assignment is correct
## however: what about long-range assignments ?
## here we have no data and return a neutral 0
##
## TODO: should maybe have some smooth weighting scheme that takes into account that information on 1-9 distances is less accurate
## than information on 1-1 and 1-2 distances
class FragDistanceScore(ScoreMethod):
	dist_lib=None
	#end DistanceMatcher

	def __init__(self,options=None):
		if not options:
			options=FragDistOptions()
		#obtain raw data, either from cmd-line or as list of stuff
		super(FragDistanceScore,self).__init__(SCORE_TYPE,options)
		if not FragDistanceScore.dist_lib:#read fragment data if necessary
			FragDistanceScore.load_library()

	@staticmethod
	def load_library():
		import time
		start = time.clock()
		args=parser.parse_args()
		if not args.dist:
			raise library.MissingInput('A fragment distance table is required. Load with -frag_dist')
		tr.Info('reading fragment distance data from %s...'%args.dist)
		FragDistanceScore.dist_lib=libraries.FragDistanceData.get_from_fragdist_table(open(args.dist,'r'))
		elapsed = (time.clock() - start)
		tr.Info('storing of histograms took %5.1f seconds.'%elapsed)

	def score(self,assignments):
		super(FragDistanceScore,self).score(assignments)
		#initialize q0 from norm
		q0=abs(math.log(1.0-self.options.prob_norm))
		total_score=0
		for assignment in assignments:
			try:
				total_score+=assignments.scores.assignments[assignment][self.name].score
			except KeyError:
				score=self.score_assignment( assignment, q0 )
				assignments.scores.assignments.setdefault(assignment,{})[self.name]=ScoreItem(score)
				total_score+=score
		return total_score

	def score_assignment( self, assignment, q0=None ):
		atom_pairs=assignment.distance_atoms()
		score=0
		if not q0: q0=abs(math.log(1.0-self.options.prob_norm))
		for pair in atom_pairs:
			if abs(pair[0].resid-pair[1].resid)<9:
				score+=self._score_atom_pair(q0,assignment.rule.distance_cutoff,pair[0],pair[1])
			#		print 'frag_score=%8.3f for %s'%(score,assignment)
		return score

	def _score_atom_pair(self,q0,noesy_dist_cutoff,atom1,atom2):
		if atom1==atom2:
			return 1.0

		distribution=self.dist_lib.by_atom_pair(atom1,atom2)
		#compute the excluded probability weight
		bad_prob=distribution.integral(noesy_dist_cutoff,100)
		try:
			q=math.log(1-bad_prob)
		except ValueError:
			q=score_definitions.VERY_LOW_LOG
#		print 'p(no_contact)=%5.2f score=%7.2f'%(bad_prob,1.0+q/self.q0)
		return 1.0+q/q0

	def notify_add(self, state, new_assignment):
		score=self.score_assignment( new_assignment )
		state.scores.assignments.setdefault(new_assignment,{})[self.name]=ScoreItem(score)

	def notify_remove(self, scores, old_assignment):
		#cannot nuke the whole assignment entry, that could mess-up other ScoreMethods, e.g., the SymmetryScore
		scores.assignments[old_assignment].pop(self.name)

	def copy_cache( self, new_cache, old_cache ):
		new_cache.simple_cache_copy_from( old_cache, self.name, ['assignments'] )

	#this fills the list of expected peaks in the PeakCollection
	# peaks are expected for every experiment with distance-based rule, where frag_dist is below a certain cutoff
	def fill_expected_peaks(self, peak_collection, molecule):
		from assignment.rules import ScoreDistanceMatcher
		expected_peaks=PeakList('expected')
		peak_collection.expected=expected_peaks
		for name,peak_list in peak_collection.experiments.iteritems():
			tr.Info('compute expected peaks for peak list %s...'%name)
			try:
				exp_peak=peak_list[0]
				peak_dim=exp_peak.dim
				rule=exp_peak.rule
				peak=Peak((None,)*peak_dim,rule,id=SCORE_TYPE,pl_name='expected_'+name)
				distance_matcher=ScoreDistanceMatcher(self, abs(math.log(1.0-self.options.prob_norm)), 0 )
				expected_peaks+=[ m for m in peak.matches( molecule, random_samples=None, frequency_matcher=None, distance_matcher=distance_matcher ) ]
			except IndexError:
				continue
#			break
#		tr.Warning('[WARNING] for testing early bail-out from peak-generation')
		tr.Info('done!')
import ScoreMethodManager
ScoreMethodManager.register(SCORE_TYPE,FragDistanceScore)


import unittest

#### unit testing
class FragDistanceScoreDistributionTestCase(unittest.TestCase):
	def setUp(self):
		s='''
DATA SEQUENCE MNLTVNGKPS TVDGAESLNV TELLSALKVA QAEYVTVELN GEVLEREAFD
DATA SEQUENCE ATTVKDGDAV EFLYFM

VARS   RESID1 ATOMNAME1 RESID2 ATOMNAME2 LOW HIGH

    1    H    1     HB2     2.171     3.609       9       4      31      34      30      38      14       9       7       4       0       0       0       0       0       1       0       0       0      19
    1    H    1     HE2     1.567     6.976       1       0       0       0       0       8       1       0       9       4      24       8      11       4      25      54      24       0      12      15
    3 HD22    4       H     3.207     9.754       6       9      18      12       9      36      40      49      21       3      10       7      16      23      28      29      29      32      10      10
    3 HD22    4     HG1     5.037     9.847       2       2      10      13      17      24      20      18      25      37      27      36      35      32      16      26      23      19      11       4
    4 HG23    5    HG12     5.416     9.715       1       0       3       2       1      14      14      13      30      58      78      82      59       8       1       0       2       3       9      21
    4 HG23    5    HG21     4.859     10.17       3       5      10       9       5       1       4      16      59      47      16      20      89      58      23      30       3       0       0       1
    8  HE2    9     HD3     1.155     8.204       1       0       0       0       1      12      16       8      17       7      20      48      55      61      48      61      29      10       4       1
    8  HE2    9     HG2     3.195     10.71       1       1       2       7       7      18      15      10      23      28      24      20      28      19      30      46      42      35      35       8
   11  HG1   13     HB3     6.272     8.808       1       2       2       4       5       9       7      19      11      13       9      10      29      14      20      13       9       9       5       5
   11  HG1   13      QG     5.386     7.105       1       1       6       5       7       6      13      22      20      12      19      18      20      14      12       9       5       2       1       3
'''
		from StringIO import StringIO
		lines=StringIO(s)
		self._frag_distance_distribution_score_method=FragDistanceScoreDistribution(lines)
		self._atom_list1=[Atom('H',1),
											Atom('H',1),
											Atom('HD22',3),
											Atom('HD22',3),
											Atom('HG23',4),
											Atom('HG23',4),
											Atom('HE2',8),
											Atom('HE2',8),
											Atom('HG1',11),
											Atom('HG1',11)]

		self._atom_list2=[Atom('HB2',1),
											Atom('HE2',1),
											Atom('H',4),
											Atom('HG1',4),
											Atom('HG12',5),
											Atom('HG21',5),
											Atom('HD3',9),
											Atom('HG2',9),
											Atom('HB3',13),
											Atom('QG',13)]

	def test_high_and_low(self):
		real_high=[3.609,6.976,9.754,9.847, 9.715, 10.17,8.204,10.71,8.808,7.105]
		real_low =[2.171,1.567,3.207,5.037,5.416,4.859,1.155, 3.195,6.272, 5.386]
		high=[]
		low=[]
		for atom1,atom2 in zip(self._atom_list1,self._atom_list2):
			high.append(self._frag_distance_distribution_score_method.distribution_library().by_atom_pair(atom1,atom2).high())
			low.append(self._frag_distance_distribution_score_method.distribution_library().by_atom_pair(atom1,atom2).low())
		for h,rh,a1,a2 in zip(high,real_high,self._atom_list1,self._atom_list2):
			self.assertAlmostEqual(h,rh,3,'highest distance of atom %s and %s should be %5.3f but we get %5.3f (which is wrong)'%(a1,a2,rh,h))
		for l,rl,a1,a2 in zip(low,real_low,self._atom_list1,self._atom_list2):
			self.assertAlmostEqual(l,rl,3,'lowest distance of atom %s and %s should be %5.3f but we get %5.3f (which is wrong)'%(a1,a2,rl,l))
	def test_seq(self):
		real_seq='MNLTVNGKPSTVDGAESLNVTELLSALKVAQAEYVTVELNGEVLEREAFDATTVKDGDAVEFLYFM'
		seq=self._frag_distance_distribution_score_method.distribution_library().seq()
		self.assertEqual(seq,real_seq,'The sequence of this frag2dist file should be %s \n but we get %s'%(real_seq,seq))

	def test_hist_and_edges(self):
		real_hist=[0.0450,  0.0200,  0.1550,  0.1700,  0.1500,  0.1900,  0.0700,  0.0450,  0.0350,  0.0200,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0050,  0.0000,  0.0000,  0.0000,  0.0950]
		real_edges=[2.1710,  2.2429,  2.3148,  2.3867,  2.4586,  2.5305,  2.6024,  2.6743,  2.7462,  2.8181,  2.8900,  2.9619,  3.0338,  3.1057,  3.1776,  3.2495,  3.3214,  3.3933,  3.4652,  3.5371,  3.6090]
		hist=self._frag_distance_distribution_score_method.distribution_library().by_atom_pair(Atom('H',1),Atom('HB2',1)).hist()
		edges=self._frag_distance_distribution_score_method.distribution_library().by_atom_pair(Atom('H',1),Atom('HB2',1)).bin_edges()
		for h,rh in zip(hist,real_hist):
			self.assertAlmostEqual(h,rh,3,'one distance hist of should be %5.3f but we get %5.3f (which is wrong)'%(rh,h))
		for e,re in zip(edges,real_edges):
			self.assertAlmostEqual(e,re,3,'one distance hist edge of should be %5.3f but we get %5.3f (which is wrong)'%(re,e))

def FragDistanceScoreDistributionTestSuite():
	suite = unittest.TestLoader().loadTestsFromTestCase(FragDistanceScoreDistributionTestCase)
	return suite


if __name__ == '__main__':
	unittest.main()
