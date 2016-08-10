#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from basic import options
from basic import Tracer
from ScoreMethod import ScoreMethod,ScoreOptions,ScoreItem
from utility.SchmidtQ import quantity_function,schmidt_final_score
import numpy
from BmrbShiftScore import BmrbShiftScore
SCORE_TYPE='consistency'
parser = options.ModuleArgumentParser("%s-score"%SCORE_TYPE, prefix=SCORE_TYPE,description='options for the BMRB shift score',add_help=False )
#parser.add_argument("-bmrb_shifts",help="instead of database file read shifts from here", default=None )
parser.add_argument("-x0", help="schmidt_final_score [default 0.5]", default=0.5)

tr=Tracer('scoring.%s'%SCORE_TYPE)

class ConsistencyScoreOptions(ScoreOptions):
	def __init__(self,x0=None):
		if not x0:
			args=parser.parse_args()
			self.x0=args.x0
		else:
			self.x0=x0

######################
##
## class ConsistencyScore
##
## a score-method to check consistency of frequency assignments to a single atom
## the deviation from the mean-assigned frequency is computed and normalized by the BMRB shift distribution of this atom
##
## TODO: maybe different options for the norm should exist:
## in a pre-assigned resonance case, for instance one would maybe take fixed tolerances for C,H,N type ?
## if more stringent frequency distributions are present (i.e., SpartaShiftDistrubtion) these could be better for the norm here...
## in such a case the selection of norming mode could be part of the options
##
class ConsistencyScore(ScoreMethod):

	def __init__(self,options=None):
		if not options:
			options=ConsistencyScoreOptions()
		super(ConsistencyScore,self).__init__(SCORE_TYPE,options)
		args=parser.parse_args()
		#configureing from cmd-line
		#need the bmrb lib for the normalization by std(bmrb_shift_distribution)
		self._bmrb_library=BmrbShiftScore()

	def score(self,assignments):
		super(ConsistencyScore,self).score(assignments)
		total_score=0
		for atom in assignments.by_atom.iterkeys():
			try: #do we have score on the atom-level already ?
				total_score+=assignments.scores.atomic[atom][self.name].score
			except KeyError:
				freqs=[w.freq for w in assignments.by_atom[atom]]
				atom_score=0
				if len(freqs)>1:
					N=len(freqs)
					mean=float(sum(freqs))/N
					norm=self._bmrb_library[atom].std()
					for freq,atomic in zip(freqs,assignments.by_atom[atom]):
						val=schmidt_final_score(quantity_function((freq-mean)/norm),0.5)/N
						atom_score+=val
						assignments.scores.atomic_matches.setdefault(atomic,{})[self.name]=ScoreItem(val)
							#finish except block: no cached atomic_assignment score
				assignments.scores.atomic.setdefault(atom,{})[self.name]=ScoreItem(atom_score)
				total_score+=atom_score
				#finish except block: no cached atom score
		return total_score
#	tr.Debug('TOT: %8.3f norm=%5.2f'%(sum(values)/len(values),norm),
#								 ', '.join(['%6.2f : %6.2f'%(freq,val) for freq,val in zip(freqs,values)]))

	def notify_add(self, assignment_collection, new_assignment):
		scores=assignment_collection.scores
		for atomic in new_assignment:
			try:
				scores.atomic[atomic.atom].pop(self.name,None)
			except KeyError:
				pass

	def notify_remove(self, scores, old_assignment):
		for atomic in old_assignment:
			try:
				scores.atomic[atomic.atom].pop(self.name,None)
			except KeyError:
				pass

	def copy_cache( self, new_cache, old_cache ):
		new_cache.simple_cache_copy_from( old_cache, self.name, ['atomic_matches','atomic'] )

import ScoreMethodManager
ScoreMethodManager.register(SCORE_TYPE,ConsistencyScore)

import unittest
# class ConsistencyScoreTestCase(unittest.TestCase):
# 	def setUp(self):
# 		self.BMRB_shift_distribution_score_method=ConsistencyScore('AHMTRIL')
# 		self._seq=self.BMRB_shift_distribution_score_method.shift_distribution_library().seq()
# 	def test_load_shift_library_from_BMRB_database(self):
# 		real_mean=[4.250,4.62,4.43,4.46,4.30,4.17,4.32]
# 		real_std=[0.46,0.63,2.92,0.49,0.48,0.58,0.72]
# 		mean=[]
# 		std=[]
# 		for i in range(1,8):
# 			mean.append(self.BMRB_shift_distribution_score_method.shift_distribution_library().by_atom(Atom('HA',i)).mean())
# 			std.append(self.BMRB_shift_distribution_score_method.shift_distribution_library().by_atom(Atom('HA',i)).std())
# 		for m,r,a in zip(mean,real_mean,self._seq):
# 			self.assertAlmostEqual(r,m,2,'shift mean of proton HA of aomino acid %c should be %5.3f but we get %5.3f (which is wrong)'%(a,r,m))
# 		for s,rs,a in zip(std,real_std,self._seq):
# 			self.assertAlmostEqual(rs,s,2,'shift std of proton HA of aomino acid %c should be %5.3f but we get %5.3f (which is wrong)'%(a,rs,s))

# def ConsistencyScoreTestSuite():
# 	suite = unittest.TestLoader().loadTestsFromTestCase(ConsistencyScoreTestCase)
# 	return suite

# if __name__ == '__main__':
# 	unittest.main()
