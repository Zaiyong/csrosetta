#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from utility import GaussianDistribution
from basic import options
from basic import database
from noe_tools import unpack_NMR_atom_group
from PDB.Polypeptide import three_to_one,one_to_three
from utility.SchmidtQ import schmidt_final_score
import unittest
from ScoreMethod import ScoreMethod,ScoreOptions,ScoreItem
from basic import Tracer

SCORE_TYPE='bmrb'
parser = options.ModuleArgumentParser("%s-score"%SCORE_TYPE,
																			prefix=SCORE_TYPE,
																			description='options for the BMRB shift score',
																			add_help=False )
parser.add_argument("-shifts",help="instead of database file read shifts from here", default=None )
parser.add_argument("-x0", help="schmidt_final_score [default 1.5]", default=1.5)
####
## simple Atom class that contains a name and residue-type but is not bound to a particular sequence position.
## these Atoms are used to specify generic knowledge about atoms of a certain kind:
## for instance, the chemical shift distribution of LEU HN vs the distribution of SER HN
##

tr=Tracer('scoring.%s'%SCORE_TYPE)
class BmrbAtomType:
	def __init__(self,aa,name):
		if len(aa)>1: aa=three_to_one(aa)
		new_atom_names=unpack_NMR_atom_group(name,aa)
		self._aa=aa
		self._name=new_atom_names[0]

	def __str__(self):
		return '%s %s'%(one_to_three(self._aa),self._name)

	def __repr__(self):
		return '%s %s'%(one_to_three(self._aa),self._name)

	def __hash__(self):
		return hash((self._aa,self._name))

	def __eq__(self,other):
		return (self._aa,self._name)==(other._aa,other._name)

	@classmethod
	def from_atom(obj,atom):
		aa=three_to_one(atom.res3_type)
		new_atom_names=unpack_NMR_atom_group(atom.name,aa)
		obj=BmrbAtomType(aa,new_atom_names[0])
		return obj

class BmrbShiftScoreOptions(ScoreOptions):
	def __init__(self,x0=None):
		if not x0:
			args=parser.parse_args()
			self.x0=args.x0
		else:
			self.x0=x0

class BmrbShiftScore(ScoreMethod):
#	standard_options=None
	_library=None
	def __init__(self,options=None):
		if not options:
			options=BmrbShiftScoreOptions()
		super(BmrbShiftScore,self).__init__(SCORE_TYPE,options)
		#
		#configureing from cmd-line
		if not BmrbShiftScore._library:
			args=parser.parse_args()
			BmrbShiftScore.load_library(args)

	@staticmethod
	def load_library(args):
		#open database file,
		if not args.shifts:
			file=database.open('cs_distribution.txt')
		else:
			file=open(args.shifts,'r')
		#read
		BmrbShiftScore._library={}
		for line in file:
			tags=line.split()
			if tags[0]=='Res': continue #skip the header line
			#add distribution for an atom
			atom_type=BmrbAtomType( tags[0], tags[1] )
			BmrbShiftScore._library[atom_type]=GaussianDistribution( float(tags[6]), float(tags[7]) )

  #this just gives the BmrbLibrary for an atom, can also be used outside of this score
	def __getitem__(self,atom):
		return BmrbShiftScore._library[BmrbAtomType.from_atom(atom)]

	def score(self,assignments):
		super(BmrbShiftScore,self).score(assignments)
		score=0
		for atom in assignments.by_atom.iterkeys():
			try:
				score+=assignments.scores.atomic[atom][self.name].score
			except KeyError:
				atom_sum=0
				distribution=None
				for atomic_assignment in assignments.by_atom[atom]:
					try:
						atom_sum+=assignments.scores.atomic_matches[atomic_assignment][self.name].score
					except KeyError:
						try:
							val, distribution = self._score_atomic_assignment( atomic_assignment, distribution )
						except KeyError:
							break #no distribution for this atom -- go to next atom
						assignments.scores.atomic_matches.setdefault(atomic_assignment,{})[self.name]=ScoreItem(val)
						atom_sum+=val
				assignments.scores.atomic.setdefault(atom,{})[self.name]=ScoreItem(atom_sum)
				score+=atom_sum
		return score

	def _score_atomic_assignment( self, atomic_assignment, distribution):
		if not distribution:
			try:
				distribution=self[atomic_assignment.atom]
			except KeyError:
				tr.Warning('[WARNING] cannot find atom %s in database'%(atomic_assignment.atom))
				raise
		freq=atomic_assignment.freq
		val=schmidt_final_score(distribution.schmidt_q(freq),self.options.x0)
		return val, distribution

	def score_assignment( self, assignment ):
		val=0
		for atomic in assignment:
			val+=self._score_atomic_assignment( atomic, None )[0]
		return val

	def notify_add(self, assignments, new_assignment):
		for atomic in new_assignment:
			#val = 0
			try:
				val, dist = self._score_atomic_assignment( atomic, None )
			except KeyError:
				continue
			assignments.scores.atomic_matches.setdefault(atomic,{})[self.name]=ScoreItem(val)
			try:
				assignments.scores.atomic[atomic.atom][self.name].score+=val
			except KeyError:
				assignments.scores.atomic.setdefault(atomic.atom,{})[self.name]=ScoreItem(val)

#			del assignments.scores.atomic_matches[atomic][self.name]

	def notify_remove(self, scores, new_assignment):
		for atomic_assignment in new_assignment:
			assign_cache=scores.atomic_matches[atomic_assignment]
			try:
				atom_correction=assign_cache[self.name].score
				atom=atomic_assignment.atom
				scores.atomic[atom][self.name].score-=atom_correction
			except KeyError:
				del scores.atomic[atom][self.name]
				tr.Info('nuked cache in BmrbShiftScore ... ')
		#	pass #in this case we might have to nuke the atom-cache ?
		#cannot delete caches before done with all atomic_assignments due to
		#peak-matches that have the same atom twice
		for atomic_assignment in new_assignment:
			assign_cache=scores.atomic_matches[atomic_assignment]
			try:
				del assign_cache[self.name]
			except KeyError:
				pass

	def copy_cache( self, new_cache, old_cache ):
		new_cache.simple_cache_copy_from( old_cache, self.name, ['atomic_matches','atomic'] )

import ScoreMethodManager
ScoreMethodManager.register(SCORE_TYPE,BmrbShiftScore)

# class BmrbShiftScoreTestCase(unittest.TestCase):
# 	def setUp(self):
# 		self.BMRB_shift_distribution_score_method=BmrbShiftScore('AHMTRIL')
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

# def BmrbShiftScoreTestSuite():
# 	suite = unittest.TestLoader().loadTestsFromTestCase(BmrbShiftScoreTestCase)
# 	return suite

# if __name__ == '__main__':
# 	unittest.main()
