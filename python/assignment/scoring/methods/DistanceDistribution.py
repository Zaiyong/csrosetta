#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


from AtomPairWise import AtomPairWise
import library

class DistanceDistribution(AtomPairWise):

	def __init__(self,name,lib):
		AtomPairWise.__init__(self,name)
		self._distribution_library=lib
		self._x0=3.0
	def score_atom_pair(self,dist,atom1,atom2):
		distribution=self._distribution_library.by_atom_pair(atom1,atom2)
		return distribution.schmidt_final_score(dist,self._x0)

	def distribution_library(self):
		return self._distribution_library
