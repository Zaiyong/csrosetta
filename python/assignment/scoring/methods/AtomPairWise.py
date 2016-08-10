#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#similar with AtomWiseScoreMethod

from ScoreMethod import ScoreMethod
import library

class AtomPairWiseScoreMethod(ScoreMethod):
	def __init__(self,name):
		ScoreMethod.__init__(self,name)

	def score_atom_pair(self,input,atom1,atom2):
		raise library.StubbedOut()

	def __call__(self,cross_peak_list):
		sum=0
		for peak in cross_peak_list.iter_peaks():
			dist=peak.volume()**(1.0/6)
			for assign in peak.iter_assignments():
				protons=assign.protons()
				try:
					score=self.score_atom_pair(dist,protons[0],protons[1])
					sum+=score
				except KeyError:#not all atom pairs have the distribution
					continue
		return sum
