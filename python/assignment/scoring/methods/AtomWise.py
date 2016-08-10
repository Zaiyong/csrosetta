#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from ScoreMethod import ScoreMethod
import library
class AtomWise(ScoreMethod):
	def __init__(self,name):
		ScoreMethod.__init__(self,name)

	def score_atom(self,input,atom):
		raise library.StubbedOut()

	def __call__(self,cross_peak_list):
		sum=0
		for peak in cross_peak_list.iter_peaks():
			freq=peak.freq()
			for assign in peak.iter_assignments():
				for i,atom in enumerate(assign.iter_atoms()):
					score=self.score_atom(freq[i],atom)
					sum+=score
		return sum
