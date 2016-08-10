#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from ScoreMethod import ScoreMethod
import library
class AssignmentWise(ScoreMethod):
	def __init__(self,name):
		ScoreMethod.__init__(self,name)

	def score_assignment(self,assignment,input):
		raise library.StubbedOut()

	def __call__(self,cross_peak_list):
		sum=0
		for peak in cross_peak_list.iter_peaks():
			for assignment in peak.iter_assignments():
				score=self.score_assignment(assignment,cross_peak_list)
				sum+=score
		return sum
