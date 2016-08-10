#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from ScoreFunction import ScoreFunction

class ScoreFunctionFactory:
	def __init__(self):
		self._score_function=ScoreFunction()
	def read_score_weights(self,file):
		lines=open(file).readlines()
		for line in lines:
			tags=line.split()
			if tags[0]=='bmrb_score':
				pass
#				self._score_function.append(BMRBShiftDistributionScoreMethod.load_shift_library_from_BMRB_database('bmrb_score'),float(tags[1]))
			elif tags[0]=='shift_score':
				pass
	#			self._score_function.append(ShiftDistributionScoreMethod('shift_score'),float(tags[1]))
		return self._score_function
