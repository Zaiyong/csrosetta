#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


#############
##
## class ScoreDistanceMatcher
##
## a generic distance matcher score-method adaptor
##
## given a score-method which implements the _score_atom_pair(q0, noesy_cutoff, atom1, atom2 ) method
## this adaptor creates a distance-matcher for peak-matching
##
class ScoreDistanceMatcher(object):
	def __init__(self, score_method, q0, threshold ):
		self.method=score_method
		self.q0=q0
		self.threshold=threshold
		self.max_sequence_separation=9
	#when is a distance of protons small enough to expect a peak
	def __call__(self, atom1, atom2, noesy_cutoff ):
		try:
			if atom1==atom2: return True
			if abs(atom1.resid-atom2.resid)>=self.max_sequence_separation: return False
			score=self.method._score_atom_pair(self.q0, noesy_cutoff,atom1,atom2)
			#if score>self.threshold: print atom1, atom2, '-->', score, self.noesy_cutoff
			return score > self.threshold
		except KeyError: #in this case, if we don't have it we don't expect a peak
			return False
