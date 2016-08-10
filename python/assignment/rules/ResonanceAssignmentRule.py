#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module

from BasicRules import *

##################
#
# class ResonanceAssignmentRule
#
# this rule is used to hijack the peak-rule system to provide pre-existing resonance assignments to the system
# in a transparent manner.
# the resonances are modelled as 1D peaks with pre-existing assignments (which is usually available)
# if multiple assignments are possible, just add more assignments to the 1D peak
# if multiple assignments are mutually exclusive use MutexPeaks
class ResonanceAssignmentRule(ExistingAssignmentRule):
	def __init__(self):
		ExistingAssignmentRule.__init__(self,self.RuleNode(1,None))
	def _iterate(self, peak, atom_tree, frequency_matcher=None, distance_matcher=None, sub_match=None ):
		if not sub_match:
			raise TypeError('This rule should not be used for peaks without pre-existing assignments')
		yield sub_match

