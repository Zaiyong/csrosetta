#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module


#provide code to convert a ResonanceList
#into a 1D PeakList with ResonanceAssignmentRule
from rules.ResonanceAssignmentRule import ResonanceAssignmentRule
from chemical import Atom
def peak_list_from_resonances( resonance_list ):
	peaks=PeakList('resonance_assignments')
	rule=ResonanceAssignmentRule()
	for resonance in resonance_list:
		assignment=Atom(resonance.name(), resonance.resid())
		id=resonance._id
		if id<0: id=None
		peak=Peak((resonance.freq(),),rule,id)
		peak.hard_assignment([(assignment,)])
		self.append(peak)

