#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from CrossPeakInfo import CrossPeakInfo
from Atom import Atom
from CrossPeak import CrossPeak
from Resonance import Resonance
from NoeStrip import NoeStrip

class SpinSystem:
	def __init__(self,resid):
		self._resid=resid
		self._strips={}

	def add_strip(self,strip):
		assert strip.resid()==self._resid
		atom=strip.proton_name()
		self._strips[atom]=strip

	def __str__(self):
		s='SpinSystem @ residue %d:\n'%self._resid
		for name, strip in self._strips.iteritems():
			s+='%(name)4s %(strip)30s\n'%locals()
		return s

	def score(self):
#		print self
		pmatch=1
		for name1, strip1 in self._strips.iteritems():
			for name2, strip2 in self._strips.iteritems():
				if name1==name2: continue
				match = strip1.pmatch(strip2)
#				print name1, name2, match
				pmatch*=max(match,1e-10)
		return pmatch




