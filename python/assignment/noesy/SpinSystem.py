#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from CrossPeakInfo import CrossPeakInfo
from Atom import Atom
from CrossPeak import CrossPeak
from Resonance import Resonance
from NoeStrip import NoeStrip
import math
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

	def score(self,type,proton_name):
 		namelist=[]
		pmatch=1
		#match={}
		match=[]
		for name1, strip1 in self._strips.iteritems():
			namelist.append(name1)
			for name2, strip2 in self._strips.iteritems():
				if name1==name2: continue
				if proton_name not in [name1,name2]: continue
				#match=strip1.pmatch(strip2)
				#print match
				#match[name1,name2]=strip1.pmatch(strip2)
				match.append(strip1.pmatch(strip2))
		score=math.sqrt(math.fsum([r*r for r in match])/len(match))
		return score
		#print match
# 		for r in namelist:
# 			for g in namelist:
# 				if r==g: continue
# 				if proton_name not in [r,g]: continue
# 				pmatch=pmatch*math.sqrt((match[r,g]*match[r,g]+match[g,r]*match[g,r])/2)
# 		return math.sqrt(pmatch)




