#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module

from collections import defaultdict

#light-weight class that returns the results of matching
class PeakMatch(object):
	__slots__=('rule_match', 'peak_match','peak')
	def __init__(self, rule_match, peak_match, peak ):
		self.rule_match=rule_match
		self.peak_match=peak_match
		self.peak=peak

	def __str__(self):
		return str(self.peak)+': '+self.match_str()

	def match_str(self):
		return '('+', '.join(('%s'%atom for atom in self.peak_match))+')'
	#convenience access to distance measurement for peak match
	def distance(self,distance_matcher):
		return self.peak.rule.distance(self,distance_matcher)

	def distance_atoms(self):
		return self.peak.rule.distance_atoms(self)

	#convenience access to frequencies in peak_match order
	@property
	def freq(self):
		return self.peak._freq

	@property
	def rule(self):
		return self.peak.rule

	def __iter__(self):
		for dim in range(0,len(self.peak_match)):
			yield AtomicPeakMatch(self,dim+1)

	def __key__(self):
		return (self.peak.__key__(),)+tuple([atom.__key__() if atom else None for atom in self.peak_match ])

	def __eq__(self,other):
		if self.peak!=other.peak: return False
		return self.peak_match==other.peak_match

	def __hash__(self):
		return hash(self.__key__())
#use this class to refer to a single atom's assignment resulting from a PeakMatch
class AtomicPeakMatch(object):
	__slots__=('freq','dim','peak_match')
	def __init__(self, peak_match, dim ):
		self.freq=peak_match.peak.freq_of_dim(dim) #direct storage saves ca. 30% of lookup time
		self.dim=dim
		self.peak_match=peak_match

	def get_atom(self):
		return self.peak_match.peak_match[self.dim-1]
	atom=property(get_atom)

	def get_peak(self):
		return self.peak_match.peak
	peak=property(get_peak)

	def __str__(self):
		return 'Atomic( d=%d, w=%5.2f, %s )'%(self.dim,self.freq,self.atom)

	def __key__(self):
		return (self.peak.__key__(), self.dim, self.atom.__key__())

	def __eq__(self,other):
		if self.peak!=other.peak: return False
		if self.dim!=other.dim: return False
		return self.atom==other.atom

	def __hash__(self):
		return hash(self.__key__())

