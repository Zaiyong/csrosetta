#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module

###############
### class Peak
# a peak consists of a collection of frequencies and a Rule (replacing former CrossPeakInfo)
# the peak itself is pretty dumb, the Rule knows how to interprete peaks and how to match peaks
# given covalent structure, frequencies and distances
#
# however, peaks can be pre-assigned (reflecting manual inputs)
from basic import Tracer
tr=Tracer('assignment.Peak')
class Peak(object):
	__slots__=('dim','freq','hard_assignments','rule','id','peak_list_name')
	def __init__( self, freq, rule, id=None, pl_name='' ):
		self.dim=len(freq)
		self.freq=tuple(freq)
		self.hard_assignments=[]
		self.rule=rule
		self.id=id
		self.peak_list_name=pl_name

	def freq_of_dim(self, dim):
		return self.freq[dim-1]

	def matches(self, protein, random_samples=None, frequency_matcher=None, distance_matcher=None, partial_match=None, match_mask=None ):
		tr.Trace('Peak::matches')
		for match in self.rule.matches( self, protein, random_samples, frequency_matcher, distance_matcher, partial_match, match_mask ):
			yield match

	def __str__(self):
		str=''
		if self.id:
			str+='%10s'%self.peak_list_name
			try:
				str+='%5d '%self.id
			except TypeError:
				str+='%5s '%self.id
		try:
			for freq in self.freq:
				str+='%8.3f '%freq
		except TypeError:
			pass
		return str[:-1]

	def __key__(self):
		#we have no fool-proof mechanism to make sure that ID is unique. use frequencies as key, too
		return (self.peak_list_name,self.id)#,self.freq)

	def __cmp__(self,other):
		if self.peak_list_name != other.peak_list_name: return cmp(self.peak_list_name,other.peak_list_name)
		return cmp(self.id,other.id)
#		return cmp( (self.peak_list_name, self.id, self.freq), (other.peak_list_name, other.id, other.freq ) )
#		return cmp( (self.peak_list_name, self.id), (other.peak_list_name, other.id ) )
#


	def __hash__(self):
		return hash(self.__key__())

	@classmethod
	def from_string( obj, s ):
		tags=s.split()
		pl_name=tags[0]
		id=int(tags[1])
		freq=[]
		for t in tags[2:]:
			freq.append(float(t))
		obj=Peak(freq,None,id,pl_name)
		return obj

###########################
#
#   MutualExclusivePeak
#
#  this class is used to reflect a (small) group of peaks with pre-existing assignments (usual identical)
#  where each assignment can exactly be made once (mutual exclusive between peaks).
#  example: we know w1 and w2 are the CA resonances of either HIS 73 or HIS 26
#     create two 1D Peak with w1 and w2 as frequency and both with CA 73 and CA 26 as pre-existing
#     assignments. Then create a MutualExclusivePeak from them.
#
#  The MEP will automatically generate a fitting MutualExclusiveRule from the Rules of the input peaks
#  input rules are asserted to be identical.
#
class MutualExclusivePeak(Peak):
	__slots__=('_peaks')
	def __init__(self,peaks):
		#check input
		if len(peaks)<1:
			raise TypeError('At least one peak must be contained in a MutualExclusivePeak, only 2 and more makes real sense...')
		#assign members
		self._peaks=tuple(peaks)
			#unpack frequencies
		self.freq=tuple(w for peak in self._peaks for w in peak.freq)
		self.dim=len(self.freq)
		#double check:
		assert self._dim == sum(peak._dim for peak in self._peaks)
		#make MutexRule
		from rules.BasicRules import MutualExclusiveRule
		dimensionalities=tuple(peak.dim for peak in self._peaks)
		self.rule=MutualExclusiveRule(peaks[0].rule,dimensionalities)
		for p in peaks:
			assert id(peaks[0].rule) == id(p.rule)

		#now we have all data together to initialize the Parent class
		Peak.__init__(self,self.freq,self._rule)

    #flatten the assignments
		self.hard_assignments=[x for x in self.enumerate_assignments(self._peaks)] # for x in group]

		print 'mutex: '
		for ass in self.hard_assignments:
			print ass

	#helper method to build the enumerated list of assignments
  # out of two assignments [A,B] for peak 1
	# this generates two assignments [A B],[B A] for the mex peak.
			# out of [A, B, C] of a triple peak-group we would get
			# [A B C],[A C B],[B A C],[B C A], [C A B], [C B A]
	@staticmethod
	def enumerate_assignments(peaks,excludes=[]):
		if len(peaks)==0:
			yield ()
			return #and finish recursion here
		head=peaks[0]
		tail=peaks[1:]
		for assign in head.hard_assignments:
			if assign in excludes: continue
			for x in MutualExclusivePeak.enumerate_assignments(tail,excludes+[assign]):
				yield assign+x


