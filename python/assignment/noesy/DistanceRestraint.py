#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from Atom import Atom
import library

class DistanceRestraint:
	def __init__( self, atom1=Atom(), atom2=Atom(), lb=1.5, ub=5.0, sd=1.0, volume=None, peakid=None ):
		self._atom1=atom1
		self._atom2=atom2
		self._lb=lb
		self._ub=ub
		self._sd=sd
		self._volume=volume
		self._peakid=peakid

	def __str__(self):
		s='AmbiguousNMRDistance %(_atom1)s %(_atom2)s BOUNDED %(_lb)8.3f %(_ub)8.3f %(_sd)8.3f NOE'%self.__dict__
		if self._peakid: s+=' Peak %(_peakid)8d'%self.__dict__
		if self._volume: s+=' Volume %(_volume)8.3e'%self.__dict__
		return s

	def __eq__(self,other):
		res1=self._atom1==other._atom1 and self._atom2==other._atom2
		res2=self._atom2==other._atom1 and self._atom1==other._atom2
		return res1 or res2

	def __hash__(self):
		return hash('%s_%s'%(self._atom1,self._atom2))

	def match_fuzzy(self, other):
		res1=self._atom1.match_fuzzy( other._atom1 ) and self._atom2.match_fuzzy( other._atom2 )
		res2=self._atom2.match_fuzzy( other._atom1 ) and self._atom1.match_fuzzy( other._atom2 )
		return res1 or res2

	@classmethod
	def read_from_line(obj,data):
		tags=data.split()
		try:
			if tags[0]!='AmbiguousNMRDistance' and tags[0]!='AtomPair':
				raise library.InconsistentInput("expected AmbiguousNMRDistance or AtomPair as first tag in line '%s'"%data)
			if tags[5]!='BOUNDED':
				raise library.InconsistentInput("expected BOUNDED as function tag in line '%s'"%data)
			atom1=Atom(tags[1],int(tags[2]))
			atom2=Atom(tags[3],int(tags[4]))
			lb=float(tags[6])
			ub=float(tags[7])
			sd=float(tags[8])
			volume=None
			peakid=None
			try: peakid=int(tags[tags.index('Peak')+1])
			except: pass
			try: volume=float(tags[tags.index('Volume')+1])
			except:	pass
			return DistanceRestraint(atom1,atom2,lb,ub,sd,volume,peakid)
		except library.InconsistentInput:
			raise
		except Exception as exc:
			raise library.InconsistentInput("problem reading DistanceRestraint from line '%s'. Caused exception: %s"%(data,exc))



def unit_test():
	s1='''
AmbiguousNMRDistance  HB3    5    H    5 BOUNDED 1.5    0.970    1.000 automatic NOE Peak 1471 Volume 995000.000
AmbiguousNMRDistance   HA    5    H    5 BOUNDED 1.5    5.661    1.000 automatic NOE Peak 1579 Volume 1335000.000
AmbiguousNMRDistance   HA    4    H    5 BOUNDED 1.5  602.246    1.000 automatic NOE Peak 1582 Volume 2906000.000
'''
	s2='''
AmbiguousNMRDistance     H    12   HB3    12 BOUNDED 1.5 4.140 0.3 NOE; rawdata 4.140
AmbiguousNMRDistance     H    13     H    14 BOUNDED 1.5 4.760 0.3 NOE; rawdata 4.760
AmbiguousNMRDistance     H    13    QB    13 BOUNDED 1.5 3.210 0.3 NOE; rawdata 3.210
AmbiguousNMRDistance     H    16    QB    16 BOUNDED 1.5 3.460 0.3 NOE; rawdata 3.460
'''
	csts=[]
	for s in (s1+s2).split('\n'):
		if len(s.strip()):
			csts.append( DistanceRestraint.read_from_line(s) )

	print s1+s2
	for cst in csts:
		print cst
