 #!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from SpinSystem import SpinSystem
from Atom import Atom
import math
class resid_cs_valid:
	def __init__(self,resid):
		self._resid=resid
		self._spinsystems=[]
		self._true_protons={}
		self._true_heavies={}
	def amount_ss(self):
		return len(self._spinsystems)

	def add_ss(self,spinsystem):
		assert spinsystem.resid()==self._resid, "residue are not matched\n"
		self._spinsystems.append(spinsystem)

	def add_true_protons(self,aa,resonancelist):
		#assert resonance.atom().resid==self._resid, "residue are not matched\n"
		#name=resonance.atom().name()
		#self._true_protons[name]=resonancelist.by_atom(Atom(name,self._resid))
		self._true_protons['N']=[resonancelist.by_atom(Atom('H',self._resid))]
		self._true_protons['CA']=[resonancelist.by_atom(Atom('HA',self._resid))]
		self._true_protons['CB']=[resonancelist.by_atom(Atom('HB',self._resid))]
		if self._true_protons['CB']==[None]:
			self._true_protons['CB']=[resonancelist.by_atom(Atom('HB2',self._resid)),resonancelist.by_atom(Atom('HB3',self._resid))]
		if self._true_protons['CB']==[None] or self._true_protons['CB']==[None,None]:
			self._true_protons['CB']=[resonancelist.by_atom(Atom('QB',self._resid))]
		if aa in 'RQEKMP':
			self._true_protons['CG']=[resonancelist.by_atom(Atom('HG2',self._resid)),resonancelist.by_atom(Atom('HG3',self._resid))]
			if self._true_protons['CG']==[None,None]:
				self._true_protons['CG']=[resonancelist.by_atom(Atom('QG',self._resid))]
		elif aa=='I':
			self._true_protons['CG1']=[resonancelist.by_atom(Atom('HG12',self._resid)),resonancelist.by_atom(Atom('HG13',self._resid))]
			if self._true_protons['CG1']==[None,None]:
				self._true_protons['CG1']=[resonancelist.by_atom(Atom('QG1',self._resid))]
			self._true_protons['CG2']=[resonancelist.by_atom(Atom('HG21',self._resid)),resonancelist.by_atom(Atom('HG22',self._resid)),resonancelist.by_atom(Atom('HG23',self._resid))]
			if self._true_protons['CG2']==[None,None,None]:
				self._true_protons['CG2']=[resonancelist.by_atom(Atom('QG2',self._resid))]
		elif aa=='L':
			self._true_protons['CG']=[resonancelist.by_atom(Atom('HG',self._resid))]
		elif aa=='T':
			self._true_protons['CG2']=[resonancelist.by_atom(Atom('HG21',self._resid)),resonancelist.by_atom(Atom('HG22',self._resid)),resonancelist.by_atom(Atom('HG23',self._resid))]
			if self._true_protons['CG2']==[None,None,None]:
				self._true_protons['CG2']=[resonancelist.by_atom(Atom('QG2',self._resid))]
		elif aa=='V':
			self._true_protons['CG1']=[resonancelist.by_atom(Atom('HG11',self._resid)),resonancelist.by_atom(Atom('HG12',self._resid)),resonancelist.by_atom(Atom('HG13',self._resid))]
			self._true_protons['CG2']=[resonancelist.by_atom(Atom('HG21',self._resid)),resonancelist.by_atom(Atom('HG22',self._resid)),resonancelist.by_atom(Atom('HG23',self._resid))]
			if self._true_protons['CG1']==[None,None,None]:
				self._true_protons['CG1']=[resonancelist.by_atom(Atom('QG1',self._resid))]
			if self._true_protons['CG2']==[None,None,None]:
				self._true_protons['CG2']=[resonancelist.by_atom(Atom('QG2',self._resid))]
		if aa in 'RKP':
			self._true_protons['CD']=[resonancelist.by_atom(Atom('HD2',self._resid)),resonancelist.by_atom(Atom('HD3',self._resid))]
			if self._true_protons['CD']==[None,None]:
				self._true_protons['CD']=[resonancelist.by_atom(Atom('QD',self._resid))]
		elif aa in 'FY':
			self._true_protons['CD1']=[resonancelist.by_atom(Atom('HD1',self._resid))]
			self._true_protons['CD2']=[resonancelist.by_atom(Atom('HD2',self._resid))]
			if self._true_protons['CD1']==[None]:
				self._true_protons['CD1']=[resonancelist.by_atom(Atom('QD',self._resid))]
			if self._true_protons['CD2']==[None]:
				self._true_protons['CD2']=[resonancelist.by_atom(Atom('QD',self._resid))]
		elif aa=='W':
			self._true_protons['CD1']=[resonancelist.by_atom(Atom('HD1',self._resid))]
		elif aa=='I':
			self._true_protons['CD1']=[resonancelist.by_atom(Atom('HD11',self._resid)),resonancelist.by_atom(Atom('HD12',self._resid)),resonancelist.by_atom(Atom('HD13',self._resid))]
			if self._true_protons['CD1']==[None,None,None]:
				self._true_protons['CD1']=[resonancelist.by_atom(Atom('QD1',self._resid))]
		elif aa=='L':
			self._true_protons['CD1']=[resonancelist.by_atom(Atom('HD11',self._resid)),resonancelist.by_atom(Atom('HD12',self._resid)),resonancelist.by_atom(Atom('HD13',self._resid))]
			self._true_protons['CD2']=[resonancelist.by_atom(Atom('HD21',self._resid)),resonancelist.by_atom(Atom('HD22',self._resid)),resonancelist.by_atom(Atom('HD23',self._resid))]
			if self._true_protons['CD1']==[None,None,None]:
				self._true_protons['CD1']=[resonancelist.by_atom(Atom('QD1',self._resid))]
			if self._true_protons['CD2']==[None,None,None]:
				self._true_protons['CD2']=[resonancelist.by_atom(Atom('QD2',self._resid))]
		elif aa=='H':
			self._true_protons['CD2']=[resonancelist.by_atom(Atom('HD2',self._resid))]
	def add_true_heavies(self,aa,resonancelist):
		#assert resonance.atom().resid==self._resid, "residue are not matched\n"
		#name=resonance.atom().name()
		#self._true_heavies[name]=resonancelist.by_atom(Atom(name,self._resid))
		self._true_heavies['N']=resonancelist.by_atom(Atom('N',self._resid))
		self._true_heavies['CA']=resonancelist.by_atom(Atom('CA',self._resid))
		self._true_heavies['CB']=resonancelist.by_atom(Atom('CB',self._resid))
		if aa in 'RQEKMP':
			self._true_heavies['CG']=resonancelist.by_atom(Atom('CG',self._resid))
		elif aa=='I':
			self._true_heavies['CG1']=resonancelist.by_atom(Atom('CG1',self._resid))
			self._true_heavies['CG2']=resonancelist.by_atom(Atom('CG2',self._resid))
		elif aa=='L':
			self._true_heavies['CG']=resonancelist.by_atom(Atom('CG',self._resid))
		elif aa=='T':
			self._true_heavies['CG2']=resonancelist.by_atom(Atom('CG2',self._resid))
		elif aa=='V':
			self._true_heavies['CG1']=resonancelist.by_atom(Atom('CG1',self._resid))
			self._true_heavies['CG2']=resonancelist.by_atom(Atom('CG2',self._resid))
		if aa in 'RKP':
			self._true_heavies['CD']=resonancelist.by_atom(Atom('CD',self._resid))
		elif aa in 'FY':
			self._true_heavies['CD1']=resonancelist.by_atom(Atom('CD1',self._resid))
			self._true_heavies['CD2']=resonancelist.by_atom(Atom('CD2',self._resid))
		elif aa=='W':
			self._true_heavies['CD1']=resonancelist.by_atom(Atom('CD1',self._resid))
		elif aa=='I':
			self._true_heavies['CD1']=resonancelist.by_atom(Atom('CD1',self._resid))
		elif aa=='L':
			self._true_heavies['CD1']=resonancelist.by_atom(Atom('CD1',self._resid))
			self._true_heavies['CD2']=resonancelist.by_atom(Atom('CD2',self._resid))
		elif aa=='H':
			self._true_heavies['CD2']=resonancelist.by_atom(Atom('CD2',self._resid))

	def resid(self):
		return self._resid

	def iter_spinsystems(self):
		for ss in self._spinsystems:
			yield ss

	def true_protons(self):
		return self._true_protons

	def true_heavies(self):
		return self._true_heavies

	def true_proton(self,name):
		try:
			return self._true_protons[name]
		except:
			assert False, "%s is not in this true protons list\n"%name


	def __str__(self):
		s='@ %d residue: \n'%self._resid
		s+='the true cs are:\n'
		for name,resonance in self._true_protons.iteritems():
			s+='%(name)4s %(freq)8.3f\n'%{'name':name,'freq':resonance.freq()}
		s+='it has %5d spinsystems\n'%self.amount_ss()
		return s

