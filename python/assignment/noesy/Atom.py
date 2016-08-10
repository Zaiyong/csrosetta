#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#-------------------------------------------
class Atom:
	#private data members
	#name, resid

	def __init__(self, name='X', resid=-1 ):
			self.set(name,resid)

	def __str__(self):
		return '%-4s %4d'%(self._name, self._resid)

	def set(self, name, resid ):
		self._name=name.strip()
		self._resid=resid
		return self

	def __eq__(self, other):
		if self._resid != other._resid: return False
		return self._name == other._name

	def match_fuzzy( self, other ):
		if self._resid != other._resid: return False
		if self._name[0:2]=='QQ' or	other._name[0:2]=='QQ':
			return self.methyl_name(2)==other.methyl_name(2)
		if self._name[0]=='Q' or other._name[0]=='Q':
			return self.methyl_name(1)==other.methyl_name(1)
		return self==other
	def name(self):
		return self._name

	#for a HB2 return QB, ...
	def methyl_name(self,depth=1):
		name=self.name()
		if self.elem()!='H' or depth==0: return name
		if name=='H': return name
		if name[0:depth]=='Q'*depth: return name
		if depth>1:
			s=self.methyl_name(depth-1)
		else:
			s=name[1:]
		return 'Q'+s[0:-1]

	def resid(self):
		return self._resid

	def set_resid(self,resid):
		self._resid=resid

	def set_name(self,name):
		self._name=name

	def element(self):
		return self.elem()

	def elem(self):
		if self._name[0]=='Q': return 'H'
		if self._name[0]=='C': return 'C'
		if self._name[0]=='N': return 'N'
		if self._name[0]=='H': return 'H'
		if self._name[0]=='O': return 'O'
		if self._name[0]=='S': return 'S'
		assert False, 'need to improve this implementation'

	def __hash__(self):
		return hash('%s'%self)


def unit_test():
	a1=Atom('HB3',1)
	a2=Atom('QB',1)
	a3=Atom('HD12',1)
	a4=Atom('QQD',1)
	a5=Atom('H',1)

	assert a1.methyl_name()=='QB'
	assert a3.methyl_name()=='QD1'
	assert a3.methyl_name(2)=='QQD', 'wrong result: %s'%a3.methyl_name(2)
	assert a4.methyl_name(1)=='QQD'
	assert a4.methyl_name(2)=='QQD'
	assert a5.methyl_name()=='H'

	assert a1.match_fuzzy(a2)
	assert a3.match_fuzzy(a4)
	assert a3.match_fuzzy(a3)
	assert a2.match_fuzzy(a2)
	assert a5.match_fuzzy(a5)
	assert not a1.match_fuzzy(a4)
	assert not a5.match_fuzzy(a1)
	assert not a5.match_fuzzy(a2)

