#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#class Atom as used in the AtomTree
import itertools
import amino_acids

class Atom:
	#assign a unique ID to every Atom that is generated.
	__slots__=('element','name','resid','res3_type','_id','_hash','type')
	def __init__(self, name='X', resid=-1, elem=None, res_type=None, types=None ):
		self.set(name,resid)
		self.element=elem
		if types:
			self.types=frozenset(types)
		else:
			self.types=None
		try:
			#store_restype
			if res_type:
				self.res3_type=amino_acids.short2long(res_type)
			else:
				self.res3_type=None
		except KeyError:
			self.res3_type=res_type
		self._hash=None
		# we compute a uniqe ID based on ASCII_VALUE of atom name
		# the maximum length of atom-names is 5 characters
		# 0..255+(0..255)*256+(0..255)*256^2+...not a problem in 64 bit long-int.
	def __str__(self):
		format='%(name)4s %(resid)4d'
#		if self.element: format+=' %(element)1s'
#		if self.res3_type: format+=' %(res3_type)3s'
		return format%self.__dict__

	def long_str(self):
		format='%(name)4s %(resid)4d'
		if self.element: format+=' %(element)1s'
		if self.res3_type: format+=' %(res3_type)3s'
		return format%self.__dict__

	def __repr__(self):
		return self.__str__()

	def __cmp__(self, other ):
		if other == None:
			return cmp( 1, None )
		assert isinstance( other, Atom )
		return cmp(self._id, other._id)

	def __key__(self):
		return self._id

	def __hash__(self):
#		if not self._hash:
#			self._hash=
		return hash(self._id)
#self._hash

	def set(self, name, resid ):
		self.name=name.strip()
		self.resid=resid
		self._compute_id()
		return self

	def _compute_id(self):
		assert len(self.name)<=5, 'Atom-Names with more than 5 characters are not supported due to ID encoding scheme'
		self._id=0
		factor=1
		for c in self.name:
			self._id+=ord(c)*factor
			factor*=256
		self._id+=self.resid*(256**5)
		self._hash=None

	def match_fuzzy( self, other ):
		if self.resid != other.resid: return False
		if self.name[0:2]=='QQ' or	other_name[0:2]=='QQ':
			return self.methyl_name(2)==other.methyl_name(2)
		if self.name[0]=='Q' or other.name[0]=='Q':
			return self.methyl_name(1)==other.methyl_name(1)
		return self==other

		#for a HB2 return QB, ...
	def methyl_name(self,depth=1):
		if self.res3_type:
			if self.res3_type=='ILE':
				if self.name[0:3]=='HD1':
					return 'QD1'
				if self.name[0:3]=='HG2':
					return 'QG2'
			elif self.res3_type=='LEU':
				if self.name[0:2]=='HD':
					return 'QD'+self.name[2]
			elif self.res3_type=='VAL':
				if self.name[0:2]=='HG':
					return 'QG'+self.name[2]
			elif self.res3_type=='ALA':
				if self.name[0:2]=='HB':
					return 'QB'
			elif self.res3_type=='MET':
				if self.name[0:2]=="HE":
					return 'QE'
			elif self.res3_type=='THR':
				if self.name[0:3]=='HG2':
					return 'QG2'
			elif self.res3_type=='LYS':
				if self.name[0:2]=='HZ':
					return 'QZ'
			return self.name
		name=self.name
		if self.element!='H' or depth==0: return name
		if name=='H': return name
		if name=='HA': return name
		if name[0:depth]=='Q'*depth: return name
		if depth>1:
			s=self.methyl_name(depth-1)
		else:
			s=name[1:]
		return 'Q'+s[0:-1]

	def methyl_atom( self ):
		return Atom( self.methyl_name(), self.resid, res_type=self.res3_type, elem=self.element, types=self.types )

	def set_resid(self,resid):
		self.resid=resid
		self._compute_id()

	def set_name(self,name):
		self.name=name
		self._compute_id()

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

