#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from Resonance import UNFOLDED
from Resonance import FoldWindow

#contains shared information about multiple cross-peaks
#i.e., which dimension is which type of atom, which is the label, etc.
class CrossPeakInfo:

	#sub-class for a single Spin-Group
	class SpinInfo:
		def __init__(self):
			self.atom=None
			self.label=None
			self.atom_col=-1
			self.label_col=-1
			self.proton_tolerance=None
			self.label_tolerance=None
			self.proton_folder=UNFOLDED
			self.label_folder=UNFOLDED

		def folder( self, id ):
			if id==1: return self.proton_folder
			if id==2: return self.label_folder

		def __str__(self):
			s='%s (column %d) label: %s (column %d)'%(self.atom, self.atom_col, self.label, self.label_col)
			if self.proton_tolerance:
				s+=' tol: %5.3f'%self.proton_tolerance
			if self.label_tolerance:
				s+=' label_tol: %5.3f'%self.label_tolerance
			return s
#			tol: %8.3f'%(, self.tolerance)

		#find the label atom connected to a given proton
		# i.e., for HN --> N
		# for HE in aa=R --> NE
		def label_atom( self, proton, aa ):
			if not self.label: return None
			if 'N' in self.label:
				if aa=='R':
					if proton=="HE": return "NE"
					if proton[0:2]=="HH": return "N"+proton[1:3]
				if aa=='K' and proton[0:2]=="HZ": return "NZ"
				if aa=='Q' and proton[0:3]=='HE2': return 'NE2'
				if aa=='N' and proton[0:3]=='HD2': return 'ND2'
				if aa=='W' and proton=='HE1': return 'NE1'
				if proton=='H' or proton=='HN': return 'N'
			if 'C' in self.label:
				if proton[0:2]=='QQ': return 'C'+proton[2:]
				if proton[0]=='Q': return 'C'+proton[1:]
				if 'HB' in proton: return 'CB'
				if 'HA' in proton: return 'CA'
				if aa=='W':
					if proton=='HH2': return 'CH2'
					if proton=='HZ2': return 'CZ2'
					if proton=='HZ3': return 'CZ3'
					if proton=='HE3': return 'CE3'
					if proton=='HD1': return 'CD1'
				if aa=='F' or aa=='Y':
					if proton=='HZ': return 'CZ'
					if proton[0:2]=='HD' or proton[0:2]=='HE': return 'C'+proton[1:3]
				if not aa=='N':
					s=len(proton)-2
					if s<1: s=1
					if proton[0:2]=='HG' or proton[0:2]=='HD': return 'C'+proton[1:s+1]
			return None
		# END SPIN GROUP

	#find the label atom connected to a given proton (select spin-group 1,2)
	def label_atom( self, dim, atom, aa ):
		if dim==1:
			return self._spin1.label_atom( atom, aa )
		elif dim==2:
			return self._spin2.label_atom( atom, aa )
		else: assert False, 'dim must be 1 or 2'

	def experiment_id( self ):
		return self._experiment_id

	def __init__(self):
		self._all={}
		self._spin1=None
		self._spin2=None
		self._dim=-1
		self._header_lines=None

	def __str__(self):
		s=''
		s+='Spin 1: %s\n'%self._spin1
		s+='Spin 2: %s\n'%self._spin2
		return s[:-1]

	def write( self, fd ):
		fd.write('# Number of dimensions %d\n'%self.dim() )
		for d in range(0,2):
			fd.write('#INAME %3d %5s\n'%(self.spin(d+1).atom_col,self._spin1.atom))
			if self.spin(d+1).label:
				fd.write('#INAME %3d %5s\n'%(self.spin(d+1).label_col,self._spin1.label))
		fd.write('#TOLERANCE')
		for d in range(0,2):
			fd.write(' %5.3f'%(self.spin(d+1).proton_tolerance))
			if self.spin(d+1).label:
				fd.write(' %5.3f'%(self.spin(d+1).label_tolerance))
		fd.write('\n')

	def dim( self ):
		return self._dim

	def spin( self, dim ):
		if dim==1:
			return self._spin1
		elif dim==2:
			return self._spin2
		else: assert False, 'dim must be 1 or 2'

	@classmethod
	def read_from_lines( obj, lines):
		obj=CrossPeakInfo()
		obj._header_lines=lines
		atom_names={}
		tolerances=None
		fold_info={}
		for l in lines:
			tags=l.lstrip('#').split()
			if len(tags)<1: continue
			obj._all[tags[0]]=" ".join(tags[1:])
			k=tags[0]
			if k=='FILENAME':
				obj._experiment_id=tags[1]
				print 'read peak list %s ...'%tags[1]
			elif k=='INAME':
				d=int(tags[1])
				name=tags[2]
				if name=='H1': name='H'
				if name=='1H': name='HN'
				if name=='C13' or  name=='13C' : name='C'
				if name=='N15' or name=='15N' : name='N'
				atom_names[d]=name
			elif k=='TOLERANCE':
				tolerances=[float(x) for x in tags[1:]]
			elif k=='FOLD':
				fold_info[int(tags[1])]=FoldWindow(float(tags[2]),float(tags[3]))
		dim=len(atom_names)
		obj._dim=dim
		if tolerances: assert len(tolerances)==dim, 'number of values for TOLERANCE entry should be consistent with the dimension'
		else:	tolerances=[None for x in range(0,dim)]
		obj._spin1=CrossPeakInfo.SpinInfo()
		obj._spin2=CrossPeakInfo.SpinInfo()
		for d, n in atom_names.items():
			if n=='h':
#				obj._col2proton[ d ]=2
#				obj._col2islabel[ d ] = False
				obj._spin2.atom='H'
				obj._spin2.atom_col=d
				obj._spin2.proton_tolerance=tolerances[d-1]
				try:
					obj._spin2.proton_folder=fold_info[d]
				except KeyError:
					pass
			elif n=='H':
#				obj._col2proton[ d ]=1
#				obj._col2islabel[ d ] = False
				obj._spin1.atom='H'
				obj._spin1.atom_col=d
				obj._spin1.proton_tolerance=tolerances[d-1]
				try:
					obj._spin1.proton_folder=fold_info[d]
				except KeyError:
					pass

			elif n=='c' or n=='n' or n=='nc' or n=='cn':
#				obj._col2proton[ d ]=2
#				obj._col2islabel[ d ] = True
				obj._spin2.label=n.upper()
				obj._spin2.label_col=d
				obj._spin2.label_tolerance=tolerances[d-1]
				try:
					obj._spin2.label_folder=fold_info[d]
				except KeyError:
					pass

			elif n=='C' or n=='N' or n=='NC' or n=='CN':
#				obj._col2proton[ d ]=1
#				obj._col2islabel[ d ] = True
				obj._spin1.label=n.upper()
				obj._spin1.label_col=d
				obj._spin1.label_tolerance=tolerances[d-1]
				try:
					obj._spin1.label_folder=fold_info[d]
				except KeyError:
					pass

		#if 3D spectrum swap them such that the labelled spin is first.
		if not obj._spin1.label and obj._spin2.label:
			spin3=obj._spin1
			obj._spin1 = obj._spin2
			obj._spin2 = spin3

		return obj



def unit_test():
	s1='''
# Number of dimensions 3
#FILENAME resort_c-ali
#FORMAT xeasy3D
#INAME 1 c
#INAME 2 H
#INAME 3 h
#TOLERANCE      0.3    0.04    0.03
'''

	s2='''
# Number of dimensions 3
#FILENAME resort_c-ali
#FORMAT xeasy3D
#INAME 1 H
#INAME 2 h
#INAME 3 N
#TOLERANCE      0.3    0.04    0.03
'''

	s3='''
# Number of dimensions 3
#FILENAME resort_c-ali
#FORMAT xeasy3D
#INAME 1 nc
#INAME 2 H
#INAME 3 h
#TOLERANCE      0.3    0.04    0.03
'''

	info1=CrossPeakInfo.read_from_lines(s1.split('\n'))
	info2=CrossPeakInfo.read_from_lines(s2.split('\n'))
	info3=CrossPeakInfo.read_from_lines(s3.split('\n'))

	assert info2.label_atom(1, 'HE1','W')=='NE1'
	assert info2.label_atom(1, 'H','W')=='N'
	assert not info1.label_atom(1, 'H','W')
	assert info3.label_atom(1, 'H','X')=='N'  #simNoesy should get the N or C label

	print info1.label_atom(1, 'HB','W')
	print info1.label_atom(1, 'HA','W')

	print info1
	print info2
	print info3

