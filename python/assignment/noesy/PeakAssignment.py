#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


from Atom import Atom
import string

#class that stands for one of multiple assignments in CrossPeaks
# this means a list of 2,  3, 4 atoms whose resonances match (not checked by class) to the CrossPeak frequencies
# use 'assign_atom' to set assigned atoms
# this class is pretty dumb, it only keeps the 2,3,4 atoms assigned to it, and has no logic on its own
# it can keep a selection of weights (also set by others)
# --> make weights a dictionary ?
class PeakAssignment:
	#private data:
	# dim, atoms
	# weights (assignment weights)
	# vc (volume contribution)
	# atom_string (atoms and residues combined in a single string -- why ? )


  def __init__(self, dim):
		self._atoms=[ Atom('CA', 0) for x in range(0,dim) ]
		self._dim=dim
		self._weights=[]
		self._vc=0.0
#		self._atom_string=[]
		#self._store_line

	def __str__(self):
		str=self.atom_str()+'    #VC %6.3f'%self._vc
    if len(self._weights):
       str=str+' #W'
			 for w in self._weights:
				 str=str+' %6.3f'%w
		return str.lstrip()

	def atom_str(self):
		str=''
		for a in self._atoms:
			str=str+' %s'%a
		return str
	#setters ----
	def assign_atom( self, i, atom ):
		self._atoms[i-1]=atom

	def set_weights( self, weights ):
		self._weights=weights

	def set_volume_contribution( self, vc ):
		self._vc=vc;

	def volume_contribution( self ):
		return self._vc

	def weight( self, idx ):
		return self._weights[ idx ]

	def dim(self):
		return self._dim

	#getters ----------
	def atom( self, i ):
		assert( i>0 and i <= self._dim ), "indices are 1, 2, 3, not 0, 1, 2"
		return self._atoms[ i-1 ]

	def has_atom( self, atom ):
		for a in self._atoms:
			if a==atom: return True
		return False

	#get the protons in the assigment
	def protons(self):
		protons=[]
		for a in self._atoms:
			if a.elem()=='H':
				protons.append(a)
		assert (len(protons)==2),"  there should be 2 protons "
		return protons

#yield atoms
	def iter_atoms(self):
		for atom in self._atoms:
			yield atom

	def has_resid( self, resid ):
		for a in self._atoms:
			if a.resid() == resid: return True
		return False

	#return the (two) residues involved in this crosspeak assignment
	def resid_pair( self ):
		pair=sorted(list( set([ a.resid() for a in self._atoms ]) ))
		return tuple(pair)

	@classmethod
	def read_from_str(obj,line,resonances=None):
		tags=line.split()

		#figure out dimension
		i=0

		#figure out whether assignments as indices (CYANA) or named atoms
		named_atoms=True
		try:
			int( tags[0] )
			named_atoms=False
		except:
			pass

		if not named_atoms and int(tags[0])==0:
			return None

		for i,t in enumerate(tags):
			if '#VC' in t or '#QU' in t or '#' in t: break
		if named_atoms:
			dim=i/2
		else:
			dim=i

		weight_start=i+1;

		assert ( dim> 0 and dim<=4 ), "cannot get fix on dimension in line: %s"%line
		obj=PeakAssignment( dim )

		for i in range(0, dim):
			if named_atoms:
				obj.assign_atom( i+1, Atom(tags[i*2],int(tags[i*2+1])) )
			else:
				assert resonances, "cannot read indexed assignments without having the resonance-list present"
				try:
					obj.assign_atom( i+1, resonances[int(tags[i])].atom() )
				except KeyError as exc:
					print "Cannot find resonance with id %d. Assignment %s ignored!"%(int(tags[i]),line), exc
					return None


		if len(tags)<=weight_start: return obj
		obj.set_volume_contribution( float( tags[weight_start] ) )

		if len(tags)<=weight_start+2: return obj
		for t in tags[weight_start+2:]:
			if t[0]=='#': break
			obj._weights.append( float( t ) )

#		obj._form_atom_string()
		return obj

def unit_test():
 	asg=PeakAssignment(3)
 	asg.assign_atom( 1, Atom('N', 15 ) )
 	asg.assign_atom( 2, Atom('H', 15 ) )
 	asg.assign_atom( 3, Atom('HA', 32 ) )
 	asg.set_weights([0.5, 0.3, 0.2, 0.1 ,0.0])
	print asg
	assert not asg.has_atom( Atom('N', 20 ) ), "has_atom: an atom is believed to be contained in the assignment which is not"
	assert asg.has_atom( Atom('HA ', 32 ) ),"has_atom: an atom which is contained in the assignment fails to appear"
 	asg2=PeakAssignment.read_from_str('N 15 H   15   HG   42    #VC 0.940 #W 0.500 0.300 0.200 0.100 0.000')
 	str_asg2=str(asg2)
	print str_asg2
#	assert str_asg2 == 'N   15    H   15   HG   42    #VC 0.940 #W 0.500 0.300 0.200 0.100 0.000','string output has changed'

	s='''
VARS INDEX SHIFT_LOW SHIFT_HIGH SIGMA ATOMNAME RESID RESNAME3 RESNAME
FORMAT %8d %5.3f %5.3f %5.3f %5s %5d %5s %3s
   1    2  4.394      0.040    HA        4 MET M
   2    1  1.905      0.040   HB2        4 MET M
   3    1  1.839      0.040   HB3        4 MET M
	4    1  1.919      0.040    QE        4 MET M
	7    1  2.415      0.040    QG        4 MET M
	9    100 175.642      0.400     C       5 MET M
	10   50  54.943      0.400    CA        5 MET M
	11   10  33.407      0.400    CB        4 MET M
	12   1  16.931      0.400    CE        4 MET M
	13   1  31.943      0.400    CG        4 MET M
'''
	from StringIO import StringIO
	from ResonanceList import ResonanceList
	rl=ResonanceList.read_from_stream(StringIO(s))
 	asg3=PeakAssignment.read_from_str('3 7 9 #VC 0.940 #W 0.500 0.300 0.200 0.100 0.000',rl)
 	asg4=PeakAssignment.read_from_str('1 4 13 #QU 0.940 #SUP 0.200',rl)
	print asg3
	print asg4

