#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
import cs
from Atom import Atom
from Resonance import *
import amino_acids

class ResonanceList(dict):
	#private data members:
   #
	# the primary container:
	# 1) a dictionary of resonances id - > resonance
	#
	# a derived lookup table
	# 2) resid -- > list of resonances
	#
	# 3) the (fasta) sequence
	def __init__(self):
		dict.__init__(self)
		self._residue_map = {}
		self._sequence = ''

#	def __len__(self):
#		return len(self._resonances)
	def __str__(self):
		s=''
		for r in self.itervalues():
			line='%5d %s %s   %s'%(r.id(),r._freq_str(),r._tail_str(),r.atom())
			if len(self._sequence)>(r.resid()-1):
				aa = self._sequence[r.resid()-1];
				aa3=amino_acids.short_to_long[aa]
				line+='  %3s %1s'%(aa3,aa)
				if r._intensity:
					line+='%3d'%r._intensity
				if r._floats:
					line+=' ['+','.join(['%d'%i for i in r._floats])+']'
			s+=line+'\n'
		return s

	def _generate_residue_map(self):
		self._residue_map = {}
		for r in self.itervalues():
			resid=r.resid()
			self._residue_map.setdefault( resid, []).append( r )

	#returns a dictionary with the two keys 'H','L' for proton and label resonance
	def coupled_atoms_resonances(self,proton_atom,label_element):
		import util
		resid=proton_atom.resid()
		proton_name=proton_atom.name()
		aa=self.sequence()[resid-1]
		label_name=util.label_atom(label_element,proton_name,aa)
		#print label_element,label_name,proton_name
		proton_resonance=self.by_atom(proton_atom)
		try:
			label_atom=Atom(label_name,resid)
			label_resonance=self.by_atom(label_atom)
		except:
			label_resonance=None
		coupled_res={'H':proton_resonance,'L':label_resonance}
		return coupled_res

	def by_residue(self, resid ):
		return self._residue_map[ resid ]

	def by_atom( self, atom ):
		for r in self._residue_map[ atom.resid() ]:
			if r.name()==atom.name(): return r
		raise KeyError( "no resonance known for atom %s"%atom )

	#we have changed how atom attributes are accessed with chemical.Atom
	def by_chemical_atom( self, atom ):
		for r in self._residue_map[ atom.resid ]:
			if r.name()==atom.name: return r
		raise KeyError

	#we have changed how atom attributes are accessed with chemical.Atom
	def has_chemical_atom( self, atom ):
		for r in self._residue_map[ atom.resid ]:
			if r.name()==atom.name: return True
		return False

	def set_by_residue(self,resid,resid_res):
		self._residue_map[resid]=resid_res
		return self

	def set_by_atom(self,atom,res):
		for r in self._residue_map[ atom.resid() ]:
			if r.name()==atom.name():
				r.set_freq(res.freq())
		return self

	def delete_by_residue(self,resid):
		for r in self.by_residue( resid ):
			del self[r.id()]
		del self._residue_map[resid]


	def delete_by_atom(self,atom):
		remove_list=[]
		for i,r in self.iteritems():
			if r.name()==atom.name() and r.resid()==atom.resid():
				remove_list.append(i)
		for i in remove_list:
			del self[i]

	def iter_by_atom_type( self, atom_type ):
		for r in self.itervalues():
			if r.name()==atom_type:
				yield r

	def iter_residues( self ):
		for i,r in self._residue_map.iteritems():
			yield i,r

	def matches( self, atom_type, freq, tol, folder=UNFOLDED, threshold=1 ):
		for r in self.iter_by_atom_type( atom_type ):
			if r.match( freq, tol, folder, threshold ):
				yield r

	def set_sequence( self, sequence ):
		self._sequence=sequence

	def sequence(self) :
		return self._sequence

	def add_resonance( self,newr ):
		new_id = newr.id()
		while new_id in self or new_id <= 0:
			new_id+=1
		newr.set_id(new_id)
		self[new_id]=newr
		self._residue_map.setdefault( newr.resid(), [] ).append( newr )
		return newr

	def remove_resonance( self, oldr ):
		del self[oldr.id()]
		self._residue_map[oldr.resid()].remove(oldr)

	def aa_from_atom( self, atom ):
		return self._sequence[ atom.resid()-1 ]

	def aa_from_resonance( self, resonance ):
		return self._sequence[ resonance.resid()-1 ]

	@classmethod
	def read_from_prot(obj, prot):
		obj=ResonanceList()
		obj._sequence = prot.sequence
		has_range = 'SHIFT_LOW' in prot
		ind_col=prot.get_col('INDEX')
		ind_err=prot.get_col('SIGMA')
		has_intensity = 'INTENSITY' in prot
		has_floats = 'FLOATS' in prot
		if has_range:
			ind_shift=prot.get_col('SHIFT_LOW')
		else:
			ind_shift=prot.get_col('SHIFT')
		if has_intensity:
			ind_intens=prot.get_col('INTENSITY')

		if has_floats:
			ind_floats=prot.get_col('FLOATS')
		has_ambiguity = 'AMBIGUITY' in prot
		if has_ambiguity:
			ind_ambiguity=prot.get_col('AMBIGUITY')

		for i, val in prot.iteritems():
			atom=Atom(i[1], i[0])
			if has_range and val[ind_shift+1]-val[ind_shift]>0.01:
				res=RangeResonance( val[ind_col], atom, val[ind_shift], val[ind_shift+1], val[ind_err] )
			else:
				res=Resonance( val[ind_col], atom, val[ind_shift], val[ind_err] )
				if has_ambiguity:
					res.ambiguity = val[ind_ambiguity]

			if has_floats and len(val)>=ind_floats+1:
				res.set_floats( val[ind_floats], obj )

			if has_intensity:
				res.set_intensity(val[ind_intens])

			obj[ val[ind_col] ]=res

		obj._generate_residue_map()
		return obj

	#generate a dictionary with all the information as used in NIH_Table.from_dict()
	#RESID--> list of residue numbers
	#ATOMNMAE --> list of atom-names
	#RESNAME --> list of residue names
	# all lists have same length
	def generate_dict( self ):
		resids=[]
		atoms=[]
		shifts=[]
		ids=[]
		err=[]
		resnames=[]
		resnames3=[]
		for r in self.itervalues():
			try:
				shifts.append( r.freq() )
				resids.append( r.resid() )
				atoms.append( r.name() )
				ids.append( r.id() )
				err.append( r.error() )
				try:
					aa=self._sequence[ r.resid()-1 ]
					resnames.append( aa )
					resnames3.append( amino_acids.short_to_long[aa] )
				except KeyError:
					print 'cannot translate %s'%aa
					pass
			except:
				print 'cannot use %s'%r
		data={}
		data['RESID']=resids
		data['INDEX']=ids
		data['ATOMNAME']=atoms
		data['SHIFT']=shifts
		data['SIGMA']=err
		if len(resnames)==len(resids):
			data['RESNAME']=resnames
		if len(resnames3)==len(resids):
			data['RESNAME3']=resnames3
		return data

	def prot_file(self):
		import cs
		prot_data=self.generate_dict()
		prot_file = cs.ProtCSFile().from_table( cs.NIH_table().from_dict( prot_data ) )
		return prot_file

	def talos_file(self):
		import cs
		prot_data=self.generate_dict()
		talos_file = cs.TalosCSFile().from_table( cs.NIH_table().from_dict( prot_data ) )
		return talos_file

	@classmethod
	def read_from_stream(obj, file):
		prot = cs.ProtCSFile()
		prot.read_stream( file, None )
		return ResonanceList.read_from_prot(prot)

def unit_test():
	s1='''
   1      4.394      0.040    HA        4 MET M
   2      1.905      0.040   HB2        4 MET M
   3      1.839      0.040   HB3        4 MET M
	4      1.919      0.040    QE        4 MET M
	7      2.415      0.040    QG        4 MET M
	9    175.642      0.400     C        5 MET M
	10     54.943      0.400    CA        5 MET M
	11     33.407      0.400    CB        4 MET M
	12     16.931      0.400    CE        4 MET M
	13     31.943      0.400    CG        4 MET M
'''

	s2='''
VARS INDEX SHIFT_LOW SHIFT_HIGH SIGMA ATOMNAME RESID RESNAME3 RESNAME
FORMAT %8d %5.3f %5.3f %5.3f %5s %5d %5s %3s
   1    2  4.394      0.040    HA        4 MET M
   2    1  1.905      0.040   HB2        4 MET M
   3    1  1.839      0.040   HB3        4 MET M
	4    1  1.919      0.040    QE        4 MET M
	7    1  2.415      0.040    QG        4 MET M
   8    175.642 -1   0.400     CG       5 MET M
	9    100 175.642   0.400     C        5 MET M
	10   50  54.943      0.400    CA        5 MET M
	11   10  33.407      0.400    CB        4 MET M
	12   1  16.931      0.400    CE        4 MET M
	13   1  31.943      0.400    CG        4 MET M
'''
	print 'START UNIT-TEST: ResonanceList'
	from StringIO import StringIO
	rl1=ResonanceList.read_from_stream(StringIO(s1))
	rl2=ResonanceList.read_from_stream(StringIO(s2))
	print 'dump list 1...'
	for r in rl1.itervalues():
		 print r

	print 'dump list 2...'
	for r in rl2.itervalues():
		print r

	print 'only output residue 5'
	for r in rl2.by_residue( 5 ):
		print r
	pass
	print 'END UNIT-TEST: ResonanceList'
