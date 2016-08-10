#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#################
##  class Conformation
##
##  stores multiple models of a protein conformation
##  coords are stored as list
##  if a particular frame does not contain a particular coordinate (e.g., changing protonation states HIS-HD vs. HIS-HE )
##  the particular number is stored as None
##
##  the coords are stored in a dictionary that uses chemical.Atom as key
##
import numpy
from numpy.linalg import norm
import chemical
from basic import Tracer
import library

tr = Tracer('conformation')

#little helper function -- allows to correct some minor problems with PDB atom names
#
#body of def from_conformation_loader()
#
def get_molecule_atom( atom, molecule ):
	while True:
		try:
			atom = molecule.atom( atom.methyl_atom() )
			return atom
		except KeyError:
			#correct N-terminus problem: if H1 at N-term make it H
			if atom.resid==molecule.first_residue and atom.name=='H1':
				atom.name='H' ## fix atom in place
			elif atom.resid==molecule.last_residue and atom.name=='OXT':
				atom.name='O'
			else: raise

class ConformationLoader:
	def __init__(self):
		self.frames=[]
		self.sequence=None

	@classmethod
	def from_pdb_file( obj, pdb_file, chains=['_','A'], molecule=None ):
		model=0
		coords={}
		obj=ConformationLoader()
		for line in open(pdb_file,'r'):
			tags=line.split()
			if tags[0]=='MODEL':
				model+=1
				if len(coords)>0:
					obj.frames.append(coords)
					coords={}
			elif tags[0]=='ATOM':
				name=tags[2]
				resname3=tags[3]
				try:
					resid=int(tags[5])
					x=(float(tags[6]),float(tags[7]),float(tags[8]))
					chain=tags[4]
					try:
						elem=tags[11]
					except IndexError:
						elem=name[0]
				except ValueError: #no chainID
					resid=int(tags[4])
					x=(float(tags[5]),float(tags[6]),float(tags[7]))
					try:
						elem=tags[10]
					except IndexError:
						elem=name[0]
					chain='_'
				if chain not in chains: continue
				pdb_atom=chemical.Atom(name,resid,res_type=resname3,elem=elem)
				if molecule:
					real_atom = get_molecule_atom( pdb_atom, molecule )# molecule.atom(pdb_atom.methyl_atom())
					coords[real_atom]=x
				else:
					coords[pdb_atom]=x
		if len(coords)>0:
			obj.frames.append(coords)
		if molecule:
			obj.sequence=molecule.sequence
			obj.start_residue=molecule.first_residue
		else:
			obj._figure_out_sequence()
		return obj

  def _figure_out_sequence(self):
		coords=self.frames[0]
		min_res=min( atom.resid for atom in coords.iterkeys() )
		max_res=max( atom.resid for atom in coords.iterkeys() )
		self.sequence=['X']*(max_res-min_res+1)
		from amino_acids import long2short
		for atom in coords.iterkeys():
			if not atom.name=='CA': continue
			self.sequence[atom.resid-min_res]=long2short(atom.res3_type)
		self.sequence=''.join(self.sequence)
		self.start_residue=min_res

	def set_of_atoms(self):
		atoms=set()
		for frame in self.frames:
			atoms=atoms.union(set(atom for atom in frame.iterkeys() ) )
		return atoms


class Conformation:
	def __init__(self,molecule=None):
		self._coords={}
		self.molecule=molecule

	def dist(self,atom1,atom2):
		pass

	@classmethod
	def from_pdb_file( obj, pdb, chains=['_','A'], molecule=None ):
		loader=ConformationLoader.from_pdb_file( pdb, chains, molecule )
		obj=Conformation.from_conformation_loader( loader, molecule )
		return obj

	@classmethod
	def from_conformation_loader( obj, loader, chains, molecule=None ):
		tr_loading=Tracer('loader',tr)
		if not molecule:		#make molecule from the loaded PDB sequence
			if 'X' in loader.sequence:
				exc=library.MissingInput('cannot read PDB file with missing density without a pre-defined molecule')
				exc.sequence=loader.sequence
				raise exc
			molecule=chemical.AtomTree.from_sequence( loader.sequence )
			molecule.offset_residue_numbers( loader.start_residue-1 )
		#get empty conformation object
		obj=Conformation(molecule)
		obj._coords={}
		# 1) initialize coords with [None,None,None,...] for every atom
		# 2) enumerate frames and store coord of atom at appropriate position: coords[atom]-->[..., x, ... ] at position i
		#initialize coord lists
		for atom in  loader.set_of_atoms():
			try:
				atom = get_molecule_atom( atom, molecule )
				obj._coords[ atom ] = [ None ] * len(loader.frames)
			except KeyError:
				tr_loading.Warning('Atom %s is not recognized and will be ignored.'%atom)
		#store coords
		for ct, frame in enumerate(loader.frames):
			for atom, x in frame.iteritems():
				try:
					atom = get_molecule_atom( atom, molecule )
					assert atom in obj._coords
					obj._coords[ atom ][ ct ] = x
				except KeyError:
					#already warned in previous loop
					pass #silently
		#now transfer to numpy
		for key,value in obj._coords.iteritems():
			obj._coords[ key ] = numpy.array( value )
		return obj

	def coords( self, atom ):
		return self._coords[ atom ].tolist()

	def numpy_coords( self, atom ):
		return self._coords[ atom ]

	def dist( self, atom1, atom2 ):
		x1=self._coords[ atom1 ]
		x2=self._coords[ atom2 ]
		return numpy.array( map( norm, x1-x2 ) )

