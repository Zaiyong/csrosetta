#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


#author: Oliver Lange

######
#  class AtomTree
#   -define a single connected graph that defines the molecular system
#   -edges reflect chemical bonds or jumps (to other molecules)
#
#
### string-output:
#   __str__() -- print name and residue range
#   print_tree() -- print tree detailed
#
### property-accessors:
#
#   first_residue
#   last_residue
#   sequence
#
### generators:
#
#  from_residue(aa_type)
#  from_sequence(sequence)
#       -make tree for a single residue or a full sequence
#
### adding atoms:
#
#  add_atom( new_atom, attachment_point, bond )
#       -add new_atom to tree and make a bond to attachment_point (which has to be an existing atom)
#
#  add_subtree( new_tree, target_atom, attachment_point, bond )
#       -add the new_tree to tree and make bond between target_atom (of new_tree) and attachment_point in tree
#
### manipulators:
#  offset_residue_numbers( self, offset )
#       -change the residue numbers in the tree
#
#
# TODO: if we have jumps, and add more residues on the chain before the jump, then residue numbers of all downstream of jump residues need to be
# changed
#
# TODO: accessors for is_bonded, counting of bonds (shortest_path in graph)

COVALENT = 1
SS_BOND  = 2
JUMP     = 4

import sys
from chemical.Atom import Atom
import library
from chemical import AllAminoAcidLib

class AtomTree:

	def __init__( self, root_atom ):
		self._residue_range=(root_atom.resid,root_atom.resid)
		self._atoms={root_atom:root_atom}   #the keys of bonds are the nodes in the graph, so also if singleton atoms are added by jumps, they need to be in bonds as 'atom:[]' pair
		#since it is slow to find Atoms if they are keys, we may want to maintain another array called
#		self._bonds={root_atom:[]}
		self._bonds={}
		self._jumps={}
		self._resid2atom={root_atom.resid:[root_atom]}
		pass

	def __str__( self ):
		return 'AtomTree from %d to %d'%(self._residue_range)

	def __iter__( self ):
		for atom in self._bonds:
			yield atom

	def print_tree( self ):
		print self
		for atom,bonded_atoms in self._bonds.iteritems():
			print '(%s) --> [%s]'%( atom, ', '.join(['(%s)'%atom for atom in bonded_atoms ]) )

	@property
  def last_residue(self):
		return self._residue_range[1]

	@property
	def first_residue(self):
		return self._residue_range[0]

	@property
	def sequence(self):
		sequence=''
		for resid in range(self._residue_range[0], self._residue_range[1]+1):
			from amino_acids import long2short
			try:
				sequence+=long2short(	self._bonds[Atom('CA',resid)][0].res3_type )
			except KeyError:
				sequence+='X'
		return sequence

	def partners( self, atom ):
		return self._bonds[atom]

	#return all bonded hydrogens for this atom
	def proton_partners( self, atom ):
		if atom.element=='H': raise ValueError('Atom %s has no proton_partners because it is a proton'%atom )
		return [ proton for proton in self._bonds[atom] if proton.element=='H' ]

	def atom( self, atom ):
		try:
			return self._atoms[atom]
		except KeyError:
			orig_exc=sys.exc_info()
      #if we look for a QX atom and its not there, we will get the first proton that fits the description
			#is this a good idea ?
			if atom.name[0]=="Q" and not 'QQ' in atom.name:
				query='C'+atom.name[1:]
				try:
					c_atom=self._atoms[Atom(query,atom.resid)]
				except KeyError:
					query=query+'1'
					try:
						c_atom=self._atoms[Atom(query,atom.resid)]
					except KeyError:
						raise orig_exc[0],orig_exc[1],orig_exc[2]
				for bonded in self.partners(c_atom):
					if bonded.element=='H': return bonded
		#couldn't handle this KeyError, re-raise
			raise orig_exc[0],orig_exc[1],orig_exc[2]

	def iter_residue( self, resid ):
		for atom in self._resid2atom[resid]:
			yield atom

	def _add_atom( self, new_atom ):
		self._atoms[new_atom]=new_atom
		self._resid2atom.setdefault(new_atom.resid,[]).append(new_atom)

	def add_atom( self, new_atom, attachment_point, bond=COVALENT ):
		#determine residue number for new_atom (if the bond is a peptide bond it is +1 or -1 from that of attachment_point
		# if bond is a 'JUMP' it is the last_residue
		resid=self._determine_residue_offset_from_attachment( new_atom, attachment_point, bond )
		new_atom.set_resid( resid )
		#new_atom should not be in the TREE yet: make sure:
		if new_atom in self._atoms:
			raise library.InconsistentInput('Trying to add atom %s to AtomTree which is already present'%new_atom)
		# add new_atom
		self._add_atom( new_atom )
		# add bond on both sides
		if bond==COVALENT:
			self.add_bond( attachment_point, new_atom )
		elif bond==JUMP:
			self.add_jump( new_atom, attachment_point )
		#update our residue ranges
		if self._residue_range[0]>resid:
			self._residue_range=(resid,self._residue_range[1])
		if self._residue_range[1]<resid:
			self._residue_range=(self._residue_range[0],resid)


	def add_subtree( self, tree, target_atom, attachment_point, bond=COVALENT ):
		#determine residue number for the target_atom at which new tree will be attached
		#if the new bond is a peptide bond it is +1 or -1 from that of attachment_point
		#if bond is a 'JUMP' it is the last_residue
		resid=self._determine_residue_offset_from_attachment( target_atom, attachment_point, bond )
		#offset the tree to the new residue number of the starting atom
		tree.offset_residue_numbers(resid-target_atom.resid)
		target_atom.set_resid(resid)
		#add the trees together, and add the bond on both sides
		self._bonds=dict(self._bonds,**tree._bonds)
		self._atoms=dict(self._atoms,**tree._atoms)
		self._jumps=dict(self._jumps,**tree._jumps)
		self._combine_resid_dictionaries(tree)
		#add the connection
		if bond==COVALENT:
			self.add_bond( target_atom, attachment_point )
		elif bond==JUMP:
			self.add_jump( attachment_point, target_atom )
		else:
			assert False
		#update residue range
		if self._residue_range[0]>tree._residue_range[0]:
			self._residue_range=(tree._residue_range[0],self._residue_range[1])
		if self._residue_range[1]<tree._residue_range[1]:
			self._residue_range=(self._residue_range[0],tree._residue_range[1])


	#add another bond between two already existing atoms, to make circles (HIS, TYR, TRP... ) or to make SS-BONDS
	def add_bond( self, from_atom, to_atom ):
		self._add_bond( from_atom, to_atom )
		self._add_bond( to_atom, from_atom )

	#not sure this should be a public method
	# offset residue numbers
	def offset_residue_numbers( self, offset ):
		#update the residue range
		self._residue_range=(self._residue_range[0]+offset,self._residue_range[1]+offset)
		#make a new tree (otherwise keys don't get properly updated)
		#first offset all atoms
		new_atoms={}
		for atom in self._atoms.iterkeys():
			atom.set_resid( atom.resid+offset )
			new_atoms[atom]=atom
		#now also regenerate the bond dictionary
		new_bonds={}
		for atom,bonds in self._bonds.iteritems():
			new_bonds[atom]=bonds
		#now also regenerate the jump dictionary
		new_jumps={}
		for atom,jumps in self._jumps.iteritems():
			new_jumps[atom]=jumps

		new_resid2atom={}
		for resid, list in self._resid2atom.iteritems():
			new_resid2atom[resid+offset]=list

		#assign new dictionaries
		self._bonds=new_bonds
		self._atoms=new_atoms
		self._jumps=new_jumps
		self._resid2atom=new_resid2atom

		return self
	#generator method, make tree from a fasta-sequence
	@classmethod
	def from_sequence( obj, sequence ):
		import sys
		obj=None
		for aa in sequence:
			try:
				new=AtomTree.from_aatype(aa)
			except Exception as exc:
				full_exc=sys.exc_info()
				msg='%s'%str(exc)+'   Failure when making AtomTree for aminoacid type %s'%aa
				raise type(exc),type(exc)(msg),full_exc[2]
			if not obj:
				obj=new
			else:
				obj.add_subtree( new, Atom('N',1), Atom('C',obj.last_residue ) )
		return obj

	#generator method, make a residue-tree for a given aa_type
	#this method uses SingleAminoAcidLib which specifies for a given aa_type,
	#which heavy-atoms and protons are present, and where the bonds are
	@classmethod
	def from_aatype( obj, aa_type ):
		aa_lib=AllAminoAcidLib()[aa_type]
		root=Atom('CA',1,'C',aa_type,aa_lib._atom_types['CA'])
		obj=AtomTree(root)
		for name in aa_lib._heavy_atoms:
			if name != 'CA':
				elem=name[0]
				atom=Atom(name,1,elem,aa_type,aa_lib._atom_types[name])
				obj._add_atom( atom )
				obj._bonds[atom]=[]
			else:
				atom=root

			#add protons
			try:
				protons=aa_lib.proton_partners(name)
				proton_atoms=[]
				for proton_name in protons:
					proton_atoms.append( Atom(proton_name,1,'H',aa_type,aa_lib._atom_types[proton_name]) )
				#correct this if there are methyl groups --> translate to QX
				if len(protons)==3 and len(aa_lib.heavy_partners(name))==1:
					methyl_name=Atom(protons[0],1,'H',aa_type).methyl_name()
					proton_atoms=[Atom(methyl_name,1,'H',aa_type,aa_lib._atom_types[protons[0]])]
				#now add the new atoms to the Tree
				for proton in proton_atoms:
					obj._add_atom( proton )
					obj._bonds[proton]=[atom]
					obj._bonds.setdefault(atom,[]).append(proton)
			except KeyError: #no protons, such as 'C' or 'O'
				pass
		#Add bonds between heavy-atoms
		for name in aa_lib.heavy_atoms():
			from_atom = Atom(name, 1)
			for partner in aa_lib.heavy_partners(name):
				to_atom = Atom( partner, 1 )
				obj._add_bond( from_atom, to_atom )
		return obj

	#add a bond (only one-direction) from from_atom to to_atom
	#use this method rather than adding it yourself, as here we make sure that the bond is made to the
	#correct instance of Atom (the one which is also used as key in the dictionary).
	def _add_bond( self, from_atom, to_atom ):
		# need to get the actual reference to the attachment_point Atom instance
		self._bonds.setdefault(from_atom,[]).append( self._atoms[to_atom] )

	def _add_jump( self, from_atom, to_atom ):
		# need to get the actual reference to the attachment_point Atom instance
		self._jumps.setdefault(from_atom,[]).append( self._atoms[to_atom] )

	def add_jump( self, atom1, atom2 ):
		self._add_jump( atom1, atom2 )
		self._add_jump( atom2, atom1 )
		# need to get the actual reference to the attachment_point Atom instance

	#resid is that of attachment_point, unless there was a peptide bond
	#there are some issues we haven't dealt with yet:
			# JUMPS if one adds another chain with a JUMP, but later adds more residues on the C-terminus of the first chain,
			# all residues connected by the JUMP need to be offset, obviously we need to keep book where the JUMPs are and stuff
			# also SS-BOND is not done, correctly, SS-BOND should only be added on-top of existing atoms, I suppose
	def _determine_residue_offset_from_attachment( self, new_atom, attachment_point, bond ):
		resid=attachment_point.resid
		if bond==COVALENT:
			if new_atom.name=="N" and attachment_point.name=="C":
				resid+=1
			elif new_atom.name=="C" and attachment_point.name=="N":
				resid-=1 #might yield negative residue numbers, YAY !
		#if it is a jump give it a new free residue number
		elif bond==JUMP:
			raise library.ProgramError('JUMP is not properly implemented. need to make sure when attaching more ',
																 'residues that numbers of JUMPed subtree are pushed forward accordingly')
			resid=self._residue_range[1]
		return resid

	def _combine_resid_dictionaries(self,other):
		for key,alist in other._resid2atom.iteritems():
			self._resid2atom.setdefault(key,[]).extend(alist)

def examples():
	root=Atom('CA',5)
	my_tree=AtomTree(root)
	my_tree.add_atom(Atom('N',5),root)
	my_tree.add_atom(Atom('C',5),root)
	my_tree.add_atom(Atom('CB',1),Atom('CA',5))
	my_tree._bonds[Atom('CA',5)]
	#my_tree.print_tree()
	print 'offset things...'
	my_tree.offset_residue_numbers(-4)
	my_tree.print_tree()
	my_tree._bonds=dict(my_tree._bonds)
	my_tree._bonds[Atom('CA',1)]

	import copy
	new_residue=copy.deepcopy(my_tree)
	new_residue.add_atom(Atom('CG',1),Atom('CB',1))
	#new_residue.print_tree()

	my_tree.add_subtree(new_residue,Atom('N',1),Atom('C',1))
	my_tree.print_tree()

	print 'final offseting'
	my_tree.offset_residue_numbers(+5)
	my_tree.print_tree()
	my_tree._bonds[Atom('N',6)]

	my_tree=AtomTree.from_aatype( 'ALA' )
	print 'ALA...'
	my_tree.print_tree()

	my_tree=AtomTree.from_aatype( 'TYR' )
	print 'TYR...'
	my_tree.print_tree()

	my_tree=AtomTree.from_sequence( 'HALLYGALLY' )
	my_tree.print_tree()
	print my_tree.sequence

#examples()
