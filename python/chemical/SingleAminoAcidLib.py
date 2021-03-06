#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os import environ
from PDB.Polypeptide import one_to_three,three_to_one
from basic import database
class SingleAminoAcidLib:
	def __init__(self,aa):
		self._aa=''
		self._atoms=[]
		self._heavy_atoms=[]
#		self._atom_types_old=[]
		self._atom_types={}
		self._protons=[]
		self._proton_partners={}
		self._heavy_partners={}
		self._partners={}
		self._nbonds={}
		self._stereo={}
		self._read_from_database(aa)
		self._check_consistency()

	def _check_consistency(self):
		for atom in self._atoms:
			self._atom_types[atom]
			self._partners[atom]
			self._heavy_partners[atom]

	def aa(self):
		return self._aa
	def atoms(self):
		return self._atoms

	def stereo(self,name):
		try: return self._stereo[name]
		except KeyError:
			#print 'the atom %s does not have stereo at %s'%(name,self._aa)
			return None

	def iter_atoms(self):
		for atom in self._atoms:
			yield atom

	def heavy_atoms(self):
		return self._heavy_atoms
	def protons(self):
		return self._protons
	def proton_partners(self,atom):
		return self._proton_partners[atom]
	def heavy_partners(self,atom):
		return self._heavy_partners[atom]
	def partners(self,atom):
		return self._partners[atom]
	def nbonds(self,atom):
		return self._nbonds[atom]

	def write_to_stream( self, fd ):
		fd.write('AMINO ACID %s\n'%self._aa)
		fd.write('HEAVY ATOMS '+' '.join(self._heavy_atoms)+'\n')
#		fd.write('ATOM_TYPES '+' '.join(self._atom_types_old)+'\n')
		fd.write('PROTONS '+' '.join(self._protons)+'\n')

		for atom in self._atoms:
			atom_types=self._atom_types[atom]
			if atom == "H":
				atom_types.add("ilvH")
			fd.write('ATOM_TYPE %s '%atom+' '.join(atom_types)+'\n')
#		for atom,atom_type in zip(self._heavy_atoms,self._atom_types_old):
		#+self._protons:
#			fd.write('ATOM_TYPE %s %s'%(atom,atom_type))
#			if 'C' in atom_type and len(atom_type)>1:
#				fd.write(' C')
#			fd.write('\n')
#		for proton in self._protons:
#			fd.write('ATOM_TYPE %s %s\n'%(proton,'H'))

		for atom, partners in self._proton_partners.iteritems():
			fd.write('PROTON PARTNERS %s '%atom+' '.join(partners)+'\n')
		visited=set()
		for atom, partners in self._heavy_partners.iteritems():
			if atom in self._protons: continue
			if atom in visited: continue
			visited.add(atom)
			for atom2 in partners:
				if atom2 in visited: continue
				fd.write('BOND %s %s\n'%(atom,atom2))
		for (atom1,atom2),length in self._nbonds.iteritems():
			if atom1<atom2: continue
			fd.write('BOND LENGTH %s %s %d\n'%(atom1,atom2,length))


	def _read_from_database(self,aa):
		if aa in 'ARNDCQEGHILKMFPSTWYV':
			aa_three=one_to_three(aa)
		elif aa in ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]:
			aa_three=aa
		else:
			print "the amino acid name is wrong, the wrong name is %s"%aa
			return None
		for line in database.open(aa_three+"_structure.txt"):
			tags=line.split()
			if tags[0]=='AMINO':
				self._aa=tags[2]
			if tags[0]=='HEAVY':
				for heavy in tags[2:]:
					self._heavy_atoms.append(heavy)
#			if tags[0]=='ATOM_TYPES':
#				for type in tags[1:]:
#					self._atom_types_old.append(type)
#				if len(self._atom_types_old) != len(self._heavy_atoms):
#					print 'ERROR: wrong number of ATOM_TYPES in database file for '+aa_three
#					assert False
			if tags[0]=='ATOM_TYPE':
				atom=tags[1]
				types=tags[2:]
				self._atom_types[atom]=set(types)
			if tags[0]=='PROTONS':
				for proton in tags[1:]:
					self._protons.append(proton)
			if tags[0]=='PROTON' and tags[1] == 'PARTNERS':
				for proton in tags[3:]:
					self._heavy_partners.setdefault(proton,[]).append(tags[2])
					self._partners.setdefault(proton,[]).append(tags[2])
					self._proton_partners.setdefault(tags[2],[]).append(proton)
					self._partners.setdefault(tags[2],[]).append(proton)
			if tags[0]=='BOND' and tags[1] != 'LENGTH':
				self._heavy_partners.setdefault(tags[1],[]).append(tags[2])
				self._partners.setdefault(tags[1],[]).append(tags[2])
				self._heavy_partners.setdefault(tags[2],[]).append(tags[1])
				self._partners.setdefault(tags[2],[]).append(tags[1])
			if tags[0]=='BOND' and tags[1] == 'LENGTH':
				self._nbonds[tags[2],tags[3]]=int(tags[4])
				self._nbonds[tags[3],tags[2]]=int(tags[4])
			if tags[0]=='STEREO':
				self._stereo[tags[1]]=tags[2]
				self._stereo[tags[2]]=tags[1]
		self._atoms=self._heavy_atoms+self._protons
		return self

def run_test():
		singlelib=SingleAminoAcidLib("TRP")
		print 'amino acid',singlelib._aa
		print 'heaby atoms',singlelib._heavy_atoms
		print 'protons',singlelib._protons
		print 'proton partners',singlelib._proton_partners
		print 'heavy partners',singlelib._heavy_partners
		print 'partners',singlelib._partners
		print 'nbonds',singlelib._nbonds
		singlelib.write_to_stream(open('test','w'))

def rewrite_all():
	for aa in 'ARNDCQEGHILKMFPSTWYV':
		lib=SingleAminoAcidLib(aa)
		aa_three=one_to_three(aa)
		out=open(aa_three+"_structure.txt",'w')
		lib.write_to_stream(out)

#rewrite_all()
