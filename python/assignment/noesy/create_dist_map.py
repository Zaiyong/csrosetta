#!/home/zak/bin/python
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


from PDB.PDBParser import PDBParser
from PDB.PDBIO import PDBIO, Select
from PDB.Structure import Structure

from assignment.noesy import Atom

def create_dist_map(pdbfile):
	pdb=PDBParser()
	structure=pdb.get_structure('pdb',pdbfile)
	atoms_coord={}
	for atom in structure.get_atoms():
		atoms_coord[atom.name, atom.parent.get_id()[1]]=atom.coord
	return atoms_coord

