#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os import environ
from PDB.Polypeptide import one_to_three,three_to_one
from SingleAminoAcidLib import SingleAminoAcidLib
if 'csrosettaDir' not in environ:
	print 'Please setup csrosettaDir to your environment'
	exit()

class AllAminoAcidLib:

	lib=None
	_aa_list=["ALA","ARG","ASN","ASP","CYS",
						"GLN","GLU","GLY","HIS","ILE",
						"LEU","LYS","MET","PHE","PRO",
						"SER","THR","TRP","TYR","VAL"]

	@classmethod
	def initialize(obj):
		obj.lib={}
		for aa in obj._aa_list:
			obj.lib[aa]=SingleAminoAcidLib(aa)

	def __init__(self):
		if self.lib==None:
			self.initialize()

  def __getitem__(self,key):
		if len(key)==1:
			key=one_to_three(key)
		return self.lib[key]

#	def aa_lib(self,aa):
#		if aa in 'ARNDCQEGHILKMFPSTWYV':
#			aa_three=one_to_three(aa)
#		elif aa in ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]:
#			aa_three=aa
#		else:
#			print "the amino acid name is wrong, the wrong name is %s"%aa
#			return None
#		return self._aa_lib[aa_three]

