#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-


global __nomenclature_table
__nomenclature_table=None

import library

#private class for use by BmrbAtomNames
# reflects a single type of aminoacid
# maps BMRB-atomnames -->
#    a) rosetta(pdb) atom names  -- via _pdb_name
#    b) methyl- or ProR, ProS names  -- via _functions
class _BmrbAtomNamesResidue:
	def __init__(self, aa):
		self._aa = aa
		self._function = {}
		self._pdb_name = {}
		self._source_atom = {}
		self._rejects = {}
		self._nmr_name = {}

	def add( self, nmr_name, function, source_atom, pdb_name, rejects ):
		fs = None
		if function == "E":
				fs = "R"
		elif function == "Z":
				fs = "S"
		else:
			fs = function
		self._function[ nmr_name ] = fs
		self._pdb_name[ nmr_name ] = pdb_name
		self._nmr_name[ pdb_name ] = nmr_name
		self._source_atom[ nmr_name ] = source_atom
		if rejects:
			self._rejects[ nmr_name ] = rejects
	@classmethod
	def read_from_lines( obj, aa_in, lines ):
		obj = _BmrbAtomNamesResidue( aa_in )
		for line in lines:
			tags=line.split()
			aa=tags[0]
			assert aa == aa_in
			nmr_name=tags[1]
			function=tags[2]
			source_atom=tags[3] #e.g., HG11 is bound to CG1, CD1 is bound to CG
			pdb_name=tags[4]
			rejects=None
			if len(tags)==6:
				rejects=tags[5].split(',')
			obj.add( nmr_name, function, source_atom, pdb_name, rejects )
		return obj

# exported class
# provides translation from BMRB-NMR atomnames to methyl-, ProR/ProS and Rosetta names
# the information is automatically read from a database file
#   currently hardcoded to: os.environ['csrosettaDir']+"/database/bmrb_nomenclature.txt"
#
def __init__():
	__read_database_file()

# returns R or S for ProR/ProS protons
def get_proR_proS( aa, name ):
	if 'R' == __nomenclature_table[aa]._function[ name ]:
		return 'R'
	if 'S' == __nomenclature_table[aa]._function[ name ]:
		return 'S'
	return None

def anti_degenerate(aa,in_name):
	nmr_names = []
	for proton, qgx in __nomenclature_table[aa]._function.iteritems():
		if qgx==in_name:
			nmr_names.append( proton )
	if len(nmr_names):
		return nmr_names
	if in_name[0]=='Q':
		in_name = 'H'+in_name[1:]
	for carbons in __nomenclature_table[aa]._nmr_name.itervalues():
		if len(in_name)>=2 and in_name in carbons:
			nmr_names.append( carbons )
	if len(nmr_names)==0:
		return [in_name]
	return nmr_names


def translate( aa, name, reverse=False ):
	if name=='1H' and reverse:
		return 'H1'
	if name=='H1' and not reverse:
		return '1H'

	if name=='2H' and reverse:
		return 'H2'
	if name=='H2' and not reverse:
		return '2H'

	if name=='3H' and reverse:
		return 'H3'
	if name=='H3' and not reverse:
		return '3H'

	if name=='OXT' and reverse:
		return 'O2'

	if 'Q' in name:
		return name

	try:
		if reverse:
			return __nomenclature_table[aa]._nmr_name[ name ]
		else:
			return __nomenclature_table[aa]._pdb_name[ name ]
	except KeyError as exc:
		raise KeyError( "For aa %s %s"%(aa,exc) )

# returns PHXXX for proR/proS protons
# returns QXXX for methyl protons
def get_combine( aa, name ):
	if get_proR_proS( aa, name ):
		return 'P'+name[:-1], __nomenclature_table[aa]._source_atom[ name ]
	f = __nomenclature_table[aa]._function[ name ]
	s = __nomenclature_table[aa]._source_atom[ name ]
	if f=='.':
		f = None
	if s =='.':
		s = None
	return f, s

def get_rejects( aa, name ):
	return __nomenclature_table[aa]._rejects[ name ]

def get_source_atom(aa,name):
	return __nomenclature_table[aa]._source_atom[ name ]

def __read_database_file():
	import os
	database_file=os.environ['csrosettaDir']+"/database/bmrb_nomenclature.txt"
	file = open( database_file, 'r' )
#	print 'open file...',database_file
	aa = None
	aa_lines = []
	global __nomenclature_table
	__nomenclature_table = {}
	for line in file:
		tags=line.split()
		if len(tags) == 0:
			continue
		if tags[0][0]=="#":
			continue
		cols=5
		if len(tags)>2 and 'QQ' in tags[2]:
			cols=6
		if len(tags) != cols:
			raise library.InconsistentInput("Error in File %s. Expected %d columns in line: %s"%( database_file, cols, line ))
		if aa and aa != tags[0]:
			residue=_BmrbAtomNamesResidue.read_from_lines( aa, aa_lines )
			__nomenclature_table[aa]=residue
			aa_lines = []
		aa = tags[0]
		aa_lines.append( line )
	if len(aa_lines):
		residue=_BmrbAtomNamesResidue.read_from_lines( aa, aa_lines )
		__nomenclature_table[aa]=residue


import unittest
class TestBmrbAtomNames( unittest.TestCase ):
	def setUp(self):
		pass

	def test_combine(self):
		cases = [('A','HB2'), ('R','HH11'),('H','HB2'),('I','HG12'),('I','HG21'),('I','QG2'),('I','QD2'),('V','CG1')]
		answers1 = ['QB','PHH1','PHB','PHG1','QG2','QQG','Exception','PCG']
		answers2 = ['CB','NH1','CB','CG1','CG2','CG2','Exception','CB']

#		translator = BmrbAtomNames
		for (aa, proton), answer, answer2 in zip(cases,answers1,answers2):
			if answer == 'Exception':
				self.assertRaises( KeyError, get_combine, aa, proton )
			else:
				result, s = get_combine( aa, proton )
				self.assertEqual( answer, result )
				self.assertEqual( answer2, s )

	def test_get_proR_proS( self ):
		cases = [('A','HB2'), ('R','HH11'),('H','HB2'),('I','HG12'),('I','HG21')]
		answers = [None,'S','S','R',None]
		for (aa, proton), answer in zip( cases, answers):
			result = get_proR_proS( aa, proton )
			self.assertEqual(answer, result )

def BmrbAtomNamesTestSuite():
	suite = unittest.TestLoader().loadTestsFromTestCase(TestBmrbAtomNames)
	return suite

__init__()
if __name__ == '__main__':
	unittest.main()

