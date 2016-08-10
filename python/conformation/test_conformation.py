#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


from os import path
import unittest
import library

#import fasta
from chemical import Atom, AtomTree
from conformation import Conformation
#########
##  general unit-test setups
#get data path
import basic
data_path = basic.get_unittest_data()
#basic.fix_unittest_args()
#
#fix Tracers
from basic.Tracer import init_tracer_from_cmdline
init_tracer_from_cmdline([])

######
## The Test
class ConformationTestCase(unittest.TestCase):
	def __init__(self,name):
		super(ConformationTestCase,self).__init__(name)

	@classmethod
	def setUpClass(cls):
	#setup code that really should only run once
		print 'initialize ConformationTestCase...'
		cls.file1=data_path+'/2lxt.pdb'
		cls.file2=data_path+'/sgr145.pdb'
		cls.file3=data_path+'/4eis.pdb'

	def test_load_pdb(self):
		test_data=[(-25.316, 3.425, 6.092),
							 (-21.552, 6.067, 9.72),
							 (-24.569, 5.275, 6.209),
							 (-20.056, 3.014, 7.935),
							 (-24.233, 3.639, 8.412),
							 (-20.647, 7.116, 8.985),
							 (-22.492, 2.639, 7.417),
							 (-23.093, 1.623, 7.553),
							 (-21.593, 3.73, 6.706),
							 (-22.281, 3.118, 7.511),
							 (-20.738, 9.594, 5.568),
							 (-21.538, 3.981, 10.233),
							 (-22.81, 3.467, 6.333),
							 (-25.934, 1.111, 7.379),
							 (-22.61, 3.478, 7.801),
							 (-25.918, -3.471, 6.048),
							 (-22.257, 2.49, 6.075),
							 (-26.592, -0.108, 7.846),
							 (-22.562, 1.757, 7.613),
							 (-27.257, -2.981, -9.321)]
		from pprint import pprint as pp
		conf = Conformation.from_pdb_file( self.file1 )
		def retrieve( atom ):
			return conf.coords( atom )

		self.assertEqual( len(retrieve( Atom( 'CA', 589 ) ) ), 20 )
		self.assertRaises( KeyError, retrieve, Atom('CA',580) )
		self.assertRaises( KeyError, retrieve, Atom('H1',589) )
#		pp( conf.numpy_coords( Atom( 'CA', 589 ) ) )
		for x1,x2 in zip(test_data,retrieve( Atom( 'CA', 589 ) ) ):
			self.assertEqual( x1, tuple(x2) )

		import numpy
		test_data=numpy.matrix('''[27.46466748  26.76340225  27.1466444   23.29558529  27.24297746
  26.55670322  25.33201058  25.44459823  24.19545616  25.27204703
  26.31901687  26.38405363  25.28402589  27.45901927  25.73467097
  26.15315901  24.56910068  27.95507832  25.44748811  25.94197807]''').tolist()[0]

		for dist1,dist2 in zip(conf.dist( Atom( 'CA', 589 ), Atom( 'CA', 606 ) ),test_data):
			self.assertAlmostEqual(dist1,dist2)

	def test_load_more_pdb( self ):
		conf = Conformation.from_pdb_file( self.file2 )
		def read_missing_input(): #read a pdb with missing density without a predefined molecule is not possible
			conf = Conformation.from_pdb_file( self.file3 )
		self.assertRaises( library.MissingInput, read_missing_input )
		try:
			read_missing_input()
		except library.MissingInput as exc:
			seq=exc.sequence
		full_sequence=seq.replace('X','A')
		mol=AtomTree.from_sequence(full_sequence).offset_residue_numbers( 1 )
		conf = Conformation.from_pdb_file( self.file3, molecule=mol )

