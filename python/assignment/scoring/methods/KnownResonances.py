#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from assignment.noesy import ResonanceList
from basic import options
from assignment.knowledge.IndividualShiftDistributionLibrary import IndividualShiftDistributionLibrary
from ShiftDistributionScoreMethod import ShiftDistributionScoreMethod
from assignment.noesy.Atom import Atom
import unittest
from os import environ

if 'csrosettaDir' not in environ:
	print 'Please setup csrosettaDir to your environment'
	exit()

parser = options.ModuleArgumentParser("prot-score",description='options for scoring against known resonance assignments',add_help=False )
parser.add_argument("-prot",help="resonance assignments in prot file format", default=None )#the default is used for unittest

class KnownResonanceScoreMethod(ShiftDistributionScoreMethod):
	def __init__(self,resonance_list=None):
		try:
			args=parser.parse_args()
			self._resonance_list=ResonanceList.read_from_stream( open(args.prot,'r') )
		except:
			self._resonance_list=resonance_list
		individual_shift_lib=IndividualShiftDistributionLibrary.load_from_resonance_list(self._resonance_list)
		ShiftDistributionScoreMethod.__init__(self, 'prot_score_method',individual_shift_lib)


class KnownResonanceScoreMethodTestCase(unittest.TestCase):
	def setUp(self):
		s='''
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
		from StringIO import StringIO
		rl=ResonanceList.read_from_stream(StringIO(s))
		self.known_resonance_score_method=KnownResonanceScoreMethod(rl)

	def test_load_resoance_from_prot(self):
		real_freq=[4.394,1.905,1.839,1.919,2.415,175.642,54.943,33.407,16.931,31.943]
		real_tol=[0.04,0.04,0.04,0.04,0.04,0.40,0.40,0.40,0.40,0.40]
		atom_list=[Atom('HA',4),#first one
							 Atom('HB2',4),
							 Atom('HB3',4),
							 Atom('QE',4),
							 Atom('QG',4),
							 Atom('C',5),
							 Atom('CA',5),
							 Atom('CB',4),
							 Atom('CE',4),
							 Atom('CG',4)]#last one
		freq=[]
		tol=[]
		for atom in atom_list:
			freq.append(self.known_resonance_score_method.shift_distribution_library().by_atom(atom).mean())
			tol.append(self.known_resonance_score_method.shift_distribution_library().by_atom(atom).std())
		for f,r,a in zip(freq,real_freq,atom_list):
			self.assertAlmostEqual(r,f,3,'resonance freq of proton %s of residue %c should be %5.3f but we get %5.3f (which is wrong)'%(a.name(),a.resid(),r,f))
		for t,rt,a in zip(tol,real_tol,atom_list):
			self.assertAlmostEqual(rt,t,3,'resonance of proton %s of residue %c should be %5.3f but we get %5.3f (which is wrong)'%(a.name(),a.resid(),rt,t))

def KnownResonanceScoreMethodTestSuite():
	suite = unittest.TestLoader().loadTestsFromTestCase(KnownResonanceScoreMethodTestCase)
	return suite


if __name__ == '__main__':
	unittest.main()
