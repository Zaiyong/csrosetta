#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os import path
import unittest
#import assignment
#import libraries
#import fasta
from assignment import AssignmentCollection
from assignment import Peak, PeakList, PeakCollection
#Collection, PeakList, AssignmentCollection
from chemical import AtomTree,Atom
from rules.BasicRules import PeakRule
from rules.NoesyRule import NoesyRule
from strips import collect_strips_from_peak_list
#########
##  general unit-test setups
#get data path
#import basic
#data_path = basic.get_unittest_data()
#basic.fix_unittest_args()
#
#fix Tracers
from basic.Tracer import init_tracer_from_cmdline
init_tracer_from_cmdline([])

######
## The Test
class StripTestCase(unittest.TestCase):

	def __init__(self,name):
		super(StripTestCase,self).__init__(name)

	@classmethod
	def setUpClass(cls):
	#setup code that really should only run once
		print 'initialize StripTestCase...'

	def setUp(self):
		self.protein=AtomTree.from_sequence('HALLYGALLY')
		RN=PeakRule.RuleNode
		self.nch=NoesyRule(((RN(None,'H'),RN(1,'N',tol=0.01)),(RN(3,'H',tol=0.03),RN(2,'C',tol=0.03))))
		self.nhh=NoesyRule(((RN(3,'H',tol=0.03),None),(RN(2,'H',tol=0.03),RN(1,'N',tol=0.01))))
		self.peak_list1=PeakList('nnh')
		self.peak_list1.append( Peak([1.0,0.5,0.3],self.nhh ) )
		self.peak_list1.append( Peak([2.0,0.2,0.15],self.nhh ) )
		self.peak_list1.append( Peak([2.0,0.2,0.35],self.nhh ) )
		self.peak_list1.append( Peak([2.0,0.2,0.55],self.nhh ) )
		self.peak_list2=PeakList('nch')
		self.peak_list2.append( Peak([1.0,0.5,0.3],self.nch ) )
		self.peak_list2.append( Peak([2.0,0.2,0.35],self.nch ) )
		self.peaks=PeakCollection()
		self.peaks.add_experiment( self.peak_list1 )
		self.peaks.add_experiment( self.peak_list2 )

	def test_extract_freqs(self):
		collect_strips_from_peak_list( self.peak_list1 )
		collect_strips_from_peak_list( self.peak_list2 )

	def test_strip_match(self):
		strips=collect_strips_from_peak_list( self.peak_list1 )[0]
		ct_strips=0
		ct_full=0
		for strip in strips:
			for strip_match in strip.matches( self.protein ):
				ct_strips+=1
				for match in strip_match.matches( self.protein ):
					ct_full+=1
		self.assertEqual(ct_full, 3024 )
		self.assertEqual(ct_strips, 24 )



