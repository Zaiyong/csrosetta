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
class CollectionTestCase(unittest.TestCase):

	def __init__(self,name):
		super(CollectionTestCase,self).__init__(name)

	@classmethod
	def setUpClass(cls):
	#setup code that really should only run once
		print 'initialize CollectionTestCase...'

	def setUp(self):
		self.protein=AtomTree.from_sequence('HALLYGALLY')
		RN=PeakRule.RuleNode
		self.nch=NoesyRule(((RN(None,'H'),RN(1,'N',tol=0.01)),(RN(3,'H',tol=0.03),RN(2,'C',tol=0.03))))
		self.nhh=NoesyRule(((RN(3,'H',tol=0.03),None),(RN(2,'H',tol=0.03),RN(1,'N',tol=0.01))))
		self.nhhn=NoesyRule(((RN(3,'H',tol=0.03),RN(4,'N',tol=0.3)),(RN(2,'H',tol=0.03),RN(1,'N',tol=0.01))))
		peak_list1=PeakList('nnh')
		peak_list1.append( Peak([1.0,0.5,0.3],self.nhh ) )
		peak_list1.append( Peak([2.0,0.2,0.35],self.nhh ) )
		peak_list2=PeakList('nch')
		peak_list2.append( Peak([1.0,0.5,0.3],self.nch ) )
		peak_list2.append( Peak([2.0,0.2,0.35],self.nch ) )
		self.peaks=PeakCollection()
		self.peaks.add_experiment( peak_list1 )
		self.peaks.add_experiment( peak_list2 )

	def test_partial_matches( self ):
		state=AssignmentCollection( self.peaks, self.protein )
		ct=0
		tuples=set()
		for peak in self.peaks.iterpeaks():
			for mask in peak.rule.spinsystem_match_masks():
				for match in peak.matches( self.protein, match_mask=mask ):
#					print match
					ct+=1
					state.add( match )
					tuples.add(match.peak_match)
		ct_again=0
		state.commit()
		for match in state:
			ct_again+=1
		self.assertEqual(ct,ct_again) #all there ?
		#can I retrieve them all individually ?
		ct_again=0
		for tuple in tuples:
#			print " <==> ".join(str( m ) for m in state[tuple])
			ct_again+=1
			self.assertEqual( len( state[tuple] ), 2) # should be two, since each list had two peaks and matching was without any filter
		self.assertEqual( ct_again, ct/2 ) #should be half  since we have two peaks in each list

	def test_unassigned_peaks( self ):
		state=AssignmentCollection( self.peaks, self.protein )
		ct_peaks=0
		for peak in self.peaks:
			ct_peaks+=1
		self.assertEqual(len(state.unassigned_peaks),ct_peaks)
		for peak in self.peaks:
			ct_match_per_peak=0
			for match in peak.matches( self.protein, 1 ):
				state.add( match )
				ct_match_per_peak+=1
			self.assertEqual(ct_match_per_peak, 1)
		self.assertEqual(len(state.unassigned_peaks),0)
		for match in [ match for match in state ]:
			state.remove( match )
		self.assertEqual(len(state.unassigned_peaks),ct_peaks)

