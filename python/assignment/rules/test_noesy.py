#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os import path
import unittest
#import assignment
#import libraries
#import fasta
from assignment import Peak
#Collection, PeakList, AssignmentCollection
from chemical import AtomTree,Atom
from BasicRules import PeakRule
from NoesyRule import NoesyRule
#########
##  general unit-test setups
#get data path
import basic
data_path = basic.get_unittest_data()
basic.fix_unittest_args()
#
#fix Tracers
from basic.Tracer import init_tracer_from_cmdline
init_tracer_from_cmdline(['assignment.collection:error','assignment.peaks:warning','scoring.methods:warning'])

######
## The Test
class NoesyRuleTestCase(unittest.TestCase):

	def __init__(self,name):
		super(NoesyRuleTestCase,self).__init__(name)

	@classmethod
	def setUpClass(cls):
	#setup code that really should only run once
		print 'initialize NoesyRuleTestCase...'

	def setUp(self):
		self.protein=AtomTree.from_sequence('HALLYGALLY')
		self.dialanin=AtomTree.from_sequence('AA')
		self.leucin=AtomTree.from_sequence('L')

		RN=PeakRule.RuleNode
		self.nch=NoesyRule(((RN(None,'H'),RN(1,'N',tol=0.01)),(RN(3,'H',tol=0.03),RN(2,'C',tol=0.03))))
		self.nhh=NoesyRule(((RN(3,'H',tol=0.03),None),(RN(2,'H',tol=0.03),RN(1,'N',tol=0.01))))
		self.nhhn=NoesyRule(((RN(3,'H',tol=0.03),RN(4,'N',tol=0.3)),(RN(2,'H',tol=0.03),RN(1,'N',tol=0.01))))
#		from assignment.Peak import Peak, MutualExclusivePeak
#		peak=Peak([1.0,0.5,0.3],noesy_nch)
#		peak2=Peak([1.3,1.5,2.0],noesy_nch)

	def test_free_nch_matching(self):
		peak=Peak([1.0,0.5,0.3],self.nch)
		#check HA,CA and QB,CB
		matches = [ match.peak_match for match in peak.matches( self.dialanin ) ]
		self.assertEqual(len(matches),8)
		for match in matches:
			self.assertEqual(match[0].name,'N')
			self.assertEqual(match[1].name[1],match[2].name[1]) # should be CB-QB and CA-HA, i.e., .B-.B and .A-.A
			self.assertTrue(match[2].name=='QB' or match[2].name=='HA')
		#check longer sidechains
		matches = [ match.peak_match for match in peak.matches( self.leucin ) ]
		true_matches=set([('CA','HA'),('CD1','QD1'),('CD2','QD2'),('CB','HB2'),('CB','HB3'),('CG','HG')])
		self.assertEqual(len(matches),len(true_matches))
		for match in matches:
			try:
				self.assertTrue((match[1].name,match[2].name) in true_matches)
			except:
				print 'wrong match: ', match
				raise
		#check distance atoms
		matches = [ match.rule_match for match in peak.matches( self.leucin ) ]
		for match in matches:
			self.assertEqual(len(match),4)  #should now be 4D
			self.assertEqual(match[0].name,'H')
			self.assertEqual(match[0].resid,match[1].resid)
			self.assertEqual(match[2].resid,match[3].resid)
			true_matches=set([('H','HA'),('H','QD2'),('H','QD1'),('H','HB2'),('H','HG'),('H','HB3')])
		for match in peak.matches( self.leucin ):
			self.assertEqual(len( match.distance_atoms()),1)
			self.assertEqual(len( match.distance_atoms()[0]),2)
			try:
				pair=match.distance_atoms()[0]
				self.assertTrue( (pair[0].name,pair[1].name) in true_matches )
			except:
				print 'wrong match: ', match.distance_atoms()[0]
				raise

	def test_free_nhh_matching(self):
#		print
		peak=Peak([1.0,0.5,0.3],self.nhh)
		#check HA,CA and QB,CB
		matches = [ match.peak_match for match in peak.matches( self.dialanin ) ]
		self.assertEqual(len(matches),12)
		for match in matches:
			self.assertEqual(match[0].resid, match[1].resid )
			self.assertEqual(match[0].name,'N')
			self.assertEqual(match[1].name,'H')
		#check distance atoms
		true_matches=set([('H','HA'),('H','H'),('H','QD2'),('H','QD1'),('H','HB2'),('H','HG'),('H','HB3')])
		for match in peak.matches( self.leucin ):
			try:
				pair=match.distance_atoms()[0]
				self.assertTrue( (pair[1].name,pair[0].name) in true_matches )
			except:
				print 'wrong match: ', match.distance_atoms()[0]
				raise
		#check side-chain amide protons in LYS residues
		self.amide_rich=AtomTree.from_sequence('K')
		matches=[ match.peak_match for match in peak.matches( self.amide_rich ) ]
		self.assertEqual( len(matches), 22 )
		true_matches=set([('N','H'),('NZ','QZ')])
		for match in matches:
			try:
				self.assertTrue( (match[0].name, match[1].name) in true_matches )
			except:
				print 'wrong match: ', match
				raise
		#check side-chain amide protons in ARG residues
		self.amide_rich=AtomTree.from_sequence('R')
		matches=[ match.peak_match for match in peak.matches( self.amide_rich ) ]
		self.assertEqual(len(matches), 78  )

		#ilvH atom-types show us which hydrogens are present when ILV samples are used. Only match
		# N-H and certain methyl H of ILV.
	def test_ilv_matching(self):
		RN=PeakRule.RuleNode
		nhh_ilv=NoesyRule(((RN(3,'H',tol=0.03,types=['ilvH']),None),(RN(2,'H',tol=0.03,types=['ilvH']),RN(1,'N'))))
		peak=Peak([1.0,0.5,0.3],nhh_ilv)
		matches = [ match.peak_match for match in peak.matches(self.dialanin) ]
		self.assertEqual(len(matches),4)
		for match in matches:
			self.assertEqual(match[0].name,'N')
			self.assertEqual(match[1].name,'H')
			self.assertEqual(match[2].name,'H')
		#check LEU which also has ilv protonated methyls
		matches = [ match.peak_match for match in peak.matches(self.leucin) ]
		self.assertEqual(len(matches),3)
		true_matches=set(['QD2','QD1','H'])
		for match in matches:
			self.assertEqual(match[0].name,'N')
			self.assertEqual(match[1].name,'H')
			try:
				self.assertTrue(match[2].name in true_matches)
			except:
				print 'wrong match: ', match
				raise

	def test_existing_assignments(self):
		peak=Peak([1.0,0.5,0.3],self.nhh)
		peak.hard_assignments=[
			(Atom('N',1),Atom('H',1),Atom('QD2',3)),
			(Atom('N',2),Atom('H',2),Atom('HA',3))]
		matches = [ match.peak_match for match in peak.matches( self.protein ) ]
		self.assertEqual(len(matches),2)
		match=matches[0]
		self.assertEqual(match[0].resid,1)
		self.assertEqual(match[2].resid,3)
		self.assertEqual(match[2].res3_type,'LEU')
		match=matches[1]
		self.assertEqual(match[0].resid,2)
		self.assertEqual(match[2].resid,3)
		self.assertEqual(match[2].res3_type,'LEU')
#now a bad pre-assignment, which cannot be fullfilled
		peak=Peak([1.0,0.5,0.3],self.nhh)
		peak.hard_assignments=[
			(Atom('C',1),Atom('H',1),Atom('QD2',3))]
		import library
		def do_bad():
			for match in peak.matches ( self.protein ):
				print match
		self.assertRaises(library.InconsistentInput, do_bad )
		peak.hard_assignments=[
			(Atom('N',1),Atom('H',1),Atom('QD2',2))]
		self.assertRaises(KeyError,do_bad)

		peak=Peak([1.0,0.5,0.3],self.nch) #NCH
		peak.hard_assignments=[
			(Atom('N',1),Atom('CD2',3),Atom('QD2',3))]
		matches = [ match.peak_match for match in peak.matches( self.protein ) ]
		self.assertEqual(len(matches),1)
		match=matches[0]
		self.assertEqual(match[0].resid,1)
		self.assertEqual(match[2].resid,3)
		self.assertEqual(match[2].res3_type,'LEU')

		peak=Peak([1.0,0.5,0.3],self.nch) #NCH
		peak.hard_assignments=[
			(Atom('N',1),Atom('CD2',3),Atom('QD2',3))]
		matches = [ match.peak_match for match in peak.matches( self.protein ) ]
		self.assertEqual(len(matches),1)
		match=matches[0]
		self.assertEqual(match[0].resid,1)
		self.assertEqual(match[2].resid,3)
		self.assertEqual(match[2].res3_type,'LEU')
		matches = [ match for match in peak.matches( self.protein ) ]
		match=matches[0].rule_match
		self.assertEqual(match[0].name,'H')
		self.assertEqual(match[1].name,'N')
		self.assertEqual(matches[0].distance_atoms()[0],(Atom('H',1),Atom('QD2',3)))
		#check that blanks are filled in with matches correctly
		peak.hard_assignments=[
			(None,Atom('CB',2),Atom('QB',2))]
		matches = [ match.peak_match for match in peak.matches( self.dialanin ) ]
		self.assertEqual(len(matches),2)
		match=matches[0]
		self.assertEqual(match[0].name,'N')
		self.assertEqual(match[1].name,'CB')
		self.assertEqual(match[0].res3_type,'ALA')
		self.assertEqual(match[0].resid,1)
		match=matches[1]
		self.assertEqual(match[0].resid,2)
		#and also the dependent blanks
		peak.hard_assignments=[
			(None,None,Atom('QB',2))]
		matches = [ match.peak_match for match in peak.matches( self.dialanin ) ]
		self.assertEqual(len(matches),2)
		match=matches[0]
		self.assertEqual(match[0].name,'N')
		self.assertEqual(match[1].name,'CB')

	def test_frequency_match(self):
		class FreqMatch:
			from chemical.Atom import Atom
			known={ Atom('QB',2):0.3, Atom('CB',2):0.51, Atom('N',2):1.0, Atom('H',7):0.5, Atom('H',2):0.2,Atom('N',1):1.0,
							Atom('CA',1):0.5,Atom('CA',5):0.5, Atom('HA',1):0.3, Atom('H',1):0.3,  Atom('H',5):0.1,Atom('N',7):1.005,
							Atom('N',3):1.3,Atom('HA',3):2.0,Atom('CA',3):1.5}
			def __init__(self):
				self.ct=0
			def __call__(self, atom, freq, tol, folder):
				if atom in self.known:
					return abs(self.known[atom]-freq)<tol
				else:
					return False
		freq_match=FreqMatch()
		peak=Peak([1.0,0.5,0.3],self.nch)
#		peak2=Peak([1.3,1.5,2.0],noesy_nch)
		matches = [ match.peak_match for match in peak.matches( self.protein, frequency_matcher=freq_match ) ]
		self.assertEqual(len(matches),6)
		for match in matches:
			self.assertLess(abs(freq_match.known[match[0]]-1.0),0.01)
			self.assertLess(abs(freq_match.known[match[1]]-0.5),0.03)
			self.assertLess(abs(freq_match.known[match[2]]-0.3),0.03)

	def test_distance_match(self):
		def switch_if_bigger(a1,a2):
			if a1.resid>a2.resid:
				return a2,a1
			else:
				return a1,a2
		class DistMatch:
			from chemical.Atom import Atom
			known={ (Atom('QB',2),Atom('H',5)):0.3 }
			def __init__(self):
				self.ct=0
			def __call__(self, atom1, atom2, cutoff):
				atom1,atom2 = switch_if_bigger(atom1,atom2)
				self.ct+=1
				try:
					return self.dist(atom1,atom2)<cutoff
				except KeyError:
					return False
			def dist(self,atom1,atom2):
				return self.known[(atom1,atom2)]

		#direct calls of DistanceMatcher -- testing the test
		dist_match = DistMatch()
		self.assertTrue(dist_match(Atom('QB',2),Atom('H',5),1))
		self.assertFalse(dist_match(Atom('HA',2),Atom('H',5),1))
		self.assertEqual(dist_match.ct,2)

		#now use it in peak-matching
		#this should give one positive match -- from 588 calls
		dist_match.ct=0
		dist_match.visited={}
		peak=Peak([1.0,0.5,0.3],self.nch)
		matches = [ match.peak_match for match in peak.matches( self.protein, distance_matcher=dist_match ) ]
		self.assertEqual(dist_match.ct, 588)
		self.assertEqual(len(matches),1)
		match=matches[0]
		self.assertEqual(match[0].resid,5)
		self.assertEqual(match[1].resid,2)
		self.assertEqual(match[2].name,'QB')

		#now set distance above cutoff of 5, then we shouldn't have a match anymore
		#this should give zero matches -- from 588 calls
		dist_match.known={ (Atom('QB',2),Atom('H',5)):5.1 }
		dist_match.ct=0
		dist_match.visited={}
		matches = [ match.peak_match for match in peak.matches( self.protein, distance_matcher=dist_match ) ]
		self.assertEqual(len(matches),0)
		self.assertEqual(dist_match.ct, 588)

		#this should get two matches, as we can have NH-H and the symmetric peak H-NH
		dist_match.known={ (Atom('H',2),Atom('H',5)):4.8 }
		dist_match.ct=0
		dist_match.visited={}
		peak=Peak([1.0,0.5,0.3],self.nhh)
		matches = [ match.peak_match for match in peak.matches( self.protein, distance_matcher=dist_match ) ]
		self.assertEqual(len(matches),2)
		self.assertEqual(dist_match.ct, 756)

	def test_masked_matching( self ):
		peak=Peak([1.0,0.5,0.3],self.nch)
		ct_matches=0
		for match in peak.matches( self.protein, match_mask=[False, True, True] ):
			for full_match in peak.matches( self.protein, partial_match=match ):
				ct_matches+=1
		#do we get just as many matches if we don't split matching into two steps ?
		for match in peak.matches( self.protein ):
			ct_matches-=1
		self.assertEqual(ct_matches,0)
			#print match, full_match.peak_match
		#repeat with a different peak-rule
		peak=Peak([1.0,0.5,0.3],self.nhh)
		ct_matches=0
		for match in peak.matches( self.protein, match_mask=[True, True, False] ):
			for full_match in peak.matches( self.protein, partial_match=match.peak_match ):
				ct_matches+=1
		#do we get just as many matches if we don't split matching into two steps ?
		for match in peak.matches( self.protein ):
			ct_matches-=1
		self.assertEqual(ct_matches,0)
		#repeat with a different peak-rule
		peak=Peak([1.0,2.0,0.5,0.3],self.nhhn)
		ct_matches=0
		for match in peak.matches( self.protein, match_mask=[True, True, False, False] ):
			for full_match in peak.matches( self.protein, partial_match=match.peak_match ):
				ct_matches+=1
		#do we get just as many matches if we don't split matching into two steps ?
		for match in peak.matches( self.protein ):
			ct_matches-=1
		for match in peak.matches( self.protein, match_mask=[False, False, True, True] ):
			for full_match in peak.matches( self.protein, partial_match=match.peak_match ):
				ct_matches+=1
		#do we get just as many matches if we don't split matching into two steps ?
		for match in peak.matches( self.protein ):
			ct_matches-=1
		self.assertEqual(ct_matches,0)
		#
		self.assertEqual( self.nhhn.spinsystem_match_masks(), [(False, False, True, True), (True, True, False, False)] )
		self.assertEqual( self.nhh.spinsystem_match_masks(), [(True, True, False)] )
		self.assertEqual( self.nch.spinsystem_match_masks(), [(False, True, True)] )
		#
		#check that exception is thrown if method is not implemented by Rule
		from BasicRules	import ExistingAssignmentRule
		class SomeRule( ExistingAssignmentRule ):
			pass
		some_rule=SomeRule([])
		def call_it():
			some_rule.spinsystem_match_masks()
		self.assertRaises( NotImplementedError, call_it )
		#
		#check that exception is thrown if partial match is attempted for rule that does not implement this
		from BasicRules	import PeakRule
		class SomeRule( PeakRule ):
			pass
		some_rule=SomeRule([])
		peak=Peak([1.0,2.0,0.5,0.3],some_rule)
		def do_match():
			for full_match in peak.matches( self.protein, partial_match=[Atom('N',1)] ):
				pass
		self.assertRaises( NotImplementedError, do_match )





def example():
	from chemical import AtomTree,Atom
	protein=AtomTree.from_sequence('HALLYGALLY')
	dialanin=AtomTree.from_sequence('AA')
	leucin=AtomTree.from_sequence('L')

	RN=PeakRule.RuleNode
	noesy_nch=NoesyRule(((RN(None,'H'),RN(1,'N')),(RN(3,'H'),RN(2,'C'))))
	noesy_nhh=NoesyRule(((RN(3,'H',tol=0.03),None),(RN(2,'H',tol=0.03),RN(1,'N'))))
	from assignment.Peak import Peak, MutualExclusivePeak

	peak=Peak([1.0,0.5,0.3],noesy_nch)
	peak2=Peak([1.3,1.5,2.0],noesy_nch)
	peak._hard_assignments=[
		(Atom('N',1),None,Atom('HA',1)),
		(Atom('N',3),None,Atom('HA',3)),
		]
	peak2._hard_assignments=[
		(Atom('N',1),None,Atom('HA',1)),
		(Atom('N',2),None,Atom('HA',2)),
		(Atom('N',3),None,Atom('HA',3)),
	]

	mutex_peak=MutualExclusivePeak([peak,peak2])
	print mutex_peak
	for match in mutex_peak.matches( protein, freq_match ):
		print 'Mutex: ',match
	exit(1)
	peak=Peak([0.3,0.5,1.0],noesy_nch)
	peak2=Peak([0.3,0.5,1.0],noesy_nch)
	peak2._hard_assignments=[
		(Atom('HA',1),None,Atom('N',2)),
		(Atom('QB',2),None,Atom('N',2)),
		(Atom('HA',5),None,Atom('N',2)),
		(None,None,Atom('N',4)),
		(Atom('HA',7),Atom('H',5),Atom('N',3))
]

	free_nhh_peak=Peak((1.0,0.5,0.3),noesy_nhh)
	nhh_peak=Peak((0.3,0.5,1.0),noesy_nhh)
	nhh_peak._hard_assignments=[
#		(Atom('HA',1),None,Atom('N',2)),
#		(Atom('QB',2),None,Atom('N',2)),
#		(Atom('HA',5),None,Atom('N',2)),
#		(None,None,Atom('N',7)),
#		(Atom('HA',1),Atom('H',1),None),
		(Atom('HA',2),Atom('HA',5),Atom('CA',5))
]
#	print 'get raw matches for nhh_peak'
#	for match in noesy_nhh.matches( free_nhh_peak, protein ):
#		print 'free match: ', match



	print 'get raw matches for nhh_peak'
	for match in free_nhh_peak.matches( protein, freq_match  ):
		print 'full match: ', match


	exit(1)
	print 'get matches for pre-assigned peak'
	for match in noesy_nch.matches( peak2, protein):
		print 'full match: ', match

	all_matches=[match for match in noesy_nch.matches( peak, protein) ]
	print 'get %d matches for generic peak'%(len(all_matches))
	print 'first 5 matches:\nfull_match: ','\nfull_match:  '.join(['%s,%s,%s,%s'%m for m in all_matches[0:5]])

#		print noesy_nch.translate_full_match_to_peak_match(match)
	print 'print matched atoms...'
	for match in noesy_nch.matches(peak,protein,freq_match):
		print 'FIN: ',noesy_nch.translate_full_match_to_peak_match(match)
	print freq_match

	print 'pre-asigned stuff without frequency comparison...'
	for match in noesy_nch.matches(peak2,protein):
		print 'FREE: ',noesy_nch.translate_full_match_to_peak_match(match)

	freq_match.reset()
	print 'print matched atoms for pre-assigned peak2...'
	for match in noesy_nch.matches(peak2,protein,freq_match):
		print 'FIN: ',noesy_nch.translate_full_match_to_peak_match(match)
	print freq_match

	exit(0)

	print 'print hnca-rule matched atoms...'
	#N-H-CA
	hnca_rule=HNCA_Rule([3,1,2],[0.03]*3)
	peak3=Peak([0.3,0.5,0.5],hnca_rule)
	for match in hnca_rule.matches( peak3, protein ):
		print hnca_rule.translate_full_match_to_peak_match(match)

	import copy
	peak4=copy.copy(peak3)
#	peak4._hard_assignments=[
	peak4._hard_assignments=[
		(None,Atom('CA',8), Atom('N',8)),
		(Atom('H',6),Atom('CA',6), Atom('N',6)),
		(None,Atom('CA',3), None),
		(None, None,Atom('N',6)),
		(None, Atom('CA',5),Atom('N',5)),
		(Atom('H',5),None,Atom('N',5))]

	print 'print hnca-rule pre-assigned matched atoms...'
	for match in hnca_rule.matches( peak4, protein ):
		print hnca_rule.translate_full_match_to_peak_match(match)

	print 'print hnca-rule freq_matched atoms...'
	freq_match.reset()
	for match in hnca_rule.matches( peak3, protein, freq_match ):
		print hnca_rule.translate_full_match_to_peak_match(match)
	print freq_match

	print 'print hnca-rule freq_matched from pre-assigned atoms...'
	freq_match.reset()
	for match in hnca_rule.matches( peak4, protein, freq_match ):
		print hnca_rule.translate_full_match_to_peak_match(match)
	print freq_match
	#hnca_peak=Peak(3,[0.3,0.5,0.5],hnca_rule)
	#for match in hnca_rule.iterate_matched( hnca_peak, protein, freq_match ):
#		print match
#example()
