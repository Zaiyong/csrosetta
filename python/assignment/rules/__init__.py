#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module
from BasicRules import PeakRule, ExistingAssignmentRule
from HN_Rules import HNCA_Rule
from NoesyRule import NoesyRule
from ScoreDistanceMatcher import ScoreDistanceMatcher
#import test_noesy
#keep this private... from ResonanceAssignmentRule import ResonanceAssignmentRule

####-----------------------------------
# example code for rules module below
#
#

def example():
	from chemical import AtomTree,Atom
	protein=AtomTree.from_sequence('HALLYGALLY')
	RN=PeakRule.RuleNode
	noesy_hcn=NoesyRule(((RN(None,'H'),RN(1,'N')),(RN(3,'H'),RN(2,'C'))))
	noesy_hhn=NoesyRule(((RN(3,'H',tol=0.03),None),(RN(2,'H',tol=0.03),RN(1,'N'))))
	from assignment.Peak import Peak, MutualExclusivePeak

	peak=Peak([1.0,0.5,0.3],noesy_hcn)
	peak2=Peak([1.3,1.5,2.0],noesy_hcn)
	peak._hard_assignments=[
		(Atom('N',1),None,Atom('HA',1)),
		(Atom('N',3),None,Atom('HA',3)),
		]
	peak2._hard_assignments=[
		(Atom('N',1),None,Atom('HA',1)),
		(Atom('N',2),None,Atom('HA',2)),
		(Atom('N',3),None,Atom('HA',3)),
	]

	class FreqMatch:
		from chemical.Atom import Atom
		known={ Atom('QB',2):0.3, Atom('CB',2):0.5, Atom('N',2):1.0, Atom('H',7):0.5, Atom('H',2):0.2,Atom('N',1):1.0,
						Atom('CA',1):0.5,Atom('CA',5):0.5, Atom('HA',1):0.3, Atom('H',1):0.3,  Atom('H',5):0.1,Atom('N',7):1,
						Atom('N',3):1.3,Atom('HA',3):2.0,Atom('CA',3):1.5}
		def __init__(self):
			self.ct=0
		def __call__(self, atom, freq, tol, folder):
			self.ct+=1
			if atom in self.known:
				print 'match atom %s with %8.3f against its %8.3f'%(atom,freq,self.known[atom])
				return abs(self.known[atom]-freq)<tol
			else:
#				print 'match unknown atom %s with %8.3f'%(atom,freq)
				return False
		def __str__(self):
			return 'FreqMatch has been used %d times'%self.ct
		def reset(self):
			self.ct=0
#	def freq_match(atom,freq,tol,folder):

	freq_match=FreqMatch()

	mutex_peak=MutualExclusivePeak([peak,peak2])
	print mutex_peak
	for match in mutex_peak.matches( protein, freq_match ):
		print 'Mutex: ',match
	exit(1)
	peak=Peak([0.3,0.5,1.0],noesy_hcn)
	peak2=Peak([0.3,0.5,1.0],noesy_hcn)
	peak2._hard_assignments=[
		(Atom('HA',1),None,Atom('N',2)),
		(Atom('QB',2),None,Atom('N',2)),
		(Atom('HA',5),None,Atom('N',2)),
		(None,None,Atom('N',4)),
		(Atom('HA',7),Atom('H',5),Atom('N',3))
]

	free_hhn_peak=Peak((1.0,0.5,0.3),noesy_hhn)
	hhn_peak=Peak((0.3,0.5,1.0),noesy_hhn)
	hhn_peak._hard_assignments=[
#		(Atom('HA',1),None,Atom('N',2)),
#		(Atom('QB',2),None,Atom('N',2)),
#		(Atom('HA',5),None,Atom('N',2)),
#		(None,None,Atom('N',7)),
#		(Atom('HA',1),Atom('H',1),None),
		(Atom('HA',2),Atom('HA',5),Atom('CA',5))
]
#	print 'get raw matches for hhn_peak'
#	for match in noesy_hhn.matches( free_hhn_peak, protein ):
#		print 'free match: ', match



	print 'get raw matches for hhn_peak'
	for match in free_hhn_peak.matches( protein, freq_match  ):
		print 'full match: ', match


	exit(1)
	print 'get matches for pre-assigned peak'
	for match in noesy_hcn.matches( peak2, protein):
		print 'full match: ', match

	all_matches=[match for match in noesy_hcn.matches( peak, protein) ]
	print 'get %d matches for generic peak'%(len(all_matches))
	print 'first 5 matches:\nfull_match: ','\nfull_match:  '.join(['%s,%s,%s,%s'%m for m in all_matches[0:5]])

#		print noesy_hcn.translate_full_match_to_peak_match(match)
	print 'print matched atoms...'
	for match in noesy_hcn.matches(peak,protein,freq_match):
		print 'FIN: ',noesy_hcn.translate_full_match_to_peak_match(match)
	print freq_match

	print 'pre-asigned stuff without frequency comparison...'
	for match in noesy_hcn.matches(peak2,protein):
		print 'FREE: ',noesy_hcn.translate_full_match_to_peak_match(match)

	freq_match.reset()
	print 'print matched atoms for pre-assigned peak2...'
	for match in noesy_hcn.matches(peak2,protein,freq_match):
		print 'FIN: ',noesy_hcn.translate_full_match_to_peak_match(match)
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

