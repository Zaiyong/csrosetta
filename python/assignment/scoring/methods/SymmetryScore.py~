#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#score method based on fragment distance distribution

import math
import itertools
import copy
from basic import options
from basic import Tracer
from assignment.scoring import score_definitions
from assignment.scoring import libraries
from ScoreMethod import ScoreMethod,ScoreOptions

SCORE_TYPE='symmetry'
#parser = options.ModuleArgumentParser("%s-score"%SCORE_TYPE,
#																			prefix=SCORE_TYPE,
#																			description='options for scoring against distance of structures',
#																			add_help=False )
#parser.add_argument("-dist",help="distance table generated from ", default=None )
#parser.add_argument("-norm",help='expect x-percent of frag_dist below noesy_cut_off [default 0.3]', type=float, default=0.3 )


tr=Tracer('scoring.%s'%SCORE_TYPE)

class SymmetryScoreOptions(ScoreOptions):
	def __init__(self):
		pass
#######################
##
## class SymmetryScore
##
## reads a dist-table generated with FragsToAtomDist.default.linuxgccrelease
## NOESY assignments to atom-pairs with distances small according to fragments
## get favorable score
## we evaluate the integral of distance distribution above the noesy cutoff
## if almost all distances are above the noesy cutoff it is unlikely that this assignment is correct
## however: what about long-range assignments ?
## here we have no data and return a neutral 0
##
## TODO: should maybe have some smooth weighting scheme that takes into account that information on 1-9 distances is less accurate
## than information on 1-1 and 1-2 distances
class SymmetryScore(ScoreMethod):
	class SymmetricPeaks(object): #or ScoreItem
#		__slots__='ct_expected','ct_sym','sym_assignments','pair'
		def __init__(self,assignment=None):
			if not assignment: return
			self.ct_expected=1 #because we already have number 1 of the group
			self.ct_sym=1
#			print 'INIT: ',assignment
			self.sym_assignments=set() #in constructor it would assign elements from iterable, we want to avoid that here
			self.sym_assignments.add(assignment)
#			self.pair=pair

		def __copy__(self):
#			raise TypeError('cannot copy SymmetricPeaks')

		# if we copy these the multiple reference network gets lost
			newone=type(self)()
			newone.ct_expected=self.ct_expected
			newone.ct_sym=self.ct_sym
			newone.sym_assignments=copy.copy(self.sym_assignments)
			return newone

		def __str__(self):
			s=''
			self.my_score=self.score
			s+='SP: %(ct_expected)5d %(ct_sym)5d %(my_score)5.3f'%self.__dict__
			s+='\n'+'\n'.join([str((w.peak.peak_list_name,w.peak.id,w.peak_match)) for w in self.sym_assignments])
			return s

		@property
		def score(self):
#			print self.ct_expected, self.ct_sym, float(self.ct_sym)/self.ct_expected
			if not self.ct_expected>1: return 0
			return min(float(self.ct_sym-1)/(self.ct_expected-1),1)

		def add(self,match):
#			from assignment.PeakMatches import PeakMatch
#			if not isinstance(match,PeakMatch):
#				raise TypeError(match)
			self.sym_assignments.add(match)
			self.ct_sym+=1

		def remove(self,match):
			self.ct_sym-=1
			try:
				self.sym_assignments.remove(match)
			except KeyError as exc:
				print 'cannot remove match: ', match
				print 'from sym-group', self
				exit(1)

	def __init__(self,options=None):
		if not options:
			options=SymmetryScoreOptions()
		#obtain raw data, either from cmd-line or as list of stuff
		super(SymmetryScore,self).__init__(SCORE_TYPE,options)

	def score(self,assignments):
		super(SymmetryScore,self).score(assignments)
		total_score=0
		total_cache_nr=0
		total_new_nr=0
		for assignment in assignments:
			try:
				total_score+=assignments.scores.assignments[assignment][self.name].score
				total_cache_nr+=1
			except KeyError:
				symmetric_peaks,correction=self.find_symmetric_assignments(assignments,assignment)
				assignments.scores.assignments.setdefault(assignment,{})[self.name]=symmetric_peaks
				total_new_nr+=1
				total_score+=symmetric_peaks.score+correction
				#end KeyError: assignment was not found in score-cache
		tr.Debug('total_cached: %5d total_new: %5d'%(total_cache_nr, total_new_nr))
		return total_score

	def notify_add(self, assignments, new_assignment):
		symmetric_peaks, correction=self.find_symmetric_assignments(assignments,new_assignment)
		assignments.scores.assignments.setdefault(new_assignment,{})[self.name]=symmetric_peaks

	def notify_remove(self, scores, old_assignment):
		#remove assignment from group of symmetric assignments and delete the entry for this one
		cache=scores.assignments[old_assignment]
		cache.pop(self.name).remove(old_assignment)

	def copy_cache( self, new_cache, old_cache ):
		for key, old_values in old_cache.assignments.iteritems():
			current=new_cache.assignments.setdefault(key,{})
			symmetric_peaks=current.get(self.name,None)
			if not symmetric_peaks:
				symmetric_peaks=copy.copy(old_values[self.name])
				for assignment in symmetric_peaks.sym_assignments:
					new_cache.assignments.setdefault(assignment,{})[self.name]=symmetric_peaks
				current[self.name]=symmetric_peaks

	def find_symmetric_assignments(self,assignments,assignment):
		pair=assignment.distance_atoms()
		if len(pair)>1: raise NotImplementedError('Currently only support for rules with a single distance connection')
		pair=pair[0]
#		print 'try to find symmetric peaks for ', pair, assignment
		symmetric_peaks=SymmetryScore.SymmetricPeaks(assignment)
		tr.Debug('')
		tr.Debug('for:   %75s'%str(assignment).strip())
		tr.Debug('distance atoms: %s'%str(pair))
		for experiment in assignments.peak_collection.experiments.itervalues():
			rule=experiment[0].rule #get rule for this peak-list, assuming rules are the same for all peaks in list
			#rule.expected_symmetric_peaks yields both (atom0, atom1) and (atom1, atom0) matches where pair=(atom0,atom1)
			for expected in rule.expected_symmetric_peaks(assignments.molecule, pair):
				if assignment.peak_match==expected: continue
				tr.Debug('expect:%75s'%str(expected).strip())
				symmetric_peaks.ct_expected+=1
				try:
					matches=assignments[expected]
					for match in matches:
						if match.peak.peak_list_name==experiment.name: #make sure we count only the matches that belong to the current experiment
							tr.Debug('found: %75s'%str(match).strip())
							try: #okay that was indeed a symmetric match, let's see if it is already a known group
								already_symm_peaks=assignments.scores.assignments[match][self.name]
								tr.Debug('found an existing symmetry group ')
								tr.Debug('for:   %75s'%str(assignment).strip())
								tr.Debug('expect:%75s'%str(expected).strip())
								tr.Debug('found: %75s'%str(match).strip())
								tr.Debug('already: %75s'%str(already_symm_peaks))
								old_score=already_symm_peaks.score
								already_symm_peaks.add(assignment)
								score_correction=already_symm_peaks.score-old_score
								return already_symm_peaks, score_correction
							except KeyError:#this match didn't belong to the group yet, add it now
								symmetric_peaks.add(match)
								assignments.scores.assignments.setdefault(match,{})[self.name]=symmetric_peaks
				except KeyError:	#expected assignment not found
					pass #do not put code here, as we don't get here, if matches were obtained which didn't match the peak_list name
		return symmetric_peaks, 0 #no correction




import ScoreMethodManager
ScoreMethodManager.register(SCORE_TYPE,SymmetryScore)


import unittest
