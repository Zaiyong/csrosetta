#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module

### module rules -- maybe split-out into files later in a /rules subdirectory
# Rules are subclasses of the BaseClass Rule
# the rule responsibilty is to know how a Peak is interpreted, i.e., what individual frequency dimensions mean,
# and how these map to covalent structure of the biomolecule: i.e., HNCA, NOESY, COSY, TOCSY...
#
# any subclass of Rule should generate a RuleNode for each matched atom (also invisble atoms). such as the proton in a CCH noesy.
from assignment.noesy.Resonance import UNFOLDED
from assignment.noesy.Resonance import FoldWindow
from chemical.Atom import Atom

from assignment import random_items
from assignment.PeakMatches import PeakMatch,AtomicPeakMatch

from basic import Tracer
tr=Tracer('assignment.rules.BasicRules')
# lookup of frequency property
#	def get_freq(self):
#		return self.peak_match.peak.freq(self.dim)
#freq=property(get_freq)


#Rules are considered const after instance creation
#A tuple of RuleNodes are kept to capture the matched atoms
#these define a sequence of matching atoms


def check_element_and_type( atom, node ):
	if atom.element != node.element: return False
	if not node.types: return True
	if not atom.types.isdisjoint(node.types): return True
	return False

class PeakRule(object):

	#RuleNode captures the information of a single matching node in a rule
	#most important is the dim field, which maps to the peak dimension or None if the matching node is invisible
	#for instance pseudo4D Noesy would have an invisible proton dimension. The proton dimension is still used for distance
	#computation but is not matched to a peak dimension
	class RuleNode:
		def __init__(self,dim,elem='H',tol=0.3,types=None,name=None,folder=UNFOLDED):
			self.element=elem
			self.dim=dim
			self.folder=folder
			self.tolerance=tol
			self.name=name
			self.types=types
		def __str__(self):
			return 'RuleNode(%(element)1s %(dim)4s)'%self.__dict__
  ### End RuleNode

	## class PeakRule ---
	def __init__(self, nodes):
		self.rule_dimension=len(nodes)
		self.nodes=tuple(nodes)
		self._update_dims()

	def __str__(self):
		str=self.__class__.__name__+'\n'
		str+='\n'.join('%s'%n for n in self.nodes)
		str+='\ndims: ('+', '.join('%4s'%d for d in self.dims)+')'
		return str

	#return a tuple of the peak's frequencies in rule-order
	def unpack_frequencies(self, peak):
		freq=[None]*self.rule_dimension
		for i,d in enumerate(self.dims):
			if d:
				freq[i]=peak.freq_of_dim(d)
		return tuple(freq)

	#update cached stuff from RuleNodes
	def _update_dims(self):
		dims=[]
		tols=[]
		folder=[]
		for i,node in enumerate(self.nodes):
			dims.append(node.dim)
			tols.append(node.tolerance)
			folder.append(node.folder)
		self.dims=tuple(dims)
		self.tolerances=tuple(tols)
		self.folder=tuple(folder)

	#a tuple of atoms in rule-order is re-ordered to peak-order
	#some columns might get lost, if the rule has invisible matches
	def translate_full_match_to_peak_match(self, match ):
		sorted_match=[None]*max(self.dims)
		for d,match in zip(self.dims,match):
			if d:
				sorted_match[d-1]=match
		return tuple(sorted_match)

	#a tupel of atoms in peak-format is translated to internal rule-format
	def translate_peak_match_to_full_match(self, match, atom_tree ):
		rule_match=tuple( match[d-1] if d else None for d in self.dims )
				#get the matching atom instances from the atom-tree
		pre_match=tuple( atom_tree.atom( atom ) if atom else None for atom in rule_match )
		for atom, node in zip(pre_match,self.nodes):
			if node and atom:
				if not check_element_and_type(atom,node):
					import library
					msg='atom %s in pre-existing assignment %s has wrong element or atom-type for %s'%(atom,match,node)
					raise library.InconsistentInput(msg)
		return pre_match

	#iterate over all tuples of atoms that match the general peak rule
	#frequency or distance matching is not happening
	#if hard-assignment exist these limit the results of iterate_raw
	def matches(self, peak, atom_tree, random_samples=None, frequency_matcher=None, distance_matcher=None, partial_match=None, match_mask=None ):
		tr.Trace('BasicRule::matches')
		if partial_match:
			raise NotImplementedError('partial matching is not implemented for PeakRule %s'%self.__class__.__name__)
		existing=peak.hard_assignments
		if len(existing):
			raise TypeError('Peak %s has pre-existing assignments but the Rule %s does not know how to treat these'%(peak,self))
	 	for match in random_items(
			self._iterate( peak, atom_tree, random_samples, frequency_matcher, distance_matcher, None )
			, random_samples):
			yield PeakMatch(match,self.translate_full_match_to_peak_match(match),peak)

	#macro-like method used to match a single rule-match position to a single frequeny
	def _match(self, freq, frequency_matcher, raw_match, offset=0, sub_match=None ):
		if not frequency_matcher: return True
		for i,match in enumerate(raw_match):
			if sub_match and sub_match[i+offset]: continue #if this is already matched it is always true
			if freq[i+offset] and not frequency_matcher(match,freq[i+offset],self.tolerances[i+offset],self.folder[i+offset]): return False
		return True

	def distance( peak_match, distance_matcher ):
		return None

	def distance_atoms(self, peak_match ):
		return ()
	#given tuples of distance-match atoms
	#how many and which peak-matches would I expect from this rule ?
	def expected_symmetric_peaks( atoms ):
		return None


class ExistingAssignmentRule(PeakRule):
	#returns tuples like [(True,True,False)] to reflect that a partial match of the 1st and 2nd dim of peak would be
	#a correlated spinsystem as e.g., (CA,HA)
	def spinsystem_match_masks(self):
		raise NotImplementedError('spin system information is not provided by PeakRule %s'%self.__class__.__name__)
	#a mask coming in peak-format is translated to internal rule-format
	def translate_mask_to_full_rule(self, mask ):
		#e.g., mask = [True, False, False ] --> [None, True, False, False]
		#for rule dimensions that are invisible we put None...
		return tuple( mask[d-1] if d else None for d in self.dims )

	def _combine(self, a, b):
		if not a: return b
		if not b: return a
		if a==b: return a
		raise ValueError('Incompatible assignments %s vs. %s'%(a,b) )

	def _combine_partial_matches(self, existing, partial_match ):
		try:
			partial_match=partial_match.peak_match
		except:
			pass #or not
		if len(existing)==0:
			if partial_match: yield partial_match
			else: return
		for existing_match in existing:
			if not partial_match:
				yield existing_match
			else:
				combined = ( self._combine(a,b) for a,b in zip(existing_match, partial_match ) )
				yield combined

	#iterate over all tuples of atoms that match the general peak rule
	#if hard-assignment exist, re-order these from peak-order into rule-order and use them to limit matches
	def matches(self, peak, atom_tree, random_samples=None, frequency_matcher=None, distance_matcher=None, partial_match=None, match_mask=None ):
		existing=peak.hard_assignments
		tr.Trace('ExistingAssignmentRule::matches')
		if match_mask: match_mask=self.translate_mask_to_full_rule( match_mask )
		if not len(existing) and not partial_match:
			for match in random_items(
				self._iterate( peak, atom_tree, random_samples, frequency_matcher, distance_matcher, None, match_mask )
				, random_samples):
				tr.Trace('ExistingAssignmentRule::matches::yield   ',PeakMatch(match,self.translate_full_match_to_peak_match(match),peak))
				yield PeakMatch(match,self.translate_full_match_to_peak_match(match),peak)
		else:
			for existing_match in self._combine_partial_matches( existing, partial_match ):
				#map from peak-dimension to rule dimension, for hidden rule-nodes the dim is None
				pre_match=self.translate_peak_match_to_full_match(existing_match,atom_tree)
				for match in random_items(
					self._iterate( peak, atom_tree, random_samples, frequency_matcher, distance_matcher, pre_match, match_mask ),
					random_samples ):
					tr.Trace('ExistingAssignmentRule::matches::yield  ',PeakMatch(match,self.translate_full_match_to_peak_match(match),peak))
					yield PeakMatch(match,self.translate_full_match_to_peak_match(match),peak)

	#helper function to split a submatch into a tuple of tuples
					#split_submatch((A,B,C,D), (3,1) ) splits to (A,B,C),(D,)
					#split_submatch((A,B,C,D), (2,2) ) splits to (A,B),(C,D)
					# if a tuple with just None is generated [e.g., (None,None)] the tuple is replaced by None
	@staticmethod
	def _split_submatch( match, splitting ):
		if not match: return [None]*len(splitting)
		matches=[]
		start=0
		for s in splitting:
			sp=match[start:(s+start)]
			if len(ExistingAssignmentRule._get_free_elements_in_match( sp ))==len(sp):
				sp=None
			matches.append( sp )
			start+=s
		return matches

	#given a submatch tubple which are the None entries ? (A,None) --> (2)
	#(A,None,None)-->(2,3)
	@staticmethod
	def _get_free_elements_in_match( match ):
		return tuple( i for i,m in enumerate(match) if m is None )

###########################
#
#   MutualExclusiveRule
#
#  this class is instantiated from a MutualExclusivePeak
#  the system is used for (small) groups of Peaks with pre-existing assignments,
#  but where each assignment choice can only be assigned to one peak (mutual exclusive)
#  example: we know w1 and w2 are the CA resonances of either HIS 73 or HIS 26
#     create two 1D Peak with w1 and w2 as frequency and both with CA 73 and CA 26 as pre-existing
#     assignments. Then create a MutualExclusivePeak from them.
#
class MutualExclusiveRule(ExistingAssignmentRule):
	def __init__(self, input_rule, dims ):
		#make copies of rule-nodes and update their dim field
		cum_peak_dim=0
		all_nodes=[]
		self.sub_rules=[]
		for peak_dim in dims:
			import copy
			new_rule=copy.deepcopy(input_rule)
			self.sub_rules.append(new_rule)
			nodes=new_rule.nodes
			new_dims=[]
			for n in nodes:
				if n.dim:
					n.dim+=cum_peak_dim
				new_dims.append(n.dim)
			new_rule.dims=new_dims
			all_nodes+=nodes
			cum_peak_dim+=peak_dim
		#use this list of nodes to initialize Super class
		super(MutualExclusiveRule,self).__init__(all_nodes)
		self.sub_match_splitter=(input_rule.rule_dimension,)*len(dims)
		self._update_dims()
#		print 'MutexRule Created..'
#		for i,rule in enumerate(self.sub_rules):
#			print 'Mutex %d\n'%(i+1),rule


	def _update_dims(self):
		for rule in self.sub_rules:
			rule._update_dims()
		super(MutualExclusiveRule,self)._update_dims()
	#iterate over all tuples of atoms that match the general peak rule
	#frequency or distance matching is not happening
	#if hard-assignment exist these limit the results of iterate_raw
	def matches(self, peak, atom_tree, random_samples, frequency_matcher=None, distance_matcher=None ):
		existing=peak.hard_assignments
		if not len(existing):
			raise TypeError('Peak %s has no pre-existing assignments, which is weird for a MutexRule'%(peak,self))

		for existing_match in existing:
			rule_match=tuple( existing_match[d-1] if d else None for d in self.dims )
			pre_match=tuple( atom_tree.atom( atom ) if atom else None for atom in rule_match )
			splitted_match=self._split_submatch( pre_match, self.sub_match_splitter )
			for match in random_items( self._recursive_matches( peak, atom_tree, frequency_matcher, distance_matcher, splitted_match, self.sub_rules), random_smaples ):
				yield PeakMatch(match,self.translate_full_match_to_peak_match(match),peak)

	@staticmethod
	def _recursive_matches( peak, atom_tree, random_samples, freq_match, dist_match, splitted_match, rules ):
		head=splitted_match[0]
		tail=splitted_match[1:]
		head_rule=rules[0]
		tail_rules=rules[1:]
#		print 'call head_rule on ', head
		for match in head_rule._iterate( peak, atom_tree, random_samples, freq_match, dist_match, head ):
			if len(tail):
				for tail_match in MutualExclusiveRule._recursive_matches( peak, atom_tree, random_samples, freq_match, dist_match, tail, tail_rules ):
					yield match+tail_match
			else: yield match
