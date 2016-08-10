#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module

from BasicRules import *

class HN_HSQC_Rule(ExistingAssignmentRule):
	def __init__(self,dims,tols):
		ExistingAssignmentRule.__init__(self,(self.RuleNode(dims[0],'N',tols[0],'N'),self.RuleNode(dims[1],'H',tols[1],'H')))
#		self._dims=dims

	def _iterate(self, peak, atom_tree, random_samples, frequency_matcher=None, distance_matcher=None, sub_match=None ):
		freq=self.unpack_frequencies(peak)
		if sub_match:
			free=self._get_free_elements_in_match( sub_match )
			full_match=[None]*self.rule_dimension
			if len(free)==1 and len(free)>0:
				non_free=1-free[0]
				fixed=sub_match[non_free]
				full_match[non_free]=fixed
				search_name=self._nodes[free[0]].name
				full_match[free[0]]=atom_tree.atom(Atom(search_name,fixed.resid()))
#				print 'full: ',full_match, '  sub: ', sub_match
				if self._match( freq, frequency_matcher, full_match, 0, sub_match ):
					yield tuple(full_match)
			elif len(free)==0:
				yield sub_match
		else:
			test_range=range( atom_tree.first_residue, atom_tree.last_residue+1)
			if random_samples:
				test_range=random_items(test_range,random_samples)
			for resid in test_range:
				try:
					nh_match=(atom_tree.atom(Atom('N',resid)), atom_tree.atom(Atom('H',resid)))
					if not self._match( freq, frequency_matcher, nh_match ): continue
					yield nh_match
				except KeyError:
					pass


class HNCA_Rule(ExistingAssignmentRule):

	def __init__(self,dims,tols):
		self._nh_rule=HN_HSQC_Rule(dims[0:2],tols[0:2])
		ExistingAssignmentRule.__init__(self,self._nh_rule._nodes+(self.RuleNode(dims[2],'C',tols[2]),))

	def _iterate(self, peak, atom_tree, random_samples, frequency_matcher=None, distance_matcher=None, sub_match=None ):
		freq=self.unpack_frequencies(peak)
		splitted_match=self._split_submatch( sub_match, (2,1) )
		if splitted_match[1] and not splitted_match[0]:
			#CA is defined but none of the NHs
			# go the inverse way
			ca_match=splitted_match[1]
			ca_resid=ca_match[0].resid()
			for resid in range(ca_resid-1,ca_resid+1):
				nh_match=(atom_tree.atom(Atom('N',resid)), atom_tree.atom(Atom('H',resid)))
				if not self._match( freq, frequency_matcher, nh_match ): continue
				yield nh_match+ca_match
		else:
			for nh_match in self._nh_rule._iterate( peak, atom_tree, random_samples, frequency_matcher, distance_matcher, splitted_match[0] ):
				resid=nh_match[0].resid()
				if splitted_match[1]:
					ca_match=splitted_match[1]
					if not ca_match[0].resid()==resid and not (ca_match[0].resid()-1)==resid:
						raise KeyError('inconsistent pre-assignment, no match possible %d-%d'%(resid,ca_match[0].resid()))
					yield nh_match+ca_match
				else: #no pre-determined match for the CA part
					try:
						ca_match=(atom_tree.atom(Atom('CA',resid)),)
						if self._match( freq, frequency_matcher, ca_match, offset=self._nh_rule.rule_dimension ): yield nh_match+ca_match
					except KeyError:
						pass
					try:
						ca_match=(atom_tree.atom(Atom('CA',resid+1)),)
						if self._match( freq, frequency_matcher, ca_match, offset=self._nh_rule.rule_dimension ): yield nh_match+ca_match
					except KeyError:
						pass
				pass
