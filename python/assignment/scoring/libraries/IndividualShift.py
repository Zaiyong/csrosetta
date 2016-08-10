#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from utility.GaussianDistribution import GaussianDistribution
from noe_tools import unpack_NMR_atom_group
from assignment.noesy.Atom import Atom

#COMMENT: what is the role of this class. What are the responsibilities of this class ? What is its scope?
class IndividualShift():
	def __init__(self):
		self._data={}
		Shift.__init__(self)
		self._seq=[]

	def by_atom(self,atom):
		aa=self._seq[atom.resid()-1]
		new_atom_names=unpack_NMR_atom_group(atom.name(),aa)
		new_atom=Atom(new_atom_names[0],atom.resid())
		return self._data[new_atom]

	def set_seq(self,seq):
		self._seq=seq

	def seq(self):
		return self._seq

	def add_distribution(self,atom,distribution):
    #use of private data not belonging to IndividualShift!!!
		self._data[atom]=distribution

	@classmethod
	def load_from_resonance_list(obj,resonance_list):
		obj=IndividualShift()
		seq=resonance_list.sequence()
		obj.set_seq(seq)
		for r in resonance_list.itervalues():
			aa=seq[r.resid()-1]
			new_atom_names=unpack_NMR_atom_group(r.name(),aa)
			atom=Atom(new_atom_names[0],r.resid())
			obj.add_distribution(atom,GaussianDistribution(r.freq(),r.error()))
		return obj
