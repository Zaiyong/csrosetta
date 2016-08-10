#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from AtomWise import AtomWise
from noe_tools import unpack_NMR_atom_group
import library
from utility.SchmidtQ import schmidt_final_score

class Shift(AtomWise):

	def __init__(self,name,lib):
		AtomWise.__init__(self,name)
		self._shift_distribution_library=lib
		self._x0=1.5
	def score_atom(self,freq,atom):
		distribution=self._shift_distribution_library.by_atom(atom)
		try:
			return schmidt_final_score(distribution.schmidt_q(freq),self._x0)

		except ValueError:

 			mean=distribution.mean()
 			std=distribution.std()
 			raise library.RunException('shift of atom %s should be %f wtih std %f but it is assigned to %f '%(atom,mean,std,freq))

	def shift_distribution_library(self):
		return self._shift_distribution_library
