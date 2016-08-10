#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

class FixCalibrator:
	def __init__( self, factor=1e-36 ):
		self._factor=factor

	def eval( self, atom1, atom2, volume ):
		return self._factor*volume**6

