#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#import options
from ExampleArgumentParser import ExampleArgumentParser
from _module_variables import _module_parsers

#instantiate a parser from this class for you main-program
# any options from imported modules will be present in the help as well
class ApplicationArgumentParser(ExampleArgumentParser):

	def __init__(self, **kwargs):
		parents=_module_parsers
		kwargs['parents']=parents
		ExampleArgumentParser.__init__(self, conflict_handler='resolve',**kwargs)
		self._args = None

	def parse_args(self, redo=False ):
		#for speed reasons
		if self._args and not redo:
			return self._args
		self._args=super(ApplicationArgumentParser,self).parse_args()

#		for m in options._module_parsers:
		for m in _module_parsers:
			m._args=self._args
		return self._args
