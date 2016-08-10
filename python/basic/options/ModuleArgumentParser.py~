#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import argparse
from _module_variables import _module_parsers, _unit_testing

#a parser that is responsible for a module should be instantiated from this class
# all options will automatically appear as an option group in the options of ApplicationArgumentParser
# careful: if options conflict the later one will overwrite the earlier one.
# calling parse_args() on a moduls parser should only be done after the global options have been parsed
# however, parse_args() output is cashed, such that a modules parse_args() can be called even in fast code
# Limitation: argparse cannot nest argument groups. So you cannot add groups of arguments to a ModuleParser
# since it is already a sub-group... I could make the automatic sub-grouping optional...
class ModuleArgumentParser(argparse.ArgumentParser):
	def __init__(self,module,**kwargs):
		self.prefix=kwargs.get('prefix')
		if self.prefix:
			del kwargs['prefix']
		argparse.ArgumentParser.__init__(self, **kwargs)

		#want to group the options of each module
		self._groups = argparse.ArgumentParser.add_argument_group(self,title=module,description=self.description)

		#make sure ApplicationArgumentParser will know about us:
		_module_parsers.append(self)

		#for cacheing the parsed options later
		self._my_args = None

	def add_argument(self, *args, **kwargs):
		if self.prefix:
			if args[0][0]=='-':
				args=('-'+self.prefix+'_'+args[0][1:],)+args[1:]
			else:
				args=(self.prefix+'_'+args[0],)+args[1:]
		return self._groups.add_argument(*args,**kwargs)

	def set_testing_args( self, **kwargs ):
		self._my_args=argparse.Namespace()
		for key, value in self._groups._option_string_actions.iteritems():
			key=key.lstrip('-')
			if self.prefix:
				my_key=key.replace(self.prefix+'_','')
			else:
				my_key=key
			if my_key in kwargs:
				setattr(self._my_args,my_key,kwargs[my_key])
			else:
				setattr(self._my_args,my_key,value.default)

	def parse_args(self, redo=False ):
		if not self._my_args and 'UNIT_TEST' in _unit_testing:
			self.set_testing_args()
		if self._my_args and not redo:
			return self._my_args
		if not self._args: super(ModuleArgumentParser,self).parse_args()
		self._my_args=argparse.Namespace()
		for key in self._groups._option_string_actions.iterkeys():
			key=key.lstrip('-')
			if self.prefix:
				my_key=key.replace(self.prefix+'_','')
			else:
				my_key=key
			setattr(self._my_args,my_key,getattr(self._args,key))
		return self._my_args

