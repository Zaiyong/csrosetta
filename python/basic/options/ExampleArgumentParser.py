## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-


import argparse
import textwrap
import library

class ExampleArgumentParser(argparse.ArgumentParser):
	def __init__(self, **kwargs):
		if 'examples' in kwargs:
			self.examples=kwargs['examples']
			del kwargs['examples']
		else:
			self.examples=None
#		print kwargs
		if 'aliases' in kwargs:
			self.aliases=kwargs['aliases']
			del kwargs['aliases']
		else:
			self.aliases=None
		argparse.ArgumentParser.__init__(self, **kwargs)

	def format_help(self):
#		return argparse.ArgumentParser.format_help(self)
		help_str=argparse.ArgumentParser.format_help(self)
		if self.examples:
			help_str=help_str+"\nexamples:"
			for ex in self.examples:
				if library.obj_is_list(ex):
					help_str=help_str+"\n\n  "+ex[0]
					help_str=help_str+"\n        "
					help_str=help_str+"\n        ".join(textwrap.wrap(ex[1],100))
				else:
					help_str=help_str+"\n  "+ex+""
			help_str=help_str+"\n\n"

		if self.aliases:
			self.aliases.remove(self.prog.replace('.py',''))
		if self.aliases and len(self.aliases)>0:
			if len(self.aliases)>1:
				s='s'
				es='es'
			else:
				s=''
				es=''
			help_str=help_str+("\nalias%s:\nThe same application can be called using the name%s: ")%(es,s)+" ".join(self.aliases)+"\n"

		return help_str%dict(prog=self.prog)



