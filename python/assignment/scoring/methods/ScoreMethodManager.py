#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#from AtomWise import AtomWise
#from AtomPairWise import AtomPairWise
#from AssignmentWise import AssignmentWise

#module as singleton class ScoreMethodManager:
score_methods={}
#	@staticmethod
def register(name,score_method):
	score_methods[name]=score_method

#		print ScoreMethodManager.score_methods
#	@staticmethod
automatic_methods={}
def generate_method(name):
	return automatic_methods.setdefault(name,score_methods[name]())

