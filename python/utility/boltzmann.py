#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from math import exp
from random import random
def boltzmann(delta,T=1):
	#print 'delta is %f'%delta
	kb=43
	rand_value=random()
	prob=exp(-1.0*delta/(kb*T))
	#print 'prob is %f'%prob
	if rand_value<=prob:
		return True
	else:
		return False

