#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from NoeStrip import NoeStrip
from Atom import Atom
from Resonance import RangeResonance, Resonance

from get_strips import read_lib
import math

def prob_strip_noise(label_freq,proton_freq,label_tol,proton_tol,label_bound,proton_bound):
	if label_freq<=label_bound[1] and label_freq>=label_bound[0]:
		label_prob=0.0
	elif label_freq>label_bound[1]:
		sigma2=label_tol*label_tol
		x=label_freq-label_bound[1]
		match_diff2=x*x/sigma2
		label_prob=1-math.exp(-match_diff2*0.5)
	elif label_freq<label_bound[0]:
		sigma2=label_tol*label_tol
		x=label_freq-label_bound[0]
		match_diff2=x*x/sigma2
		label_prob=1-math.exp(-match_diff2*5)

	if proton_freq<=proton_bound[1] and proton_freq>=proton_bound[0]:
		proton_prob=0.0
	elif proton_freq>proton_bound[1]:
		sigma2=proton_tol*proton_tol
		x=proton_freq-proton_bound[1]
		match_diff2=x*x/sigma2
		proton_prob=1-math.exp(-match_diff2*0.5)
	elif proton_freq<proton_bound[0]:
		sigma2=proton_tol*proton_tol
		x=proton_freq-proton_bound[0]
		match_diff2=x*x/sigma2
		proton_prob=1-math.exp(-match_diff2*5)

	prob_noise=math.sqrt((proton_prob*proton_prob+label_prob*label_prob)*0.5)
	return prob_noise

def prob_strip_miss(num):
	return math.exp(-num*num*0.1)

def prob_strip_off():
	return 0.0

def spinsystem_evaluation(ss,lib,aa):
	prob=[]
	for name, strip in ss._strips.iteritems():
		label_name=strip.label_name()
		proton_name=strip.proton_name()
		label_bound=read_lib(lib,aa,label_name)
		proton_bound=read_lib(lib,aa,proton_name)
		label_freq=strip.label().freq()
		proton_freq=strip.proton().freq()
		label_tol=strip.label().error()
		proton_tol=strip.proton().error()
		indirect_proton_num=len(strip.indirect_protons())
		prob_noise=prob_strip_noise(label_freq,proton_freq,label_tol,proton_tol,label_bound,proton_bound)
		prob_miss=prob_strip_miss(indirect_proton_num)
		prob_off=prob_strip_off()
		prob.append(math.sqrt(prob_noise*prob_noise+prob_miss*prob_miss+prob_off*prob_off))
	return math.fsum(prob)
