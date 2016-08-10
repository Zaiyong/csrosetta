#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from math import erf,fabs,sqrt,log

ERF_CUTOFF=1e-5
def quantity_function(input):
	return log(max(ERF_CUTOFF,1-erf(fabs(input)/sqrt(2))))

def schmidt_final_score(q,x0):
	return 1+q/fabs(quantity_function(x0))
