#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

def output():
	for i in range(0,10):
		yield i

def output2():
	for y in output():
		yield y

for x in output2():
	print x
