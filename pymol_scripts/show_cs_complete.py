#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from pymol import cmd
import pymol

def color_bar(num):
	assert num<=1.0,'the percentage is more than 100%'
	if num<=0.1:
		return 'red'
	elif num<=0.2:
		return 'orange'
	elif num<=0.3:
		return 'brightorange'
	elif num<=0.4:
		return 'yellow'
	elif num<=0.5:
		return 'limon'
	elif num<=0.6:
		return 'green'
	elif num<=0.7:
		return 'cyan'
	elif num<=0.8:
		return 'marine'
	elif num<=0.9:
		return 'blue'
	elif num<=1.0:
		return 'magenta'

def show_cs_complete( cs_complete_file ):

	cmd.show_as('cartoon')
	cmd.color('gray','all')
	for line in open(cs_complete_file,'r'):
		tags=line.split()
		level=color_bar(float(tags[3]))
		cmd.color(level,'resi '+tags[1])

cmd.extend("show_cs_complete",show_cs_complete)
