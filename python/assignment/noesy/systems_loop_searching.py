#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os.path import basename
import os
import string
import sys
from assignment.noesy import SpinSystem
import fasta
def aa_strip_couple(aa):
	if aa in 'RQEKMPLT':
		couple=3
	elif aa in 'IV':
		couple=4
	elif aa=='G':
		couple=1
	else:
		couple=2
	return couple

def systems_loop_searching(resid,seq,resonance_list,crosspeaks,lib,systems):
	final_systems=[]

	if len(systems)==0:
		return final_systems
	elif len(systems)==1:
		final_systems=systems[0]
	elif len(systems)>1:
		max_long_connection=max([ r.long_connection() for r in systems])
		mid_systems=[]
		for r in systems:
			if r.long_connection()==max_long_connection:
				mid_systems.append(r)
		if len(mid_systems)==1:
			final_systems=mid_systems[0]
		elif len(mid_systems)>1:
			dim=len(mid_systems[0].strips())*2
			length=len(mid_systems)
			centroid={}
			for name,st in mid_systems[0].strips().iteritems():
				proton_name=name
				label_name=st.label().name()
				centroid[label_name+'H']=[]
				centroid[label_name]=[]
			for sp in mid_systems:
				for name,st in sp.strips().iteritems():
					proton_name=name
					label_name=st.label().name()
					centroid[label_name+'H'].append(st.proton().freq())
					centroid[label_name].append(st.label().freq())
			initial_distance=0.0
			initial_ss_index=0
			for name,value in centroid.iteritems():
				initial_distance+=value[0]**2
			for i in range(0,length):
				mid_distance=0.0
				for name,value in centroid.iteritems():
					mid_distance+=value[i]**2
				if mid_distance<initial_distance:
					initial_ss_index=i
					initial_distance=mid_distance
			final_systems=mid_systems[initial_ss_index]

	return final_systems
