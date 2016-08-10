#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#Distance distribution library generated from fragments

from chemical import Atom
from math import sqrt,fabs
from utility import UniformHistogram
from noe_tools import unpack_NMR_atom_group
import library
from scipy import linspace

def switch(a,b):
	return b,a

class FragDistanceData:
	def __init__(self):
		self._data={}

	def by_atom_pair(self,atom1,atom2):
		if atom1>atom2:
			atom1,atom2=switch(atom1,atom2)
#		new_atom1=Atom(unpack_NMR_atom_group(atom1.name(),self._seq[atom1.resid()-1])[0],atom1.resid())
#		new_atom2=Atom(unpack_NMR_atom_group(atom2.name(),self._seq[atom2.resid()-1])[0],atom2.resid())
		return self._data[atom1,atom2]

	def add_distribution(self,atom1,atom2,distribution):
		self._data[atom1,atom2]=distribution

	@classmethod
	def get_from_fragdist_table(obj,file):
		#table is the distance distribution file like:
		#     1    H    1     HB2     2.171     3.609       9       4      31      34      30      38      14       9       7       4       0       0       0       0       0       1       0       0       0      19
		#the first 6 variables are atom1 name, atom1 resid, atom2 name,atom2 resid, lowest value,highest value
		#the rest are hists.
		obj=FragDistanceData()
		sequence=''
		vars={}
		atoms_dist_list={}
		for line in file:
			tags=line.split()
			if len(tags)<6 and ( 'DATA' not in line ): continue
			if tags[0]=='DATA' and tags[1]=='SEQUENCE':
				for tag in tags[2:]:
					sequence=sequence+tag
			elif tags[0]=='VARS':
				for i,tag in enumerate(tags[1:]):
					vars[tag]=i
			else:
				atom1=Atom(tags[1],int(tags[0]))
				atom2=Atom(tags[3],int(tags[2]))
				if atom1>atom2: atom1,atom2 = switch(atom1,atom2)
				low=float(tags[4])
				high=float(tags[5])
				num_bin=len(tags)-6
				#try allocating memory at once
				table=range(0,num_bin)
				sum_bin=0
				for i,bin in enumerate(tags[6:]):
					table[i]=float(bin)

				if (high-low<0.001):
					high=low+0.1
				dist_hist=UniformHistogram(table,low,high)
				obj.add_distribution(atom1, atom2, dist_hist )
#				atoms1=unpack_NMR_atom_group(atom1.name(),sequence[atom1.resid()-1])
#				atoms2=unpack_NMR_atom_group(atom2.name(),sequence[atom2.resid()-1])
#				for a1 in atoms1:
#					for a2 in atoms2:
						#in the distance distribution file, there are atom combination like QQG, QG2, etc.
#						obj.add_distribution(Atom(a1,atom1.resid()),Atom(a2,atom2.resid()),dist_hist)
		obj._sequence=sequence
		return obj

