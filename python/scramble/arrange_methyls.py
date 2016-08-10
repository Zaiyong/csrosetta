#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os.path import basename
import os
import string
import StringIO
from assignment.noesy import Atom, Resonance

def arrange_methyls(resonancelist):
	seq=resonancelist.sequence()
	deleted_id=[]
	for i,aa in enumerate(seq):
		resid=i+1
		protons=[]
		try:
			res_list=resonancelist.by_residue(resid)
		except KeyError:
			continue
		for r in res_list:
			if r.atom().elem()=='H': protons.append(r.name())
		if aa=='A':
			if 'QB' in protons:
				continue
			elif 'HB1' in protons:
				id1=resonancelist.by_atom(Atom('HB1',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HB2',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HB3',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
#				new_res=Resonance(id=-1,atom=Atom('QB',resid),freq=resonancelist.by_atom(Atom('HB1',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
				resonancelist[id1].atom().set_name('QB')
			elif 'HB' in protons:
				id=resonancelist.by_atom(Atom('HB',resid)).id()
				resonancelist[id].atom().set_name('QB')
#				new_res=Resonance(id=-1,atom=Atom('QB',resid),freq=resonancelist.by_atom(Atom('HB',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
		elif aa=='I':
			if 'QG2' in protons:
				continue
			elif 'HG21' in protons:
				id1=resonancelist.by_atom(Atom('HG21',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HG22',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HG23',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
#				new_res=Resonance(id=-1,atom=Atom('QG2',resid),freq=resonancelist.by_atom(Atom('HG21',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
				resonancelist[id1].atom().set_name('QG2')

			elif 'HG2' in protons:
				id=resonancelist.by_atom(Atom('HG2',resid)).id()
				resonancelist[id].atom().set_name('QG2')
#				new_res=Resonance(id=-1,atom=Atom('QG2',resid),freq=resonancelist.by_atom(Atom('HG2',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
			if 'QD1' in protons:
				continue
			elif 'HD11' in protons:
				id1=resonancelist.by_atom(Atom('HD11',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HD12',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HD13',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
#				new_res=Resonance(id=-1,atom=Atom('QD1',resid),freq=resonancelist.by_atom(Atom('HD11',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
				resonancelist[id1].atom().set_name('QD1')
			elif 'HD1' in protons:
				id=resonancelist.by_atom(Atom('HD1',resid)).id()
				resonancelist[id].atom().set_name('QD1')
#				new_res=Resonance(id=-1,atom=Atom('QD1',resid),freq=resonancelist.by_atom(Atom('HD1',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
		elif aa=='L':
			if 'QD1' in protons:
				continue
			elif 'HD11' in protons:
				id1=resonancelist.by_atom(Atom('HD11',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HD12',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HD13',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
#				new_res=Resonance(id=-1,atom=Atom('QD1',resid),freq=resonancelist.by_atom(Atom('HD11',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
				resonancelist[id1].atom().set_name('QD1')
			elif 'HD1' in protons:
				id=resonancelist.by_atom(Atom('HD1',resid)).id()
				resonancelist[id].atom().set_name('QD1')
#				new_res=Resonance(id=-1,atom=Atom('QD1',resid),freq=resonancelist.by_atom(Atom('HD1',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
			if 'QD2' in protons:
				continue
			elif 'HD21' in protons:
				id1=resonancelist.by_atom(Atom('HD21',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HD22',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HD23',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
				resonancelist[id1].atom().set_name('QD2')
#				new_res=Resonance(id=-1,atom=Atom('QD2',resid),freq=resonancelist.by_atom(Atom('HD21',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
			elif 'HD2' in protons:
				id=resonancelist.by_atom(Atom('HD2',resid)).id()
				resonancelist[id].atom().set_name('QD2')
#				new_res=Resonance(id=-1,atom=Atom('QD2',resid),freq=resonancelist.by_atom(Atom('HD2',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]

		elif aa=='V':
			if 'QG1' in protons:
				continue
			elif 'HG11' in protons:
				id1=resonancelist.by_atom(Atom('HG11',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HG12',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HG13',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
				resonancelist[id1].atom().set_name('QG1')
#				new_res=Resonance(id=-1,atom=Atom('QG1',resid),freq=resonancelist.by_atom(Atom('HG11',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
			elif 'HG1' in protons:
				id=resonancelist.by_atom(Atom('HG1',resid)).id()
				resonancelist[id].atom().set_name('QG1')
#				new_res=Resonance(id=-1,atom=Atom('QG1',resid),freq=resonancelist.by_atom(Atom('HG1',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
			if 'QG2' in protons:
				continue
			elif 'HG21' in protons:
				id1=resonancelist.by_atom(Atom('HG21',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HG22',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HG23',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
				resonancelist[id1].atom().set_name('QG2')
#				new_res=Resonance(id=-1,atom=Atom('QG2',resid),freq=resonancelist.by_atom(Atom('HG21',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
			elif 'HG2' in protons:
				id=resonancelist.by_atom(Atom('HG2',resid)).id()
				resonancelist[id].atom().set_name('QG2')
#				new_res=Resonance(id=-1,atom=Atom('QG2',resid),freq=resonancelist.by_atom(Atom('HG2',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
		elif aa=='M':
			if 'QE' in protons:
				continue
			elif 'HE1' in protons:
				id1=resonancelist.by_atom(Atom('HE1',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HE2',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HE3',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
				resonancelist[id1].atom().set_name('QE')
#				new_res=Resonance(id=-1,atom=Atom('QE',resid),freq=resonancelist.by_atom(Atom('HE1',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
			elif 'HE' in protons:
				id=resonancelist.by_atom(Atom('HE',resid)).id()
				resonancelist[id].atom().set_name('QE')
#				new_res=Resonance(id=-1,atom=Atom('QE',resid),freq=resonancelist.by_atom(Atom('HE',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
		elif aa=='T':
			if 'QG2' in protons:
				continue
			elif 'HG21' in protons:
				id1=resonancelist.by_atom(Atom('HG21',resid)).id()
				try:
					id2=resonancelist.by_atom(Atom('HG22',resid)).id()
					deleted_id.append(id2)
				except KeyError:
					pass
				try:
					id3=resonancelist.by_atom(Atom('HG23',resid)).id()
					deleted_id.append(id3)
				except KeyError:
					pass
				resonancelist[id1].atom().set_name('QG2')
#				new_res=Resonance(id=-1,atom=Atom('QG2',resid),freq=resonancelist.by_atom(Atom('HG21',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id1]
			elif 'HG2' in protons:
				id=resonancelist.by_atom(Atom('HG2',resid)).id()
				resonancelist[id].atom().set_name('QG2')
#				new_res=Resonance(id=-1,atom=Atom('QG2',resid),freq=resonancelist.by_atom(Atom('HG2',resid)).freq(),error=0.02)
#				resonancelist.add_resonance(new_res)
#				del resonancelist[id]
	for id in deleted_id:
		del resonancelist[id]
	return resonancelist
