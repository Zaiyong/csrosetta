#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from NoeStrip import NoeStrip
from Atom import Atom
from Resonance import RangeResonance, Resonance
from PDB.Polypeptide import one_to_three

def read_lib(lib,aa,name):
	bound=[]
	aa_three=one_to_three(aa)
	for r in lib:
		tags=r.split()
		if aa_three==tags[0] and name==tags[1]:
			bound=[float(tags[8]),float(tags[9])]
	assert bound, "aa name or atom name is wrong"
	return bound


def create_single_strip_family_2D(resid,aa,lib,L_name,H_name,crosspeaks):
	L_atom=Atom(L_name,resid)
	H_atom=Atom(H_name,resid)
	L_cs_bound=read_lib(lib,aa,L_name)
	H_cs_bound=read_lib(lib,aa,H_name)
	L_range=RangeResonance(-1,L_atom,L_cs_bound[0],L_cs_bound[1],0.4)
	H_range=RangeResonance(-1,H_atom,H_cs_bound[0],H_cs_bound[1],0.04)
	strip=NoeStrip.generate_strip_family_2D(crosspeaks,1,H_range,L_range)
	return strip

def heavy_atom_based_strips(heavy_atoms,aa,resid,lib,resonance_list,crosspeaks):
	strips=[]
	for r in heavy_atoms:
		if r=='N':
			N_atom=Atom('N',resid)
			H_atom=Atom('H',resid)
			N_resonance=resonance_list.by_atom( N_atom )
			H_resonance=resonance_list.by_atom( H_atom )
			strips.append(NoeStrip.generate_strip(crosspeaks,1,H_resonance,N_resonance))
		elif r=='CA':
			CA_atom=Atom('CA',resid)
			HA_atom=Atom('HA',resid)
			CA_resonance=resonance_list.by_atom( CA_atom )
			HA_range=RangeResonance(-1,HA_atom,3,6,0.03)
			new_ss=[]
			for g in NoeStrip.generate_strip_family(crosspeaks,1,HA_range,CA_resonance):
				new_ss.append([strips[0],g])
			strips=new_ss
		elif r=='CB':
			CB_atom=Atom('CB',resid)
			HB_atom=Atom('HB',resid)
			CB_resonance=resonance_list.by_atom( CB_atom )
			HB_range=RangeResonance(-1,HB_atom,0,4,0.04)
			if not resonance_list.by_atom( HB_atom ):## sometimes, in prot file, it is HB1 instead of HB
				HB_atom=Atom('HB2',resid)
				HB_range=RangeResonance(-1,HB_atom,0,4,0.04)
			new_ss=[]
			for r in strips:
				for g in NoeStrip.generate_strip_family(crosspeaks,1,HB_range,CB_resonance):
					new_ss.append(r+[g])
			strips=new_ss
		elif r=='CG':
			if aa in 'ANDCGHFSWY':
				continue
			elif aa in 'RQEKMP':
				strip_list=create_single_strip_family_2D(resid,aa,lib,'CG','HG2',crosspeaks)
			elif aa=='I':
				strip_list=create_single_strip_family_2D(resid,aa,lib,'CG1','HG12',crosspeaks)
				strip_list=create_single_strip_family_2D(resid,aa,lib,'CG2','HG2',crosspeaks)
			elif aa=='L':
				strip_list=create_single_strip_family_2D(resid,aa,lib,'CG','HG',crosspeaks)
			elif aa=='T':
				strip_list=create_single_strip_family_2D(resid,aa,lib,'CG2','HG2',crosspeaks)
			elif aa=='V':
				strip_list=create_single_strip_family_2D(resid,aa,lib,'CG1','HG1',crosspeaks)
			new_ss=[]
			for r in strips:
				for g in strip_list:
					new_ss.append(r+[g])
			strips=new_ss
		else: continue
	return strips

def get_strips( resid,seq,lib,resonance_list,crosspeaks ):
	aa=seq[resid-1]
	aa_map='ARNDCQEGHILKMFPSTWYV'
	aa_heavy_map=[
		['N','CA','CB'],
		['N','CA','CB','CG','CD','CZ','NE'],
		['N','CA','CB','CG','ND2'],
		['N','CA','CB','CG'],
		['N','CA','CB','CG'],
		['N','CA','CB','CG','CD','NE2'],
		['N','CA','CB','CG','CD'],
		['N','CA'],
		['N','CA','CB','CG','CD2','CE1','ND1','NE2'],
		['N','CA','CB','CG1','CG2','CD1'],
		['N','CA','CB','CG','CD1','CD2'],
		['N','CA','CB','CG','CD','CE','NZ'],
		['N','CA','CB','CG','CD','CE','NZ'],
		['N','CA','CB','CG','CD1','CD2','CE1','CE2','CZ'],
		['N','CA','CB','CG1','CD'],
		['N','CA','CB'],
		['N','CA','CB','CG2'],
		['N','CA','CB','CG','CD1','CD2','CE2','CE3','CZ2','CZ3','CH2','NE1'],
		['N','CA','CB','CG','CD1','CD2','CE1','CE2','CZ'],
		['N','CA','CB','CG1','CG2']]
	strips=heavy_atom_based_strips(aa_heavy_map[aa_map.find(aa)],aa,resid,lib,resonance_list,crosspeaks)
	return strips



