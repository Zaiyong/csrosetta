#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import argparse
import sys
import os
import structure_analyse
from scipy import linspace
from math import fabs
def Ramachandran_plot_score(phi_psi_file,rigid_file):
	phi_psi_pair=[]
	file_list=open(phi_psi_file,'r').readlines()
	rigid_residues=[]
	sel_2_list=open(rigid_file,'r').readlines()
	for line in sel_2_list:
		tags=line.split()
		if tags[0]=='RIGID':
			for  r in range(int(tags[1]),int(tags[2])+1):
				rigid_residues.append(r)

	for line in file_list:
		tags=line.split()
		if int(tags[1]) in rigid_residues:
			phi_psi_pair.append([float(tags[2]),float(tags[3])])

	Ramachandran_plot_file=os.environ['csrosettaDir']+"/database/Ramachandran_plot.txt"
	phi_dir=linspace(-180,180,361)
	psi_dir=linspace(-180,180,361)
	map={}
	R_plot_matrix=[]
	Ramachandran_plot=open(Ramachandran_plot_file,'r').readlines()
	for line in Ramachandran_plot:
		tags=line.split()
		R_plot_matrix.append(tags)
	for phi in phi_dir:
		for psi in psi_dir:
			map[phi,psi]=int( R_plot_matrix[int(phi)+180][int(psi)+180] )
	score=0
	for pair in phi_psi_pair:
		phi=round(pair[0])
		psi=round(pair[1])
		if fabs(phi)>180: phi=phi-360*phi/fabs(phi)
		if fabs(psi)>180: psi=psi-360*psi/fabs(psi)
		score+=map[phi,psi]
	score=1.0*score/(2*len(phi_psi_pair))
	return score

def beta_sheet_content(top_file):
	file_list=open(top_file,'r').readlines()
	sheet_content=0
	for line in file_list:
		if 'PAIRSTAT_ENTRY' not in line: continue
		tags=line.split()
		tags.reverse()
		for r in tags:
			if r != 'pleating:':
				sheet_content+=2
			else:
				break
	return sheet_content

def backbone_contact_order(top_file):
	file_list=open(top_file,'r').readlines()
	max_order=0
	for line in file_list:
		if 'PAIRSTAT_ENTRY' not in line: continue
		tags=line.split()
		for r in tags:
			if '-' in r:
				pair=r.split('-')
				max_order=max(max_order,int(pair[1])-int(pair[0]))
	return max_order

def sidechain_contact_order(contact_sidechain_file):
	file_list=open(contact_sidechain_file,'r').readlines()
	max_order=0
	for line in file_list:
		tags=line.split()
		max_order=max(max_order,int(tags[1])-int(tags[0]))
	return max_order

parser = argparse.ArgumentParser( description='structure evaluation', add_help=True)
parser.add_argument("-torsion", help="phi_psi.txt",default=None);
parser.add_argument("-rigid", help="rigid file",default=None);
parser.add_argument("-top", help="top file",default=None);
parser.add_argument("-sc_contact", help="sidechain contact order file",default=None);
args = parser.parse_args()

ramachandran_plot_score=structure_analyse.Ramachandran_plot_score(args.torsion,args.rigid)
beta_sheet_content=structure_analyse.beta_sheet_content(args.top)
backbone_contact_order=structure_analyse.backbone_contact_order(args.top)
sidechain_contact_order=structure_analyse.sidechain_contact_order(args.sc_contact)

print ramachandran_plot_score,beta_sheet_content,backbone_contact_order,sidechain_contact_order
