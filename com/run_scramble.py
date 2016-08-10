#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os.path import basename, splitext
import string
import argparse
import sys
import StringIO
import math
import fasta
import scramble
from assignment.noesy import Atom, ResonanceList
import BmrbAtomNames
#from PDB.Polypeptide import one_to_three
import random
from PDB.Polypeptide import one_to_three
from math import fabs
from chemical import AllAminoAcidLib
#parser0 =argparse.ArgumentParser(description="",add_help=False)
#args0 = parser0.parse_args()
def pool_name( name, aa ):
	import BmrbAtomNames
	try:
		pn, s = BmrbAtomNames.get_combine( aa, name )
		if not pn:
			pn = name
		rejects = None
		if 'QQ' in pn:
			rejects = BmrbAtomNames.get_rejects( aa, name )
		return pn, s, rejects
	except KeyError:
		return name, None, None

class FileArgumentParser(argparse.ArgumentParser):
	def convert_arg_line_to_args(self, arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        yield arg

parser = FileArgumentParser(description="assign noesy peaks",fromfile_prefix_chars='@',add_help=True)
parser.add_argument('-swap_methyl',help='NO. of residue pairs  to be swapped methyls',type=int,default=0)
parser.add_argument('-swap_sidechain',help='NO. of sidechain pairs to be swapped',type=int,default=0)
parser.add_argument('-combine_methyl',help='percentage to combine methyls',type=float,default=0)
parser.add_argument('-miss_methyl',help='percentage to miss methyl',type=float,default=0)
parser.add_argument('-miss_sidechain',help='percentage to miss residue',type=float,default=0)
parser.add_argument('-miss_proton',help='percentage to miss proton',type=float,default=0)
parser.add_argument('-swap_proton',help='percentage to swap proton',type=float,default=0)
parser.add_argument('-swap_carbon',help='percentage to swap carbon',type=float,default=0)
parser.add_argument('-swap_coupled',help='NO. of pairs of C-H to be swapped',type=int,default=0)
parser.add_argument('-combine_stereo',help='percentage to combine stereo',type=float,default=0)
parser.add_argument('-swap_stereo',help='percentage to swap stereo',type=float,default=0)
parser.add_argument("-in", dest='input', help="chemical shift file",default=None)
parser.add_argument("-out", help="scrambled chemical shift file",default=None)
parser.add_argument("-random_seed", help="set a random seed, i.e,. for testing", default=None, type=int)

args = parser.parse_args()
if args.random_seed:
	random.seed(args.random_seed)

def swap_coupled(res_list,num):
	import BmrbAtomNames
	swapped_protons=[]
	seq=res_list.sequence()
	for i in range(0,len(seq)):
		resid=i+1
		try:
			resid_res_list=res_list.by_residue(resid)
		except KeyError:
			continue
		for r in resid_res_list:
			if r.atom().elem()=='H' and r.name() !='H' and r.name() !='HA':
				#rand_prob=random.uniform(0,1)
				#if rand_prob<prob:
				swapped_protons.append(r)
	#swapped_protons=random.sample(all_protons,int(len(all_protons)*PCT+0.5))
	random.shuffle(swapped_protons)
	midpoint=len(swapped_protons)/2
	for (m1,m2) in zip(swapped_protons[0:midpoint],swapped_protons[midpoint:]):
		resid1=m1.resid()
		resid2=m2.resid()
		try:
			heavy_res1=res_list.by_atom(Atom(BmrbAtomNames.get_source_atom(seq[resid1-1],m1.name()),resid1))
			heavy_res2=res_list.by_atom(Atom(BmrbAtomNames.get_source_atom(seq[resid2-1],m2.name()),resid2))
			if (fabs(m1.freq()-m2.freq()) <= 2) and (fabs(heavy_res1.freq()-heavy_res2.freq()) <= 2) and num >0:
				scramble.swap_frequencies(m1,m2)
				scramble.swap_frequencies(heavy_res1,heavy_res2)
				print 'swapped proton and its source: %15s    <-> %15s'%(m1.atom(), m2.atom() )
				num-=1
		except KeyError:
			continue


def swap_stereo(res_list,PCT):
	all_aa_lib=AllAminoAcidLib()
	seq=res_list.sequence()
	stereo_couples=[]
	for r in res_list.itervalues():
		resid=r.resid()
		if r.atom().elem()=='H':
			aa=res_list.sequence()[r.resid()-1]
			if aa=='-': continue
			stereo_name=all_aa_lib[aa].stereo(r.name())
			if stereo_name and ([res_list.by_atom(Atom(stereo_name,resid)),r] not in stereo_couples):
				try:
					stereo_couples.append([r,res_list.by_atom(Atom(stereo_name,resid))])
				except KeyError:
					continue
	swap_stereos=random.sample(stereo_couples,int(len(stereo_couples)*PCT+0.5))
	for stereo_couple in swap_stereos:
		freq0=stereo_couple[0].freq()
		freq1=stereo_couple[1].freq()
		stereo_couple[0].set_freq(freq1)
		stereo_couple[1].set_freq(freq0)

# 		for r in resid_res_list:
# 			if r.atom().elem() !='H': continue
# 			pn,s,rejects=pool_name(r.name(),seq[i])
# 			if ('Q'+pn[2:] in combined_stereo):
# 				print 'delete %s'%r.atom()
# 				res_list.delete_by_atom(Atom(r.name(),resid))


def combine_stereo(res_list,PCT):
	seq=res_list.sequence()
	for i in range(0,len(seq)):
		resid=i+1
		all_stereo_names=[]
		try:
			resid_res_list=res_list.by_residue(resid)
		except KeyError:
			continue
		for r in resid_res_list:
			#rand_prob=random.uniform(0,1)
			try:
				if BmrbAtomNames.get_proR_proS(seq[i],r.name()):
					pn,s,rejects=pool_name(r.name(),seq[i])
					if 'Q'+pn[2:] not in all_stereo_names:
						#if rand_prob < prob:
						#r.atom().set_name('Q'+pn[2:])
						all_stereo_names.append('Q'+pn[2:])
						#print r.atom()
			except KeyError:
				continue
		combined_stereo=random.sample(all_stereo_names,int(len(all_stereo_names)*PCT+0.5))
		finished_stereo=[]
		for r in resid_res_list:
			if r.atom().elem() !='H': continue
			pn,s,rejects=pool_name(r.name(),seq[i])
			if ('Q'+pn[2:] in combined_stereo) and ('Q'+pn[2:] not in finished_stereo):
				r.atom().set_name('Q'+pn[2:])
				finished_stereo.append('Q'+pn[2:])
				print 'combine stereo on %15s of residue %d '%(s,r.resid())
		for r in resid_res_list:
			if r.atom().elem() !='H': continue
			pn,s,rejects=pool_name(r.name(),seq[i])
			if ('Q'+pn[2:] in combined_stereo):
				print 'delete %s'%r.atom()
				res_list.delete_by_atom(Atom(r.name(),resid))

# def swap_protons( res_list, resid, aa):
#   #which protons to work on
# 	aromatic_protons={ 'F': ['HD1','HE1','HZ','HE2','HD2'],
# 										 'Y': ['HD1','HD2','HE1','HE2'],
# 										 'W': ['HE3','HZ3','HZ2'] }

# 	if aa not in 'FYW': return
# 	midpoint=len(aromatic_protons[aa])/2
# 	random.shuffle(aromatic_protons[aa])
# 	for (m1,m2) in zip(aromatic_protons[aa][0:midpoint],aromatic_protons[aa][midpoint:]):
# 		try:
# 			r1=res_list.by_atom( Atom(m1,resid) )
# 			r2=res_list.by_atom( Atom(m2,resid) )
# 			scramble.swap_frequencies(r1,r2)
# 			print 'swapped aromatic proton: %15s    <-> %15s'%(r1.atom(), r2.atom() )
# 		except KeyError:
#           #only swap aromatic protons
# 			continue


def swap_carbons(res_list,PCT):
	seq=res_list.sequence()
	all_carbons={}
	carbon_list=['CB','CG','CD','CZ','CE','CG1','CG2','CD1','CD2','CE1','CE2','CE3','CZ2','CZ3','CH2']
	for carbon in carbon_list:
		all_carbons[carbon]=[]
		for i in range(0, len(seq) ):
			resid=i+1
			#rand_prob=random.uniform(0,1)
			#if rand_prob<prob:
			try:
				all_carbons[carbon].append(res_list.by_atom(Atom(carbon,resid)))
			except KeyError:
				continue
		swapped_carbons=random.sample(all_carbons[carbon],int(len(all_carbons[carbon])*PCT+0.5))
		midpoint=len(swapped_carbons)/2
		for (m1,m2) in zip(swapped_carbons[0:midpoint],swapped_carbons[midpoint:]):
			print 'swap carbon %15s    <-> %15s'%(m1.atom(),m2.atom())
			scramble.swap_frequencies(m1,m2)

def swap_methyls( res_list, num ):
	####------swap methyl--------####
	all_methyls=[]
	seq=res_list.sequence()
	for i in range(0, len(res_list.sequence()) ):
		resid=i+1
		#rand_prob=random.uniform(0,1)
		#if rand_prob<prob:
		if seq[i]  in 'AILVMT':
			try:
				all_methyls+=scramble.search_methyls(resid,seq[i],res_list)
			except KeyError as exc:
				print exc
				print 'cannot find methyl for resid %s amino acid %s, probably got lost already'%(resid,seq[i])
				continue
	swap_methyls=random.sample(all_methyls,min(num,len(all_methyls)))

	## switch methyls. Switch frequencies of resonances in place
	midpoint=len(swap_methyls)/2
	for (m1,m2) in zip(swap_methyls[0:midpoint],swap_methyls[midpoint:]):
		print 'swap methyls %15s    <-> %15s'%(m1['H'].atom(),m2['H'].atom())
		scramble.switch_coupled_res(m1,m2)


def swap_sidechains( res_list, num ):
	####------swap residues--------####
	#swap_aa=['KR','M',]#only these amino acids can be swapped
	seq=res_list.sequence()
	N=0
	residue_pair=[]
	for i in range(0,len(res_list.sequence())-1):
		aai=seq[i]
		for j in range(i+1,len(res_list.sequence())):
			aaj=seq[j]
			#for aa in swap_aa:
			if aai == aaj:
				N+=1
				residue_pair.append([i+1,j+1])
				#rand_prob=random.uniform(0,1)
				#if rand_prob<prob:
	selected_pair=random.sample(residue_pair,num)
	for p in selected_pair:
		print 'swap sidechains: %3d %3s  <-> %3d %3s'%(p[0],seq[p[0]-1],p[1],seq[p[1]-1])
		scramble.swap_sidechains(p[0],p[1],res_list)
	print 'there are totally %d same amino acid pairs'%N

def miss_protons(res_list,resid,PCT):
	protons=[]
	for r in res_list.by_residue(resid):
		if r.atom().elem()=='H' and (r.name() not in ['H','HN','HA','HA2','HA3','QA']): protons.append(r.atom())
	missed_protons=random.sample(protons,int(len(protons)*PCT+0.5))
	for proton in missed_protons:
		#rand_prob=random.uniform(0,1)
		#if rand_prob < prob:
		try:
			res_list.delete_by_atom(proton)
			print 'delete  resonance of  proton %s'%proton
		except IndexError:
        #didn't find any protons in this residue
			continue

resonances=ResonanceList.read_from_stream( open(args.input,'r') )

# type_prob={'combine_methyl':args.combine_methyl,
# 					 'swap_methyl':args.swap_methyl,
# 					 'miss_methyl':args.miss_methyl,
# 					 'swap_residue':args.swap_residue,
# 					 'miss_residue':args.miss_residue,
# 					 'miss_proton':args.miss_proton,
# 					 'swap_carbon':args.swap_carbon,
# 					 'combine_stereo':args.combine_stereo,
# 					 'swap_coupled':args.swap_coupled,
# 					 'swap_proton':args.swap_proton}

swap_methyls( resonances, 2*args.swap_methyl )
swap_sidechains( resonances, args.swap_sidechain )
swap_carbons(resonances, args.swap_carbon )
swap_coupled(resonances, args.swap_coupled)
combine_stereo(resonances,args.combine_stereo)
swap_stereo(resonances,args.swap_stereo)
####------other types--------####


seq=resonances.sequence()
size=len(seq)

two_methyl_residue_list=[]
methyl_residue_list=[]
for i in range(0,size):
	if seq[i] in 'ILV':
		two_methyl_residue_list.append(i+1)
		methyl_residue_list.append(i+1)
	elif seq[i] in 'AMT':
		methyl_residue_list.append(i+1)
combine_methyl_residues=random.sample(two_methyl_residue_list,int(len(two_methyl_residue_list)*args.combine_methyl+0.5))
for r in combine_methyl_residues:
	resid=r
	scramble.combine_methyls( resonances, resid, seq[resid-1] )

miss_methyl_residues=random.sample(methyl_residue_list,int(len(methyl_residue_list)*args.miss_methyl+0.5))
for r in miss_methyl_residues:
	resid=int(r)
	scramble.miss_methyl( resonances, resid, seq[resid-1] )

import numpy
residue_list=numpy.linspace(1,size,size)
miss_sidechain_residues=random.sample(residue_list,int(size*args.miss_sidechain+0.5))
for r in miss_sidechain_residues:
	resid=int(r)
	scramble.miss_sidechains( resonances, resid )

for i in range(0,size):
	resid=i+1
	try:
		miss_protons(resonances,resid,args.miss_proton)
	except KeyError:
		continue
# for i in range(0,size):
# 	resid=i+1
# 	aa = seq[i]

# 	for t in type_prob.iteritems():
# 		if t !='miss_proton':
# 			rand_prob=random.uniform(0,1)
# 			if rand_prob>p: continue
# 			else:
# 				elif t=='combine_methyl': scramble.combine_methyls( resonances, resid, aa )
# 				elif t=='miss_methyl':
# 					try:
# 						scramble.miss_methyl( resonances, resid, aa )
# 					except KeyError:
# 						continue
# 				elif t=='miss_residue': scramble.miss_residue( resonances, resid )
# #			elif t=='miss_proton': scramble.miss_proton( resonances, resid )
# 				else: continue
# 		else:
# 			try:
# 				miss_protons(resonances,resid,p)
# 			except KeyError:
# 				continue
prot_file=resonances.prot_file()
prot_file.write( open(args.out,'w') )




