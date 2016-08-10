#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import scramble
import numpy
from assignment.noesy import Atom, ResonanceList
import random
import argparse

class FileArgumentParser(argparse.ArgumentParser):
	def convert_arg_line_to_args(self, arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        yield arg

parser = FileArgumentParser(description="swap and miss chemical shifts of side chains",fromfile_prefix_chars='@',add_help=True)
parser.add_argument('-swap',help='percentage of sidechain pairs to be swapped',type=float,default=0)
parser.add_argument('-miss',help='percentage to miss residue',type=float,default=0)
parser.add_argument("-in", dest='input', help="chemical shift file",default=None)
parser.add_argument("-out", help="scrambled chemical shift file",default=None)
args = parser.parse_args()

res_list=ResonanceList.read_from_stream( open(args.input,'r') )

def swap_sidechains( res_list, PCT ):
	####------swap residues--------####
	#swap_aa=['KR','M',]#only these amino acids can be swapped
	seq=res_list.sequence()
	seq_len=len(seq)
	seq_index=range(0,seq_len,1)
	residue_pair=[]
	N=0
	random.shuffle(seq_index)
	for i,r1 in enumerate(seq_index):
		aai=seq[r1]
		for j,r2 in enumerate(seq_index):
			if j<=i: continue
			aaj=seq[r2]
			#for aa in swap_aa:
			if aai == aaj:
				N+=1
				residue_pair.append([r1+1,r2+1])
				break
				#rand_prob=random.uniform(0,1)
				#if rand_prob<prob:
	selected_pair=random.sample(residue_pair,int(seq_len*PCT+0.5))
	for p in selected_pair:
		print 'swap sidechains: %3d %3s  <-> %3d %3s'%(p[0],seq[p[0]-1],p[1],seq[p[1]-1])
		scramble.swap_sidechains(p[0],p[1],res_list)
	print 'there are totally %d same amino acid pairs'%N

def miss_sidechains(res_list, PCT):
	seq=res_list.sequence()
	seq_len=len(seq)
	seq_index=range(1,seq_len+1,1)
	miss_sidechain_residues=random.sample(seq_index,int(seq_len*PCT+0.5))
	for r in miss_sidechain_residues:
		scramble.miss_sidechains( res_list,r )

swap_sidechains(res_list,args.swap)
miss_sidechains(res_list,args.miss)

prot_file=res_list.prot_file()
prot_file.write( open(args.out,'w') )

