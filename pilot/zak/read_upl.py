#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from sys import argv

assert( len(argv)>2)
upl_file = argv[1]
threshold=int(argv[2])

upl_line=open(upl_file,'r').readlines()
lr_nh_nh=0
lr_nh_methyl=0
lr_methyl_methyl=0

def methyl_detect(aa,name):
	if aa=='ALA' and ('B' in name):
		return True
	elif aa=='MET' and ('E' in name):
		return True
	elif aa=='THR' and ('G2' in name):
		return True
	elif aa=='VAL' and ('G' in name):
		return True
	elif aa=='LEU' and ('D' in name):
		return True
	elif aa=='ILE' and ( ('D1' in name) or ('G2' in name)):
		return True
	return False

for line in upl_line:
	tags=line.split()
	if abs(int(tags[0])-int(tags[3]))>threshold:
		if tags[2]=='H':
			if tags[5]=='H':
				lr_nh_nh+=1
			if methyl_detect(tags[4],tags[5]):
#				print line
				lr_nh_methyl+=1
		if tags[5]=='H':
			if methyl_detect(tags[1],tags[2]):
				lr_nh_methyl+=1
		if  methyl_detect(tags[4],tags[5]) and  methyl_detect(tags[1],tags[2]):
			print line
			lr_methyl_methyl+=1

print 'LR_NH_NH  %d'%lr_nh_nh
print 'LR_NH_Methyl  %d'%lr_nh_methyl
print 'LR_Methyl_Methyl  %d'%lr_methyl_methyl



