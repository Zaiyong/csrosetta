#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from sys import argv
import fasta
from PDB.Polypeptide import one_to_three
assert( len(argv)>3)
rdc_file = argv[1]
fasta_file = argv[2]
orientation = argv[3]

rdc_line=open(rdc_file,'r').readlines()
seq=fasta.read_fasta(fasta_file)

if int(orientation)==1:
	error=3.2
elif int(orientation)==2:
	error=4.5

print "#  First atom      Second atom           RDC   Error  Weight Orientation"
for line in rdc_line:
	tags=line.split()
	resid1=int(tags[0])
	atom1=tags[1]
	resid2=int(tags[2])
	atom2=tags[3]
	rdc_value=float(tags[4])
	aa1=one_to_three(seq[resid1-1])
	aa2=one_to_three(seq[resid2-1])
	print'%5d %4s %3s   %5d %4s %3s      %8.3f  %5.3f  1.000 %5d'%(resid1,aa1,atom1,resid2,aa2,atom2,rdc_value,error,int(orientation))
