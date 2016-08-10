#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


import string

from operator import add
from os.path import basename
import argparse
import sys
import traceback
from assignment import noesy
### toolbox library
import library
import fasta
from assignment.noesy.util import miss_res_atoms, deviant_res_atoms


parser = argparse.ArgumentParser(prog=basename(__file__), description="renumber residues in chemical shift file")
parser.add_argument("-prot", help="chemical shift file");
parser.add_argument("-fasta",help="sequence in fasta format");
parser.add_argument("-t",help="threshold of deviant resonances",type=float,default=0.1)
library.add_standard_args( parser )
args = parser.parse_args()

from cs import ProtCSFile
tab=ProtCSFile()
tab.read_file( args.prot )
sequence=tab.sequence
if not sequence and args.fasta:
	sequence=fasta.read_fasta(args.fasta)
	tab.set_sequence(sequence)

res_in=noesy.ResonanceList.read_from_prot( tab )

#print res_in.sequence()
# o_res_atoms=miss_res_atoms(res_in)
# for atom in no_res_atoms:
# 	print '%s     has no res'%atom

wired_res_atoms=deviant_res_atoms(res_in,args.t)
for atom in wired_res_atoms:
 	print 'res of %s   is wired'%atom
