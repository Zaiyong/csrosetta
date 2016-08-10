#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from sys import argv
import fasta

number_C_H={'A':1,
            'R':3,
            'N':1,
            'D':1,
            'C':1,
            'Q':2,
            'E':2,
            'G':0,
            'H':4,
            'I':4,
            'L':4,
            'K':4,
            'M':3,
            'F':6,
            'P':3,
            'S':1,
            'T':2,
            'W':6,
            'Y':5,
            'V':3}

assert( len(argv)>1)
fasta_file = argv[1]
sequence=fasta.read_fasta(fasta_file)
all_num=0
for aa in sequence:
   all_num+=number_C_H[aa]

print all_num
