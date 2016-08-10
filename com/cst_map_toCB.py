#!/usr/bin/python
## make mammoth structure alignments

import string
from glob import glob
from sys import argv,stderr,exit,stdout
from os import popen,system
from os.path import exists
from operator import add
from math import sqrt
import argparse
import traceback
#from pymol.vfont import plain
import sys
from noe_tools import cst_map_to_CB
import library
import os

#############################
def Help():
    print '\n'
    print '-'*75
    print 'USAGE: %s fa_cst'%argv[0]
    print '\n will map all QX to CB and add 1.5 padding per methyl .. result printed to screen'
    print '\n\n\n\n'
    exit()

#if len(argv) <=1:
#    Help()
if 'csrosettaDir' in os.environ:
    CB_to_proton_distance_library=os.environ['csrosettaDir']+"/database/C_to_proton_distances.txt"
else:
    CB_to_proton_distance_library=None

#file = argv[1]
parser =argparse.ArgumentParser(description="mapt the cst between H-H to the cst between C-C",
                                add_help=True)
parser.add_argument("-fasta",help="fasta file",required=True)
parser.add_argument("-library", help="library of CB-H dist",default=CB_to_proton_distance_library )
parser.add_argument("-c_type",help="which C atom is mapped",default="CB")
parser.add_argument("infile",help="(fullatom) constraint file")
parser.add_argument("outfile", help="(centroid) constraint file",nargs="?",default="stdout");
parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )

args = parser.parse_args()



#main program
try:
#output:
    verbose=1
    if args.outfile=="stdout":
        outfile=sys.stdout
        verbose=0
    else:
        outfile=open(args.outfile,'w');

    seq=library.read_fasta(args.fasta)

    librarylist=open(args.library).readlines()
    dist_map_library=[]
    for m in librarylist:
        dist_map_library.append(m.split())

    #print seq
    cst_map_to_C(args.infile, outfile,seq,dist_map_library, args.c_type, verbose=verbose )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)

