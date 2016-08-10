#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


import string
from glob import glob
#from sys import argv,stderr,exit
#import sys
from os import popen,system,fdopen
from os import dup2
from os.path import exists
from operator import add
from math import sqrt
from os.path import basename
import argparse
import sys
from basic.options import ExampleArgumentParser
### toolbox library
import library
import traceback
import fasta


parser = ExampleArgumentParser(prog=basename(__file__),
															 description="renumber residues in rdc file with headers",
                               examples=[('(%(prog)s in.tab out.tab -s 5 -e 76','keep residue 5 to 76, start counting in output at 1'),
																				 ('(%(prog)s in.tab -fasta target.trim.fasta > out.tab',
																					'figure out correct trimming from target.trim.fasta, trim output, start counting at 1 for output file')])

parser.add_argument("infile", help="rdc file");
parser.add_argument("outfile", help="rdc file",nargs="?",default="stdout");
parser.add_argument("-s","--start",dest="start",default="1",type=int, help="starting residue");
parser.add_argument("-e","--end",dest="end",default="0",type=int,help="ending residue");
parser.add_argument("-rosetta", help="write in rosetta 3.4 compatible format", default=False, action='store_true' );
mutex=parser.add_mutually_exclusive_group()
mutex.add_argument("-fasta",help="figure out trimming from given sequence");
mutex.add_argument("-rigid",help="use first and last rigid residue from .rigid file as written by pred2rigid",default=None)
library.add_standard_args( parser )
args = parser.parse_args()

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


#output:
verbose=1
if args.outfile=="stdout":
    outfile=sys.stdout
    verbose=0
else:
    outfile=open(args.outfile,'w');

####### program start
if verbose:
    library.hello( __file__ )

try:
	target=0

	sequence="";
	end=args.end

	from rdc_library import RDCFile
	tab=RDCFile()
	tab.read_file( args.infile )
	sequence=tab.sequence
	start=args.start
	end=args.end

	if args.fasta:
    target=fasta.read_fasta(args.fasta)
		start=-fasta.find_fasta_offset(target,sequence)+1
		end=start+len(target)-1;

	if args.rigid:
		start,end = library.read_rigid_file( args.rigid )


	if start>1 or end!=0:
		if verbose>0: print 'Will trim from %d to %d'%(start,end)
		if sequence:
			sequence, end=fasta.cut_sequence(sequence,start,end,verbose)
		tab.renumber( start, end )
	if not args.rosetta:
		tab.write( outfile )
	else:
		atom2=tab.get_slice( 'ATOMNAME2' )
		rdc=tab.get_slice( 'RDC' )
		resid2=tab.get_slice( 'RESID2' )
		keys=rdc.keys();
		keys.sort();
		for key in keys:
			#first residue RDCs often make problem due to atom-names... throw them out !
			if key[0]==1 or resid2[key]==1: continue
			outfile.write( "%8d %5s %8d %5s %8.3f\n"%( key[0],key[1], resid2[key], atom2[key], rdc[key] ) )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)


