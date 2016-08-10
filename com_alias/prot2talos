#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import string
import argparse
import sys
import traceback

from os import path
### toolbox library
from cs import ProtCSFile
from cs import TalosCSFile

from library import read_aa3_sequence
from fasta import read_fasta
import library

parser = argparse.ArgumentParser(prog=path.basename(__file__), description="renumber residues in chemical shift file")
parser.add_argument("infile", help="chemical shift file");
parser.add_argument("outfile", help="chemical shift file",nargs="?",default="stdout");
parser.add_argument("-seq", help="sequence file");
parser.add_argument("-fasta", help="fasta file");
parser.add_argument("-noheader", help="suppress printing of the header", action='store_false', dest='header', default=True )
parser.add_argument("-header", help="print the header", action='store_true', dest='header', default=True )
library.add_standard_args( parser )

args = parser.parse_args()

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
	sequence=None
	if args.fasta:
		sequence=read_fasta( args.fasta )
	if args.seq:
		sequence=read_aa3_sequence( args.seq)

	prot = ProtCSFile()
	prot.read_file( args.infile, sequence )

	if not sequence:
		sequence = prot.sequence

	talos = TalosCSFile()
	talos.from_table( prot, sequence=sequence )
	talos.write( outfile )


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
