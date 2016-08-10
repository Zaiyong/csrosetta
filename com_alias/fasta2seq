#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-


import string
from os.path import exists
from os.path import basename
import argparse
from basic.options import ExampleArgumentParser
import sys
import library
import fasta


parser = ExampleArgumentParser(prog=basename(__file__), description="obtain sequence file (CYANA) from provided fasta-file",
examples=['%(prog)s target.fasta target.seq',
          '%(prog)s target.fasta > target.seq'])

parser.add_argument("infile", metavar='fasta',nargs="?", help="fasta file", default='stdin' );
parser.add_argument("outfile", metavar='seq',nargs="?", help="seq file", default='stdout' );
parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )

args = parser.parse_args()

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello( __file__ )

if args.infile=='stdin':
	infile=sys.stdin
else:
	infile=open(args.infile,'r')


try:
	sequence=fasta.read_fasta_stream( infile )
	library.write_aa3_sequence_stream( outfile, sequence )


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)

