#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-

#application specific headers

#default headers
import argparse
from basic.options import ExampleArgumentParser
from os.path import basename
import traceback, sys
import fasta

#toolbox headers
import library

parser = ExampleArgumentParser(prog=basename(__file__),
															 description="extract sequence information from upl restraint file")
parser.add_argument("infile", metavar='upl', help="upl restraint file");
parser.add_argument("outfile", metavar='fasta', help="fasta file", nargs="?",default="stdout");
library.add_standard_args( parser )
args = parser.parse_args()

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello( __file__ )

try:
	lines = open( args.infile,'r').readlines();
	outfile.write('>UPL-FASTA %s\n%s\n'%(args.infile,fasta.upl2fasta(lines)))

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)




