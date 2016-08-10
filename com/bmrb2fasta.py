#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

from os.path import basename
import argparse
from basic.options import ExampleArgumentParser
import sys
import bmrb
### toolbox library
import library
import traceback
import textwrap
#############################

parser = ExampleArgumentParser(prog=basename(__file__), description="extract fasta sequence from BMRB file",
examples=['%(prog)s 2jrm.bmrb 2jrm.fasta',
					'%(prog)s 2jrm.brmb > 2jrm.fasta'])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("bmrb", help="A bmrb file or the chemical shift section of an bmrb file");
parser.add_argument("outfile", metavar="fasta", help="chemical shift file",nargs="?",default="stdout")
parser.add_argument("-header", help="write a header into the file", action='store_true', default=False )
library.add_standard_args( parser )

args = parser.parse_args()
library.init(args)

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello(__file__)

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

try:
	if verbose:
		print "reading bmrb file %s..."%args.bmrb

	bmrb_file = bmrb.BmrbFile(args.bmrb)

	try:
		seq_frame=bmrb_file._frames['monomeric_polymer']
		if verbose:
			print "Molecule(s) found in BMRB:\n",seq_frame
	except KeyError:
			loops, fields = bmrb_file.capture_loops(open(args.bmrb,'r'))
			bmrb_file.add_frame('generic',loops,fields,'monomeric_polymer');

  (fasta, fastas, nmol)=bmrb.get_sequence( bmrb_file )
	for name, seq in fastas.iteritems():
			outfile.write("\n".join(textwrap.wrap(">%s\n%s\n"%(name,seq)))+"\n")

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)



