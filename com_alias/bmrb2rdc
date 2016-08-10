#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-
from os.path import basename
import sys
import argparse
from basic.options import ExampleArgumentParser
import bmrb
### toolbox library
import library
from os.path import exists
import cs
import traceback
#############################

parser = ExampleArgumentParser(prog=basename(__file__), description="extract rdcs from BMRB file if present",
examples=['%(prog)s 2jrm.bmrb 2jrm.tab',
					'%(prog)s 2jrm.brmb > 2jrm.tab'])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("bmrb", help="A bmrb file or just the relevant section of an bmrb file");
parser.add_argument("outfile", metavar="rdcfile", help="rosetta-format RDC file",nargs="?",default="stdout");

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
	frame_categories=['residual_dipolar_couplings','RDCs']
	bmrb_file,fasta=bmrb.read_fullbmrb_or_trunkfile(args.bmrb,frame_categories,verbose)
	rdc_data=bmrb_file.process_frame(frame_categories, bmrb.rdc_loop_cols )
	if verbose:
		extra_msg=''
#		if args.header: extra_msg='with header to '
		print "writing rdcs %sfile in tab-format to %s..."%(extra_msg,args.outfile)
#	cs_table=cs.NIH_table()
#	cs_table.from_dict( cs_data )
#	print rdc_data
  for i,r in enumerate(rdc_data):
		raw_table=cs.NIH_table()
		raw_table.from_dict( r[2] )
		rdc_table=rdc_library.RDCFile()
		rdc_table.from_table( raw_table, sequence=fasta )
		rdc_table.write_file( 'med_%s.rdc'%r[1])


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)



