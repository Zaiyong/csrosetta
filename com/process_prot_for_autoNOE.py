#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from os.path import basename
import library
import tools
import traceback
import sys
parser=tools.ProcessProtForAutoNOE.get_parser()
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
	pass
#	library.hello( __file__ )

try:
	doit=tools.ProcessProtForAutoNOE(args)
	doit( open(args.infile,'r'), outfile )
except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_info()[0], sys.exc_info()[1], None)


