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
import rdc_library

parser = ExampleArgumentParser(prog=basename(__file__),
															 description="simple analysis of RDC values using the histogram method",
                       )

parser.add_argument("infile", help="rdc file");

library.add_standard_args( parser )
args = parser.parse_args()

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)



try:
	rdcs=rdc_library.RDC_Data()
	rdcs.read_file( args.infile )
	Dxx, R,range, inv_range = rdcs.estimate_Da_and_R_hist()
	print 'Estimated:  Da %8.3f\nRhombicity: R  %8.3f'%(Dxx,R)


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)


