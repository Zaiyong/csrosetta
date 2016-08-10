#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


#from sys import argv,stderr,exit
#import sys
from os.path import exists
from os.path import basename
#from os import path
#from glob import glob
import traceback
#import argparse
import sys
#import os
import subprocess
from basic.options import ExampleArgumentParser
### toolbox library
try:
	import automatic_setup
	from library import LibException, MissingInput
	import library

except ImportError as exc:
	traceback.print_exc(exc)
	print "\ncall 'source %s/init'"%path.dirname(__file__).replace('com_alias','com')
	print "before using the toolbox or put this call into your .bashrc"
	exit()

parser =ExampleArgumentParser(prog=basename(__file__),
										description="setup a swarm of autoNOE-Rosetta runs with different cst-weights (for more options see setup_run -method autoNOE -h)",
										add_help=True)

#add special options that are only relevant for setup_run
parser.add_argument("-run_label", help="add this to the name of the directory", default="" );
parser.add_argument("-cst", nargs='*', help="a list of cst weights to run [5,10,25,50]", type=int, default=[5,10,25,50]);
library.add_standard_args( parser )
#parse cmd-line and check usage
args,argv = parser.parse_known_args()
library.init(args)
print argv

try:

	for cst in args.cst:
		run_label=(args.run_label+'/cst_%02d'%cst).lstrip('/')
		cmdstr='setup_run -method autoNOE -run_label %(run_label)s -noesy_cst_strength %(cst)d '%locals()+" ".join(argv)
		print '***********************************************************'*3
		print 'EXECUTE META SCRIPT: %s'%cmdstr
		print '***********************************************************'*3
		try:
			pipe=subprocess.check_call(cmdstr, shell=True, stdout=sys.stdout )#subprocess.PIPE)
		except subprocess.CalledProcessError:
			print 'Error occurred in command "%s"'%cmdstr
			exit(1)
#		for line in pipe.stdout:
#			if 'substitute flag file' in line: continue
#			if 'integrate job ' in line: continue
#			print line[:-1]

		print '***********************************************************'*3
		print 'SUCCESS SETTING UP CST-WEIGHT: %5d'%cst
		print '***********************************************************'*3

except LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
	exit(1)
