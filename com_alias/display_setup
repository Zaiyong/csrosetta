#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

#from sys import argv,stderr,exit
#import sys

from os.path import exists
from os.path import basename
from os import path
from glob import glob

import argparse
import sys
import traceback

### toolbox library
try:
	import automatic_setup
	from library import LibException, MissingInput
	import library
except ImportError as exc:
	traceback.print_exc(exc)
	print "\ncall 'source %s/init'"%path.dirname(__file__)
	print "before using the toolbox or put this call into your .bashrc"
	exit()

name=sys.argv[0]
if 'display_setup' in path.basename( name ):
	alias='setup_display.py'
else:
	alias='display_setup.py'


desc= 'This tool allows to view setups stored for the selected target. Alias: %s'%alias
application = automatic_setup.SetupApplication( desc, __file__, sys.argv, method_integration=False )
parser=application.parser
parser.add_argument("-target", help="target_dir is target_prefix/target", default="t000_" );
library.add_standard_args( parser )

application.exit_if_no_method_selected()
application.check_absence_method_options()

try:
	args=parser.parse_args()
	library.init(args)
	application.exit_if_no_method_selected()
	target=application.get_target( args.target )
	setup = automatic_setup.Setup( target, application.method.name, args.label )
	if not setup.exists():
		raise MissingInput("%s does not exist!\nCreate one with setup_target.py"%setup)
	print '\n',setup
	setup.print_file_list()
	application.load_setup(setup)


except LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
#        print sys.exc_type
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
 #       print inst

