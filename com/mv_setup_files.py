#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

import string
import argparse
from get_setup_files import get_setup_files
from os import system
import sys
from os import path
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
if 'list_setup' in path.basename( name ):
	alias='setup_list.py'
else:
	alias='list_setup.py'

#setup common-application cmdline options
application = automatic_setup.SetupApplication( 'This tool lists all setups in the target-library. One can filter-out setups using options -method/-methods, '+
																								'-label/-labels and -target/-targets. Use flag -detail to get the full information of the listed setups. '+
																								'This application does not change any setup'+
																								'(requires choice of -method) Alias: %s'%alias,__file__, sys.argv, method_integration=False )

parser=application.parser
parser.add_argument("-target_name", help="")
parser.add_argument("-method_name", help="")
parser.add_argument("-label_name", help="")
parser.add_argument("-option_name", help="")
application.check_absence_method_options()
args=parser.parse_args()

files=get_setup_files(args.target_name,args.method_name,args.label_name,args.option_name,args)
for r in files:
	system('ln -s '+r[1]+' .')
