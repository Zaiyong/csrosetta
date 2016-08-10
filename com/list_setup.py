#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

#from sys import argv,stderr,exit
#import sys
from os.path import exists
from os.path import basename
from os import path
from glob import glob
import os
import traceback
import argparse
import sys
import re
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
if 'list_setup' in path.basename( name ):
	alias='setup_list.py'
else:
	alias='list_setup.py'


#setup common-application cmdline options
application = automatic_setup.SetupApplication( 'This tool lists all setups in the target-library. One can filter-out setups using options -method/-methods, '+
																'-label/-labels and -target/-targets. Use flag -detail to get the full information of the listed setups. '+
																'This application does not change any setup'+
																'(requires choice of -method) Alias: %s'%alias,__file__, sys.argv, method_integration=False )

#setup specific cmdline options
parser=application.parser
parser.add_argument("-targets","-target",dest='targets',nargs='*', help="you can filter the display to a subset of targets");
parser.add_argument("-labels",nargs='*', default=None)
parser.add_argument("-methods",nargs='*')
parser.add_argument("-details",help="show full details of all listed setups (requires option -method)", default=False, action='store_true')
parser.add_argument("-show", help="show details for selected options in table form (requires option -method)", nargs='*', default=None )
parser.add_argument("-parse", help="show output like this xya|eas|alkd", default=False, action='store_true')
parser.add_argument("-fix", help="show always all columns", default=False, action='store_true')
parser.add_argument("-file", help="get full filename(s) of this entry")
parser.add_argument("-exact", help="only accept exact matches of target", default=False, action='store_true' )

#check usage  -- no specific method options allowed
application.check_absence_method_options()
args=parser.parse_args()
if not args.quiet and ( (args.file or args.show) and not args.parse ):
	library.hello(__file__ )


if ( args.details or args.show ) and not application.method:
	print "Option -detail or -show only works when you choose a method with -method"
	exit()

#main program
try:
	#look in targetlib for subdirectories
	targets=os.listdir(args.target_prefix)
	if args.show:
		msg="%-15s %-15s "%("Target","Label") +" %-15s"*len(args.show)%(tuple(args.show))
		print msg
		print "-"*len(msg)
	for t in targets:
		if args.targets:
			found=False
			for pattern in args.targets:
				if re.match(pattern, t):
					if args.exact: found=pattern==t
					else:	found=True
					break
			if not found: continue
		#remove things that are just a file and not a subdirectory
		if not path.isdir( args.target_prefix+"/"+t ): continue
#		#attempt to instantiate a TargetDir object
		target=automatic_setup.TargetDir( target=t, prefix=args.target_prefix )
		if not target.exists(): continue
		#at this point target is a valid target

		#look for subdirectories that are methods
		methods=os.listdir(target.dir())
		for m in methods:
			#filter: a) it must be in method_choices,
			#        b) cmdline method selection, single method
			#        c) cmdline method selection, method list
			if not m in application.method_choices: continue
			if application.method and m!=application.method.name: continue
			if args.methods and not m in args.methods: continue
			if args.method and not m in args.method: continue
			#okay, this is a valid method -- get the variant labels
			labels=os.listdir(target.dir()+"/"+m)
			for l in labels:
				#obtain clear method-instance
				application.fresh_method()
				#apply cmdline filters. a) list, b) single label
				if args.labels and not l in args.labels: continue
				if '-label' in sys.argv and l!=args.label: continue

				#instantiate Setup for this set (target, method, label )
				setup=automatic_setup.Setup(target,m,l)

				#prints a line like: Setup( target, method_label )
				set_args=None
				if args.show or args.file or args.details:
					if application.method:
						set_args=application.load_setup(setup, verbose=args.details)

				if args.show:
					vals = [ getattr(set_args,x) for x in args.show ]
					if args.parse:
						str="%s|%s|"%(t,setup.label)
						for val in vals:
							if library.obj_is_list(val):
								str+="@".join(["%s"%s for s in val ])+"|"
							else:
								str+=val+"|"
						print str
					else:
						print "%-15s %-15s "%(t,setup.label) + " %-15s"*len(vals)%(tuple(vals))

				elif args.file:
					val = None
					val = getattr(set_args,args.file)

#					if args.parse:
					if not '-label' in sys.argv and len(labels)>1 or not '-target' in sys.argv or args.fix:
						str="%-15s %-15s "%(t,setup.label)
					else:
						str=""
					if library.obj_is_list(val):
						str+=" ".join(["%s/%s"%(target.dir(),s) for s in val ])
					else:
						if val:
							str+='%s/%s'%(target.dir(),val)
						else:
							str+='None'
					if args.parse:
						str="|".join(str.split())
					print str
#				else:
#						print "%-15s %-15s "%(t,setup.label) + " %-15s/%s"%(target.dir(),val)

				else:
					print setup


				#print details, (need to load them)

#handle exceptions
except LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
#        print sys.exc_type
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
 #       print inst
