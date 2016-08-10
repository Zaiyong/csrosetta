#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


#from sys import argv,stderr,exit
#import sys
from os.path import exists
from os.path import basename
from os import path
from glob import glob
import traceback
import argparse
import sys
import os
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

import rosetta_environment as rosetta


application = automatic_setup.SetupApplication( 'generate run-director from a set of targets with existing setup',__file__, sys.argv, examples=[
		'%(prog)s -dir ~/structure_calculations/ -target t000_ -method rasrec -label with_rdc -job slurm',
		'%(prog)s -dir run -target t000_ -method abrelax -extras default -job interactive',
		'%(prog)s -dir ~/structure_calculations/ -target t000_ -method rasrec -run_label try_again -regenerate_derived_files'],setup_run=True )
parser=application.parser

#now include the flag_lib/methods/XXX/options.py files
job_choices=[]
for jobf in glob (application.flag_lib+"/jobtemplates/production.*.job" ):
	jobf=jobf.replace("production.","").replace(".job","")
	print "integrate job: %s from library"%basename(jobf)
	job_choices.append(basename(jobf))

#add special options that are only relevant for setup_run
parser.add_argument("-dir", help="setup a run-directory here",required=True );
parser.add_argument("-run_label", help="add this to the name of the directory" );
parser.add_argument("-job", help="run-scripts for which job-type? -- default all", choices=job_choices );
parser.add_argument("-database", help="rosetta database location", default=rosetta.db);
parser.add_argument("-bin", help="rosetta binaries location", default=rosetta.bin);
parser.add_argument("-platform", help="which platform", choices=['linux','macos'], default=rosetta.platform );
parser.add_argument("-comp", help="which compiler", choices=['gcc','icc'], default='gcc' );
parser.add_argument("-extras", help="which extras", choices=['default','static','mpi'], default=rosetta.build );
parser.add_argument("-overwrite", help="overwrite run-directory", default=False, action='store_true' )
parser.add_argument("-target","-targets", dest='targets',nargs='*',help="target_dir is target_prefix/target", default=['t000_'] );
parser.add_argument("-regenerate_derived_files", action='store_true', default=False, help='files that are generated during processing of a method (e.g., *cst.centroid ) are regenerated' )
parser.add_argument("-input", help="look for previous results in this directory" );
parser.add_argument("-label2dir", help='instead of attaching target-label to target-name like this: target_label, create another subdirectory: target/label/run_label', default=False, action='store_true' )
library.add_standard_args( parser )
#parse cmd-line and check usage
args = parser.parse_args()
library.init(args)
application.exit_if_no_method_selected()
application.check_absence_method_options()

#start the main program
try:
	if not args.dir:
		raise MissingInput("need to specify a run-directory")

	#figure out job-scripts to build
	jobs=glob( args.flag_lib+"/jobtemplates/*job" )
	selected_jobs=[]
	if args.job:
		for p in jobs:
			if args.job+".job" in p:
				selected_jobs.append(p)
	else:
		selected_jobs=jobs
	from sets import Set
	print "Build job-files "," ".join( Set([x.split('.')[1] for x in  map( basename,selected_jobs )] ))

	#figure out binary postfix
	binary_postfix=args.extras+"."+args.platform+args.comp+"release";
	print binary_postfix

	#setup runs for each target given in cmdline
	for target_name in args.targets:

    #obtain RunGenerator
		subdir_level = 0
		run_dir = args.dir
		if not (target_name=="t000_" and len(args.targets)==1):
			run_dir = args.dir+"/"+target_name
			subdir_level = subdir_level + 1
		if not args.label=="standard":
			if not args.label2dir:
				run_dir = run_dir+"_"+args.label
			else:
				run_dir = run_dir+"/"+args.label
				subdir_level += 1
		if args.run_label:
			run_dir = run_dir+"/"+args.run_label
			import string
			subdir_level = subdir_level + 1 + args.run_label.strip('/').count('/')

		my_run = automatic_setup.RunGenerator( run_dir,
														application.flag_lib,
														selected_jobs,
														args.bin,
														database=args.database,
														extension=binary_postfix, subdir_level = subdir_level )

    #obtain clear method-instance
		application.fresh_method()

    #get the target and setup
		target=application.get_target( target_name, rundir=args.dir )
		setup = automatic_setup.Setup( target, application.method.name, args.label, args.regenerate_derived_files )

		if args.input:
			if exists( args.input+"/"+target_name ):
				application.method.set_input_dir( args.input+"/"+target_name, args.dir )
			else:
				application.method.set_input_dir( args.input )

		#check existence
		if not setup.exists():
			raise MissingInput("%s does not exist!\nCreate one with setup_target.py"%setup)
		print '\n',setup

		#load setup
		args=application.load_setup(setup)

		#setup method-flag-files etc. and apply patches
		application.method.setup_file_library()

		if args.patch:
			exec open( args.patch, 'r' )
			patch_flags( application.method.file_library )

		#setup run
		my_run.generate_dirs( args.overwrite )
		my_run.generate_flags( setup, application.method )

		#final message with possible user instructions
		application.method.motd(my_run.rundir())


except LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
	exit(1)
