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

parser = ExampleArgumentParser(prog=basename(__file__), description="extract chemical section from BMRB file",
examples=['%(prog)s 2jrm.bmrb 2jrm.tab',
					'%(prog)s 2jrm.brmb > 2jrm.tab'])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("bmrb", help="A bmrb file or the chemical shift section of an bmrb file");
parser.add_argument("outfile", metavar="tab", help="chemical shift file",nargs="?",default="stdout");
parser.add_argument("-backbone", action='store_true',
										help="remove all atoms that are not used by fragment picker: keep CA, CB, C, N, HN, HA", default=False )
parser.add_argument("-ignore_errors", action='store_true', help='some bmrb files have corrupted data entries that do not affect the chemical shift data. These errors can be ignored', default=False ),
parser.add_argument("-ambiguity", action='store_true', help='read out the AmbiguityCode field of the BMRB shift table', default=False )
#parser.add_argument("-header", help="write a header into the file", action='store_true', default=False )
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

	if args.ignore_errors:
		errors=[]
	else:
		errors=None

	bmrb_file,fasta=bmrb.read_fullbmrb_or_trunkfile(args.bmrb,['assigned_chemical_shifts'],verbose,errors)

	try:
		cs_data=bmrb_file.process_frame(['assigned_chemical_shifts'], bmrb.cs_loop_cols )
	except KeyError:
		raise library.InconsistentInput("File %s does not contain a frame of the category 'assigned_chemical_shifts'. Maybe not a proper BMRB file?"%args.bmrb)
	if verbose:
		extra_msg=''
#		if args.header: extra_msg='with header to '
		print "writing chemical shifts %sfile in tab-format to %s..."%(extra_msg,args.outfile)
	cs_table=cs.NIH_table()
#	print cs_data
	if len(cs_data)<1:
		raise library.InconsistentInput("No chemical shift data found in file %s"%args.bmrb)
	if args.ambiguity:
		if not 'AMBIGUITY' in cs_data[0][2]:
			raise library.InconsistentInput("No AmbiguityCode found in bmrb file %s"%args.bmrb)
		cs_table.from_dict( cs_data[0][2] )
	else:
		cs_table.from_dict( cs_data[0][2], exclude=['AMBIGUITY'] )
	#cs_table.add_column( 'AMBIGUITY', cs_data[0][2]['AMBIGUITY'], format='%2d' )

#	print fasta
	if args.backbone:
		import sets
		keep=sets.Set(['HN','HA','CO','CA','CB','N','H','HA1','HA2','HA3','C','HA'])
		filtered_table={}
		for key,entry in cs_table.table.iteritems():
			if key[1] in keep:
				filtered_table[key]=entry
		cs_table.table=filtered_table

	talos=cs.TalosCSFile()
	talos.from_table( cs_table, sequence=fasta )

	talos.write( outfile )#, header=args.header )


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
		print '\nsuggestion: run with -ignore_errors '



