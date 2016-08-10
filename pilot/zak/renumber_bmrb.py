#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


import string
import bmrb
from ExampleArgumentParser import ExampleArgumentParser
from os.path import basename
import library
from os.path import exists
import cs
import sys
parser = ExampleArgumentParser(prog=basename(__file__), description="Extract chemical shifts from BMRB file and renumber it",
examples=['%(prog)s 2jrm.bmrb 2jrm_trim.bmrb',
					'%(prog)s 2jrm.brmb > 2jrm_trim.bmrb'])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("bmrb", help="A bmrb file or the chemical shift section of an bmrb file");
parser.add_argument("outfile", metavar="bmrb", help="chemical shift file",nargs="?",default="stdout")
parser.add_argument("-s","--start",dest="start",default=1,type=int, help="starting residue");
parser.add_argument("-e","--end",dest="end",default=1000,type=int,help="ending residue");
parser.add_argument("-header", help="write a header into the file", action='store_true', default=False )
library.add_standard_args( parser )

args = parser.parse_args()

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
	bmrb_file, fasta=bmrb.read_cs_bmrb(args.bmrb,verbose)
	cs_data=bmrb_file.process_frame('assigned_chemical_shifts', bmrb.cs_loop_cols )
except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
length=len(cs_data['ID'])
cs_keys=['ID','RESID','RESNAME','ATOMNAME','ATOMTYPE','SHIFT','SIGMA','AMBIGUITY']
bmrb_lines=[]
max_front_trimmed_id=0
for i in range(0,length):
	if cs_data['RESID'][i]<args.start:
		max_front_trimmed_id=max(cs_data['ID'][i],max_front_trimmed_id)
	elif cs_data['RESID'][i]<=args.end and cs_data['RESID'][i]>=args.start:
		try:
			print "%5d %5d %5s %5s %3s %10.3f %8.3f %5d"%(cs_data['ID'][i]-max_front_trimmed_id,cs_data['RESID'][i]-args.start+1,cs_data['RESNAME'][i],cs_data['ATOMNAME'][i],cs_data['ATOMTYPE'][i],cs_data['SHIFT'][i],cs_data['SIGMA'][i],cs_data['AMBIGUITY'][i])
		except:
			print "%5d %5d %5s %5s %3s %10.3f   0.00  %5d"%(cs_data['ID'][i]-max_front_trimmed_id,cs_data['RESID'][i]-args.start+1,cs_data['RESNAME'][i],cs_data['ATOMNAME'][i],cs_data['ATOMTYPE'][i],cs_data['SHIFT'][i],cs_data['AMBIGUITY'][i])
