#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

import string
from glob import glob
from os import popen,system,fdopen
from os import dup2
from os.path import exists
from math import sqrt
from os.path import basename
import argparse
from basic.options import ExampleArgumentParser
import sys

### toolbox library
import library
import fasta

parser = ExampleArgumentParser(prog=basename(__file__), description="renumber and remove residues in pdb-file",
examples=['%(prog)s -s 5 -e 25 full.pdb trim.pdb',
          '%(prog)s -fasta new_sequence.fasta full.pdb trim.pdb',
          '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])

parser.add_argument("infile", help="pdb file");
parser.add_argument("outfile", help="pdb file",nargs="?",default="stdout");
parser.add_argument("-s","--start",dest="start",default="1",type=int, help="starting residue");
parser.add_argument("-e","--end",dest="end",default="0",type=int,help="ending residue");
mutex=parser.add_mutually_exclusive_group()
mutex.add_argument("-fasta",help="figure out trimming from given sequence");
mutex.add_argument("-rigid",help="use first and last rigid residue from .rigid file as written by pred2rigid",default=None)
library.add_standard_args( parser )

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
	library.hello( __file__ )


try:
	target=0
	start=args.start
	end=args.end
  pdb_fasta=fasta.pdb2fasta(args.infile)
	if args.fasta:
		target=fasta.read_fasta(args.fasta)
		start=-fasta.find_fasta_offset(target,pdb_fasta)+1
		end=start+len(target)-1;
    print pdb_fasta
    print '-'*(start-1)+target
    if verbose:
			print "worked out trimming from fasta-sequences: start: %d end: %d"%(start,end)

	if args.rigid:
		start,end=library.read_rigid_file( args.rigid )

	if verbose>0: print 'Will trim from %d to %d'%(start,end)

	if pdb_fasta:
    pdb_fasta, end=fasta.cut_sequence(pdb_fasta,start,end,verbose)
#input:
	lines = open( args.infile,'r').readlines();
	oldresnum = '   '
	count = 0;
	outid  = outfile
	for line in lines:
    line_edit = line
    if line[0:3] == 'TER':
        continue

    if line_edit[0:4] == 'ATOM' or line_edit[0:6] == 'HETATM':
        if not (line[16]==' ' or line[16]=='A'): continue
        resnum = line_edit[23:26]
        if not resnum == oldresnum:
            count = count + 1
        oldresnum = resnum
        if ( count>=start and count<=end ):
            newnum = '%3d' % (count-start+1)
            line_edit = line_edit[0:23] + newnum + line_edit[26:]
            outid.write(line_edit)
	outid.close()

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)


