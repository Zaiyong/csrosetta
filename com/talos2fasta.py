#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-

#application specific headers
import string

#default headers
import argparse
from basic.options import ExampleArgumentParser
from os.path import basename
import traceback, sys

#toolbox headers
import library


parser = ExampleArgumentParser(prog=basename(__file__),
                               description="extract fasta sequence from chemical shift file (Talos Format)",
															 examples=['%(prog)s 2jrm.tab 2jrm.fasta','%(prog)s 2jrm.tab > 2jrm.fasta'])
parser.add_argument("infile", help="chemical shift file");
parser.add_argument("outfile", help="chemical shift file",nargs="?",default="stdout");
library.add_standard_args( parser )
args = parser.parse_args()

#input:
lines = open( args.infile,'r').readlines();

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello( __file__ )

try:
	sequence="";
	for line in lines:
    tags=string.split(line);
    if (len(tags)<2):
			continue
    if tags[0]=="DATA" and tags[1]=="SEQUENCE":
			for i in range(2, len(tags) ):
				sequence=sequence+tags[i];
    elif tags[0]=="VARS":
			vars=line
    elif tags[0]=="FORMAT":
			format=line[(line.find("FORMAT")+6):]

	outfile.write(">t000\n")
	outfile.write(sequence+"\n")

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)






# def cut_sequence(sequence,start,end,verbose):
#     if ( end==0 ):
#         end=len(sequence)
#     if ( end > len(sequence) ):
#         sys.stderr.write("your end value is out of range\n");
#     new_seq=sequence[(start-1):end];
#     if verbose:
#         print "\nold fasta:\n %s"%sequence
#         print "\nnew fasta:\n %s%s%s"%('-'*(start-1),new_seq,'-'*(len(sequence)-end))
#     return new_seq, end

# def write_sequence(sequence,outfile):
#     for i in range(0,len(sequence),10):
#         if i%50==0: outfile.write("\nDATA SEQUENCE");
#         outfile.write(" "+sequence[i:(i+10)]);
#     outfile.write("\n\n");


