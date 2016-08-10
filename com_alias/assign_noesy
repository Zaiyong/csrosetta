#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
## make mammoth structure alignments


import string
import argparse
import sys
import library
from assignment import noesy
import fasta
#############################
#if len(argv) <=1:
#    Help()

#file = argv[1]
parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("peaks",nargs='*',help="files with peak-lists in xeasy format")
parser.add_argument("-resonances", help="a file with resonances for assignment in prot format",required='True');
parser.add_argument("-fasta", help='sequence information, if not already contained in resonance file', default=None );
parser.add_argument("-o", help='output file with assigned peaks', default='NOE.dat' );
parser.add_argument("-cst", help='output file with generated rosetta restraints', default=None );
parser.add_argument("-skip", help='skip assignment step', action='store_true', default=False );
parser.add_argument("-ignore", help='ignore assignments in peak file', action='store_true', dest='ignore', default=True )
parser.add_argument("-noignore", help='ignore assignments in peak file', action='store_false', dest='ignore', default=True )
library.add_standard_args( parser )
args = parser.parse_args()

def dump_list( cl, fd ):
	curr_header=cl._headers[0]
	curr_header.write(fd)
	for c in cl._peaks:
		if not c._info is curr_header:
			curr_header=c._info
			curr_header.write(fd)
		fd.write('%s\n'%c)

try:
#output:
	library.hello( __file__ )
	res_file = open( args.resonances, 'r')
	resonances=noesy.ResonanceList.read_from_stream( res_file )

	if args.fasta:
		sequence=fasta.read_fasta( args.fasta )

	if args.fasta and resonances.sequence():
		offset=fasta.find_fasta_offset( resonances.sequence(), sequence )
		if offset:
			raise library.InconsistentInput('Sequence in %s has offset of %d from sequence in %s'%(args.fasta, offset, args.resonances))
	if args.fasta:
		resonances.set_sequence(sequence)


	cl=noesy.read_peak_files(args.peaks, args.ignore, resonances )

	if not args.skip:
		print "assign peaks"
		cl.assign_resonances( resonances )
	else:
		print "skip assignment..."

	if args.cst:
		cst_file=open(args.cst,'w')
		for c in cl._peaks:
			cst=c.generate_rosetta_cst()
			if cst: cst_file.write('%s\n'%cst )

#	print resonances
	if args.o:
		dump_list( cl, open(args.o,'w') )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
