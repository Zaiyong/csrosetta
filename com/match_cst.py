#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'


import string
import argparse
import sys
import library
from assignment import noesy
import traceback

parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)

parser.add_argument("cst1",help="restraints with distance")
parser.add_argument("cst2",help="restraints with volume")

library.add_standard_args( parser )
args = parser.parse_args()

def read_restraints(fd):
	csts=[]
	for line in fd:
		if len(line):
			csts.append(noesy.DistanceRestraint.read_from_line(line))
	return csts

try:
	library.hello( __file__ )

	cst1=read_restraints(open(args.cst1,'r'))
	cst2=read_restraints(open(args.cst2,'r'))
	for cst in cst1:
		try:
#			vols = [ c._volume for c in cst2 if c.match_fuzzy( cst ) ]
			matches=[ c for c in cst2 if c.match_fuzzy( cst ) ]
			print 'matching: %s'%cst
			if len(matches):
				print '\n'.join(['match %d:  %s'%(i,m) for i,m in enumerate(matches) ])
				vol=0
				cst._volume=sum([ c._volume for c in matches ])
		except ValueError:
			pass
		print cst



except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
