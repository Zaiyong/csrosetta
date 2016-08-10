#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-

#application specific headers

#default headers
import argparse
from basic.options import ExampleArgumentParser
from os.path import basename
import traceback, sys
import amino_acids
import fasta
#toolbox headers
import library
import BmrbAtomNames

parser = ExampleArgumentParser(prog=basename(__file__),
															 description="Convert rosetta distance restraint file to CYANA upl format. ",
															 examples=[('%(prog)s noe_autoassign.cst final.upl','convert constraints in noe_autoassign.cst into final.upl')]
															 )
parser.add_argument("infile", metavar='cst', help="ROSETTA cst file");
parser.add_argument("outfile", metavar='upl', help="CYANA upl file",nargs="?",default="stdout");
parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )
parser.add_argument("-fasta", help="fasta file to fill in sequence information");
parser.add_argument("-seq", help="sequence file to fill in sequence information");
parser.add_argument("-verbose", default=1, type=int )

args = parser.parse_args()

#output:
verbose=args.verbose
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	if verbose: library.hello( __file__ )


sequence=None
if args.fasta:
	sequence=fasta.read_fasta( args.fasta )
if args.seq:
	sequence=library.read_aa3_sequence( args.seq)


normal=['AmbiguousNMRDistance','AtomPair']
ambiguous=['AmbiguousNMRConstraint']

class Restraint:
	def __init__(self,name1,resid1,name2,resid2,lb,ub,peakid=None):
		self.resid1=resid1
		self.resid2=resid2
		self.lb=lb
		self.ub=ub
		aa1=sequence[resid1-1]
		aa2=sequence[resid2-1]
		self.resname1=amino_acids.short_to_long[aa1]
		self.resname2=amino_acids.short_to_long[aa2]

		try:
			self.name1=BmrbAtomNames.translate(aa1,name1,reverse=True)
		except KeyError:
			self.name1=name1

		try:
			self.name2=BmrbAtomNames.translate(aa2,name2,reverse=True)
		except KeyError:
			self.name2=name2

		self.peakid=peakid

	def __str__(self):
		s = '%(resid1)4d %(resname1)-4s %(name1)5s %(resid2)4d %(resname2)-4s %(name2)5s %(ub)8.2f'%self.__dict__
		if self.peakid:
			s+='    #peak %d'%self.peakid
		return s

def read_peak(tags,offset):
	try:
		if tags[offset-1]=="Peak":
			peakid=int(tags[offset])
			return peakid
		return None
	except IndexError:
		pass
	except ValueError:
		pass
	return None

def parse_single_cst(line):
	tags=line.split()
	if 'End_' == tags[0][0:4]:
		return None
	if tags[0] in normal:
		name1=tags[1]
		name2=tags[3]
		resid1=int(tags[2])
		resid2=int(tags[4])
		lb=float(tags[6])
		ub=float(tags[7])
		peakid=read_peak(tags,12)
		return Restraint(name1,resid1,name2,resid2,lb,ub,peakid)
	else:
		raise 'cannot read line %s'%line


try:
	file = open(args.infile,'r')
	for line in file:
		tags=line.split()
		if tags[0] in normal:
			restraint=parse_single_cst(line)
			outfile.write('%s\n'%restraint)
		elif tags[0] in ambiguous:
			lb=float(tags[2])
			ub=float(tags[3])
			peakid=read_peak(tags,8)
			for line in file:
				restraint=parse_single_cst(line)
				if not restraint:
					break
				restraint.ub=ub
				restraint.lb=lb
				restraint.peakid=peakid
				if ub:
					ub=0
					lb=0
				outfile.write('%s\n'%restraint)



except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
