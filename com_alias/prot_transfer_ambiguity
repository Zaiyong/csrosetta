#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os import dup2
from os.path import exists
from os.path import basename
import argparse
import sys
import traceback
import cs
from assignment import noesy
import library
import BmrbAtomNames
from basic.options import ExampleArgumentParser

parser = ExampleArgumentParser(prog=basename(__file__), description="extract ambiguity information from one prot-file and add to second file where atoms match", examples=[]   )
parser.add_argument("-ref", help="A prot file that provides that reference ambiguity info");
parser.add_argument("infile", metavar="prot", help="input prot file ")
parser.add_argument("outfile", metavar="prot", help="output prot file");
parser.add_argument("-noheader", dest='header', help='output header or not', action='store_false', default=True)

#library.add_standard_args( parser )
library.add_standard_args( parser )


args = parser.parse_args()
library.init(args)

def fix_name( reso, aa ):
	if reso.name()=="QD" and aa in 'L':
		reso._atom._name="QQD"
	if reso.name()=="QG" and aa in 'V':
		reso._atom._name="QQG"

def unpack(atom,aa):
	if atom.elem() not in 'HC': return []
	val = BmrbAtomNames.anti_degenerate(aa,atom.name())
#	print atom, aa, val
	return val

from cs import ProtCSFile
tab=ProtCSFile()
tab.read_file( args.ref )
sequence=tab.sequence
reference=noesy.ResonanceList.read_from_prot( tab )

from cs import ProtCSFile
tab=ProtCSFile()
tab.read_file( args.infile )
sequence_in=tab.sequence
payload=noesy.ResonanceList.read_from_prot( tab )

if sequence!=sequence_in:
	raise library.InconsistentInput("Sequence Mismatch")


for reso in payload.itervalues():
	aa = sequence[reso.resid()-1]
	fix_name( reso, aa )
	try:
		ref_reso = reference.by_atom( reso.atom() )
		reso.ambiguity = ref_reso.ambiguity
	except KeyError:
		reso.ambiguity = 4
		unpacked = unpack( reso.atom(), aa )
		for name in unpacked:
			try:
				ref_reso = reference.by_atom( noesy.Atom( name, reso.resid() ) )
				reso.ambiguity = ref_reso.ambiguity
				continue
			except KeyError:
				pass
		try:
			pn, s = BmrbAtomNames.get_combine( aa, reso.name() )
			pn = pn.replace('P','')
			pn = pn.replace('H','Q')
#			print pn, reso
			try:
				ref_reso = reference.by_atom( noesy.Atom( pn, reso.resid() ) )
				reso.ambiguity = ref_reso.ambiguity
			except KeyError:
				pass
		except:
			pass


prot_data=payload.generate_dict()
ambiguity=[]
for r in payload.itervalues():
	ambiguity.append( r.ambiguity )
	prot_data['AMBIGUITY']=ambiguity

nih_table = cs.NIH_table().from_dict( prot_data )
#print 'convert to ProtCS-File'
prot_file = cs.ProtCSFile().from_table( nih_table )
prot_file.write_file( args.outfile, header=args.header )
