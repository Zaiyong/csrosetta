#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import string
import argparse
import sys
import traceback

from os import path
### toolbox library
from cs import ProtCSFile
from cs import TalosCSFile

from library import read_aa3_sequence
from fasta import read_fasta
import library
from assignment import noesy

parser = argparse.ArgumentParser(prog=path.basename(__file__), description="check completeness of chemical shift assignments")
parser.add_argument("infile", help="chemical shift file");
#parser.add_argument("-labeled", help="residues that are methyl-labelled: e.g., ILV or ILVAMT");
parser.add_argument("-fasta", help="take sequence from fasta file, not from prot file");
parser.add_argument("-seq", help="sequence in three letter format");
parser.add_argument("-verbose", help='how verbose?' );

library.add_standard_args( parser )

args = parser.parse_args()

#output:
verbose=args.verbose

####### program start
if verbose:
	library.hello( __file__ )

expected_atoms={ 'L':['CD1','CD2','HD1','HD2'],'V':['CG1','CG2','HG1','HG2'],'I':['CD1','HD1'] }

class Counter:
	def __init__(self, res_list, verbose=0 ):
		self.res_list = res_list
		self.should_have_ct = 0
		self.actual_has_ct =0
		self.verbose = verbose

	def test_atom(self, atom ):
		try:
			self.should_have_ct+=1
			self.res_list.by_atom( atom )
			self.actual_has_ct+=1
		except KeyError as ex:
			if atom.name()[0]=='H':
				atom._name=atom._name.replace('H','Q')
				try:
					self.res_list.by_atom( atom )
					self.actual_has_ct+=1
					return
				except KeyError as ex:
					pass
			if self.verbose>1: print 'no shifts for Atom %s'%atom

	def completeness(self):
		return 1.0*self.actual_has_ct/self.should_have_ct


try:
	sequence=None
	if args.fasta:
		sequence=read_fasta( args.fasta )
	elif args.seq:
		sequence=read_aa3_sequence( args.seq )

	resonance_list=noesy.ResonanceList.read_from_stream( open(args.infile,'r') )

	if not sequence:
		prot = ProtCSFile()
		prot.read_file( args.infile, sequence )
		sequence = prot.sequence

	if not sequence:
		raise library.MissingInput("require sequence information to check for completeness")

	counts = Counter( resonance_list, verbose )
	counts_backbone = Counter( resonance_list )
	counts_ilv = Counter( resonance_list )

	for pos,aa in enumerate(sequence):
		#check backbone
		#C, CA, CB, N, H
		for atom_name in ['C','CA','CB','N','H']:
			atom = noesy.Atom(atom_name,pos+1)
			counts.test_atom( atom )
			counts_backbone.test_atom( atom )
		#if args.labeled and aa in args.labeled:
		if aa in 'ILV':
			for atom_name in expected_atoms[aa]:
				atom = noesy.Atom(atom_name, pos+1 )
				counts.test_atom( atom )
				counts_ilv.test_atom( atom )


  print 'overall completeness:  %5.1f%%'%(counts.completeness()*100)
	print 'backbone completeness: %5.1f%%'%(counts_backbone.completeness()*100)
	print 'ilv completeness:      %5.1f%%'%(counts_ilv.completeness()*100)




except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
