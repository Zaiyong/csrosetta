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

parser = argparse.ArgumentParser(prog=path.basename(__file__), description="renumber residues in chemical shift file")
parser.add_argument("infile", help=".prot chemical shift file");
parser.add_argument("outfile", help=".prot chemical shift file",nargs="?",default="stdout");
parser.add_argument("-stereo", nargs='*', choices=['beta','aromatics','combine'], help='define all beta-protons of a residue as floating assignments', default=None );

library.add_standard_args( parser )

args = parser.parse_args()

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

def beta_proton( atom, letter = 'B', elem='HQ' ):
	return (atom.name()[0] in elem) and letter in atom.name()

try:
	res_list = noesy.ResonanceList.read_from_stream( open(args.infile,'r') )
	for resid,resonances in res_list.iter_residues():
		aa = res_list.sequence()[ resid - 1]
		if not aa in ['G','A']:
			floats = []
			if args.stereo and 'beta' in args.stereo:
				for betas in res_list.by_residue( resid ):
					if not beta_proton( betas.atom() ): continue
					floats.append( betas.id() )
				if len(floats)>1:
					for betas in res_list.by_residue( resid ):
						if not beta_proton( betas.atom() ): continue
						betas.set_floats( floats, res_list )
			if args.stereo and 'aromatics' in args.stereo and aa in ['F','W','Y','H']:
				for rank in ['D','E']:
					floats = []
					for betas in res_list.by_residue( resid ):
					#	print betas.atom(), rank, beta_proton( betas.atom(), rank )
						if not beta_proton( betas.atom(),letter=rank ): continue
						floats.append( betas.id() )
					if len(floats)>1:
						for betas in res_list.by_residue( resid ):
							if not beta_proton( betas.atom(), letter=rank ): continue
							betas.set_floats( floats, res_list )
			if args.stereo and 'combine' in args.stereo:
				for rank in ['B','G','D','E','Z']:
					by_freq={}
					labels={}
					for r in res_list.by_residue( resid ):
						if beta_proton( r.atom(), rank, 'C' ):
							labels[r.atom()]=r
						if not beta_proton( r.atom(), rank ): continue
						by_freq.setdefault(r.freq(),[]).append(r)
					for f,rs in by_freq.iteritems():
						if len(rs)>1:
							for r in rs:
								label=r.atom().name().replace('H','C')
								label2=label[:-1]
								print r.atom(), label, label2
							catom = noesy.Atom( rs[0].atom().name().replace('H','Q')[:-1],rs[0].atom().resid() )
							combined = noesy.Resonance(rs[0].id(), catom, f, rs[0].error() )
							for r in rs:
								res_list.remove_resonance( r )
							res_list.add_resonance( combined )


		#some have the same frequency -- these should be grouped into QXXX

	print res_list






except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
