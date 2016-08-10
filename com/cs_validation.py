#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os.path import basename
import string
import argparse
import sys
#sys.path.append('/home/zak/Downloads/numpy-1.6.2')
import library
import StringIO
import noesy
#from noesy import CrossPeakInfo
#from noesy import CrossPeakList
#from noesy import Resonance
#from noesy import ResonanceList
from noesy import SpinSystem
import fasta

def get_strips( resid ):
	from noesy import NoeStrip
	from noesy import Atom
	from noesy import RangeResonance
	from noesy import Resonance
#get HN_strip
	N_atom=Atom('N',resid)
	H_atom=Atom('H',resid)
	N_resonance=resonance_list.by_atom( N_atom )
	H_resonance=resonance_list.by_atom( H_atom )
	HN_strip=NoeStrip.generate_strip(crosspeaks,1,H_resonance,N_resonance)

#get HA_strip
	CA_atom=Atom('CA',resid)
	HA_atom=Atom('HA',resid)
	CA_resonance=resonance_list.by_atom( CA_atom )
	HA_range=RangeResonance(-1,HA_atom,3,6,0.03)
	HA_strips=NoeStrip.generate_strip_family(crosspeaks,1,HA_range,CA_resonance)

#get HB_strip
	CB_atom=Atom('CB',resid)
	HB_atom=Atom('HB',resid)
	CB_resonance=resonance_list.by_atom( CB_atom )
	HB_range=RangeResonance(-1,HB_atom,0,2,0.02)
	HB_strips=NoeStrip.generate_strip_family(crosspeaks,1,HB_range,CB_resonance)

#	print HN_strip
#	print 'Family HA: '
	print 'true HA %d: %s'%(args.resid, resonance_list.by_atom( HA_atom ).freq() )
#	print '\n Strip: '.join(['%s'%s for s in HA_strips])

#	print 'Family HB: '
	print 'true HB: %d: %s'%(args.resid, resonance_list.by_atom( HB_atom ).freq() )
#	print '\n Strip: '.join(['%s'%s for s in HB_strips])
	return HN_strip, HA_strips, HB_strips


parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-resid",help='whichr residue is validated',type=int, default=1);
parser.add_argument("-putative_num", help='how many putative stripes will be generated for one proton',type=int, default=50 )
parser.add_argument("-new_label", help='Besides HA-CA and HB-CB, another label atom', default='CG' )
parser.add_argument("-new_H", help='Besides HA-CA and HB-CB, another label atom', default='HG' )
parser.add_argument("-library", help="library of CS distribution bounds",default=None )
parser.add_argument("-fasta", help="fasta file",default=None )
library.add_standard_args( parser )
args = parser.parse_args()

#librarylist=open(args.library,'r').readlines()
resonance_list=noesy.ResonanceList.read_from_stream( open(args.prot,'r') )
crosspeaks=noesy.read_peak_files(args.peaks)

if args.fasta:
	seq=fasta.read_fasta(args.fasta)

systems=[]
HN, HAs, HBs = get_strips(args.resid)
for HA in HAs:
	for HB in HBs:
		sp=SpinSystem(args.resid)
		sp.add_strip(HN)
		sp.add_strip(HA)
		sp.add_strip(HB)
		if sp.score() > 0.5:
			systems.append(sp)

systems=sorted(systems, key=SpinSystem.score, reverse=True )
for sys in systems:
	print '%5.3f %s'%(sys.score(), sys)


