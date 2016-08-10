#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from assignment.noesy.CrossPeakList import CrossPeakList
from assignment.noesy import read_peak_files
from assignment.noesy.ResonanceList import ResonanceList
from basic import options
import score

parser = options.ApplicationArgumentParser(description="assign noesy peaks",
                                add_help=True)

parser.add_argument("-peaks",nargs='*', help='file with assigned peak-list in xeasy format',default=None);
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-weights",help='weight file including score methods and weights',default='weights.txt');

args = parser.parse_args()
crosspeaks=read_peak_files(args.peaks)
resonance_list=ResonanceList.read_from_stream( open(args.prot,'r') )


def simple_test(crosspeaks,resonance_list):

	scfxn=score.ScoreFunction()
	scfxn.set_weight('bmrb',1)
	scfxn.set_weight('prot',1)
	scfxn.set_weight('frag',1)
	scfxn.print_scores(crosspeaks)

simple_test(crosspeaks,resonance_list)


