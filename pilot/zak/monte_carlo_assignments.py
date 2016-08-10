#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from assignment.noesy.CrossPeakList import CrossPeakList
from assignment.noesy import read_peak_files
from assignment.noesy.ResonanceList import ResonanceList
from basic import options
import score
import random
from utility.boltzmann import boltzmann
parser = options.ApplicationArgumentParser(description="assign noesy peaks",
                                add_help=True)

parser.add_argument("-peaks",nargs='*', help='file with assigned peak-list in xeasy format',default=None);
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-weights",help='weight file including score methods and weights',default='weights.txt');
parser.add_argument("-steps",help='',type=int,default=100);
args = parser.parse_args()
crosspeaks=read_peak_files(args.peaks)
resonance_list=ResonanceList.read_from_stream( open(args.prot,'r') )

old_score=100000.0

def simple_test(crosspeaks,resonance_list):

	scfxn=score.ScoreFunction()
	scfxn.set_weight('bmrb',1)
#	scfxn.set_weight('prot',1)
#	scfxn.set_weight('frag',1)
#	scfxn.print_scores(crosspeaks)
	return scfxn.total_score(crosspeaks)

for r in crosspeaks.iter_peaks():
	r.keep_assignments()

npeaks=crosspeaks.npeaks()

for i in range(args.steps):
	for r in range(npeaks):
		#n=random.randrange(npeaks)
		crosspeaks.peak(r).random_out_assignments()
		new_score=simple_test(crosspeaks,resonance_list)
		delta=new_score-old_score
		if delta<=0:
			print new_score
			old_score=new_score
		else:
			if boltzmann(delta):
				print new_score
				old_score=new_score
			else:
				print 'reject'
				for r in range(npeaks):
					crosspeaks.peak(r).inverse_random_out_assignments()
