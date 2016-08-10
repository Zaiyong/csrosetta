#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os.path import basename
import os
import string
import argparse
import sys
#sys.path.append('/home/zak/Downloads/numpy-1.6.2')
import library
import StringIO
from assignment import noesy

parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*', help='NOE_out.dat with assigned peak-lists in xeasy format',default=None);

library.add_standard_args( parser )
args = parser.parse_args()

crosspeaks=noesy.read_peak_files(args.peaks)
no_assign_peaks=[]
for peak in crosspeaks:
	if peak.nassign()==0:
		no_assign_peaks.append(peak)

print 'number of no assign peaks ',len(no_assign_peaks)
# for peak in no_assign_peaks:
# 	crosspeaks.remove_crosspeak(peak)

# crosspeaks.write_split_files()
