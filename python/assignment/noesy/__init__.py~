#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'

# please keep to the Coding Guidelines in this module
# Each class should have its own file
# import each class in this module
# add a unit-test to each class, add the unit test below.
# don't access private members from outside the class ( private are names starting with _ , such as crosspeak._dim )

from PeakAssignment import PeakAssignment
from CrossPeak import CrossPeak
from CrossPeakList import CrossPeakList
from ResonanceList import ResonanceList
from Resonance import ResonanceBase, Resonance, RangeResonance
from CrossPeakInfo import CrossPeakInfo
from DistanceRestraint import DistanceRestraint
from Atom import Atom
from NoeStrip import NoeStrip
from SpinSystem import SpinSystem

def read_peak_files(files, ignore_assignments = False, resonances=None):
	cp_list=CrossPeakList()
	for file in files:
		cplnew=CrossPeakList.read_from_stream(open(file,'r'), ignore_assignments, resonances )
		if not hasattr(cplnew._headers[0],'_experiment_id' ):
			import os
			from os import path
			cplnew._headers[0]._experiment_id = path.basename( file )
		cp_list.merge(cplnew)

	return cp_list

def unit_test():
	 import Atom
	 Atom.unit_test()

	 import PeakAssignment
	 PeakAssignment.unit_test()

	 import CrossPeak
	 CrossPeak.unit_test()

	 import CrossPeakList
	 CrossPeakList.unit_test()

	 import CrossPeakInfo
	 CrossPeakInfo.unit_test()

	 import Resonance
	 Resonance.unit_test()

	 import ResonanceList
	 ResonanceList.unit_test()

	 import DistanceRestraint
	 DistanceRestraint.unit_test()

	 import NoeStrip
	 NoeStrip.unit_test()

	 s='''
# Number of dimensions 3
#FILENAME resort_c-ali
#FORMAT xeasy3D
#INAME 1 c
#INAME 2 H
#INAME 3 h
#CYANAFORMAT cHh
#TOLERANCE      0.3    0.04    0.03
    1   40.932    4.805    0.522  1 U 2.114E+05  0.000E+00  e 0
    5   45.226    6.787    2.695  1 U 1.161E+05  0.000E+00  e 0
    6   40.719    1.831    1.764  1 U 1.355E+06  0.000E+00  e 0
    7   40.847    0.934    1.743  1 U 1.674E+05  0.000E+00  e 0
'''

	 r='''
1      1.794      0.040    HA        4 MET M
2      1.75       0.040   HB2        4 MET M
3      1.79      0.040   HB3        4 MET M
4      1.919      0.040    QE        4 MET M
7      0.525      0.040    QG        4 MET M
9    175.642      0.400     C        5 MET M
10     54.943      0.400    CA        5 MET M
11     40.407      0.400    CB        4 MET M
12     16.931      0.400    CE        4 MET M
13     40.943      0.400    CG        4 MET M
14     4.8 0.4 H 12 LEU L
'''

	 from StringIO import StringIO
	 rl=ResonanceList.ResonanceList.read_from_stream(StringIO(r))
	 cl=CrossPeakList.CrossPeakList.read_from_stream(StringIO(s))

	 cl.assign_resonances( rl )

	 print cl
	 for c in cl._peaks:
			print c

	 print rl
