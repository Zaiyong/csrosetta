#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

#The function of this program is load Crosspeak assignment file (NOE_out.dat)
from sys import argv,stderr,stdout
from os import popen,system
from os.path import exists
import sys
sys.path.append('/Users/zak/csrosetta3/pymol')

import argparse
from Noesy_Data_Collection import Noesy_Data_Collection

#parser = argparse.ArgumentParser(description="load Crosspeak assignment file",
#                                 add_help=True)

#parser.add_argument("-input", help="NOE_out.dat");
#parser.add_argument("-noe_type", help="certain noesy type",type=int,default=1)

#args = parser.parse_args()


#assert( len(argv)>1)
input = argv[1]

#from PeakAssignment import PeakAssignment


#--------------------------------------------

#--------------------------------------------
#from NOESY_Data_collection import NOESY_Data_collection
#from Peaklist import Peaklist
#from PeakAssignment import PeakAssignment
infile=Noesy_Data_Collection.read_from_file(input)
print infile
print infile._peaklists[2]
print infile._peaklists[1]
print infile._peaklists[2]._crosspeaks[2]
print infile._peaklists[2]._crosspeaks[0]
print infile._peaklists[2]._crosspeaks[0]._assignments[0]
print infile._peaklists[2]._crosspeaks[0]._assignments[3]

#print(r)
#print('\n')
    #for m in Peaks[0][0:10]:
#print m
#print('\n')
    #for m in peakassign[0][0:30]:
#print m
