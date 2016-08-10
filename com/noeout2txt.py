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
from math import *
#from noesy import CrossPeakInfo
#from noesy import CrossPeakList
#from noesy import Resonance
#from noesy import ResonanceList
from assignment.noesy import SpinSystem, get_strips,spinsystem_evaluation,ResonanceBase, ResonanceList
import fasta
import copy
from sets import Set

parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*', help='NOE_out.dat with assigned peak-lists in xeasy format',default=None);
parser.add_argument("-split", type=int, dest='split_level',  default=None);
parser.add_argument("-split_level", type=int, help='how deep should we split into individual file-names',default=0 );
parser.add_argument("-min_vc", default=0.5, type=float, help='in the last stage of autoNOE-Rosetta assignments with a volume contribution of less than 0.5 are usually ignored' )

library.add_standard_args( parser )
args = parser.parse_args()


class PeakCounter:
	def __init__(self,indent=0):
		self._indent=indent
		self._VALUE_COL=80
	def count_peak(self,p):
		pass
	def show_counter(self):
		pass
	def set_indent(self,val):
		self._indent=val
	def get_indent(self):
		return self._indent
	def str_indented(self,str,offset = 0,fillchar='.'):
		indent=self._indent+offset
		s=' '*indent+str
		fill=self._VALUE_COL-len(s)
		s+=fillchar*fill+'| '
		return s

class FileNameCounter(PeakCounter):
	def __init__(self,indent):
		PeakCounter.__init__(self,indent)
		self.files={}
	def count_peak(self,p):
		file=p.info().experiment_id()
		self.files.setdefault(file,0)
		self.files[file]+=1
	def sum_count(self):
		sum = 0
		for v in self.files.itervalues():
			sum+=v
		return sum

	def __str__(self):
		s=''
		if args.split_level==None or (args.split_level+1)*2>self.get_indent():
			for f,v in self.files.iteritems():
				s+='\n'+self.str_indented('in '+f,offset=2,fillchar=' ')+'%8d'%v
		return s

class SimplePeakCounter(PeakCounter,FileNameCounter):
	def __init__(self,indent,name):
		PeakCounter.__init__(self,indent)
		FileNameCounter.__init__(self,indent)
		self._name=name
		self._count=0

	def __str__(self):
		s=''
		s+=CompoundCounter.__str__(self)
		s+=FileNameCounter.__str__(self)
		return s

	def count_peak(self,p):
		if self.applies(p):
			self._count+=1
			FileNameCounter.count_peak(self,p)
			return True
		return False

	def applies(self,p):
		pass

	def __str__(self):
		s=self.str_indented(self._name)
		s+='%8d'%self._count
		s+=FileNameCounter.__str__(self)
		return s

	def count(self):
		return self._count


class CompoundCounter(PeakCounter,FileNameCounter):
	def __init__(self,indent,name):
		PeakCounter.__init__(self,indent)
		FileNameCounter.__init__(self,indent)
		self._reasons=[]
		self._name=name

	def append( self, counter ):
	  counter.set_indent(self.get_indent()+2)
		self._reasons.append( counter )

	def count_peak(self,p):
		applied=False
		for r in self._reasons:
			applied=r.count_peak(p) or applied
		if applied:
			FileNameCounter.count_peak(self,p)
		return applied

	def __str__(self):
		s=''
		ct=0
		for r in self._reasons:
			ct+=r.count()
		s+=self.str_indented(self._name)
		s+='%8d'%ct
		s+=FileNameCounter.__str__(self)
		for r in self._reasons:
			s+='\n%s'%r
		return s

class UnassignedPeaks(CompoundCounter):
	def __init__(self,indent):
		name='unassigned'
		CompoundCounter.__init__(self,indent,name)

class PickedCounter(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'picked')
	def applies(self,p):
		return True

class ZeroIntensityCounter(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'zero intensity')
	def applies(self,p):
		return p.volume()==0

class AssignedCounter(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'assigned')
	def applies(self,p):
		return p.nassign()>0 and not 'eliminated' in p.comment

class NoInitialAssignment(SimplePeakCounter):
	def __init__(self,indent=0):
		SimplePeakCounter.__init__(self,indent,'without assignment possibility')
	def applies(self,p):
		return p.nassign()==0

class Eliminated(SimplePeakCounter):
	def __init__(self,why,indent=0):
		SimplePeakCounter.__init__(self,indent,'eliminated due to '+why)
		self._why=why
	def applies(self,p):
		return p.nassign()>=0 and 'eliminated' in p.comment and self._why in p.comment

class UniqueAssignments(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'with unique assignment')
	def applies(self,p):
		return p.nassign()==1 and not 'eliminated' in p.comment

class MultipleAssignments(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'with >=5 initial assignments')
	def applies(self,p):
		return p.nassign()>=5 and not 'eliminated' in p.comment


class IntraResidueAssignments(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'with intraresidue assignment')
	def applies(self,peak):
		if 'eliminated' in peak.comment or peak.nassign()==0: return False
		min_seq_delta=100
		for a in peak:
			min_seq_delta=min(min_seq_delta,abs( a.atom(1).resid()-a.atom(2).resid() ))
		return min_seq_delta==0

class ShortRangeAssignments(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'with sequential assignment |i-j|=1')
	def applies(self,peak):
		if 'eliminated' in peak.comment or peak.nassign()==0: return False
		min_seq_delta=100
		for a in peak:
			min_seq_delta=min(min_seq_delta,abs( a.atom(1).resid()-a.atom(2).resid() ))
		return min_seq_delta==1

class MediumRangeAssignments(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'with medium-range assignment 1<|i-j|<5')
	def applies(self,peak):
		if 'eliminated' in peak.comment or peak.nassign()==0: return False
		min_seq_delta=100
		for a in peak:
			min_seq_delta=min(min_seq_delta,abs( a.atom(1).resid()-a.atom(2).resid() ))
		return min_seq_delta > 1 and min_seq_delta <5

class LongRangeAssignments(SimplePeakCounter):
	def __init__(self,indent):
		SimplePeakCounter.__init__(self,indent,'with long-range assignment |i-j|>=5')
	def applies(self,peak):
		if 'eliminated' in peak.comment or peak.nassign()==0: return False
		min_seq_delta=100
		for a in peak:
			min_seq_delta=min(min_seq_delta,abs( a.atom(1).resid()-a.atom(2).resid() ))
		return min_seq_delta >= 5

class ViolationCounter(SimplePeakCounter):
	def __init__(self,lo,hi=None,indent=0):
		if not hi:
			msg='above %3.1f A'%lo;
		else:
			msg='between %3.1f and %3.1f A'%(lo,hi);
		SimplePeakCounter.__init__(self,indent,msg)
		self._cutoff=lo
		self._ignore=hi
	def applies(self,p):
		if p.nassign()==0: return False
		start=p.comment.find('DistViol')
		if start<0: return False
		start=p.comment.find('>',start)
		end=p.comment.find('A',start+1)
		if start<0: return False
		viol=float(p.comment[start+1:end])
		return 'DistViol' in p.comment and viol >= self._cutoff and (not self._ignore or viol < self._ignore)

class AllViolationCounter(CompoundCounter):
	def __init__(self,indent):
		name='eliminated due to violations (DistViol)...'
		CompoundCounter.__init__(self,indent,name)
		self.append( ViolationCounter(0,0.5) )
		self.append( ViolationCounter(0.5,2) )
		self.append( ViolationCounter(2,5) )
		self.append( ViolationCounter(5) )

initial_peak_stats=[]
peak_stats=[]
indent=2
peak_stats.append(PickedCounter(indent))
peak_stats.append(ZeroIntensityCounter(indent))
peak_stats.append(AssignedCounter(indent))
initial_peak_stats.append(MultipleAssignments(indent))

unassigned_counter=UnassignedPeaks(indent)
unassigned_counter.append( NoInitialAssignment() )
unassigned_counter.append( Eliminated('Network') )
unassigned_counter.append( Eliminated('MinPeakVol') )
unassigned_counter.append( Eliminated('DistViol') )
unassigned_counter.append( Eliminated('MaxAssign') )
peak_stats.append( unassigned_counter )

assign_stats=[]
assign_stats.append( UniqueAssignments(indent) )
assign_stats.append( IntraResidueAssignments(indent) )
assign_stats.append( ShortRangeAssignments(indent) )
assign_stats.append( MediumRangeAssignments(indent) )
assign_stats.append( LongRangeAssignments(indent) )
assign_stats.append( AllViolationCounter(indent) )


crosspeaks=noesy.read_peak_files(args.peaks)

import itertools
for peak in crosspeaks:
#	print peak
	for stat in initial_peak_stats:
		stat.count_peak(peak)

	if args.min_vc and args.min_vc>0:
		op=copy.copy(peak)
		peak.clear_assignments()
		for a in op:
			vc = a.volume_contribution()
			if vc < args.min_vc:
				continue
			peak.add_assignment(a)

	for stat in itertools.chain(peak_stats,assign_stats):
		stat.count_peak(peak)

print 'Peaks:'
for stat in itertools.chain(peak_stats,initial_peak_stats):
	print stat
print 'Assignments:'
for stat in assign_stats:
	print stat

