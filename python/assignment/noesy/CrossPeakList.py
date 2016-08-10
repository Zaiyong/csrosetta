#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from CrossPeak import CrossPeak
from CrossPeakInfo import CrossPeakInfo

#holds all CrossPeaks of a project,
#possibly multiple files (i.e., with different headers)
#CrossPeakLists (read from single file) can be merged using CrossPeakList::merge() into a single list
class CrossPeakList:
	def __init__(self):#nfn means noesy file name
		self._peaks=[]
		self._headers=[]

	def __str__(self):
		s=''
		for info in self._headers:
		  s+='%s \n'%info
		return s[:-1]
	#	return 'I am a peaklist with %d crosspeaks' %self.npeaks()

	#return number of cross-peaks (or diagonal peaks)
	def ncrosspeak(self):
		return len(self._peaks)

	#return number of cross-peaks (or diagonal peaks)
	def npeaks(self):
		return len(self._peaks)

	#add another cross-peaks
	def add_crosspeak(self,crosspeak):
		self._peaks.append(crosspeak)

	def append(self, cp ):
		self._peaks.append(cp)

#return nth peak
	def peak(self,i):
		assert i<self.npeaks(),'the index out of peak range'
		return self._peaks[i]

	#assign all matching resonances to all cross-peaks
	def assign_resonances(self, resonances ):
		for c in self._peaks:
			try:
				c.assign_resonances( resonances )
			except KeyError:
				pass
	#merge multiple CrossPeakLists
	def merge(self, peak_list ):
		import types
		if isinstance( peak_list,CrossPeakList ):
			self._peaks.extend( peak_list._peaks )
			self._headers.extend( peak_list._headers )
		else:
			for i in peak_list:
				self.merge( i )


	def iter_peaks(self):
		for p in self._peaks:
			yield p

	def __iter__(self):
		return self.iter_peaks()

	def write_to_stream( self, fd ):
		current_info = None
		for p in self._peaks:
			if current_info != p.info():
				current_info = p.info()
				if current_info._header_lines:
					fd.write('\n'.join(current_info._header_lines)+'\n')
#				current_info.write(fd)
			fd.write('%s\n'%p)

	def write_split_files( self, prefix='', suffix='.peaks' ):
		current_info = None
		for p in self._peaks:
			if current_info != p.info():
				current_info = p.info()
				file=open(prefix+current_info._experiment_id+suffix,'w')
				if current_info._header_lines:
					file.write('\n'.join(current_info._header_lines))
					file.write('\n')
#				current_info.write(fd)
			file.write('%s\n'%p)

	#generate from FileObject
	@classmethod
	def read_from_stream(obj, file, ignore_assignments = False, resonances=None ):
		obj=CrossPeakList()
		info=None
		cp_lines=[]
		for l in file:
			l=l.rstrip()
			if len(l)==0: continue

			#no indentation -- make new CrossPeak out of previously collected lines.
			if len(cp_lines) and len(l[0:20].lstrip())>0:
#				print cp_lines
				cp=CrossPeak.read_from_lines(cp_lines, ignore_assignments, resonances )
				cp.set_info( info )
#				print 'P: ',cp
				obj.add_crosspeak(cp)
				cp_lines=[]

			#new header -- update info object
			if l[0]=='#':
				header_lines=[]
				for hl in file:
					if len(hl)==0 or hl[0]!='#': break
					header_lines.append(hl.rstrip())
				l=hl.rstrip()

#				print 'HEADER: ',header_lines
#				print l
				info=CrossPeakInfo.read_from_lines(header_lines)
#				print info
				obj._headers.append(info)

			#a line for collection
			cp_lines.append(l)
		if len(cp_lines)>0:
			cp=CrossPeak.read_from_lines(cp_lines, ignore_assignments, resonances )
			cp.set_info( info )
			obj.add_crosspeak(cp)
		return obj

# 	@classmethod
# 	def read_from_lines(obj,lines,ignore_assignments = False, resonances=None):
# 		print "get read_from_lines in CrossPeakList.py"
# 		obj=CrossPeakList()
# 		info=None
# 		head_start=0
# 		head_end=0
# 		cp_start=[]
# 		for i in range(0,len(lines)):
# 			tags=lines[i].split()
# 			if lines[i].find("Number of dimensions")>=0:
# 				head_start=i
# 			elif lines[i].find("#TOLERANCE")>=0:
# 				head_end=i
# 			elif tags[0][0]<='9' and tags[0][0]>='0':
# 				cp_start.append(i)
# 		info=CrossPeakInfo.read_from_lines(lines[head_start:head_end+1])
# 		for i in range(0,len(cp_start)-1):
# 			cp=CrossPeak.read_from_lines(lines[cp_start[i]:cp_start[i+1]], ignore_assignments, resonances )
# 			cp.set_info( info )
# 			obj.add_crosspeak(cp)
# 		cp=CrossPeak.read_from_lines(lines[cp_start[-1]:], ignore_assignments, resonances )
# 		cp.set_info( info )
# 		obj.add_crosspeak(cp)
# 		obj._headers.append(info)
# 		return obj

def unit_test():
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
    5   45.226    6.787    2.695  1 U 1.161E+05  0.000E+00  e 0  CB  110    HD2   22    HB2  110   #VC 0.000 #W 0.249 1.000 1.000 0.167 0.000 0.000 23.165  #d 5.89447 #eliminated: Network
                                                                 CB  110    HE3   36    HB2  110   #VC 0.899 #W 0.087 1.000 1.000 0.167 0.506 0.506 14.898  #d 5.89447 #eliminated: Network
                                                                 CB  110     QD   64    HB2  110   #VC 0.000 #W 0.194 10.000 1.000 0.167 0.000 0.000 22.100  #d 5.89447 #eliminated: Network
                                                                 CB  110   HE22   84    HB2  110   #VC 0.000 #W 0.477 1.000 1.000 0.167 0.000 0.000 16.346  #d 5.89447 #eliminated: Network
                                                                 CB  110    HE3   85    HB2  110   #VC 0.101 #W 0.087 1.000 1.000 0.167 0.057 0.057 9.344  #d 5.89447 #eliminated: Network
                                                                 CB  110    HZ2   85    HB2  110   #VC 0.000 #W 0.477 1.000 1.000 0.167 0.000 0.057 13.228  #d 5.89447 #eliminated: Network
   6   40.719    1.831    1.764  1 U 1.355E+06  0.000E+00  e 0
   7   40.847    0.934    1.743  1 U 1.674E+05  0.000E+00  e 0
'''

	sfold='''
# Number of dimensions 3
#FILENAME resort_c-ali
#FORMAT xeasy3D
#INAME 1 c
#INAME 2 H
#INAME 3 h
#CYANAFORMAT cHh
#TOLERANCE      1  0.04    0.03
#FOLD 1 40 50
    1   40.932    4.805    0.522  1 U 2.114E+05  0.000E+00  e 0
    5   45    6.787    2.695  1 U 1.161E+05  0.000E+00  e 0
   6   40.719    1.831    1.764  1 U 1.355E+06  0.000E+00  e 0
   7   40.847    0.934    1.743  1 U 1.674E+05  0.000E+00  e 0
'''

	from StringIO import StringIO
	pseudo_file=StringIO(s)
 	pl=CrossPeakList.read_from_stream(pseudo_file)
	print "freshly read cross-peak list..."
	print pl
	assert pl.npeaks()==4, 'wrong number of peaks in list'
	assert pl._peaks[1].id()==5, 'wrong peak-id %d'%pl._peaks[1].id()
	assert pl._peaks[1].nassign()==6, 'wrong number of assignments for peak #5 ( should be 6, found %d )'%pl._peaks[1].nassign()



	from Resonance import Resonance
	from ResonanceList import ResonanceList
	from Atom import Atom
	from StringIO import StringIO
	from sys import stdout
	pseudo_file2=StringIO(sfold)
 	cpl_fold=CrossPeakList.read_from_stream(pseudo_file2)

	my_res=Resonance(1, Atom( "CA", 1), 35, 0.3 )
	my_Hres=Resonance(1, Atom( "HA", 1), 2.695, 0.01 )
	my_hres=Resonance(1, Atom( "HZ", 2), 6.787, 0.01 )
	res_list=ResonanceList()
	res_list.add_resonance(my_res)
	res_list.add_resonance(my_Hres)
	res_list.add_resonance(my_hres)
	res_list.set_sequence('AA')
	cpl_fold.assign_resonances(res_list)
	print "cross-peak list after assigning a folded Resonance"
	print cpl_fold
	for c in cpl_fold._peaks:
		stdout.write('%s\n'%c)
