#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

from CrossPeak import CrossPeak

class PeakList:
	def __init__(self,dim,nfn):#nfn means noesy file name
		self._dim=dim
		self._nfn=nfn
		self._crosspeaks=[]
		#self._store_lines=[]
	def __str__(self):
		return 'I am a peaklist with %d crosspeaks' %self.ncrosspeak()

	def ncrosspeak(self):
		return len(self._crosspeaks)

	def add_crosspeak(self,crosspeak):
		self._crosspeaks.append(crosspeak)

	@classmethod
	def read_from_lines(obj,lines):
		if lines[0].split()[0][0]=='#' and lines[1].split()[0]=='#FILENAME':
			#self._store_lines=lines
			for line in lines:
				if line.split()[0]=='#FILENAME':
					nfn=line.split()[1]
					continue
				if line.split()[0]=='#TOLERANCE':
					dim=len(line)-1
					break
			obj=PeakList(dim,nfn)
		else:
			print('wrong lines for one peaklist')
			assert(False)
		cp_line_num=[]#the lines having crosspeaks

#find the line-numbers where cross-peak begins
		for i in range (0,len(lines)):
			if lines[i].split()[0][0]<='9' and lines[i].split()[0][0]>='0':
				cp_line_num.append(i)
		cp_line_num.append(len(lines))
#for
		for j in range (0,len(cp_line_num)-1):
			cp=CrossPeak.read_from_lines(lines[cp_line_num[j]:cp_line_num[j+1]])
			obj.add_crosspeak(cp)
		return obj


# class peaklist_test():
# 	lines=['','','','','','','','','','','','','','','','','']
# 	lines[0]='# Number of dimensions 3'
# 	lines[1]='#FILENAME resort_c-ali'
# 	lines[2]='#FORMAT xeasy3D'
# 	lines[3]='#INAME 1 c'
# 	lines[4]='#INAME 2 H'
# 	lines[5]='#INAME 3 h'
# 	lines[6]='#CYANAFORMAT cHh'
# 	lines[7]='#TOLERANCE      0.3    0.04    0.03'
# 	lines[8]='    1   40.932    4.805    0.522  1 U 2.114E+05  0.000E+00  e 0'
# 	lines[9]='   5   45.226    6.787    2.695  1 U 1.161E+05  0.000E+00  e 0  CB  110    HD2   22    HB2  110   #VC 0.000 #W 0.249 1.000 1.000 0.167 0.000 0.000 23.165  #d 5.89447 #eliminated: Network'
# 	lines[10]='                                                                  CB  110    HE3   36    HB2  110   #VC 0.899 #W 0.087 1.000 1.000 0.167 0.506 0.506 14.898  #d 5.89447 #eliminated: Network'
# 	lines[11]='                                                                  CB  110     QD   64    HB2  110   #VC 0.000 #W 0.194 10.000 1.000 0.167 0.000 0.000 22.100  #d 5.89447 #eliminated: Network'
# 	lines[12]='                                                                 CB  110   HE22   84    HB2  110   #VC 0.000 #W 0.477 1.000 1.000 0.167 0.000 0.000 16.346  #d 5.89447 #eliminated: Network'
# 	lines[13]='                                                                  CB  110    HE3   85    HB2  110   #VC 0.101 #W 0.087 1.000 1.000 0.167 0.057 0.057 9.344  #d 5.89447 #eliminated: Network'
# 	lines[14]='                                                                  CB  110    HZ2   85    HB2  110   #VC 0.000 #W 0.477 1.000 1.000 0.167 0.000 0.057 13.228  #d 5.89447 #eliminated: Network'
# 	lines[15]='     6   40.719    1.831    1.764  1 U 1.355E+06  0.000E+00  e 0'
# 	lines[16]='     7   40.847    0.934    1.743  1 U 1.674E+05  0.000E+00  e 0'
# 	pl=PeakList.read_from_lines(lines)
# 	print pl

# peaklist_test()
