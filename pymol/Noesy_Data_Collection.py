#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

from PeakList import PeakList

class Noesy_Data_Collection:
	def __init__(self, nafn):#nafn means noesy assignment file name
		self._nafn=nafn
		self._peaklists=[]
		#self._store_lines=[]
	def __str__(self):
		return 'I am a Noesy data collection file with %d peaklists' %self.npeaklist()

	def npeaklist(self):
		return len(self._peaklists)

	def add_peaklist(self,peaklist):
		self._peaklists.append(peaklist)

	@classmethod
	def read_from_file(obj,file):
		na_file=open(file)
		filelist=na_file.readlines()
		#self._store_lines=filelist
		nfn_index=[]
		for i in range (0, len(filelist)-1):
			if filelist[i].split()[0][0]=='#' and filelist[i+1].split()[0]=='#FILENAME':
				nfn_index.append(i)
		nfn_index.append(len(filelist))
		obj=Noesy_Data_Collection(file)
		for j in range(0,len(nfn_index)-1):
			pl=PeakList.read_from_lines(filelist[nfn_index[j]:nfn_index[j+1]])
			obj.add_peaklist(pl)
		return obj


# class noesy_data_test():
# 	noesy=Noesy_Data_Collection.read_from_file('/Users/zak/csrosetta3/pymol/NOE_out.dat')
# 	print noesy
# noesy_data_test()
