#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from Peak import Peak

from basic import Tracer

tr=Tracer('assignment.peaks')


#small storage class to keep the experimental peak-lists
#and expected Peaks
class PeakCollection:
	def __init__(self):
		self.experiments={}
		self.expected=PeakList('expected')

	def add_experiment(self,peak_list):
		self.experiments[peak_list.name]=peak_list

	def iterpeaks( self ):
		for peak_list in self.experiments.itervalues():
			for peak in peak_list:
				yield peak

	@classmethod
	def from_peak_files( peak_collection, files, ignore=False ):
		from os import path
		from basic import Tracer
		tr = Tracer('assignment.io')
		peak_collection=PeakCollection()
		for file in files:
			name=path.splitext(path.basename(file))[0]
			tr.Info('read peak-list %s from %s...'%(name,file))
			peak_collection.add_experiment( PeakList.read_from_stream( name, open(file,'r'), ignore ) )
			tr.Info('done reading... ')
		return peak_collection

	#same as iterpeaks
	def __iter__(self):
		for peak_list in self.experiments.itervalues():
			for peak in peak_list:
				yield peak

class PeakList(list):

	def __init__(self,name='Unknown'):
		self.name=name
		super(PeakList,self).__init__()
		self._max_id=0

	def __str__(self):
		return 'PeakList %s with %5d peaks'%(self.name,len(self))

	def append( self, peak ):
		if not peak.id:
			self._max_id+=1
			peak.id=self._max_id
		else:
			self._max_id=max(self._max_id, peak.id )
#		if peak.peak_list_name != '' and peak.peak_list_name != self.name:
#			raise ValueError('adding peak from peak-list %s to peak-list %s. replace name by "" first'%(peak.peak_list_name, self.name ))
		if not peak.peak_list_name or peak.peak_list_name == '':
			peak.peak_list_name=self.name
		super(PeakList,self).append(peak)

	def write_to_stream( self, fd, file_adaptor=None ):
		rule_dict={}
		for p in self:
			rule_dict.setdefault(id(p._rule),[]).append(p)
		for peaks in rule_dict.itervalues():
			if len(peaks)==0: continue
			rule=peaks[0]._rule
			#write rule
			if file_adaptor: file_adaptor.write_rule( fd, rule )
			else: fd.write('#RULE\n%s\n#END_RULE\n'%rule)

			#write peaks belonging to this rule
			for peak in peaks:
				if file_adaptor: file_adaptor.write_peak( fd, peak )
				else: fd.write('%s\n'%peak)

	#generate from FileObject
	#read method using the old-style CrossPeak, CrossPeakInfo FILE-IO
	#this reads old-style NOESY PeakLists
	#TODO: new FILE-IO with generic headers for all Rules
	@classmethod
	def read_from_stream(obj, list_name, file, ignore_assignments = False, resonances=None ):

#little helper code
		def crosspeak_to_peak(cp, rule, filter,list_name):
			import chemical
			freq=tuple([cp._freq[i] for i in filter])
			peak=Peak( freq, rule, cp._id,list_name )
			if cp._assignments:
				peak.hard_assignments=[]
				for assignment in cp._assignments:
					new_ass=[]
					for i in filter:
						a=assignment._atoms[i]
						new_ass.append(chemical.Atom(a._name,a._resid) )
					peak.hard_assignments.append(tuple(new_ass))
			return peak

		from assignment.noesy	import CrossPeak, CrossPeakInfo
		from rules import NoesyRule
		obj=PeakList( list_name )
		info=None
		cp_lines=[]
		for l in file:
			l=l.rstrip()
			if len(l)==0: continue
			#no indentation -- make new CrossPeak out of previously collected lines.
			if len(cp_lines) and len(l[0:20].lstrip())>0:
				cp=CrossPeak.read_from_lines(cp_lines, ignore_assignments, resonances )
				obj.append(crosspeak_to_peak(cp,rule,freq_filter,obj.name))
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
				rule,pseudo4D_column=NoesyRule.from_old_cross_peak_info(info)
				freq_filter=range(0,max(rule.dims))

				if pseudo4D_column:
					import copy
					nodes=copy.deepcopy(rule.nodes)
					for i,n in enumerate(nodes):
						if n.dim>pseudo4D_column[1]: n.dim-=1
						if i==pseudo4D_column[0]: n.dim=None
					corrected_rule=NoesyRule(((nodes[0],nodes[1]),(nodes[2],nodes[3])))
					rule=corrected_rule
					del freq_filter[freq_filter.index(pseudo4D_column[1]-1)]
			#a line for collection

			cp_lines.append(l)
		if len(cp_lines)>0:
			cp=CrossPeak.read_from_lines(cp_lines, ignore_assignments, resonances )
			obj.append(crosspeak_to_peak(cp,rule,freq_filter,obj.name))

		if info.spin(1).label=='C' or info.spin(2).label=='C':
			tr.Warning('C-dimension detected: checking frequencies to figure out type assignment: aliC or aroC or both ? ')
			for node in obj[0].rule.nodes:
				types=set()
				if node.element=='C' and node.dim:
					freqs=[p.freq[node.dim-1] for p in obj]
					low=min(freqs)
					high=max(freqs)
					tr.Info('C frequency varies between %8.3f and %8.3f'%(low,high))
					if low<100: types.add('aliC')
					if high>100: types.add('aroC')
					tr.Warning('setting '+str(types)+' for '+str(node))
					node.types=types
		return obj


def unit_test():
	print 'Test PeakList.py...'
	s='''
# Number of dimensions 3
#FILENAME resort_c-ali
#FORMAT xeasy3D
#INAME 1 c
#INAME 2 h
#INAME 3 N
#INAME 4 H
#CYANAFORMAT chN
#TOLERANCE      0.3    099    0.03 0.3
    1   40.932   1 4.805    0.522  1 U 2.114E+05  0.000E+00  e 0
    5   45.226   1 6.787    2.695  1 U 1.161E+05  0.000E+00  e 0  CB  110  H 2  HD2   22    HB2  110   #VC 0.000 #W 0.249 1.000 1.000 0.167 0.000 0.000 23.165  #d 5.89447 #eliminated: Network
                                                                 CB  110  H 2 HE3   36    HB2  110   #VC 0.899 #W 0.087 1.000 1.000 0.167 0.506 0.506 14.898  #d 5.89447 #eliminated: Network
                                                                 CB  110 H 2     QD   64    HB2  110   #VC 0.000 #W 0.194 10.000 1.000 0.167 0.000 0.000 22.100  #d 5.89447 #eliminated: Network
                                                                 CB  110  H 2 HE22   84    HB2  110   #VC 0.000 #W 0.477 1.000 1.000 0.167 0.000 0.000 16.346  #d 5.89447 #eliminated: Network
                                                                 CB  110  H 2  HE3   85    HB2  110   #VC 0.101 #W 0.087 1.000 1.000 0.167 0.057 0.057 9.344  #d 5.89447 #eliminated: Network
                                                                 CB  110 H 2    HZ2   85    HB2  110   #VC 0.000 #W 0.477 1.000 1.000 0.167 0.000 0.057 13.228  #d 5.89447 #eliminated: Network
   6   40.719   1 1.831    1.764  1 U 1.355E+06  0.000E+00  e 0
   7   40.847   1 0.934    1.743  1 U 1.674E+05  0.000E+00  e 0
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
 	pl=PeakList.read_from_stream('test',pseudo_file)
	print "freshly read cross-peak list..."
	print pl
	assert len(pl)==4, 'wrong number of peaks in list'
	assert pl[1]._id==5, 'wrong peak-id %d'%pl[1]._id
	assert len(pl[1]._hard_assignments)==6, 'wrong number of assignments for peak #5 ( should be 6, found %d )'%len(pl[1]._hard_assignments)

	print pl[0]._rule
	print '\n\nFILE OUTPUT: '
	import sys
	pl.write_to_stream(sys.stdout)

#unit_test()
