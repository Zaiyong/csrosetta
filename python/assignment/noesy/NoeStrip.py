#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#from PeakAssignment import PeakAssignment
from CrossPeakInfo import CrossPeakInfo
from Atom import Atom
from CrossPeak import CrossPeak
from Resonance import Resonance, RangeResonance
from CrossPeakList import CrossPeakList


def _group_frequencies(freqs,tolerance):
	def _compare( f1, f2, tolerance ):
		return abs(f1-f2)<=tolerance

	freqs=sorted(freqs)
	groups=[]
	for f1 in freqs:
		matched=False
		for g in groups:
			for f2 in g:
				if _compare( f1,f2,tolerance ):
					g.append(f1)
					matched=True
					break
			if matched: break
		if not matched:
			g=[f1]
			groups.append(g)
	mid_points=[ sum(g)/len(g) for g in groups ]
	return mid_points


class NoeStrip:
	class Undefined( Exception ):
		pass

	def __init__(self,dim,proton_resonance,label_resonance=None,matches=[]):
		self._dim=dim
		self._proton=proton_resonance
		self._label=label_resonance
		self._matches=matches

	def resid(self):
		if not ( self._proton.resid()>0 or self._label.resid()>0 ):
			raise Undefined
		if self._proton.resid()>0 and self._label.resid()>0:
			assert self._proton.resid() == self._label.resid()
		if self._proton.resid()>0: return self._proton.resid()
		return self._label.resid()

	def label_name(self):
		return self._label.name()

	def label(self):
		return self._label

	def proton_name(self):
		return self._proton.name()

	def proton(self):
		return self._proton

	def indirect_dim(self):
		return self._dim%2+1

	def indirect_protons(self):
		return sorted( [cp.proton_resonance(self.indirect_dim()) for cp in self._matches ] )

	def pmatch(self, strip):
		#match direct proton of self to indirect protons of strip
		indirect=self.indirect_dim()
#		print 'match self: %s'%self
#		print 'against strip: %s'%strip
#		for cp in strip._matches:
#			print self._proton.freq_str(), cp.proton_resonance(self.indirect_dim()), self._proton.peak_pmatch( cp, indirect )
		if strip._matches:#sometimes, the strip is empty
			return max([ self._proton.peak_pmatch( cp, indirect ) for cp in strip._matches ])
		else:
			return 0.0

	def __str__(self):
		s='Strip(%3s %s |%3s %s )'%(self._proton.name(), self._proton._freq_str(),self._label.name(), self._label._freq_str())
		s+=' '.join(['%8.3f '%f for f in self.indirect_protons()])
		return s

	def long_repr(self):
		s='[ %(_proton)s , %(_label)s ]:\n'%self.__dict__
		s+='matched peaks: %d\n'%len(self._matches)
		if len(self._matches):
			s+='\n'.join(['%s'%cp for cp in self._matches ])
		return s

	@classmethod
	def generate_strip(obj,cpl,dim,proton,label=None):
		#cpl: CrossPeakList
		#dim 1,2 refers to the entry in CrossPeak we want to match
		#proton, label of type ResonanceBase
		#
		assert proton
		assert proton.atom().elem() == 'H'
		matches=[]
		for cp in cpl.iter_peaks():
			if not proton.peak_match( cp, dim ): continue
			if label: assert cp.spin_info(dim).label
			if label and not label.peak_match( cp, dim, 'label_resonance' ): continue
			matches.append(cp)
		obj=NoeStrip(dim,proton,label,matches)
		return obj

	@classmethod
	def generate_strip_family(obj,cpl,dim,proton,label=None,tolerance=0.03):
		broad_strip=NoeStrip.generate_strip(cpl,dim,proton,label)
		proton_freqs=[ cp.proton_resonance( dim ) for cp in broad_strip._matches ]
		mid_points=_group_frequencies(proton_freqs,tolerance)
#		print 'matched fres: ',proton_freqs
#		print 'mid points: ', mid_points
		strips=[]
		for freq in mid_points:
			strip_proton=Resonance(freq=freq,error=0.03,atom=proton.atom())
			strips.append(NoeStrip.generate_strip(cpl,dim,strip_proton,label))
		return strips

	@classmethod
	def generate_strip_family_2D(obj,cpl,dim,proton,label=None,tolerance=[0.03,0.3]):
		broad_strip=NoeStrip.generate_strip(cpl,dim,proton,label)
		freqs_2D=[ [ cp.proton_resonance( dim ),cp.label_resonance( dim ) ] for cp in broad_strip._matches ]
		proton_mid_points=_group_frequencies(freqs_2D[0],tolerance[0])
		label_mid_points=_group_frequencies(freqs_2D[1],tolerance[1])
#		print 'matched fres: ',proton_freqs
#		print 'mid points: ', mid_points
		strips=[]
		for freq_proton in proton_mid_points:
			strip_proton=Resonance(freq=freq_proton,error=0.03,atom=proton.atom())
			for freq_label in label_mid_points:
				strip_label=Resonance(freq=freq_label,error=0.3,atom=label.atom())
				strips.append(NoeStrip.generate_strip(cpl,dim,strip_proton,strip_label))
		return strips

#test
def unit_test():
	s='''
# Number of dimensions 3
#FILENAME resort_c-ali
#FORMAT xeasy3D
#INAME 1 c
#INAME 2 H
#INAME 3 h
#CYANAFORMAT cHh
#FOLD 1 40 50
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
   8   45.226    3.787    2.695  1 U 1.161E+05  0.000E+00  e 0
# Number of dimensions 3
#INAME 1 N
#INAME 2 h
#INAME 3 H
#CYANAFORMAT NhH
#TOLERANCE 0.3 0.04 0.03
       1  126.668    2.191    8.933 3 U 3.569e+05
       2  128.843    1.729   10.218 3 U 1.347e+05
       3  121.139    8.912    8.320 3 U 1.226e+05
       4  121.153    4.957    8.328 3 U 4.684e+05
       5  122.140    4.528    8.287 3 U 3.174e+05
       6  127.663    7.060    9.615 3 U 1.133e+05
       7  128.747    4.665    9.006 3 U 2.131e+05
       8  125.034    0.551    8.995 3 U 8.296e+05
       9  126.779    1.415    8.877 3 U 2.358e+05
      10  126.849    5.313    8.874 3 U 3.267e+05
      11  126.789    5.061    8.884 3 U 1.498e+06
      12  127.589    7.305    9.613 3 U 3.501e+05
      13  127.592    9.014    9.614 3 U 1.559e+05
      14  125.343    4.345    9.645 3 U 4.988e+06
      15  125.347    2.227    9.645 3 U 5.825e+05
      16  125.349    2.070    9.645 3 U 1.234e+06
      17  125.347    9.647    9.645 3 U 3.147e+07
      18  125.345    5.434    9.645 3 U 1.111e+06
      19  125.340    6.939    9.646 3 U 1.502e+05
      20  125.342    7.947    9.646 3 U 3.795e+05
      21  125.347    2.858    9.645 3 U 3.123e+06
      22  125.336    8.771    9.645 3 U 2.717e+05
      23  125.360    0.759    9.648 3 U 2.315e+05
      24  125.377    8.100    9.643 3 U 2.091e+05
      25  125.355    4.541    9.645 3 U 1.142e+06
      26  125.337    2.694    9.645 3 U 4.860e+05
      27  125.335    9.330    9.647 3 U 4.600e+05
      28  115.444    9.648    6.938 3 U 2.061e+05
      29  115.438    9.647    7.945 3 U 4.623e+05
      30  115.457    5.431    7.947 3 U 1.844e+05
      31  115.473    6.935    7.945 3 U 2.110e+07
      32  115.469    7.946    6.936 3 U 2.060e+07
      33  115.457    2.858    6.936 3 U 1.124e+06
      34  115.464    2.860    7.944 3 U 3.063e+06
      35  115.426    1.145    6.936 3 U 2.811e+05
      36  115.539    1.141    7.943 3 U 2.189e+05
      37  121.826    9.649    8.769 3 U 2.303e+05
      38  121.844    5.432    8.767 3 U 6.358e+06
      39  121.859    2.860    8.767 3 U 1.454e+06
      40  121.835    4.953    8.768 3 U 9.074e+05
      41  121.824    1.400    8.768 3 U 5.927e+05
'''
	print 'NOE-STRIP TEST...'
	from StringIO import StringIO
	pseudo_file=StringIO(s)
 	pl=CrossPeakList.read_from_stream(pseudo_file)
	print pl
	print 'generate strips'
	noestrip_test=NoeStrip.generate_strip(pl,1,Resonance(1,Atom('H',110),2.695,0),Resonance(2,Atom('CB',110),35.226,0))
	print noestrip_test

	strips=NoeStrip.generate_strip_family(pl,1,RangeResonance(1,Atom('H',110),7,8,0.02),Resonance(2,Atom('N',110),115.444,0.3))
	print 'Family: '
	print '\nNew Strip '.join(['%s'%s for s in strips])

