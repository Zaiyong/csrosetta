#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'

from Atom import Atom
import math

class UnFolded:
	def __call__( self, freq ):
		return freq
	def is_folded(self):
		return False

	def inverse( self, freq, ref_freq ):
		return freq

	def __str__( self ):
		return "UNFOLDED"


class FoldWindow:
	def __init__( self, lb, ub ):
		self._lb=lb
		self._ub=ub
		self._ws=ub-lb
		assert self._ws>0

	def __call__( self, freq ):
		while self._lb>freq:
			freq+=self._ws
		while self._ub<freq:
			freq-=self._ws
		return freq

	def __str__( self ):
		return "FOLD: [ %8.3f %8.3f ] - sw %8.3f"%(self._lb,self._ub,self._ws)

	def inverse( self, freq, ref_freq ):
#		diff=freq-ref_freq
#		sign=math.copysign(1, diff)
#		return ref_freq+sign*(abs(diff)%self._ws)
		return ref_freq + (freq-self(ref_freq))

	def is_folded(self):
		return True

UNFOLDED=UnFolded()


#A resonance is an entry in a .prot file:
# atom, and frequency including tolerance
# through inheritance a Resonance can also be a RangeResonance that
# matches to a range of resonances
class ResonanceBase:
	#private data members
 	#id, atom
	def __init__( self, id=-1, atom=None, error=999 ):
		self._id = id  #int
		self._atom = atom
		self._error = error
		self._intensity = 1
		self._floats = None
		self._reslist = None
		self.ambiguity = None
	#getters ---

	#tolerance/error
	def error(self):
		return self._error

	#atom
	def atom(self):
		return self._atom

	#atom name
	def name(self):
		return self._atom.name()

	#atom residue number
	def resid(self):
		return self._atom.resid()

	#id in .prot file
	def id(self):
		return self._id

	def intensity(self):
		return self._intensity

	def set_intensity( self, set):
		self._intensity = set

	#setters -------------
	def set_atom(self, atom):
		self._atom=atom
		return self

	def set_id( self, label ):
		self._id = label
		return self

	def set_floats( self, setting, reslist ):
		self._floats= setting
		self._reslist = reslist

	def _unpack_numbers( self, crosspeak, dim, attr ):
		freq=getattr(crosspeak,attr)(dim)
		tol=getattr(crosspeak.spin_info(dim),attr.replace('resonance','tolerance'))
		folder=getattr(crosspeak.spin_info(dim),attr.replace('resonance','folder'))
		return freq,tol,folder

	def peak_match(self, crosspeak, dim, attr='proton_resonance', threshold=1 ):
		freq, tol, folder = self._unpack_numbers( crosspeak, dim, attr )
		return self.match( freq, tol, folder, threshold)

	def peak_pmatch(self, crosspeak, dim, attr='proton_resonance' ):
		freq, tol, folder = self._unpack_numbers( crosspeak, dim, attr )
		return self.pmatch( freq, tol, folder )

	#for probability matching: return similar to boltzmann-factor: exp( -E**2/2 ), where E is the error in units of tolerance/sigma
	def pmatch(self, freq, tol, folder=UnFolded() ):
		return math.exp(-self._float_match( freq, tol, folder)*0.5 ) # normierung nicht so gut hier: /math.sqrt(2*math.pi*sigma2)

	#return score between 0 and 1, 1 for perfect match.
	def _match_diff2( self, freq, tol, folder=UnFolded() ):
		assert False, 'BaseClass:match should never be called'
		return 0.0

	def _float_match( self, freq, tol, folder ):
		if self._floats:
			best_match = 1000000
			for f in self._floats:
				best_match = min( best_match, self._reslist[ f ]._match_diff2( freq, tol, folder ) )
			return best_match
		else:
			return self._match_diff2( freq, tol, folder)

	#returns True if squared difference (in units of tolerance/sigma) is below 1
	def match(self, freq, tol, folder=UnFolded(), threshold=1 ): #threshold is given in sigma
		return (self._float_match( freq, tol, folder) < threshold*threshold)

	def __str__(self):
		return '%(_id)8d %(_atom)s '%self.__dict__+self._freq_str()+self._tail_str()

	def _freq_str( self ):
		return 'baseclass'

	def _tail_str( self ):
		return ""

	#class for resonances with a single frequency determined up to a small error
class Resonance(ResonanceBase):
	#private data
	def __init__( self, id=-1, atom=None, freq=0, error=0 ):
		ResonanceBase.__init__( self, id, atom, error )
		self._freq = freq

	def freq(self):
		return self._freq

	def set_freq(self,freq):
		self._freq=freq

	def _match_diff2(self, freq, tol, folder=UnFolded() ):
		sigma2=tol*tol+self._error*self._error
		x=freq-folder(self._freq)
		return x*x/sigma2

	def _tail_str(self):
		return ' %(_error)8.3f'%self.__dict__+ResonanceBase._tail_str(self)

	def _freq_str( self ):
		return '%(_freq)8.3f'%self.__dict__


	#class for Atoms whose resonance frequency can be narrowed down to a range of frequencies (e.g., HB: 0-2ppm)
class RangeResonance(ResonanceBase):
	#private data members
	#freq_lb, freq_ub

	def __init__( self, id=-1, atom=None, freq_lb=0, freq_ub=0, error=0 ):
		ResonanceBase.__init__( self, id, atom, error )
		self._freq_lb = freq_lb
		self._freq_ub = freq_ub

	def _tail_str(self):
		return ' %(_error)8.3f'%self.__dict__+ResonanceBase._tail_str(self)

	def _freq_str( self ):
		return '[ %(_freq_lb)8.3f %(_freq_ub)8.3f ]'%self.__dict__

	def freq(self):
		return (self._freq_lb, self._freq_ub)


	def _match_diff2( self, freq, tol, folder=UnFolded() ):
		sigma2=tol*tol+self._error*self._error
		ub=folder(self._freq_ub)
		lb=folder(self._freq_lb)

		#interval bigger than folding window...always a match.
		if (( self._freq_ub != ub or self._freq_lb !=lb ) and lb<ub ) or ( self._freq_ub != ub and self._freq_lb !=lb ):
			return 0

		#one bound was folded and now lb>ub... interval is relatively small and crossed one of the folding window's edges.
		if ( self._freq_ub != ub or self._freq_lb !=lb ) and lb>ub:
			if freq>=lb or freq<=ub: return 0
			x=min(freq-ub,lb-freq)
		elif freq>=lb and freq<=ub: return 0
		elif freq>ub:
			x=freq-ub
		elif freq<lb:
			x=lb-freq
		return x*x/sigma2



def unit_test():
	r=RangeResonance(1, Atom('H',5), 120.0, 128.3, 0.02)
	assert r.match(122,0.3,UNFOLDED, 2)
	assert not r.match(119,0.3,UNFOLDED,2)
	assert r.match(119.8,0.3,UNFOLDED,2)
	print r

	r2=RangeResonance(2, Atom('CA',5),50, 55, 0.02)
	assert r2.match(30,0.3,FoldWindow(30,52))
	assert r2.match(51,0.3,FoldWindow(30,52))
	assert r2.match(33.1,0.3,FoldWindow(30,52))

	r3=RangeResonance(3,Atom('CA',5),20,40,0.02)
	assert r3.match(30,0.3,FoldWindow(30,52))
	assert r3.match(40,0.3,FoldWindow(30,52))
	assert r3.match(50,0.3,FoldWindow(30,52))

	r4=RangeResonance(3,Atom('CA',5),37,40,0.02)
	assert not r4.match(36,0.3,FoldWindow(30,52))
	assert not r3.match(41,0.01,FoldWindow(30,52))

