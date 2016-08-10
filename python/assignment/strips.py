#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from PeakList import PeakList
from Peak import Peak

class StripMatch( object ):
	def __init__(self, peak_match, strip ):
		self.strip=strip
		strip_match = tuple( peak_match[d] for d in strip.dims )
		self.strip_match=strip_match
		self._peak_match=peak_match

	def __str__(self):
		return str(self.strip)+': '+self.match_str()

	def match_str(self):
		return '('+', '.join(('%s'%atom for atom in self.strip_match if atom))+')'

	def __key__(self):
		#muss hier noch die Peak-Ids verhackstuecken....
		return (self.strip.__key__(), tuple([atom.__key__() if atom else None for atom in self.strip_match ]))

	def __hash__(self):
		return hash(self.__key__())

	def matches( self, protein, random_samples=None, frequency_matcher=None , distance_matcher=None ):
		for peak in self.strip:
			for match in peak.matches( protein, random_samples, frequency_matcher, distance_matcher, partial_match=self._peak_match ):
				yield match

	@property
	def rule(self):
		return self.strip.rule

	@property
	def freq(self):
		return self.strip.freq

	#for compatibility with AssignmentCollection
	@property
	def peak(self):
		return self.strip

	@property
	def peak_match(self):
		return self.strip_match

	def __iter__(self):
		for dim in range(0,len(self.strip_match)):
			yield AtomicStripMatch(self,dim+1)

	def __cmp__(self,other):
		if self.strip!=other.strip: return cmp( self.strip, other.strip )
		return cmp( self.peak_match, other.strip )

#use this class to refer to a single atom's assignment resulting from a PeakMatch
class AtomicStripMatch(object):
	__slots__=('freq','dim','peak_match')
	def __init__(self, strip_match, dim ):
		self.freq=strip_match.strip.freq[dim-1] #direct storage saves ca. 30% of lookup time
		self.dim=dim
		self.peak_match=strip_match

	@property
	def atom(self):
		return self.peak_match.strip_match[self.dim-1]

	@property
	def strip(self):
		return self.peak_match.strip

	@property
	def peak(self):
		return self.peak_match.strip

	def __str__(self):
		return 'Atomic( d=%d, w=%5.2f, %s )'%(self.dim,self.freq,self.atom)

	def __key__(self):
		return (self.strip.__key__(), self.dim, self.atom.__key__())

	def __eq__(self,other):
		if self.strip!=other.strip: return False
		return self.atom==other.atom

	def __hash__(self):
		return hash(self.__key__())

class Strip( PeakList ):
	def __init__( self, mask, dims, name='strip', id=-1 ):
		self.mask=mask
		self.dims=[ i for i,m in enumerate(mask) if m ]
		super(Strip,self).__init__( name )
		self._average_freqs=[]
		self.id=id

	def __str__(self):
		return super( Strip, self ).__str__().replace('PeakList','Strip')+' %5d'%self.id

	def matches( self, protein, random_samples=None, frequency_matcher=None ):
		if len(self)==0: return
		freq=list(self[0].freq)
		for w,d in zip(self.freq,self.dims):
			freq[d]=w
		peak=Peak(freq, self[0].rule)
		for match in peak.matches( protein, random_samples, frequency_matcher, match_mask=self.mask ):
			yield StripMatch( match.peak_match, self )

	@property
	def rule(self):
		if len(self)==0: return None
		return self[0].rule

	@property
	def freq(self):
		if len(self._average_freqs)==0 and len(self)>0:
			import numpy as np
			freqs=np.empty([len(self),len(self.dims)])
			for i,peak in enumerate(self):
				for j,d in enumerate(self.dims):
					freqs[ i, j ]=peak.freq[ d ]
			self._average_freqs=tuple(np.mean(freqs,axis=0))
		return self._average_freqs

	def __lt__(self,other):
		if self.name != other.name: return self.name < other.name
		if self.id and other.id:
			return self.id < other.id
		else:
			return tuple([peak.__key__() for peak in self]) < tuple([peak.__key__() for peak in other])

	def __cmp__(self,other):
		if self.name != other.name: return cmp(self.name,other.name)
		if self.id and other.id:
			return cmp(self.id,other.id)
		else:
			return cmp( tuple([peak.__key__() for peak in self]), tuple([peak.__key__() for peak in other]) )

  def __key__(self):
		if self.id:
			return (self.name, self.id)
#, self.freq )
		else:
			return (self.name, tuple([peak.__key__() for peak in self]))

	def __hash__(self):
		return hash( self.__key__() )

	def append( self, peak):
		self._average_freqs=[]
		super(Strip, self).append(peak)

	@classmethod
	def key_from_string( obj, line ):
		#can generate a Strip which is just a husk but good as key for lookup in a strip direcotry
		obj=Strip([],[])
		tags=line.split()
		if tags[0]!='Strip': raise ValueError
		if tags[2]!='with': raise ValueError
		if tags[4]!='peaks': raise ValueError
		obj.id=int(tags[5])
		obj.name=tags[1]
		return obj

def collect_strips_from_peak_collection( peak_collection ):
	coll_strips=[]
	id=1
	for plist in peak_collection.experiments.itervalues():
		for strips in collect_strips_from_peak_list( plist ):
			for strip in strips:
				coll_strips.append(strip)
				strip.id=id
				id+=1
	return coll_strips


def collect_strips_from_peak_list( peak_list ):
#	from scipy.cluster.hierarchy import linkage, fcluster
# in matlab clustering with Z=linkage(x,'average','seuclidean')
# and cutting into flat clusters with
# c=cluster( Z, 'cutoff', 0.02, 'criterion','distance')
# gave reasonable results.
	from scipy.spatial.distance import pdist
	from scipy.cluster.hierarchy import linkage, fcluster
	import numpy as np

	if len(peak_list)==0: return
	rule=peak_list[0].rule
	strip_name=peak_list.name+'_strip'
	strip_list=[]
	for mask in rule.spinsystem_match_masks():
		dims=[ i for i,m in enumerate(mask) if m ]
		freqs=np.empty([len(peak_list),len(dims)])
		for i,peak in enumerate(peak_list):
			for j,d in enumerate(dims):
				freqs[ i, j ]=peak.freq[d]
		Y=pdist( freqs, 'seuclidean' )
		Z=linkage( Y, method='average', metric='seuclidean' )
		cluster_indices=fcluster( Z, 0.03, criterion='distance', depth=0 )
		clusters=set(cluster_indices)
		strips=[]
		for i in range(0, len(clusters) ):
			strips.append(Strip(mask,dims,strip_name))
		for c,p in zip(cluster_indices,peak_list):
			strips[c-1].append(p)
		strip_list.append(strips)
	return strip_list
