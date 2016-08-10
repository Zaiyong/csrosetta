#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#author: Oliver Lange

#as part of general assignment module

###############
### class AssignmentCollection
# this class stores the complete state of assignments
# this is done in multiple redundant dictionaries to allow fast queries of the state of assignments
#
# the dictionaries are:
#       atom2assignment(s)	  (which assignments does this atom have ? )
# 	peak2assignment(s)	  (which assignments does this peak have ? )
# 	assignment2assignment  (do I have this assignment?)
#
# additionaly we store the peak_collection and molecule for convience
# (i.e., refer to the same old peak-lists which are probably somewhere else also
# which will not cost anything but helps us to keep track of things...
#
# the class foremost behaves like a set of assignments:
#    ass_coll[PeakMatch(None,(atom1,atom2,atom3),None)]   retrieves an assignment with atom1,atom2,atom3
# mmh, danger--> there might be multiple assignments (different peaks) to the same triple of atoms...
#
# TODO: correct copy behavior
#    the score-cache needs to be copied the moment we make add/remove changes
#    same with the assignment dictionaries...
#
from collections import defaultdict
from PeakMatches import PeakMatch
from scoring.methods.ScoreMethod import ScoreCache

from basic import Tracer
import library
import copy

tr=Tracer('assignment.collection')
def copy_dict_of_list(d):
	new_dict=defaultdict(list)
	for key, value in d.iteritems():
		new_dict[key]=copy.copy(value)
	return new_dict

class AssignmentCollection(object):
	def __init__(self, peak_collection, molecule):
		#these lists should maybe become sets ? -- faster removal
		self._atom2assignment=defaultdict(list)
		self.peak2assignment=defaultdict(list)
		self.atomtuple2assignment=defaultdict(list)
		self._unassigned_peaks=set(peak_collection)
		self.peak_collection=peak_collection
		self.molecule=molecule
		self.scores=ScoreCache()
		self._remove_set=set()
		self._add_set=set()
		self._lock_commit=False

	def __copy__(self):
		self.commit()
		newone = self.copy_without_scores()
		newone.scores=copy.copy(self.scores)
		return newone

	def copy_without_scores(self):
		self.commit()
		newone = type(self)(self.peak_collection,self.molecule)
		#copy the dictionaries
		newone._atom2assignment=copy_dict_of_list(self._atom2assignment)
		newone._unassigned_peaks=copy.copy(self._unassigned_peaks)
		newone.peak2assignment=copy_dict_of_list(self.peak2assignment)
		newone.atomtuple2assignment=copy_dict_of_list(self.atomtuple2assignment)
		return newone

	def clear_scores(self):
		self.scores=ScoreCache()

	def commit(self):
		if self._lock_commit: return
		self._lock_commit=True
		try:
			for m in self._remove_set:
				self._remove( m )
			for m in self._add_set:
				self._add( m )
			self._remove_set.clear()
			self._add_set.clear()
		finally:
			self._lock_commit=False

	def remove( self, peak_match ):
		if peak_match in self._add_set:
			self._add_set.remove( peak_match )
			return False
		else:
			if not peak_match in self.peak2assignment[peak_match.peak]: raise ValueError
			self._remove_set.add( peak_match )
			return True

	def add( self, peak_match ):
		if peak_match in self._remove_set:
			self._remove_set.remove( peak_match )
			return False
		else:
			self._add_set.add( peak_match )
			return True

	def _add( self, peak_match ):
		self.atomtuple2assignment[peak_match.peak_match].append(peak_match)
		for atomic in peak_match:
			self._atom2assignment[atomic.atom].append(atomic)
		self.peak2assignment[peak_match.peak].append(peak_match)
		self.scores.notify_add( self, peak_match )
		self._unassigned_peaks.discard(peak_match.peak)

	def _remove( self, peak_match ):
		#remove from peak-lists
		peaklist = self.peak2assignment[peak_match.peak]
		peaklist.remove(peak_match)
		if len(peaklist)==0:
			del self.peak2assignment[peak_match.peak]
			self._unassigned_peaks.add(peak_match.peak)
		#remove from tuple-lists
		tuple_list=self.atomtuple2assignment[peak_match.peak_match]
		tuple_list.remove(peak_match)
		if len(tuple_list)==0:
			del self.atomtuple2assignment[peak_match.peak_match]
		#remove from atom lists
		for atomic in peak_match:
			atomic_list=self._atom2assignment[atomic.atom]
			atomic_list.remove(atomic)
			if len(atomic_list)==0:
				del self._atom2assignment[atomic.atom]
		#update scores
		self.scores.notify_remove( peak_match )

	@property
	def by_atom(self):
		self.commit()
		return self._atom2assignment

	@property
	def by_peak(self):
		self.commit()
		return self.peak2assignment

	@property
	def unassigned_peaks(self):
		self.commit()
		return self._unassigned_peaks

	def __iter__(self):
		self.commit()
		for ass in self.atomtuple2assignment.itervalues():
			for a in ass:
				yield a

	def __getitem__(self,tuple):
		self.commit() #cannot do this here, as score-methods might use get_item
		#use get here to avoid creating empty lists where a non-existing assignment has been queried
		#get returns None if no item present. raise an exception in this case
		val=self.atomtuple2assignment.get(tuple)
		if not val: raise KeyError(tuple)
		return val

	def write_to_stream(self, fd, unassigned=True, write_scores=True ):
		self.commit()
		if unassigned:
			for peak in sorted(self.peak_collection):
				peak_str='%-45s '%(str(peak)+':')
				fd.write(peak_str)
				first=True
				for assignment in self.by_peak[ peak ]:
					if not first:
						fd.write(' '*len(peak_str))
					first=False
					fd.write('%s'%assignment.match_str())
					if write_scores: self.write_score_line( fd, assignment )
					fd.write('\n')
				if first: fd.write('\n')
		else:
			for assignment in sorted(self,key=lambda ass: ass.peak):
				fd.write('%s\n'%assignment)

	def write_score_line(self, fd, assignment ):
		fd.write(' SCORES: ')
		sum=defaultdict(int)
		#read-out atomic-assignment scores
		for atomic_assignment in assignment:
			try:
				cache=self.scores.atomic_matches[atomic_assignment]
				for type,score in cache.iteritems():
					sum[type]+=score.score
			except KeyError:
				pass
		for type,val in sum.iteritems():
			fd.write(' %s=%8.3f'%(type,val))
		#read-out assignment scores
		try:
			cache=self.scores.assignments[assignment]
			for type,score in cache.iteritems():
				fd.write(' %s=%8.3f'%(type,score.score))
		except KeyError:
			pass

	@classmethod
	def from_hard_assignments(obj,peak_collection,molecule):
		obj=AssignmentCollection(peak_collection,molecule)
		for name,peak_list in peak_collection.experiments.iteritems():
			for peak in peak_list:
				if len(peak.hard_assignments)==0: continue
				try:
					for match in peak.matches( molecule ):
						obj.add(match)
						#code below was wrong as it wouldn't fill in the hidden atoms into the assignment
#					for ass in peak.hard_assignments:
#					try:

	#					if None in ass: continue
	#					rule_match=peak.rule.translate_peak_match_to_full_match(ass,molecule)
	#					peak_match=peak.rule.translate_full_match_to_peak_match(rule_match)
	#					if None in peak_match:
	#							tr.Debug('Cannot use hard-assignment: Unknown atom in ',ass)
	#							continue
	#					obj.add(PeakMatch(rule_match,peak_match,peak))
				except Exception:
					tr.Warning('[WARNING] pre-existing assignment of Peak %s has problems and is ignored'%peak)
					tr.Warning('Details: ',library.exception2str())
		obj.commit()
		return obj

	@classmethod
	def from_saved_state(obj, peaks, molecule, fd):
		class NoMatch( Exception ):
			pass
		from Peak import Peak
		from chemical import Atom
		from strips import Strip
		def read_matched_atoms( molecule, line ):
			atoms=line.replace('SCORES','').strip().rstrip(')').lstrip('(').split(',')
			if len(atoms)==1 and len(atoms[0])==0: raise NoMatch
			if len(atoms)==0: return NoMatch
			peak_match=[]
			for atom_str in atoms:
				tags=atom_str.split()
				atom=molecule.atom(Atom(tags[0],int(tags[1])))
				peak_match.append(atom)
			return peak_match

		def read_peak_match( fast_peaks, molecule, line, last_peak ):
			line=line.replace('SCORES:','|')
			if last_peak and not ':' in line:
				peak=last_peak
				match_str=line.split('|')[0]
			else:
				tags=line.split(':')
				read_peak=Peak.from_string(tags[0])
				peak=fast_peaks[read_peak]
				match_str=tags[1].split('|')[0]
			peak_match=tuple(read_matched_atoms( molecule, match_str ) )
			rule_match=peak.rule.translate_peak_match_to_full_match( peak_match, molecule )
			return PeakMatch( rule_match, peak_match, peak ), peak

		def read_strip_match( fast_strips, molecule, line, last_strip ):
			if last_strip and not ':' in line:
				strip=last_strip
				match_str=line
			else:
				tags=line.split(':')
				dtags=tags[0].split()
				if dtags[0]!='Strip': raise ValueError('line %s is not a Strip'%line.strip() )
				strip_key=Strip.key_from_string(tags[0])
				print strip_key
				strip=fast_peaks[strip_key]
				match_str=tags[1]
			peak_match=tuple(read_matched_atoms( molecule, match_str ))
			return StripMatch( peak_match, strip ), strip
		#End helper routines

		obj=AssignmentCollection( peaks, molecule )
		#get dictionary of peaks or strips for fast lookup
		fast_peaks={}
		for peak in peaks:
			fast_peaks[peak]=peak
		#read file line by line
		last_peak=None
		for line in fd:
			try:
				try:
					match, last_peak = read_peak_match( fast_peaks, molecule, line, last_peak )
					obj._add( match )
				except ValueError:
					match, last_peak = read_strip_match( fast_peaks, molecule, line, last_peak )
					obj._add( match )
			except NoMatch:
				continue
			except:
				print 'problems reading line: %s'%line.strip()
				raise
		return obj
