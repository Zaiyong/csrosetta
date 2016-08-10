#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import library
from basic import Tracer
import copy

tr=Tracer('scoring.methods')
#############
##
## class ScoreCache
##
## used to store scores on a per-atom or per-assignment basis.
## you'll find this as attribute .scores in an AssignmentCollection after scoring.
## it also allows to 'look-at' scores in case one wants to base deciscions for future assignments on these
##
## the cache provides
##     atomic { Atom --> { score_name --> val } }
##     atomic_matches { AtomicPeakMatche --> { score_name --> val } }
##
## Note, one could probably make these dictionaries the other way round
## i.e., atomic { score_name --> { Atom --> val } }
## but i thought as a user of these one is mostly interested in getting all scores to a given Atom
## for the score-caching it probably would make more sense to have them the other way around.
##
def _copy_score_cache_dict(target_dict,input_dict,score_type):
	for key,scores in input_dict.iteritems():
		target_cache=target_dict.setdefault(key,{})
		try:
			target_cache[score_type]=copy.copy(scores[score_type])
		except KeyError:
			#nothing to copy
			pass

class ScoreCache(object):
	def __init__(self):
			#score caching -- do not manipulate out of this class - could put this into self.score_cache.notify, self.score_cache.atomic, self.score_cache.atomic_match
		self.notify={}
		self.atomic={}
		self.atomic_matches={}
		self.assignments={}
  #called by AssignmentCollection when a PeakMatch has been added
	def notify_add(self, assignment_collection, peak_match ):
		for id,score_method in self.notify.itervalues():
			score_method.notify_add( assignment_collection, peak_match )

  #called by AssignmentCollection when a PeakMatch has been removed
	def notify_remove(self, peak_match ):
		for id,score_method in self.notify.itervalues():
			score_method.notify_remove(self, peak_match )

	def __copy__(self):
		newone = type(self)()
		newone.notify=self.notify.copy()
		newone.atomic={}
		newone.atomic_matches={}
		newone.assignments={}
		for id,score_method in self.notify.itervalues():
			try:
				score_method.copy_cache( newone, self )
			except NotImplementedError:
				del newone.notify[score_method.name]
#		newone.atomic=copy_score_cache_dict(self.atomic)
#		newone.atomic_matches=copy_score_cache_dict(self.atomic_matches)
#		newone.assignments=copy_score_cache_dict(self.assignments)
		return newone

	def simple_cache_copy_from(self, old_cache, score_type, cache_types=[] ):
		for cache_type in cache_types:
			_copy_score_cache_dict( getattr(self,cache_type), getattr(old_cache,cache_type), score_type )
			# 	if cache_type=='atomic':
# 					_copy_score_cache_dict( self.atomic, old_cache.atomic, score_type )
# 				if cache_type=='atomic_matches':
# 					_copy_score_cache_dict( self.atomic_matches, old_cache.atomic_matches, score_type )
# 				if cache_type=='assignments':
# 					_copy_score_cache_dict( self.assignments, old_cache.assignments, score_type )



######################
##
##  use this class to store a score in the ScoreCache
##
##  you can overload this to store more complex things in the ScoreCache
##  in this case you should make sure that the total score is still returned via .score
##  Example: SymmetryScore
class ScoreItem(object):
	__slots__='score'
	def __init__(self,score):
		self.score=score

	def __str__(self):
		return str(self.score)
#################################
##
##  use this class for the options (parameters) of you're ScoreMethod
##
class ScoreOptions(object):
	def __setattr__(self, *args):
		try:
			frz=self._frozen
		except AttributeError:
			super(ScoreOptions,self).__setattr__( *args )
			return
		if not frz:
			super(ScoreOptions,self).__setattr__( *args )
		else:
			raise TypeError('Tried to change a ScoreFunction Option after it has been used on an Assignment.\n'
											'To resolve use options.clone() to get a fresh instance and then change the options of the fresh instance')

	__delattr__ = __setattr__

	def freeze(self):
		super(ScoreOptions,self).__setattr__('_frozen',True)

	def clone(self):
		vals=self.__dict__
		del vals['_frozen']
		obj=ScoreOptions()
		obj.__dict__=vals
		return obj

######################
##
## class ScoreMethod
##
## derive custom ScoreMethods from this class, each class should have a unique name
##
## Important! A ScoreMethod should have no parameters. All parameters should live under
## score.options.param1 ... score.options.paramN   you can choose the names .param freely
## put all options here, and nuke the instance if you want to change options
## this allows the ScoreFunction to nuke the score-cache if this is the case
##
class ScoreMethod(object):
	def __init__(self,name, options=None):
		self.name=name
		self.options=options

	def __call__(self,assignments):
		return self.score(assignments)

	def notify_add(self, assignments, new_assignment):
		pass

	def notify_remove(self, assignments, new_assignment):
		pass

	def copy_cache(self, new_cache, old_cache ):
		raise NotImplementedError

	def __key__(self):
		return self.name

	def __str__(self):
		return self.name

	def _get_id(self):
		my_id=id(self)
		if self.options:
			option_id=id(self.options)
			self.options.freeze()
		else:
			option_id=None
		#we would return (my_id,option_id) but a ScoreMethod should always do the same thing unless options change...
		return (my_id,option_id)

	##IMPORTANT: always call super class Method at beginning of scoring
	def score(self,assignments):
		self._check_cache(assignments)

	def _check_cache(self,assignments):
		try:
			if assignments.scores.notify[self.name][0]!=self._get_id():
		    tr.Info('nuking scores for %s...'%str(self))
				self.nuke_scores(assignments, 'atomic_matches' )
				self.nuke_scores(assignments, 'atomic' )
				self.nuke_scores(assignments, 'assignments' )
		except KeyError:
			#if our score was not in there, nothing to be nuked.
			pass
		except AttributeError:
			assignments.scores=ScoreCache()
	  #put call-backs into place
	  self.install_notify_handler(assignments)

	#invoked if a new instance of ScoreMethod (with the same .name attribute) or a new instance of option is used
	#as option is frozen that SHOULD be the only way how ScoreMethods can change there behaviour
  #nuking will not happen very often, as score-methods and/or score-options shouldn't change too much
	#however, we might also just nuke everything, rather than fixing all these dictionaries...
	def nuke_scores(self,assignments,attribute):
		cache=assignments.__getattribute__(attribute)
		for entry in cache.itervalues():
			try:
				del entry[key]
			except KeyError:
				pass
		pass

	def install_notify_handler(self,assignments):
		assignments.scores.notify[self.name]=(self._get_id(),self)
