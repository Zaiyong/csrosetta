#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#just for now replace with something more useful later
#import ScoreMethodManager

import methods


class ScoreFunction:

	def __init__(self,**kwargs):
		self.weights={}
		self.methods={}
		self.set_weights(kwargs)

	def __str__(self):
			output=''
			for type,weight in self.weights.iteritems():
				output+='%20s'%type+'%20.3f'%weight
			return output

	def _evaluate(self, assignments):
		#run the individual score_methods
		scores=[]
		assignments.commit()
		for method in self.methods.itervalues():
			scores.append( method( assignments ) )
		return scores

	#expect a dictionarity SCORE_TYPE --> ScoreItem
	# as found in AssignmentColleciton.scores.XXX[locator]
	def evaluate_local_scores( self, items ):
		sum=0
		for score_type,score in items.iteritems():
			try:
				sum+=self.weights[score_type]*score.score
			except KeyError: #if we don't have the score_type in weights it is 0 contribution
				pass
		return sum

	def total_score(self,scores):
		return sum([s*w for s,w in zip(scores,self.weights.itervalues())])

	def local_score_str( self, state, item ):
		self( state )
		s=''
		try:
			cache=state.scores.atomic[item]
			for type, val in cache.iteritems():
				s+=' %s=%8.3f'%(type,val.score)
		except KeyError:
			pass
		try:
			cache=state.scores.atomic_matches[item]
			for type, val in cache.iteritems():
				s+=' %s=%8.3f'%(type,val.score)
		except KeyError:
			pass
		if s=='':
			raise KeyError
		return s

	def score_str(self,assignments):
		scores=self._evaluate(assignments)
		spacing=10
		fmt_score='%%%ds'%spacing
		fmt_number='%%%d.3f'%spacing
		score_head_out='SCORE:  '+fmt_score%'total'+' '.join([fmt_score%s for s in self.weights.keys()])
		weight_out=    'WEIGHT: '+fmt_score%('')+' '.join([fmt_number%f for f in self.weights.values()])
		score_out=     'SCORE:  '+fmt_number%self.total_score(scores)+' '.join([fmt_number%f for f in scores])
		return score_head_out+'\n'+weight_out+'\n'+score_out

	def print_scores(self,assignments):
		print self.score_str(assignments)

	def __call__(self,assignments):
		scores=self._evaluate(assignments)
		return self.total_score(scores)
	#manually change the weight of a score-term
	#calls upon ScoreMethodManager if new score-method has to be generated

	def set_weight(self,type,weight):
		from methods import ScoreMethodManager
		if weight==0:
			try:
				del self.methods[type]
				del self.weights[type]
			except KeyError:
				pass
		if weight!=0 and not type in self.methods.keys():
			self.methods[type]=ScoreMethodManager.generate_method(type)
			self.weights[type]=weight
		return self

	def set_weights(self,weight_dict):
		for type,weight in weight_dict.iteritems():
			self.set_weight(type,weight)
		return self

	def setw(self, **kwargs):
		self.set_weights(kwargs)

	def append(self,method,weight):
		type=method.name
		self.methods[type]=method
		self.weights[type]=weight
		return self
