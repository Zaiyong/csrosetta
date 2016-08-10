#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from assignment.noesy import Atom
import library
from math import erf,fabs,sqrt,log
from utility.SchmidtQ import quantity_function
class Distribution:
	def __init__(self):
		pass
	def probability(self,x):
		raise library.StubbedOut()

				#stubbed out
				#return None
	def mean(self):
		raise library.StubbedOut()

	def std(self):
		raise library.StubbedOut()
        #return self._mean
	def max(self,threshold=0.95):
		raise library.StubbedOut()

	def min(self,threshold=0.95):
		raise library.StubbedOut()

	def integral(self,low,high):
		raise library.StubbedOut()


	def schmidt_q(self,value):
		mean=self.mean()
		std=self.std()
		x=float(value-mean)/std
		return quantity_function(x)

