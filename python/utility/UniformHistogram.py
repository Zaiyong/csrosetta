#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
import math

import unittest

##class for a histogram with equal sized bins
## input: hist=list(float)
class UniformHistogram():
	def __init__(self, hist, low, high):
		self.hist=hist
		self.low=low
		self.high=high
		self.nbins=len(self.hist)
		self.step=float(high-low)/float(self.nbins)
		self._normalize()

	def __str__(self):
		return 'Hist[%5.2f,%5.2f] %s'%(self.low,self.high,str(self.hist))

	def index(self,x):
		if x<self.low or x>=self.high: raise KeyError('x out of bounds [%8.3f,%8.3f)'%(self.low,self.high))
		return int(math.floor((x-self.low)/self.step))

	def bin_center(self,index):
		return self.step*0.5+index*self.step+self.low

	def _normalize(self):
		psum=sum(self.hist)*self.step
		try:
			self.hist=[x/psum for x in self.hist]
		except ZeroDivisionError:
			print 'division by zero in normalize: ',self
			raise
	def probability(self,x):
		try:
			return self.hist[self.index(x)]
		except KeyError:
			return 0.0

	def first_moment(self):
		mean=0
		for i,prob in enumerate(self.hist):
			mean+=self.bin_center(i)*prob
		return mean*self.step

	def integral(self,low,high):
		sum=0
		#get bounds (using the bin-center)
		start=int(math.floor((low-(self.low+self.step*0.5))/self.step))
		stop=int(math.ceil((high-(self.low+self.step*0.5))/self.step))
		if start<0: start=0
		if stop>self.nbins: stop=self.nbins
		#integrate over selected bins

		for i in xrange( start, stop ):
			sum+=self.hist[i]
		return sum*self.step


class UniformHistogramTestCase(unittest.TestCase):
	def setUp(self):
		data=[1.00, 4.00, 27.00, 89.00, 331.00, 688.00, 863.00, 1238.00, 1085.00, 991.00, 754.00, 764.00, 635.00, 280.00, 95.00, 84.00, 45.00, 11.00, 3.00, 12.00]
		self.hist=UniformHistogram(data,2.9,18.7)


	def test_mean(self):
		mean=self.hist.first_moment()
		should=9.96
		self.assertAlmostEqual(mean,should,2,
													 'mean value of UniformHistogram is %5.3f (true= %5.3f)'%(mean,should))
	def test_integral(self):
		limit=[(0,100),(5,10)]
		results=[1,0.537]
		for lim,res in zip(limit,results):
			val=self.hist.integral(lim[0],lim[1])
			self.assertAlmostEqual(val,res,2,'integral on [%5.2f,%5.2f] is %5.3f (true= %5.3f)'%(lim[0],lim[1],val,res))

def UniformHistogramTestSuite():
	suite = unittest.TestLoader().loadTestsFromTestCase(UniformHistogramTestCase)
	return suite

if __name__ == '__main__':
	unittest.main()

