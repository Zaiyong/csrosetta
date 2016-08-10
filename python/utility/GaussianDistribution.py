#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from assignment.noesy import Atom
from Distribution import Distribution
from os import environ
from PDB.Polypeptide import one_to_three,three_to_one
from math import pi,exp,erf,sqrt,fabs,log
import random
import unittest

#scipy might not be installed and many tools of the toolbox don't require this
try:
	from scipy.stats import norm
except ImportError:
	pass

class GaussianDistribution(Distribution):
	def __init__(self,mean,std):
		Distribution.__init__(self)
		self._std=std
		self._mean=mean

	def probability(self,freq):
		return exp(-0.5*((freq-self._mean)/self._std)**2)/(self._std*(2*pi)**0.5)

	def std(self):
		return self._std

	def mean(self):
		return self._mean

	def min(self,threshold=0.95):
		assert threshold<1,'threshold should be smaller than 1'
		min=self._mean
		while norm.cdf(2*self._mean-min,self._mean,self._std)-norm.cdf(min,self._mean,self._std)<threshold:
			min-=0.01*self._std
		return min

	def max(self,threshold=0.95):
		assert threshold<1,'threshold should be smaller than 1'
		max=self._mean
		while norm.cdf(max,self._mean,self._std)-norm.cdf(2*self._mean-max,self._mean,self._std)<threshold:
			max+=0.01*self._std
		return max

	def integral(self,low,high):
		assert high>=low,'up bound of the integration should not be smaller than down bound'
		return norm.cdf(high,self._mean,self._std)-norm.cdf(low,self._mean,self._std)



class GaussianDistributionTestCase(unittest.TestCase):
	def setUp(self):
		self.gaussian_distribution1=GaussianDistribution(1,0.2)
		self.gaussian_distribution2=GaussianDistribution(-2,1)
	def test_probability(self):
		from numpy import linspace
		x=linspace(-4,4,10)
		y1=[self.gaussian_distribution1.probability(r) for r in x ]
		y2=[self.gaussian_distribution2.probability(r) for r in x ]
		real_y1=[0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0421,  0.4974,  0.0000,  0.0000,  0.0000]
		real_y2=[0.0540,  0.2152,  0.3892,  0.3194,  0.1190,  0.0201,  0.0015,  0.0001,  0.0000,  0.0000]
		for t,s,r in zip(x,y1,real_y1):
			self.assertAlmostEqual(r,s,3,'N(1,0.2) distribution at position x=%5.3f should be %5.3f but we get %5.3f (which is wrong)'%(t,r,s))
		for t,s,r in zip(x,y2,real_y2):
			self.assertAlmostEqual(r,s,3,'N(-2,1) distribution at position x=%5.3f should be %5.3f but we get %5.3f (which is wrong)'%(t,r,s))

	def test_min_and_max(self):
		x1_bound=[self.gaussian_distribution1.min(),self.gaussian_distribution1.max()]
		x2_bound=[self.gaussian_distribution2.min(),self.gaussian_distribution2.max()]
		real_y1_bound=[0.608,1.392]
		real_y2_bound=[-3.96,-0.04]
		for s,r in zip(x1_bound,real_y1_bound):
			self.assertAlmostEqual(r,s,2,'N(1,0.2) distribution at position x=%5.3f should be bound but we get %5.3f (which is wrong)'%(r,s))
		for s,r in zip(x2_bound,real_y2_bound):
			self.assertAlmostEqual(r,s,2,'N(-2,1) distribution at position x=%5.3f should be bound but we get %5.3f (which is wrong)'%(r,s))

	def test_integral(self):
		x1_integral=[[0.5,1.2],[0.4,0.6],[1.3,1.8]]
		x2_integral=[[-3.5,1],[-4,-2.5],[-2,1]]
		y1_integral=[self.gaussian_distribution1.integral(r[0],r[1]) for r in x1_integral]
		y2_integral=[self.gaussian_distribution2.integral(r[0],r[1]) for r in x2_integral]
		real_y1_integral=[0.8351,0.0214,0.0668]
		real_y2_integral=[0.9318,0.2858,0.4987]
		for x,s,r in zip(x1_integral,y1_integral,real_y1_integral):
			self.assertAlmostEqual(r,s,3,'the integration of N(1,0.2) distribution from position %5.3f to position %5.3f should be %5.3f but we get %5.3f (which is wrong)'%(x[0],x[1],r,s))
		for x,s,r in zip(x2_integral,y2_integral,real_y2_integral):
			self.assertAlmostEqual(r,s,3,'the integration of N(-2,1) distribution from position %5.3f to position %5.3f should be %5.3f but we get %5.3f (which is wrong)'%(x[0],x[1],r,s))

def GaussianDistributionTestSuite():
	suite = unittest.TestLoader().loadTestsFromTestCase(GaussianDistributionTestCase)
	return suite

if __name__ == '__main__':
	unittest.main()

