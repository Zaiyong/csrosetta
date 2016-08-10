#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from assignment.noesy import Atom
from Distribution import Distribution
from GaussianDistribution import GaussianDistribution
from math import exp,sqrt
import unittest
from numpy import linspace

class SumGaussianDistribution(Distribution):
	def __init__(self,weights,gaussians):
		self._gaussians=gaussians #these are GaussianDistribution
		self._weights=weights

	def probability(self,x):
		prob=0
		S=sum(self._weights)
		for r,w in zip(self._gaussians,self._weights):
			prob+=w*r.probability(x)
		return prob/S

	def mean(self):
		mean=0
		S=sum(self._weights)
		for r,w in zip(self._gaussians,self._weights):
			mean+=w*r.mean()
		return mean/S

	def std(self):
		std=0
		S=sum(self._weights)
		for r,w in zip(self._gaussians,self._weights):
			std+=w*r.std()**2
		return sqrt(std)/S

	def min(self,threshold=0.95):
		assert threshold<1,'threshold should be smaller than 1'
		min=self.mean()
		mean=self.mean()
		integration=0
		S=sum(self._weights)
		while True:
			for r,w in zip(self._gaussians,self._weights):
				integration+=1.0*w*r.integral(min,mean)/S+1.0*w*r.integral(mean,2*mean-min)/S
			if integration>=threshold:
				break
			else:
				min=min-0.01*self.std()
				integration=0
		return min

	def max(self,threshold=0.95):
		assert threshold<1,'threshold should be smaller than 1'
		max=self.mean()
		mean=self.mean()
		integration=0
		S=sum(self._weights)
		while True:
			for r,w in zip(self._gaussians,self._weights):
				integration+=1.0*w*r.integral(mean,max)/S+1.0*w*r.integral(2*mean-max,mean)/S
			if integration>=threshold:
				break
			else:
				max=max+0.01*self.std()
				integration=0
		return max

	def integral(self,low,high):
		assert high>=low,'up bound of the integration should not be smaller than down bound'
		integration_value=0
		S=sum(self._weights)
		for r,w in zip(self._gaussians,self._weights):
			integration_value+=w*r.integral(low,high)
		return integration_value/S



	def append(self, weight, gaussian ):
		self._weights.append(weight)
		self._gaussian.append(gaussian)



class SumGaussianDistributionTestCase(unittest.TestCase):
	def setUp(self):
		gaussians=[ GaussianDistribution(1,0.2),GaussianDistribution(-2,1),GaussianDistribution(0.7,0.5) ]
		weights=[1,2,3]
		self.sum_gaussian_distribution=SumGaussianDistribution(weights,gaussians)

	def test_probability(self):
		x=linspace(-4,4,10)
		y=[self.sum_gaussian_distribution.probability(r) for r in x ]
		real_y=[0.0180,  0.0717,  0.1297,  0.1066,  0.0687,  0.3638,  0.2623,  0.0039,  0.0000,  0.0000]
		for t,s,r in zip(x,y,real_y):
			self.assertAlmostEqual(r,s,3,'Sum distribution of [1,0.2] [-2,1] [0.7,0.5] with weights [1,2,3]  at position x=%5.3f should be %5.3f but we get %5.3f (which is wrong)'%(t,r,s))

	def test_mean(self):
		mean=self.sum_gaussian_distribution.mean()
		real_mean=-0.150
		self.assertAlmostEqual(mean,real_mean,3,'mean value of Sum distribution of [1,0.2] [-2,1] [0.7,0.5] with weights [1,2,3] should be %5.3f but we get %5.3f (which is wrong)'%(real_mean,mean))

	def test_min_and_max(self):
		bound=[self.sum_gaussian_distribution.min(),self.sum_gaussian_distribution.max()]
		real_bound=[-3.037,2.737]
		for s,r in zip(bound,real_bound):
			self.assertAlmostEqual(r,s,2,' Sum distribution of [1,0.2] [-2,1] [0.7,0.5] with weights [1,2,3] at position x=%5.3f should be bound but we get %5.3f (which is wrong)'%(r,s))

	def test_integral(self):
		x_integral=[[-1.5,-1.2],[0.4,0.6],[0.3,1.8]]
		y_integral=[self.sum_gaussian_distribution.integral(r[0],r[1]) for r in x_integral]
		real_y_integral=[0.0323,0.0780,0.5573]
		for x,s,r in zip(x_integral,y_integral,real_y_integral):
			self.assertAlmostEqual(r,s,3,'the integration of N(1,0.2) distribution from position %5.3f to position %5.3f should be %5.3f but we get %5.3f (which is wrong)'%(x[0],x[1],r,s))

def SumGaussianDistributionTestSuite():
	suite = unittest.TestLoader().loadTestsFromTestCase(SumGaussianDistributionTestCase)
	return suite

if __name__ == '__main__':
	unittest.main()
