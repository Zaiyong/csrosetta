#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import unittest
import itertools
from numpy import linspace
from assignment.noesy import Resonance
from assignment.noesy import Atom
from utility.SchmidtQ import quantity_function
def best_match(res_list,freq_list):
	assert len(res_list)==len(freq_list),'the length of res_list and freq_list should be the same, however they are %s(res_list) and %s(freq_list), seperately'%(len(res_list),len(freq_list))
	#linspace(1,len(res_list),len(res_list))
	index=[]
	length=len(res_list)
	for i in range(0,length):
		index.append(i)
	permutes=itertools.permutations(index,length)
	max_match_prob=-1
	selected_permute=[]
	for permute in permutes:
		match_prob=1.0
		for i in range(0,length):
			match_prob=match_prob*res_list[permute[i]].pmatch(freq_list[i],0)
		if match_prob>max_match_prob:
			max_match_prob=match_prob
			selected_permute=permute
	return max_match_prob

class spectrum_util_TestCase(unittest.TestCase):
	def setUp(self):
		self.res1=Resonance(1,Atom('N',1),121.1,0.2)
		self.res2=Resonance(2,Atom('H',1),8.13,0.02)
		self.res3=Resonance(3,Atom('CA',1),59.3,0.2)
		self.res4=Resonance(4,Atom('HA',1),4.25,0.02)
		self.res5=Resonance(5,Atom('CB',1),31.6,0.2)
		self.res6=Resonance(6,Atom('HB2',1),2.33,0.02)
		self.res_list0=[self.res1,self.res2]
		self.res_list1=[self.res2,self.res3,self.res4]
		self.res_list2=[self.res3,self.res4,self.res5,self.res6]
		self.freq_list0=[8.15,121.3]
		self.freq_list1=[59.3,4.26,8.15]
		self.freq_list2=[2.31,4.25,57.3,25.4]
		self.right_result=[0.3679,0.5353,2.4526e-231]
	def test_best_match(self):
		self.assertAlmostEqual(self.right_result[0],best_match(self.res_list0,self.freq_list0),4,'the best match prob of res_list and freq_list should be %8.4f but we got %8.4f'%(self.right_result[0],best_match(self.res_list0,self.freq_list0)))
		self.assertAlmostEqual(self.right_result[1],best_match(self.res_list1,self.freq_list1),4,'the best match prob of res_list and freq_list should be %8.4f but we got %8.4f'%(self.right_result[1],best_match(self.res_list1,self.freq_list1)))
		self.assertAlmostEqual(self.right_result[2],best_match(self.res_list2,self.freq_list2),4,'the best match prob of res_list and freq_list should be %8.4f but we got %8.4f'%(self.right_result[2],best_match(self.res_list2,self.freq_list2)))

def spectrum_util_TestSuite():
	suite=unittest.TestLoader().loadTestsFromTestCase(spectrum_util_TestCase)
	return suite

if __name__ == '__main__':
	unittest.main()
