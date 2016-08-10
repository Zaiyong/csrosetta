#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from Peak import Peak
from PeakList import PeakList, PeakCollection
from AssignmentCollection import AssignmentCollection
from PeakMatches import PeakMatch, AtomicPeakMatch
#import test_collection
#import test_strips
#obtain a subset of items from an iterable randomly
def random_items(iterable, k=1):
	from random import random, shuffle
	if not k:
		for m in iterable:
			yield m
		return

	result = [None,] * k
	for i, item in enumerate(iterable):
		if i < k:
			result[i] = item
		else:
			j = int(random() * (i+1))
			if j < k:
				result[j] = item
	#remove elements with None
  result = [x for x in result if x ]
	shuffle(result)
	for m in result:
		yield m
