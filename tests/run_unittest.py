#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

####################################
###
###  UNIT_TEST  HOW-TO
###
###  run all tests: ./run_unittest.py
###  run individual tests from directory of test: python -m unittest -v test_score_cache
###  run individual tests from top-level directory: python -m unittest -v assignment.scoring.test_score_cache
###
###
###  How to make unittests :
###
###  1) don't edit this file
###  2) write test_XXX.py modules and put them into the python/...-subtree of csrosetta3
###  3) if you need data files for your unit-test in python/XXX/YYY
###       then you should store all data files in ./unit_test_data/XXX/YYY/data/
###  4) to get the data path for you're testing module use
###            data_path = basic.get_unittest_data()
###  5) see example for correct unittest in python/assignment/scoring/test_scores.py
###
###

import unittest

#have to have this here, so that ModuleArgumentParser modules can be loaded
from basic.options._module_variables import _unit_testing
_unit_testing.append('UNIT_TEST')

if __name__ == '__main__':
	suite = unittest.TestLoader().discover('../python',pattern='test*py',top_level_dir='../python')
	unittest.TextTestRunner(verbosity=2).run(suite)
