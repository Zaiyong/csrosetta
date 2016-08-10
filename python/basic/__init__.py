#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import options
from Tracer import Tracer


def fix_unittest_args():
	from basic.options._module_variables import _unit_testing
	_unit_testing.append('UNIT_TEST')

def get_unittest_data():
	import os.path
	import inspect
	frm = inspect.stack()[1]
	module=frm[1]
	if module[0]!='/':
		module=os.getcwd()+'/'+module
	return os.path.dirname(module).replace('python','tests/unit_test_data')+'/data/'
