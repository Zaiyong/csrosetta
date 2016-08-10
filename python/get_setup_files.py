#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

from os.path import exists
from os.path import basename
from os import path,listdir
import traceback
import sys

try:
	import automatic_setup
	from library import LibException, MissingInput
	import library

except ImportError as exc:
	traceback.print_exc(exc)
	print "\ncall 'source %s/init'"%path.dirname(__file__)
	print "before using the toolbox or put this call into your .bashrc"
	exit()


def get_setup_files(target_name,method_name,label_name,option,args):
	option_files=[]
### toolbox library

	try:
		targets=listdir(args.target_prefix)
		target=automatic_setup.TargetDir( target=target_name, prefix=args.target_prefix )
		file_dir=target.dir()
		setup=automatic_setup.Setup(target,method_name,label_name)
		options=setup.load_options()
		for r in options:
			if r[0]==option:
				files=r[1]
				break
		for r in files:
			file_name=r.split('/')[-1]
			option_files.append([ option,file_dir+'/'+r,file_name])
	except LibException as inst:
		if args.traceback:
			print traceback.print_exc(inst )
		else:
			print traceback.print_exception(sys.exc_type, sys.exc_value, None)
	return option_files

# class unit_test:
# 	def __init__(self):
# 		a=get_setup_files('casd_gmr137','autoNOE','standard','frags')
# 		print a
# #unit_test()
