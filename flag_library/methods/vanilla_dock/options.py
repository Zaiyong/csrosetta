##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

import string
from glob import glob
from os import popen,system,fdopen,mkdir,makedirs
from os import dup2,path
from os.path import exists
from operator import add
from math import sqrt
from os.path import basename
import argparse
import sys
import shutil

### toolbox library

import traceback


#print "ABRELAX: ", method_path, flag_lib, method_name
sub_method_code = flag_lib+"/methods/_docking_base/docking_options.py"
if exists( sub_method_code ):
    exec open(sub_method_code, 'r' )
else:
    print "CANNOT FIND METHOD CODE %s"%sub_method_code
    exit()

# definition of options for the method RASREC
# add_argument("-nstruct", help="how many decoys should be produced", default=1 );


class VanillaDockingMethod(DockingBaseMethod):

	def __init__(self,name,path):
		print "VanillaDockingMethod: ", name
		DockingBaseMethod.__init__(self,name,path)
		self.non_file_options.append('nstruct')

	def make_target_flags(self, run, setup, filename, flags, subs  ):
		DockingBaseMethod.make_target_flags( self, run, setup, filename, flags, subs )
		args=self.get_args()
#		flags.write("-nstruct %d\n"%int(args.nstruct)); # not in Rasrec !

method = VanillaDockingMethod(method_name, method_path)
