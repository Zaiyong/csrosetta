##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

import string
from glob import glob
#from sys import argv,stderr,exit
#import sys
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
#import automatic_setup


import traceback


#print "ABRELAX: ", method_path, flag_lib, method_name
sub_method_code = flag_lib+"/methods/_denovo_base/denovo_options.py"
if exists( sub_method_code ):
    exec open(sub_method_code, 'r' )
else:
    print "CANNOT FIND METHOD CODE %s"%sub_method_code
    exit()

# definition of options for the method RASREC
if 'group' in locals():
	group.add_argument('-cs',help='add chemical shift score to final output',metavar='<file>.tab')
	group.add_argument("-fix_topol", help="provide topology file (obtained with r_pdb2top from pdb-file" );
if 'run_group' in locals():
	run_group.set_defaults(cycle_factor=10.0)
	run_group.add_argument("-relax",  action='store_true', dest='relax', help="relax the structures after abinitio sampling ")
	run_group.add_argument("-norelax",  action='store_false', dest='relax', help="relax the structures after abinitio sampling ")
	run_group.add_argument("-nstruct", help="how many decoys should be produced", default=10000 );
#group.add_argument("-cycle_factor", help="the length of each abinitio stage (in terms of monte-carlo cycles) can be increase/decreased using this flag", default=10, type=float );


class AbrelaxMethod(DenovoBaseMethod):

	def __init__(self,name,path):
		print "AbrelaxMethod: ", name
		DenovoBaseMethod.__init__(self,name,path)
		self.non_file_options.append('nstruct')
		self.non_file_options.append('relax')
		self.option2dir['cs']=self.nmr_sub_dir
		self.option2dir['fix_topol']='structural_knowledge'

	def make_target_flags(self, run, setup, filename, flags, subs  ):
		DenovoBaseMethod.make_target_flags( self, run, setup, filename, flags, subs )
		args=self.get_args()
		if args.cs:
			flags.write("-evaluation:chemical_shifts %s chem_shifts\n"%setup.cm_path(args.cs));

		if args.fix_topol:
			raise StubbedOut("-fix_topol flag is not implemented yet")

		flags.write("-out:nstruct %d\n"%int(args.nstruct)); # not in Rasrec !

	def setup_file_library( self ):
		DenovoBaseMethod.setup_file_library( self )
		args=self.get_args()
		fl = self.file_library
		if args.relax:
			fl.add("flags", "flags_denovo", "@flags_fullatom\n")


method = AbrelaxMethod(method_name, method_path)
