##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
## make mammoth structure alignments

import string
from glob import glob
from os import dup2,path
from os.path import exists
from os.path import basename
import argparse
import sys
import shutil

### toolbox library
from library import Tracer
import library
from automatic_setup import BasicMethod
import traceback

# definition of options for the method RASREC
if 'group' in locals():
#group=parser.add_argument_group(title='method-options', description='extra options for method')
	 group.add_argument("-partners", help="partner-string to specify which chains are docked", default="A_B" );
	 group.add_argument("-partners_file", help="partner file to specify which chains are docked");
	 group.add_argument("-pdb", help="start-pdb with both chains placed randomly");
#group.add_argument("-decoys", nargs='*',help="bunch of starting structures in silent-file format");
	 group.add_argument("-native", help="supply a native pdb for RMSD calculation" );
#	 group.add_argument("-cst_file",help="constraint file, using given name Encounter.cst" );
#group.add_argument("-out", help="file name for output");


if 'run_group' in locals():
	 run_group.add_argument("-protocol", choices=['rep_high','rep_cen','rep_high_rigid','rep_high_frozen'], help="choose the mode", default="rep_cen");
	 run_group.add_argument("-score",help="which score weight to use?, options: interchain_cen, docking, score12");
	 run_group.add_argument("-extra_score", choices=['docking_interface_score'],help="if to use docking_interface_score as extra_score");
	 run_group.add_argument("-min_score_score",help="replica hack for low-res score fuction sucks");
	 run_group.add_argument("-n_replica", help="how many replicas" );

tr = Tracer( "replica_docking_base_method" )

# method-based code
class ReplicaDockingBaseMethod(BasicMethod):
	 def __init__(self,name,path):
			BasicMethod.__init__(self,name,path)
			self.non_file_options.append('partners')
			self.non_file_options.append('run')
#			self.non_file_options.append('out')
			self.non_file_options.append('score')
			self.non_file_options.append('extra_score')
			self.non_file_options.append('n_replica')
			self.non_file_options.append('min_score_score')
			self.option2dir['pdb']='start'
#			self.option2dir['decoys']='start'
			self.option2dir['native']='native'
			self.option2dir['partners_file']='partners_file'
#			self.option2dir['cst_file']='cst_file'

	 def make_target_flags(self, run, setup, filename, flags, subs ): # for the non-fixed options
			tr.out("make target flags...")
			args=self.get_args()

			if not args.partners_file:
				 flags.write("-docking:partners %s\n"%args.partners)
			else:
				 run.add_subst('CM_PARTNERS_FILE', setup.cm_path(args.partners_file))
				 flags.write("@$CM_PARTNERS_FILE\n")

#			flags.write("-docking:partners %s\n"%args.partners)

			if not args.pdb:
				 raise library.MissingInput("you need an input structure (pdb) for method %s"%self.name)
			run.add_subst('CM_INPUT_PDB', setup.cm_path(args.pdb))
			flags.write("-in:file:s $CM_INPUT_PDB\n")

			if not args.native:
				 raise library.MissingInput("you need an native input structure (pdb) for method %s"%self.name)
			run.add_subst('CM_NATIVE_PDB', setup.cm_path(args.native))
			flags.write("-in:file:native $CM_NATIVE_PDB\n")

			if not args.n_replica:
				 raise library.MissingInput("you need to specify n_replica for method %s"%self.name)
			flags.write("-n_replica %d\n"%int(args.n_replica))

			if args.protocol == "rep_cen":
				 if args.score == "interchain_cen":
						flags.write("-score:weights %s\n"%args.score)
				 else:
						raise library.MissInput("wrong score weights for low resolution docking")
				 if args.min_score_score:
						flags.write("-score:min_score_score %f\n"%float(args.min_score_score))
			else:
				 if not args.score or args.score=="interchain_cen":
						raise library.MissInput("you need to give score weights correctly for high res docking")
				 flags.write("-score:weights %s\n"%args.score)
				 if args.extra_score:
						flags.write("-score:docking_interface_score")


# 			if args.score:
# 				 flags.write("-score:weights %s\n"%args.score)
# 			if args.out:
# 				 flags.write("-out:file:silent %s\n"%args.out)


	 def setup_file_library( self ):
 			BasicMethod.setup_file_library( self )
			args = self.get_args()
			fl = self.file_library

			path = flag_lib+"/methods/_docking_base/"

			fl.executable  = "r_play_with_etables"
#		fl.provide_file( "flags", path, "flags_docking" )
#		fl.add_string( "commandline", "@flags_docking @$CM_FLAGFILE");
		#alternative to flags.write("-in:file:silent CM_DECOYS" ) up in make_target_flags
		#you could write here:
#		args=self.get_args()
#		if args.decoys:
#			fl.add("flags","flags_docking","-in:file:silent $CM_DECOYS")

			fl.provide_file( "flags", path, "flags_docking_replica" )
			fl.add_string( "commandline", "@flags_docking_replica @$CM_FLAGFILE ")
			if args.protocol == 'rep_cen':
				 fl.provide_file( "flags", path, "dock_cen.xml" )
				 fl.add_string( "commandline", " -parser:protocol @@dock_cen.xml ")
				 fl.provide_file( "flags", path, "hamiltonians_cen.txt")
			if args.protocol == 'rep_high':
				 fl.provide_file( "flags", path, "dock_high.xml" )
				 fl.add_string( "commandline", "-parser:protocol @@dock_high.xml ")
			if args.protocol == 'rep_high_rigid':
				 fl.provide_file( "flags", path, "dock_high_rigid.xml" )
				 fl.add_string( "commandline", "-parser:protocol @@dock_high_rigid.xml ")
			if args.protocol == 'rep_high_frozen':
				 fl.provide_file( "flags", path, "dock_high_frozen.xml" )
				 fl.add_string( "commandline", "-parser:protocol @@dock_high_frozen.xml ")


