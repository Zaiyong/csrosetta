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
	 group.add_argument("-partners", help="partner-string to specify which chains are docked", default="A_B" );
	 group.add_argument("-partners_file", help="partner file to specify which chains are docked");
	 group.add_argument("-pdb", help="start-pdb with both chains placed randomly" );
	 group.add_argument("-decoys", nargs='*',help="bunch of starting structures in silent-file format");
	 group.add_argument("-native", help="supply a native pdb for RMSD calculation" );
	 group.add_argument("-cst_file",help="constraint file")
	 group.add_argument("-cst_weight", help="constraint weight", default=5)

if 'run_group' in locals():
	 run_group.add_argument("-protocol", choices=['standard','centroid','refine'], help="choose the mode", default="standard");
	 run_group.add_argument("-out", help="file name for output");
	 run_group.add_argument("-score",help="which score weight to use?, options: interchain_cen, docking, score12...");
	 run_group.add_argument("-extra_score", choices=['docking_interface_score'], help="if to use docking_interface_score as extra_score");
	 run_group.add_argument("-nstruct",help="how many decoys should be produced", default=1 );
	 run_group.add_argument("-batches",help="how many batches to run")

tr = Tracer( "docking_base_method" )

# method-based code
class DockingBaseMethod(BasicMethod):
	 def __init__(self,name,path):
			BasicMethod.__init__(self,name,path)
			self.non_file_options.append('partners')
			self.non_file_options.append('run')
			self.non_file_options.append('out')
			self.non_file_options.append('score')
			self.non_file_options.append('extra_score')
			self.non_file_options.append('cst_weight')
			self.non_file_options.append('batches')
			self.non_file_options.append('nstruct')
			self.option2dir['pdb']='start'
			self.option2dir['decoys']='start'
			self.option2dir['native']='native'
			self.option2dir['partners_file']='partners_file'
			self.option2dir['cst_file']='cst_file'

	 def make_target_flags(self, run, setup, filename, flags, subs ):
			tr.out("make target flags...")
			args=self.get_args()
			if not args.pdb and not args.decoys:
				 raise library.MissingInput("you need an input structure (pdb or decoys) for method %s"%self.name)
			if not args.decoys: # docking low or high res
				 run.add_subst('CM_INPUT_PDB', setup.cm_path(args.pdb))
				 flags.write("-in:file:s $CM_INPUT_PDB\n")


			else: # if args.decoys: # refinement
				 if len(args.decoys)>1:
						batches=[]
						for i,file in enumerate(args.decoys):
							 flag_file=setup.create_file("flags_batch%04d"%(i+1))
							 open(setup.abspath(flag_file),'w').write("-in:file:silent %s\n-out:file:silent refine_%s\n"%(setup.abspath(file),basename(file)) )
							 batches.append(setup.abspath(flag_file))
						flags.write("-run:batches "+" ".join(batches)+"\n")  # after the loop finished, so same ident with for
						if args.pdb:
							 run.add_subst('CM_INPUT_PDB', setup.cm_path(args.pdb))
							 flags.write("-docking:recover_sidechains $CM_INPUT_PDB\n")
				 else:
						run.add_subst('CM_DECOYS',setup.cm_path(args.decoys))
						flags.write("-in:file:silent $CM_DECOYS\n")
						flags.write("-out:file:silent refine_%s\n"%basename(args.decoys));

			if args.batches:
				 batches=[]
				 for i in range( 0, int(args.batches) ):
						flag_file=setup.create_file("flags_batch%04d"%(i+1))
						open(setup.abspath(flag_file),'w').write("-out:file:silent decoys_%04d.out\n-out::suffix _%04d\n"%( (i+1), (i+1) ))
						batches.append(setup.abspath(flag_file))
				 flags.write("-run:batches "+" ".join(batches)+"\n")

			if args.out and not args.batches and not args.decoys:
				 flags.write("-out:file:silent %s\n"%args.out)

			if not args.native:
				 raise library.MissingInput("you need an native structure for rms calculation for method %s"%self.name)
			run.add_subst('CM_NATIVE_PDB', setup.cm_path(args.native))
			flags.write("-in:file:native $CM_NATIVE_PDB\n")

			if not args.partners_file:
				 flags.write("-docking:partners %s\n"%args.partners)
			else:
				 run.add_subst('CM_PARTNERS_FILE', setup.cm_path(args.partners_file))
				 flags.write("@$CM_PARTNERS_FILE\n")

			if args.score:
				 flags.write("-score:weights %s\n"%args.score)
			if args.out:
				 flags.write("-out:file:silent %s\n"%args.out)
			if args.cst_file:
				 run.add_subst('CM_CST_FILE', setup.cm_path(args.cst_file))
				 flags.write("-cst_file $CM_CST_FILE\n")
				 flags.write("-cst_weight %s\n"%float(args.cst_weight))
			flags.write("-nstruct %d\n"%int(args.nstruct) )


	 def setup_file_library( self ):
			BasicMethod.setup_file_library( self )
			args = self.get_args()
			fl = self.file_library

			path = flag_lib+"/methods/_docking_base/"

			fl.executable  = "docking_protocol"
#		fl.provide_file( "flags", path, "flags_docking" )
#		fl.add_string( "commandline", "@flags_docking @$CM_FLAGFILE");
		#alternative to flags.write("-in:file:silent CM_DECOYS" ) up in make_target_flags
		#you could write here:
#		args=self.get_args()
#		if args.decoys:
#			fl.add("flags","flags_docking","-in:file:silent $CM_DECOYS")
			if args.protocol == 'refine':
				 fl.provide_file( "flags", path, "flags_docking_refine" )
				 fl.add_string("commandline", "@flags_docking_refine @$CM_FLAGFILE");
			if args.protocol == 'standard':
				 fl.provide_file( "flags", path, "flags_docking_std" )
				 fl.add_string( "commandline", "@flags_docking_std @$CM_FLAGFILE");
			if args.protocol == 'centroid':
				 fl.provide_file( "flags", path, "flags_docking_cen" )
				 fl.add_string( "commandline", "@flags_docking_cen @$CM_FLAGFILE");


