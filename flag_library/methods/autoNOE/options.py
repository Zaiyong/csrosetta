##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'


from os import path
import argparse
import os
### toolbox library
from library import MethodException

import traceback
from silent_lib import ReadSilentData


# definition of options for the method RASREC
sub_method_code = flag_lib+"/methods/_rasrec_base/rasrec_options.py"
if path.exists( sub_method_code ):
    exec open(sub_method_code, 'r' )
else:
    print "CANNOT FIND METHOD CODE %s"%sub_method_code
    exit()

if 'group' in locals():
	group.add_argument('-peaks', nargs='*', metavar='<noe>.peaks', action='append', help='add peak files for automatic assignment')
	group.add_argument("-shifts", nargs='*', metavar='<noe>.prot', help="add 1 or N files with chemical shifts, N=number of peak-files" );
	group.add_argument('-silent', help='start assignment from the provided ensemble (silent file)')
	group.add_argument('-no_auto_intensities', help='add intensity values for protons, override whats in shift files', default=False, action='store_true' )
	group.add_argument('-local_distances', help='provide a file of local distances, e.g., predicted with FragsToAtomDist from fragment library', default=None )

if 'run_group' in locals():
	run_group.add_argument('-noesy_cst_strength',type=float, help='increase to make noesy restraints more powerful', default=10) #former default 4
	run_group.add_argument('-no_auto_clean_prot', dest='auto_clean_prot', default=True, action='store_false', help='clean up the prot file automatically')

	run_group.add_argument('-upweight_cst_relax',type=float, help='increase weight in relax from default value by this factor', default=10) #former default 1

	run_group.add_argument('-nonetwork',help='switch off network anchoring during autoNOE assignment', action='store_true', default=False );
	run_group.add_argument('-nonetwork_filt',default=False, action='store_true', help="don't use network anchoring for the filter-restraints" )

	run_group.add_argument('-calibration_target',type=float, default=3.8,help='first calibration target for structure independent calibration')
	run_group.add_argument('-calibration_convergence', type=float, default=0, help='set a max stddev for distances used for calibration' )
	run_group.add_argument('-calibration_max_dist', type=float, default=0, help='set a maximum distance for noe bounds' )
	run_group.add_argument('-calibration_cycles', default=3, type=int, help='multiple iterations with distance-violation elimination to improve calibration')

	run_group.add_argument('-no_aggressive', dest='later_aggressive', default=True, action='store_false', help='first local-distviol, later full-distviol' )

#benchmark options
	run_group.add_argument("-scramble", help="scramble configuration file", default=None );

#legacy options
	run_group.add_argument('-later_aggressive', default=True, action='store_true', help='first local-distviol, later full-distviol' )
	run_group.add_argument('-network_mode', default="clean",help='select network implementation (orig, clean[clean])')
	run_group.add_argument('-cb_mapping',default=False,action='store_true',help='switch to vintage CB-mapping')
	run_group.add_argument('-uniform_calib',default=False,action='store_true',help='switch to atom independent calibration')
	run_group.add_argument('-qual_params',default=None, help='point to flag file for the noesy:prob settings')
	run_group.add_argument('-split_classes',default=False, action='store_true',help='split classes of restraints into HI and MED')
	run_group.add_argument('-old_symmetry_treatment', default=False, action='store_true', help='treat symmetry contribution as in CANDID' )
	run_group.add_argument('-update_filter_noe', help='update filter cst when going to fullatom phase', action='store_true', default=False )
	run_group.add_argument('-continously_update_filter_noe', help='keep continsously updating the restraints used for deciding about the archived structures from newest assignments', action='store_true', default=False )
	run_group.add_argument('-dcut', help='distance cutoff for assignment with initial ensemble')
	run_group.add_argument('-combine', help='should we use combination of restraints 1==no ', choices=['auto','1','2'], default='auto');
	run_group.add_argument('-drop_randomly', default=0, type=float, help='take the constraints that could be violated by distance and sample them sometimes')
	run_group.add_argument('-padding', default=None, type=float, help='add padding to final phaseII+phaseIII csts' )
	run_group.add_argument('-ignore_local_distances', default=False, action='store_true', help='even if local-distance data is present it will be ignored' )


class AutoNOEMethod(RasrecBaseMethod):
	def __init__(self,name,path):
		RasrecBaseMethod.__init__(self,name,path)

		self.non_file_options.append('dcut')
		self.non_file_options.append('combine')
#		self.non_file_options.append('continuous')
		self.non_file_options.append('no_auto_intensities')

		self.option2dir['peaks']=self.nmr_sub_dir
		self.option2dir['shifts']=self.nmr_sub_dir
		self.option2dir['silent']='structural_knowledge'
		self.option2dir['local_distances']='structural_knowledge'
		self.double_file_options.append('shifts')

		self.phase=1
		self.silent_input_file=None


	def prepare_prot( self, setup, input_shifts ):
		new_name=path.splitext( basename( input_shifts ))[0]+"_autoprepare.prot"
		has_prot, prot = setup.has_file( new_name, subdir=self.nmr_sub_dir )
		if not has_prot:
			print "prepare prot files for autoNOE-Rosetta ..."
			prot = setup.create_file( new_name, subdir=self.nmr_sub_dir)
			tr.out("should generate files %s..."%prot)
			tr.out("will go %s --> %s"%(setup.abspath(input_shifts), setup.abspath(prot)))
			import tools
			processor=tools.ProcessProtForAutoNOE.from_options( stereo=False, header=False, v=0 )
			processor( open(setup.abspath(input_shifts),'r'), open(setup.abspath(prot),'w'), fasta=self.fasta )
		return prot

	def make_target_flags(self, run, setup, filename, flags, subs ):
		RasrecBaseMethod.make_target_flags( self, run, setup, filename, flags, subs)
		fl = self.file_library
		args=self.get_args();
		if not args.peaks:
			raise MethodException( self, "require peak-files for automatic assignment" );
		if not args.shifts:
			raise MethodException( self, "requires at least one resonance file, or as many resonance files as peak files");
		if len(args.shifts) > 1 and not len(args.shifts)==len(args.peaks):
			raise MethodException( self, "either provide a single resonance file that is used for all peak files, "+
											 "or provide the same number of file-names as peak-files such that each shift file corresponds to one peak-file"
											 )


		if len(args.shifts) == 1:
			prots=args.shifts[0]
			if args.auto_clean_prot:
				prots=self.prepare_prot( setup, prots )
			flags.write( "-noesy:in:resonances %s\n"%setup.cm_path( prots ) )

			for p in args.peaks:
				flags.write( "-noesy:in:peaks %s\n"%setup.cm_path( p ) )

			if args.scramble:
				shift_files={}
				for s in args.shifts:
					fnames=path.splitext(basename(s))
					fname=fnames[0].replace('.','_')+fnames[1]
					shift_files[s]=fname
					system("cat %s | sed s/HN/\ H/ > %s"%(setup.abspath(s),fname))
				pwd=os.getcwd()
				os.chdir(pwd)
				shifts_in = shift_files[args.shifts[0]]
				shifts_tmp_out=path.splitext( shifts_in )[0]+'.scramble'+'.tmp'+path.splitext( shifts_in )[1]
				shifts_out = path.splitext( shifts_in )[0]+'.scramble'+path.splitext( shifts_in )[1]
				shift_files[args.shifts[0]]=shifts_out
				config_scramble=basename(args.scramble)
				config_shift=basename(setup.abspath(args.shifts[0]))
				shutil.copy( args.scramble, config_scramble )
				fl.provide_file("others",path.dirname(args.scramble),config_scramble )
				fl.add("others","initialize_assignments_phaseI.sh","process_prot_for_autoNOE -nostereo %s %s"%(setup.abspath(args.shifts[0]),shifts_tmp_out))
				fl.add("others","initialize_assignments_phaseI.sh","run_scramble -in %s @%s -out %s"%(shifts_tmp_out,config_scramble,shifts_out))
				flags.write("-noesy:in:resonances %s\n"%shifts_out)

		else:
			if args.scramble:
				raise library.MethodException(self,'cannot scramble shifts if we have multiple shift files')

			for (i,p) in enumerate( args.peaks):
				prots=setup.cm_path( args.shifts[ i ] )
				if args.auto_clean_prot:
					prots=self.prepare_prot( setup, prots )
				flags.write( "-noesy:in:peak_resonance_pairs %s %s\n"%( setup.cm_path( p ), prots) )

#		fl.add("flags","flags_denovo","@flags_nmr_patches")
		self.add_restraint_scores()
		fl.add("flags","flags_denovo","@flags_noe_assign")
#		flags.write("@flags_nmr_patches\n")

		if args.local_distances:
			if not args.ignore_local_distances:
				flags.write( "-noesy:in:local_dist_table %s\n"%setup.cm_path( args.local_distances ) )

		if self.silent_input_file: silent_file = self.silent_input_file
		else: silent_file = args.silent

		if silent_file:
			if not args.dcut:
				raise MethodException( self, "if -silent is chosen a distance cutoff with -dcut has to be provided (usually between 1...7)" )
			subs["CM_AUTONOE_DECOYS"] = setup.cm_path( silent_file )
			subs["CM_AUTONOE_DCUT"] = args.dcut
			flags.write("@flags_phaseII\n")
			self.phase = 2
		else:
#			flags.write("@flags_phaseI\n")
			self.phase = 1

		if self.combine == '1' and silent_file:
			flags.write("-iterative:initial_beta_topology %s\n"%"init_phaseII/beta.top")


	def set_input_dir( self, dir, rundir ):
		RasrecBaseMethod.set_input_dir(  self, dir, rundir )
		in_file=self.input_dir()+"/run/fullatom_pool/decoys.out"
		if not path.exists( in_file ):
			raise library.MissingInput(" cannot find file run/fullatom_pool/decoys.out in directory "+self.input_dir() );

		out_file=path.dirname(in_file)+'low_score_30.out';
		tags=silent_lib.read_lowscore_tags( in_file, 'score', 30 )
		silent_lib.extract_decoys( in_file, out_file, tags )
		self.silent_input_file = out_file

	def setup_file_library( self ):
		RasrecBaseMethod.setup_file_library( self )
		args=self.get_args()
		fl = self.file_library

		path = flag_lib+"/methods/autoNOE/"
#		fl.provide_file( "patches", path, "nmr_autoNOE_pool_patch" )
		fl.provide_file( "others", path, "initialize_assignments_phaseI.sh" )
		fl.provide_file( "others", path, "final_assignment.sh" )
		fl.provide_file( "others", path, "initialize_assignments_phaseII.sh" )

		 #fl.executable  = "minirosetta"
		fl.provide_file( "flags", path, "flags_noe_assign" )
		fl.provide_file( "flags", path, "flags_phaseI" )
		fl.provide_file( "flags", path, "flags_initnoe_filter" )
		fl.provide_file( "flags", path, "flags_initnoe_sampling" )

#		 fl.provide_file( "flags", path, "flags_closeloops_relax" )

	#	 fl.add_string( "commandline", "@flags_denovo @$CM_FLAGFILE -increase_cycles 2.0");
# --- figure out whether restraint combination is still required
		if self.silent_input_file: silent_file = self.silent_input_file
		else: silent_file = args.silent
		if args.noesy_cst_strength:
			print 'replacing the restraint strenght with cmd-line controlled weight...'
			fl.add("flags", "flags_noe_assign", "-noesy_weights:cst_strength %10.3e"%args.noesy_cst_strength)
		if args.combine=="auto" and silent_file:
			try:
				print "Automatic assignment of -combine flag -- check whether input decoys have less than 50 in noesy_autoassign_cst"
				sfd=ReadSilentData(['noesy_autoassign_cst'])
				file=open( silent_file, 'r' )
				min_score = 1000000
				for l in file:
					res = sfd.read_line( l )
					if res and min_score < float( res ):
						min_score = float( res )
				if min_score < 50:
					self.combine= '1'
				else:
					self.combine= '2'
			except library.MissingInput as inp:
				inp.add("Input decoys have not been made with phaseI-rasrec protocol. Cannot figure out optimal -combine setting in 'auto' - mode")
				inp.add("choose the combine mode 1 or 2 with option -combine ")
				raise
		elif args.silent:
			self.combine = args.combine
		else:
			self.combine = '2'

		fl.override( "flags", "flags_noe_assign", "-constraints:combine %s"%self.combine )
		if args.qual_params:
			for line in open( args.qual_params, 'r' ).readlines():
				tags=line.split()
				if len(tags)>0 and len(tags[0])>1 and tags[0][0]=='-':
					fl.add("flags","flags_noe_assign", line)
		if args.split_classes or args.qual_params:
			fl.add("flags","flags_iterative","-iterative:split_autoNOE_restraints")
			fl.add("flags","flags_initnoe_sampling", "-noesy:out:split" )
			fl.add("flags","flags_initnoe_sampling", "-noesy:out:worst_prob_class 4" )

		if self.combine == '1':
			fl.override("flags", "flags_iterative", "-iterative:max_nstruct -1 -1 0 0 -1 -1 0 0" )

		fl.override( "flags", "flags_iterative", "-iterative:max_nstruct 0 0 0 0 0 0 0 0" )

		if args.continously_update_filter_noe:
			fl.add( "flags", "flags_iterative", "-iterative:update_noe_filter_cst" )
			fl.remove( "flags", "flags_iterative", "-iterative:evaluate_only_on_slaves" )
			fl.remove( "flags", "flags_noe_assign", "-iterative:never_update_noesy_filter_cst" )

		if args.update_filter_noe:
			fl.remove( "flags", "flags_noe_assign", "-iterative:never_update_noesy_filter_cst" )

		if args.upweight_cst_relax:
			setting = float(fl.get_line("patches", "nmr_relax_patch", "atom_pair_constraint" ).split()[2])
			fl.override("patches", "nmr_relax_patch", "atom_pair_constraint = %f"%(setting*args.upweight_cst_relax) )

		if args.uniform_calib:
			fl.remove( "flags", "flags_noe_assign", "-noesy:atom_dependent_calibration")

		if args.calibration_target:
			fl.add("flags","flags_noe_assign", "-noesy_weights:defaults:calibration_target %4.1f 0.15 0.15 0.1 0.1 0.1 0.1"%args.calibration_target)

		if args.cb_mapping:
			fl.remove("flags","flags_noe_assign", "-noesy:map_to_cen_atom")

		if args.nonetwork:
			fl.add("flags","flags_noe_assign","-noesy:no_network")

		if args.network_mode:
			fl.override("flags", "flags_noe_assign", "-noesy:network:mode %s"%args.network_mode)

		if args.nonetwork_filt:
			fl.add( "flags", "flags_initnoe_filter", "-noesy:no_network")

		if args.old_symmetry_treatment:
			fl.add("flags","flags_noe_assign", "-noesy_weights:min_symmetry_reinforcement 1.0")

		if args.calibration_convergence:
			fl.add("flags","flags_noe_assign", "-noesy:calibration:convergence %8.1f"%args.calibration_convergence );

		if args.calibration_max_dist:
			fl.add("flags","flags_noe_assign", "-noesy:calibration:max_noe_dist %8.1f"%args.calibration_max_dist );

		if args.no_auto_intensities:
			print 'no_auto_intensities is True'
			fl.remove("flags", "flags_noe_assign", "-noesy:ignore_resonancefile_intensities" );

#		if args.guntert_calib:
#			fl.provide_file( "flags", path, "flags_guntert" )
#			fl.add("flags", "flags_iterative", "-iterative:refine_auto_noe_flags @@flags_guntert" );
#			fl.add("flags", "flags_noe_assign", "-noesy_weights:calibration_target 4.0" );
#			fl.add("flags", "flags_noe_assign", "-noesy:calibration:use_median" );


		fl.provide_file( "flags", path, "flags_phaseII" )
		fl.provide_file( "flags", path, "flags_phaseIII" )

		#iterative calibration together with phaseIII (aggressive distviol) doesn't make sense.
		#the effect is that violations are gone and thus the target of 10% violations cannot be reached
		#a work-around is to set the targets down: -noesy_weights:calibration_target 0.05 -noesy_weights:dcalibrate 0.05
		# however, tested on casd_ar3436 this didn't help. -- should do some testing on nmr_xray_pairs...
		fl.add("flags", "flags_phaseII", "-noesy:calibration:cycles %d"%args.calibration_cycles )
		if args.calibration_cycles > 1:
			fl.add("flags", "flags_phaseII","-noesy:calibration:ignore_eliminated_peaks" )

		fl.add("flags", "flags_noe_assign", "-iterative:randomize_elimination_candidates %f"%args.drop_randomly );

		use_assign_phases = False
		if args.later_aggressive:
			use_assign_phases = True

		if use_assign_phases:
			fl.add("flags","flags_iterative", "-iterative::staged_auto_noe_flags NONE NONE NONE NONE @@flags_phaseII @@flags_phaseII @@flags_phaseIII @@flags_phaseIII" )

		if args.padding:
			fl.add("flags","flags_phaseII","-noesy:out:padding %5.2f"%args.padding)
			fl.add("flags","flags_phaseIII","-noesy:out:padding %5.2f"%args.padding)


	def motd( self, rundir ):
		 if self.phase == 1:
			 script="initialize_assignments_phaseI.sh"
		 else:
			 script="initialize_assignments_phaseII.sh"

		 print "\n\nAttention: before starting the structure calculation in %s/run or %s/test\n"%(rundir,rundir),
		 print "the initialization script %s has to be started"%script
		 print "cd %s/run;"%rundir
		 print ". %s"%script
		 RasrecBaseMethod.motd( self, rundir )


#echo "Starting from beta.top stage3 .... Make sure it is not a helical protein, eventually remove this line from flags_noe_assign"
#echo -iterative:initial_beta_topology init_phaseII/beta.top >> flags_noe_assign
#fi

method = AutoNOEMethod(method_name, method_path)
