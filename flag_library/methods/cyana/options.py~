##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'


from os import path
import argparse

### toolbox library
from library import MethodException
import automatic_setup
import traceback
from silent_lib import ReadSilentData

import shutil
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

if 'run_group' in locals():
	run_group.add_argument("-scramble", help="scramble configuration file", default=None );
	run_group.add_argument("-ignore_rdcs", help="don't use rdcs even if they are present", action='store_true', default=False );
	run_group.add_argument("-no_find_tensor", help="don't use CYANA macro find-tensor", action='store_true', default=False );
	run_group.add_argument("-zero_tensor", help="start with dummy values for tensor", action='store_true', default=False );
	run_group.add_argument("-loose_aco", help="take the recommended bounds on the dihedral restraints", action='store_true', default=False );
	run_group.add_argument("-loose_tolerances", help="use 0.5/0.03/0.03 for tolerances", default=False, action='store_true' );
	pass
#run_group.add_argument('');

class CyanaMethod(RasrecBaseMethod):
	def __init__(self,name,path):
		RasrecBaseMethod.__init__(self,name,path)

		self.option2dir['peaks']=self.nmr_sub_dir
		self.option2dir['shifts']=self.nmr_sub_dir
		self.double_file_options.append('shifts')

	def make_target_flags(self, run, setup, filename, flags, subs ):
		#RasrecBaseMethod.make_target_flags( self, run, setup, filename, flags, subs )
		args=self.get_args();
		if not args.peaks:
			raise MethodException( self, "require peak-files for automatic assignment" );
		if not args.shifts:
			raise MethodException( self, "requires at least one resonance file, or as many resonance files as peak files");
		if len(args.shifts) > 1 and not len(args.shifts)==len(args.peaks):
			raise MethodException( self, "either provide a single resonance file that is used for all peak files, "+
											 "or provide the same number of file-names as peak-files such that each shift file corresponds to one peak-file")
		cyana_peaks=[]
		cyana_shifts=[]
		fl = self.file_library
		if not args.cs:
			raise MethodException( self, "require cs.tab file for TALOS+ restraints");

		if not args.native:
			raise MethodException( self, "require native to generate seq file -- could do some programming to get this from FASTA file");

		if not args.fasta:
			raise library.MissingInput("fasta-file required for method %s"%self.name)

		run.add_subst('CM_FASTA',setup.cm_path(args.fasta))

		self.fasta=fasta.read_fasta(setup.abspath(args.fasta))

		add_tensors_to_seq=0
		if args.rdc and len(args.rdc)>0:
			print 'Have %d sets of RDCs...'%len(args.rdc)
			if not args.ignore_rdcs:
				add_tensors_to_seq=len(args.rdc)
			else:
				print '....RDCs are ignored per user request'

		run.add_file_to_rundir( setup, args.cs )
		import os
		pwd=os.getcwd()
		os.chdir( run._fn_rundir+"/run/" )
		os.mkdir( "talos" )
		os.chdir("talos")
		print 'run talos+...'
		system("talos+ -in ../%s > talos.log"%basename(args.cs))
		if not args.loose_aco:
			system("~/bin/talos+2cyana.awk pred.tab > ../talos.aco")
		else:
			system("~/bin/talos+2dyana.awk pred.tab > ../talos.aco")

		system("grep CA %s | awk 'BEGIN {started=0} started==0 {offset=$6-1; started=1} {print $4, NR}' > ../target.seq"%setup.abspath(args.native))
		os.chdir('..')


#write RDCs to a cyana-readable file
		if add_tensors_to_seq<=0:
			#rather work by removal than by adding since we require certain position for these...
			fl.remove("others","CALC.cya","weight_rdc = 0.02               # weight for RDC restraints")
			fl.remove("others","CALC.cya","cut_rdc      = 0.2                # cutoff for RDC violation output")
		else:
			fl.add("others","init.cya","#read @@rdcdistances.cya")
			fl.add("others","init.cya","rdcdistances")
			fl.override("others","CALC.cya","constraints := talos.aco,cyana.rdc    # additional (non-NOE) constraints")
			print 'process RDCs...'
			#first add pseudo-residues for tensors to sequence files
			seqfile=open('target.seq','a')
			seqfile.write('PL    3050\n')
			for i in range(0,add_tensors_to_seq):
				resid=3050+10*i
				elements=[('LL5',1),('LL5',1),('LL5',1),('LL5',1),('LL5',1),('ORI',5)]
				for e,offset in elements:
					resid+=offset
					seqfile.write('%-10s %5d\n'%(e,resid))

			seqfile.close()

			#read sequence as dictionary for converting RDC file
			seqfile=open('target.seq','r')
			seq_data={}
			for line in seqfile:
				tags=line.split()
				if len(tags)==2:
					seq_data[int(tags[1])]=tags[0]
			seqfile.close()

			rdcout=open('cyana.rdc','w')
			orientation=0
			rdcout.write('# Orientation  Magnitude  Rhombicity  ORI residue number\n')

			#get tensor estimates:
			for rdc_file in args.rdc:
				orientation+=1
				rdcs=RDC_Data()
				rdcs.read_file(setup.abspath(rdc_file))
				Dzz,R,a,b=rdcs.estimate_Da_and_R_hist()
				ori_res=3050+10*orientation
				if args.zero_tensor:
					Dzz=1
					R=0
				rdcout.write('%(orientation)10d %(Dzz)10.5f %(R)10.5f %(ori_res)d\n'%locals() )

			#now convert rdc-files
			orientation=0
			rdcout.write('#  First atom      Second atom                    RDC      Error  Weight  Orientation\n')
			for rdcin_filename in args.rdc:
				rdcin=open(setup.abspath(rdcin_filename),'r')
				orientation+=1
				for rdc_line in rdcin:
					tags=rdc_line.split()
					if len(tags)==0: continue
					if tags[0]=='#': continue
					res1=int(tags[0])
					res2=int(tags[2])
					atom1=tags[1]
					atom2=tags[3]
					try:
						resname1=seq_data[res1]
						resname2=seq_data[res2]
						rdc=float(tags[4])
						error=2.0
						weight=1.0
						rdcout.write('%(res1)10d %(resname1)10s %(atom1)5s %(res2)10d %(resname2)10s %(atom2)5s %(rdc)8.3f %(error)8.3f %(weight)8.3f %(orientation)5d\n'%locals())
					except KeyError:
						pass
			rdcout.close()
#			initcya=open('init.cya','w')
#			initcya.write('''
#cyanalib
#read seq target.seq
#rdcdistances
#''')
#			initcya.close()
			if not args.no_find_tensor and not args.zero_tensor:
				os.chdir(pwd)
				run.substitute( self, '@@FindTensor.cya' )
				run.substitute( self, '@@init.cya')
				run.substitute( self, '@@rdcdistances.cya' )
				os.chdir( run._fn_rundir+"/run/" )
				print 'Calculate Tensor with macro FindTensor.cya ...'
			  #now run cyana FindTensor.cya
  #			macro=open('FindTensor.cya','w')
  #			macro.write('''

  #read rdc cyana.rdc
  #print "    Input alignment tensor:"
  #do i 1 orientations
  #  print "    Orientation $i: magnitude = $magnitude(i) Hz, rhombicity = $rhombicity(i)."
  #end do
  #
  #rdc fittensor method=simplex       # (can take several minutes)
  #''')

  #			macro.close()
				import subprocess
				pipe=subprocess.Popen('module load cyana; cyana FindTensor.cya', shell=True, stdout=subprocess.PIPE)
				blabla=True
				axial=[]
				rhombic=[]
				for line in pipe.stdout:
					print line[:-1]
					if blabla:
						if 'Fitted axial and rhombic components' in line:
							blabla=False
						else:
							continue
					tags=line.split()
					if 'Axial' in line:
						axial.append(float(tags[-1]))
					if 'Rhombic' in line:
						rhombic.append(float(tags[-1]))
				print axial
				print rhombic
				assert len(axial)==len(rhombic)
				assert orientation==len(axial)


			  #now write RDCs again with new tensor
				rdcout=open('cyana.rdc','w')
				orientation=0
				rdcout.write('# Orientation  Magnitude  Rhombicity  ORI residue number\n')

			  #get tensor estimates:
				for rdc_file in args.rdc:
					orientation+=1
					rdcs=RDC_Data()
					rdcs.read_file(setup.abspath(rdc_file))
					ori_res=3050+10*orientation
					rdcout.write('%10d %10.5f %10.5f %d\n'%(orientation,axial[orientation-1],rhombic[orientation-1],ori_res) )

			  #now convert rdc-files
				orientation=0
				rdcout.write('#  First atom      Second atom                    RDC      Error  Weight  Orientation\n')
				for rdcin_filename in args.rdc:
					rdcin=open(setup.abspath(rdcin_filename),'r')
					orientation+=1
					for rdc_line in rdcin:
						tags=rdc_line.split()
						if len(tags)==0: continue
						if tags[0]=='#': continue
						res1=int(tags[0])
						res2=int(tags[2])
						atom1=tags[1]
						atom2=tags[3]
						try:
							resname1=seq_data[res1]
							resname2=seq_data[res2]
							rdc=float(tags[4])
							error=2.0
							weight=1.0
							rdcout.write('%(res1)10d %(resname1)10s %(atom1)5s %(res2)10d %(resname2)10s %(atom2)5s %(rdc)8.3f %(error)8.3f %(weight)8.3f %(orientation)5d\n'%locals())
						except KeyError:
							pass
				rdcout.close()


		peak_files={}
		for p in args.peaks:
			fnames=path.splitext(basename(p))
			fname=fnames[0].replace('.','_')+fnames[1]
			peak_files[p]=fname
			origin_peak_lines=open(setup.abspath(p),'r').readlines()
			peak_dim=3
			for line in origin_peak_lines:
				tags=line.split()
				if line.find('Number of dimensions')>=0:
					peak_dim=int(tags[4])
					break
			if peak_dim==3:
				cmdstr="cat %s | awk '/^#/ {print; next}"%(setup.abspath(p)) + '{printf(" %7d %8.3f %8.3f %8.3f %d %s %8.4e 0 e 0   0    0    0 0\\n",$1,$2,$3,$4,$5,$6,$7)}'+"' > %s"%fname
			elif peak_dim==4:
				cmdstr="cat %s | awk '/^#/ {print; next}"%(setup.abspath(p)) + '{printf(" %7d %8.3f %8.3f %8.3f %8.3f %d %s %8.4e  0.000E+00 e 0     0     0     0     0\\n",$1,$2,$3,$4,$5,$6,$7,$8)}'+"' > %s"%fname
			print cmdstr
			system(cmdstr)
			peak_lines=open(fname,'r').readlines()
			fd=[open(fname,'w'),None]
			for line in peak_lines:
				line_alt=None
				tags=line.split()
				if tags[0]=="#TOLERANCE" and args.loose_tolerances:
					line_alt=tags[0]+' '
					tol_dict={0.3:0.5, 0.04:0.03, 0.03:0.03, 999.0:999.0, 0.4:0.5}
					for t in tags[1:]:
						line_alt+='%5.2f '%tol_dict[float(t)]
					line_alt=[line_alt+'\n',line_alt+'\n']
				if tags[0]=="#CYANAFORMAT":
					symstr=['[NC]' ,'[CN]','[cn]','[nc]']
					for s in symstr:
						i=line.find(s)
						if i>=0: break
					if len(tags[1])==3 and (tags[1].find('c')>=0 or tags[1].find('n')>=0):#if it's 3D, the label should be upper case.
					#if str.islower(tags[1][0]):
					   tags[1]=tags[1].swapcase()
					if i>=0:
						symstr=line[i:(i+4)]

						fd[1]=open(fname+'symC','w')
						peak_files[p+'symC']=fname+'symC'
						fd[1].write('# Number of dimensions 3\n')
						line_alt=[line.replace(symstr,symstr[1]),line.replace(symstr,symstr[2])]
						assert len(args.shifts) == 1, 'cannot del with different shift-files and symmetric peaks, currently'
					line=" ".join(tags)+"\n"
				if line_alt:
					fd[0].write(line_alt[0])
				else:
					fd[0].write(line)
				if fd[1]:
					if line_alt:
						fd[1].write(line_alt[1])
					else:
						fd[1].write(line)
		fd[0].close()
		if fd[1]: fd[1].close()
		shift_files={}
		for s in args.shifts:
			fnames=path.splitext(basename(s))
			fname=fnames[0].replace('.','_')+fnames[1]
			shift_files[s]=fname
			system("cat %s | sed s/HN/\ H/ > %s"%(setup.abspath(s),fname))
		os.chdir(pwd)
		if len(args.shifts) == 1:
			if args.scramble:
				shifts_in = shift_files[args.shifts[0]]
				shifts_tmp_out=path.splitext( shifts_in )[0]+'_scramble'+'.tmp'+path.splitext( shifts_in )[1]
				shifts_out = path.splitext( shifts_in )[0]+'_scramble'+path.splitext( shifts_in )[1]
				shift_files[args.shifts[0]]=shifts_out
				config_scramble=basename(args.scramble)
				shutil.copy( args.scramble, config_scramble )
				fl.provide_file("others",path.dirname(args.scramble),config_scramble )
				fl.add("others","prepare.sh","process_prot_for_autoNOE -nostereo %s %s"%(shifts_in,shifts_tmp_out))
				fl.add("others","prepare.sh","run_scramble -in %s @%s -out %s"%(shifts_tmp_out,config_scramble,shifts_out))

			for p in peak_files.keys():
				cyana_peaks.append(peak_files[p])
				cyana_shifts.append(shift_files[args.shifts[0]])
		else:
			if args.scramble:
				raise library.MethodException(self,'cannot scramble shifts if we have multiple shift files')

			for (i,p) in enumerate( peak_files.keys() ):
				cyana_peaks.append(peak_files[p])
				cyana_shifts.append(shift_files[args.shifts[i]])

# else:
# 			for s in args.scramble_shifts:
# 				fnames=path.splitext(basename(s))
# 				fname=fnames[0].replace('.','_')+fnames[1]
# 				shift_files[s]=fname
# 				system("cat %s | sed s/HN/\ H/ > %s"%(s,fname))
# 			os.chdir(pwd)
# 			if len(args.scramble_shifts) == 1:
# 				for p in peak_files.keys():
# #					flags.write( "-noesy:in:peaks %s\n"%setup.cm_path( p ) )
# 					cyana_peaks.append(peak_files[p])
# 					cyana_shifts.append(shift_files[args.scramble_shifts[0]])
# 			else:
# 				for (i,p) in enumerate( peak_files.keys() ):
# 					cyana_peaks.append(peak_files[p])
# 					cyana_shifts.append(shift_files[args.scramble_shifts[i]])

		subs["CM_CYANA_PEAKS"]=",".join(cyana_peaks)
		subs["CM_CYANA_PROT"]=",".join(cyana_shifts)
		#system("rm "+flag_lib+"/methods/cyana/flags_cyana_rescore")
		#tmp=open(flag_lib+"/methods/cyana/flags_cyana_rescore",'r')
		#flags_cyana_rescore=open(flag_lib+"/methods/cyana/flags_cyana_rescore",'w')
		if args.native:
			#flags_cyana_rescore.write("-in:file:native $CM_NATIVE\n");
			subs["CM_NATIVE"] = setup.cm_path( args.native );
			if args.native_restrict:
				for i,rigid in enumerate(args.native_restrict):
					if i==0:
						#flags_cyana_rescore.write("-evaluation:rmsd NATIVE _full $CM_NATIVE_RIGID\n")
						subs["CM_NATIVE_RIGID"] = setup.cm_path( rigid )
						continue
					#flags_cyana_rescore.write("-evaluation:rmsd NATIVE _%s %s\n"%(rigid.split('.')[0].split('/')[-1], setup.cm_path( rigid )) )
		#flags_cyana_rescore.close()
		#RasrecBaseMethod.make_target_flags(self, run, setup, filename, flags, subs)
		fl = self.file_library

	def setup_file_library( self ):
		#RasrecBaseMethod.setup_file_library( self )
		automatic_setup.BasicMethod.setup_file_library(self)
		args=self.get_args()
		fl = self.file_library

		path = flag_lib+"/methods/cyana/"
		fl.provide_file( "others", path, "CALC.cya" )
		fl.provide_file( "others", path, "init.cya" )
		fl.provide_file( "others", path, "rdcdistances.cya" )
		fl.provide_file( "others", path, "FindTensor.cya" )

		fl.provide_file( "others", path, "flags_rescore" )
		fl.provide_file( "others", path, "rescore.sh" )
		fl.provide_file( "others", path, "prepare.sh" )


		if args.native:
			fl.add("others", "flags_rescore", "-in:file:native $CM_NATIVE")
			if args.native_restrict:
				fl.add("others", "flags_rescore", "-evaluation:rmsd NATIVE _full $CM_NATIVE_RIGID")

		fl.executable  = "cyana"
		fl.add_string( "commandline", "@@CALC.cya")
	def motd( self, rundir ):
		print "cd %s/run; sbatch cyana.slurm.job; or cyana CALC.cya > cyana.log"%rundir
		RasrecBaseMethod.motd( self, rundir )


#echo "Starting from beta.top stage3 .... Make sure it is not a helical protein, eventually remove this line from flags_noe_assign"
#echo -iterative:initial_beta_topology init_phaseII/beta.top >> flags_noe_assign
#fi

method = CyanaMethod(method_name, method_path)
