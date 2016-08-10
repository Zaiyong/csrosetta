#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

from os import path
import shutil

import argparse
from basic.options import ExampleArgumentParser
import sys
import subprocess
import os
### toolbox library
import library as csrosettalib
import traceback
import textwrap
import fasta
import xray2talos

def process_command_line(argv):
	parser = ExampleArgumentParser(
		prog=path.basename(__file__),

		description="Pick fragments based on chemical shifts and sequence information. This application is a wrapper that integrates several steps required to obtain fragments. First, BLAST is used to obtain a sequence profile from multiple sequence alignment. This step also provides the names of homologous proteins that can be excluded from fragment picking in benchmark mode using the -nohom flag. Second, TALOS+ is executed to obtain secondary structure predictions based on the chemical shift data. Finally, the ROSETTA application fragment_picker is started to assemble the fragment libraries using the VALL. Obviously, this wrapper has a great many dependencies. These dependencies are configured in csrosetta3/frag_picker/setup_paths.pl. Run install.py if dependencies have changed or to update the BLAST sequence database.",
		examples=[('%(prog)s -cs 2jrm_trim.tab -nohom','pick fragments from non-homologous proteins (for benchmarking)'),
							('%(prog)s -cs exp.tab','pick fragments for chemical shifts in exp.tab'),
							('%(prog)s -cs 1yyv.pdb -xray', 'pick fragments using fake TALOS prediction drawn from 1yyv.pdb'),
							]
		)

	parser.add_argument("-cs", help="chemical shifts in TALOS format", metavar="cs.tab", default=None, required='true' )
	parser.add_argument("-sizes", help="which sizes of fragments shall be build", type=int, nargs='*', default=[3, 9] )
	parser.add_argument("-nfrags", help="how many frags per size-class and sequence position to collect",
											type=int, default=200 )
	parser.add_argument("-hom", help="pick fragments also from homologous proteins [default]", action='store_true',
											default=True, dest='hom' )
	parser.add_argument("-nohom", help="do not pick fragments from homologous proteins", action='store_false', dest='hom')
	parser.add_argument("-outlier", help="report chemical shift outliers that are x*limit (no effect on fragments)",
											default=1.5)
	parser.add_argument("-trim", help='trim the sequence according to TALOS+ output and pred2rigid',
											action='store_true', default=False)
	parser.add_argument("-fasta", help="a target sequence", default=None )
	parser.add_argument("-nocheck", dest='check', action='store_false', default=True,
											help="don't run TALOS to check for chemical shift offsets or trimmin")
	parser.add_argument("-xray", help="Use a pdb file to fake TALOS-predicted secondary structure and phi/psi angles.",
											default=False, action="store_true")
#	parser.add_argument("-q", dest="verbose", action="store_false",
#											help="Reduce output (unimplimented).")
	csrosettalib.add_standard_args( parser )

	args = parser.parse_args()
	configs = argparse.Namespace()


	#Environmental Setup: If the cs file isn't in the current directory, copy it there as a temporary file
	if not path.split(os.path.abspath(args.cs))[0] == os.getcwd():
		if(path.basename(args.cs) in os.listdir(os.getcwd())):
			sys.stderr.write(" ".join([path.basename(args.cs), "already exists in the working directory. Remove it and try again."]))
			sys.exit(0)

		print "Copying cs file into local directory as temporary file..."
		shutil.copy(args.cs, path.basename(args.cs))

	configs.target_name = path.basename(args.cs).split('.')[0]

	#Runtime Environment Check: Verify that csrosetta's location is known
	try:
		configs.cs_root=os.environ['csrosettaDir']
	except KeyError:
		sys.exit("\n".join(["ERROR: Cannot find environment variable csrosettaDir.",
												"You probably forgot to run 'source <yourpath>/csrosetta3/com/init'.",
												"It is recommended to put this line into your .bashrc or .cshrc file.",
												]))
	#Runtime Environment Check: Verify that Rosetta's location is known
	try:
		configs.rosettabin = os.environ['ROSETTA3_BIN']
		configs.rosettadb  = os.environ['ROSETTA3_DB']
	except KeyError:
		try:
			rpath = os.environ['ROSETTA3_PATH']
			configs.rosettabin = os.path.abspath(rpath)+"/rosetta_source/bin"
			configs.rosettadb  = os.path.abspath(rpath)+"/rosetta_database"
		except KeyError:
			if(args.xray):
				sys.exit("\n".join(["ERROR: Cannot find envirionment variables ROSETTA3_PATH or.",
														"ROSETTA3_DB and ROSETTA3_BIN.",
														"You probably forgot to set this value in your .bashrc file.",
														"(The -xray flag requires Rosetta.)"
														]))
			else:
				pass


	#Environmental Setup: clean out old frag files and configure fragment file names
	configs.frag9_file=configs.target_name+'.frags9.dat'
	configs.frag3_file=configs.target_name+'.frags3.dat'

	if path.exists( configs.frag9_file+'.gz' ) and path.exists( configs.frag3_file+'.gz' ):
		sys.exit(" ".join(['fragments exist already - remove', configs.frag3_file+".gz", 'and',
											 configs.frag9_file+".gz", 'before running.']))

	#Environmental Setup: if both final frags files don't exist, intermediate frag files cause the fragment_pl script to skip; remove them.
	else:
		for r in ['frags.fsc.score.200.9mers','frags.fsc.score.200.3mers',
							configs.frag3_file+'.gz', configs.frag9_file+'.gz']:
			try:
				os.remove(r)
			except OSError:
				continue

	#Basic Options Conversion: parse 'hom' option
	if args.hom:
		configs.hom_flag='-hom'
	else:
		configs.hom_flag='-nohom'

	#Basic Options Converstion: parse 'frag_sizes' option
	configs.size_flags=" ".join([ "-frag_sizes "+str(x) for x in args.sizes])

	#Basic Options Converstion: parse 'nfrags' option
	if args.nfrags:
		configs.n_flag="-n_frags %d"%args.nfrags
	else:
		configs.n_flag=""

	return configs, args

def check_offsets(target, cs_file, outlier):
	talos_dir='%(target)s.fasta.talos'%locals()
	print 'Run TALOS+...'
	pipe=subprocess.Popen('mkdir -p %(talos_dir)s; cd %(talos_dir)s; talos+ -in ../%(cs_file)s 2>/dev/null'%locals(), shell=True, stdout=subprocess.PIPE)
	outlier_line="Checking for Chemical Shift Outliers ..."
	offset_str=""
	for line in pipe.stdout:
		line=line.strip('\n')
		tags=line.split()
		if len(tags)>4 and tags[3]=='Secondary' and abs(float(tags[5])/float(tags[7]))>outlier:
			if outlier_line:
				print outlier_line
				outlier_line=None
			print line
		if 'Checking Chemical Shift Referencing' in line:
			atom=None
			offset=None
			for line in pipe.stdout:
				if 'ANN prediction' in line: break
				line=line.strip('\n')
				print line
				tags2=line.split()
				if 'Referencing' in line:
					atom=tags2[4]
					if len(tags2)>5:
						offset=float(tags2[5])
				elif atom and len(tags2)>=1:
					offset=float(tags2[0])
				if atom and offset:
					if atom=='CA/CB:':
						offset_str+='CA: %(offset)f CB: %(offset)f '%locals()
					elif atom=='HA:':
						offset_str+='HA: %(offset)f '%locals()
					elif atom=="C':":
						offset_str+='C: %(offset)f '%locals()
					elif atom=="N:" or atom=="HN:" or atom=="H:":
						offset_str+='%(atom)s: %(offset)f '%locals()
					atom=None
					offset=None

	if len(offset_str):
		original_shifts=cs_file
		cs_file='%s.corrected.tab'%(path.splitext(cs_file)[0])
		print 'Re-referencing chemical shifts: the re-referenced shifts are stored to %(cs_file)s...'%locals()
		subprocess.call('correctCSoffset %(original_shifts)s -offsets "%(offset_str)s" > %(cs_file)s'%locals(), shell=True)
		print 'Rerun TALOS+ with %(cs_file)s...'%locals()
		pipe=subprocess.Popen('mkdir -p %(talos_dir)s; cd %(talos_dir)s; talos+ -in ../%(cs_file)s 2>/dev/null'%locals(), shell=True, stdout=subprocess.PIPE)
		for line in pipe.stdout:
			if 'Offset' in line:
				print 'WARNING still offsets in chemical shifts', line
	return cs_file,talos_dir

def write_fasta_file(cs_file, target, args, configs):

	from cs import TalosCSFile
	tab=TalosCSFile()
	tab.read_file( cs_file )
	sequence=tab.sequence
	fasta_file=target+'.fasta'

	if args.fasta:
		print 'read fasta sequence from %s'%args.fasta
		target_fasta=fasta.read_fasta(args.fasta)
		print target_fasta
		if not fasta.compare_fasta( target_fasta, sequence, strict_length=True ):
			sys.exit("\n".join([
						'fasta-sequence does not match sequence in chemical shift file!!! ',
						'FASTA: '+target_fasta,
						'CS:    '+sequence]))

		sequence=fasta.fill_gaps(target_fasta,sequence)

	if path.exists( fasta_file ):
		target_fasta=fasta.read_fasta(fasta_file)
		if not '-' in sequence and target_fasta!=sequence:
			print "inconsistent fasta sequence: between chemical shifts %(cs_file)s and fasta file %(fasta_file)s"%locals()
			print "will overwrite fasta file, and create backup of original fasta file .bak"
			shutil.copy(fasta_file,fasta_file+".bak")
			fasta.write_fasta(fasta_file, sequence, configs.target_name )

	fasta_file=target+'.fasta'
	if '-' in sequence:
		exit('require full sequence information, provide fasta file with -fasta')

	fasta.write_fasta(fasta_file, sequence, configs.target_name )

	return fasta_file

def main(argv=None):

	configs, args = process_command_line(argv)

	#the perl script called from this wrapper will otherwise produce new files in a remote directory
	# this is clearly unexpected behavior hence we copy file to local dir.
	cs_file=path.basename(args.cs)

	if args.xray:

    (sequence, ss_assignment, phi_psi) = xray2talos.run_jd2_score(
			args.cs, os.getcwd()+"/",
			#TODO: modifiy this to look for any score_jd2 (e.g. linuxgccdebug, osxgccrelease).
			#This is probably better done in the process_command_line function
			configs.rosettabin+"/score_jd2.linuxgccrelease",
			configs.rosettadb)

    predSS_outfilename = configs.target_name+".predSS.tab"
    xray2talos.write_fake_talos_predSS_file(sequence, ss_assignment, cs_file, predSS_outfilename)
    pred_outfilename = configs.target_name+".pred.tab"
    xray2talos.write_fake_talos_pred_file(sequence, phi_psi, cs_file, pred_outfilename)

		fasta_file = configs.target_name+".fasta"
		fasta.write_fasta(fasta_file, sequence, configs.target_name)

		fasta_talos_dirname = configs.target_name+".fasta.talos"
		try:
			os.mkdir(fasta_talos_dirname)
		except OSError:
			#TODO: this strikes me as dangerous. If the dir already exists, it's fine. But if it can't be created, now's the
			# time to react.
			pass

		shutil.copy(configs.target_name+".predSS.tab", fasta_talos_dirname+"/predSS.tab")
		shutil.copy(configs.target_name+".pred.tab", fasta_talos_dirname+"/pred.tab")

	  #shutil.copy(args.cs, path.basename(args.cs))

		perl_cmd = ['perl', configs.cs_root+"/frag_picker/pick_fragments_without_TALOS.pl"]

	else: #not args.xray

		if args.check:
			cs_file,talos_dir=check_offsets(configs.target_name, cs_file, args.outlier)

			if args.trim:
				original_shifts=cs_file
				cs_file='%s.autotrim.tab'%(path.splitext(cs_file)[0])
				subprocess.call('cd %(talos_dir)s; pred2rigid pred.tab > pred.rigid; cd ..; renumber_talos %(original_shifts)s %(cs_file)s -rigid %(talos_dir)s/pred.rigid'%locals(),
												shell=True, cwd=os.getcwd())
				print 'Rerun TALOS+ with %(cs_file)s...'%locals()
				pipe=subprocess.Popen('mkdir -p %(talos_dir)s; cd %(talos_dir)s; talos+ -in ../%(cs_file)s 2>/dev/null'%locals(),
															shell=True, stdout=subprocess.PIPE)

		fasta_file = write_fasta_file(cs_file, configs.target_name, args, configs)
		perl_cmd = ['perl', configs.cs_root+"/frag_picker/pick_fragments.pl"]


	import rosetta_environment as re
	perl_cmd.extend(['-talos_fn', cs_file,
									 '-fasta', fasta_file, configs.hom_flag, configs.size_flags, configs.n_flag,
									 '-picker %s/fragment_picker.%s.%s%s%s'%(re.bin,re.build,re.platform,'gcc','release'),
									 '&& mv frags.score.200.9mers', configs.frag9_file,
									 '&& mv frags.score.200.3mers', configs.frag3_file,
									 '&& gzip', configs.target_name+'.frags?.dat'])

	print " ".join(perl_cmd)
	subprocess.call(" ".join(perl_cmd), shell=True, cwd=os.getcwd())

	#subprocess.call(('perl %(cs_root)s/frag_picker/pick_fragments.pl -talos_fn %(cs_file)s -fasta %(fasta_file)s %(hom_flag)s %(size_flags)s %(n_flag)s'+
	#								'&& mv frags.score.200.9mers %(frag9_file)s && mv frags.score.200.3mers %(frag3_file)s && gzip %(target)s.frags?%(hom_name)s.dat')%locals(), shell=True, cwd=os.getcwd())


	return True

if __name__ == "__main__":
	csrosettalib.hello(__file__)

	status = main(sys.argv)
	sys.exit(status)
