#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

import argparse
import traceback
import string
import os
from os.path import basename
from os import path
import subprocess
import library
import shutil
import sys

parser = argparse.ArgumentParser(prog=basename(__file__), description="file silent file to extract decoys with lowest/highest score values")

#parser.add_argument("-select", choices=['weights','rmsd'], default='weights');
#parser.add_argument("-weights", nargs=2, type=float, default=[5,12]);
#parser.add_argument("-cst", choices=['log','sqrt'], default='log' )
parser.add_argument("-batch", action='store_true', default=False );
parser.add_argument("-force", help='force re-evaluation', action='store_true', default=False );
#parser.add_argument("-get_dir", default=False, action='store_true' );
parser.add_argument("-ensemble", help='write ensemble of 10 lowest energy structures from optimal run to this pdb file', default=None )
import glob


args=parser.parse_args()


class Run:
	def __init__(self,root,sc,prec,cst,target_fct,converged_fraction,rmsd=None):
		self.score=sc
		self.prec=prec
		self.rmsd=rmsd
		self.id=root
		self.cst=cst
		self.converged_fraction=converged_fraction
		self.target_fct=target_fct

	def __str__(self):
		s='%(id)40s cst=%(cst)5.1f score=%(score)8.1f prec=%(prec)5.1f fraction_converged=%(converged_fraction)4.2f target_fct=%(target_fct)5.1f'%self.__dict__
		if self.rmsd:
			s+=' rmsd=%(rmsd)5.1f'%self.__dict__
#		if self._key:
#			s+=' key=%(_key)5.1f'%self.__dict__
		return s


def read_noe_strength(file):
	for line in open(file,'r'):
		if '-noesy_weights:cst_strength' in line:
			tags=line.split()
			cst=float(tags[1])
			return cst
	raise KeyError('-noesy_weights:cst_strength not found in file %s'%file)
	return -1

def read_median_score( file, col ):
	import math
	from silent_lib import ReadSilentData
	data = ReadSilentData( [col], throw_exception=True ).read_column( open(file,'r' ))
	data=sorted(data)

	if len(data)%2:
		mid=len(data)/2
		median=(data[mid-1]+data[mid])/2
	else:
		median=(data[int(math.floor(len(data)/2))])

	return median


def read_lowest_score( file, col ):
	import math
	from silent_lib import ReadSilentData
	data = ReadSilentData( [col], throw_exception=True ).read_column( open(file,'r' ))
	return min(data)


def count_decoys( file ):
	import math
	from silent_lib import ReadSilentData
	data = ReadSilentData( ['score'], throw_exception=True ).read_column( open(file,'r' ))
	return len(data)

def determine_nr_flexible_talos( cache={} ):
	# get chemical shift file
	for line in open('../../flags_cs_rescore','r'):
		if '-evaluation:chemical_shifts' in line:
			talos_cs_file=line.split()[1]
			break
	if not path.isabs(talos_cs_file):
		talos_cs_file='../../../'+talos_cs_file

#if file not found try to repair under the assumption that we just have a different location of the cs_targetlib
	targetlib=os.environ['CS3_BENCH_TARGETLIB']
	targetlib_root=path.dirname(targetlib)
	targetlib_base=path.basename(targetlib)
	if not path.exists(talos_cs_file) and targetlib_base in talos_cs_file:
		start_good=string.find(talos_cs_file, targetlib_base )
		new_talos_cs_file=targetlib_root+'/'+talos_cs_file[start_good:]
		if path.exists(new_talos_cs_file):
			talos_cs_file=new_talos_cs_file

	if talos_cs_file in cache:
		return cache[talos_cs_file]

	#we don't have the data yet. run TALOS+
	if not path.exists('talos_run'):
		os.mkdir('talos_run')

	if not path.exists('talos_flexible.txt'):
		os.chdir('talos_run')

		subprocess.check_call('talos+ -in %(talos_cs_file)s'%locals(), shell=True)
		subprocess.check_call('pred2rigid pred.tab -flexible -smooth 1 > ../talos_flexible.txt', shell=True)
		os.chdir('..')

	nr_flexible=0
	for line in open('talos_flexible.txt','r'):
		tags=line.split()
		if tags[0]=='LOOP':
			nr_flexible+=int(tags[2])-int(tags[1])+1

	cache[talos_cs_file]=nr_flexible
	return nr_flexible

####
##  main program
####
start_dir=os.getcwd()
runs=[]
sequence=None

## parse local directory tree to discover all runs
## if 'fullatom_pool_stage8' is present then we consider this a finished run
for root,dirs,files in os.walk('.',followlinks=False):
	if 'fullatom_pool_stage8' in dirs:
		cst=read_noe_strength(root+'/flags_noe_assign')
		try:
			if path.exists(root+'/fullatom_pool/_auto_analysis') and args.force:
				shutil.rmtree(root+'/fullatom_pool/_auto_analysis')

			if not path.exists(root+'/fullatom_pool/_auto_analysis'):
				os.mkdir(root+'/fullatom_pool/_auto_analysis')
			if cst > 50 or cst < 1:
				print 'WARNING: cst-weights below 1 or above 50 are ignored as they did not improve results in benchmark'
				continue

			os.chdir(root+'/fullatom_pool/_auto_analysis')
			if path.exists('low_10.out') and not count_decoys('low_10.out')>=10:
				os.remove('low_10.out')

			if not path.exists('low_10.out'):
				if not args.batch: print 'extract 10 low-energy decoys for %s  - %5.1f'%(root,cst)
				subprocess.check_call('extract_decoys ../decoys.out -score 10 > low_10.out',shell=True)
			if not path.exists('ensemble_precision.txt'):
				if not args.batch: print 'execute ensemble_anslysis to obtain precision... '
				subprocess.check_call('ensemble_analysis.default.linuxgccrelease -residues:patch_selectors replonly -in:file:silent low_10.out -wRMSD 2 -calc:rmsd -mute all -unmute main > ensemble_precision.txt',shell=True)

			#read sequence from low_10.out to check that all runs are consistent...
			sequence_tags=open('low_10.out','r').readline().split()

			if 'SEQUENCE:' != sequence_tags[0]:
				raise library.InconsistentInput('Expected SEQUENCE at beginning of silent-file low_10.out in directory %s/fullatom_pool/_auto_analysis'%root)
			if not sequence:
				sequence=sequence_tags[1]
			elif sequence!=sequence_tags[1]:
				raise library.InconsistentInput('Sequences mismatch!!! call this program in a subtree that only contains runs of the same protein !')

			median_score=read_median_score('low_10.out','score')
			median_rmsd=None
			try:
				median_rmsd=read_median_score( 'low_10.out', 'rms_full' )
			except library.MissingInput:
				try:
					median_rmsd=read_median_score( 'low_10.out', 'rms' )
				except library.MissingInput:
					pass

			precs=[]
			for line in open('ensemble_precision.txt','r'):
				if not 'mean' in line: continue
				if not 'RMSD' in line: continue
				if 'number of atoms from' in line:
					tags=line.split()
					total_nr_residues=int(tags[5])
					converged_nr_residues=int(tags[-1])

				precs.append(float(line.split()[-1]))

			#determine target function equivalent for second success criterion
			mfs=read_lowest_score('../../centroid_pool_stage4/decoys.out','noesy_autoassign_cst')
			cst_file=open('../../initial_assignment/noe_auto_assign.cst','r')
			strength=None
			for line in cst_file:
				tags=line.split()
				if tags[0]=='AmbiguousNMRDistance':
					strength=tags[8]
				elif tags[0]=='AmbiguousNMRConstraint':
					strength=tags[4]
				elif tags[0]=='AtomPair':
					strength=tags[8]
				if strength: break
			target_fct=mfs*float(strength)/cst

			#determine fraction of converged residues normalized by residues predicted to be flexible
			nr_flexible_talos=determine_nr_flexible_talos()
			converged_fraction=min(1.0*converged_nr_residues/(total_nr_residues-nr_flexible_talos),1.0)

			runs.append(Run(root.strip('.').strip('/'),median_score,precs[2],cst,target_fct,converged_fraction,median_rmsd))
			os.chdir(start_dir+'/'+root)
			if ( not path.exists('final_assignment') or not path.exists('final_assignment/NOE_final.out') ) and path.exists('final_assignment.sh'):
				print 'Generate final NOE assignment...'
				subprocess.check_call('bash final_assignment.sh > final_assignment.log', shell=True )
				print 'done!'

			if path.exists('final_assignment') and path.exists('final_assignment/NOE_final.out'):
				print 'found final_assignment for %(root)s'%locals()
				pass
		except (subprocess.CalledProcessError, KeyboardInterrupt) as exc:
			full_trace=sys.exc_info()
			print 'execution failed...cleaning up'
			print exc
			os.chdir(start_dir)
			os.chdir(root+'/fullatom_pool')
			try:
				shutil.rmtree('_auto_analysis')
			except OSError: #for some reason after ^-C the directory removal fails, but directory is empty at that point after program finishes.. might be that subprocesses are not fully terminated yet.
				pass #just ignore this now.
			raise full_trace[0],full_trace[1],full_trace[2]
		os.chdir(start_dir)


runs=sorted(runs,key = lambda run: run.cst)

by_cst={}
for r in runs:
	by_cst.setdefault(r.cst,[]).append(r)

number_of_runs=len(runs)
number_of_levels=len(by_cst)

if len(runs)<1:
	print 'no finished runs found'
	exit()

import math
wcst=12
wprec=5

for run in runs:
	run._key=run.score+wprec*run.prec-wcst*math.log(run.cst)

weighted_runs=sorted(runs,key = lambda run: run._key )
if not args.batch:
		print ('------------------------------------- all runs -----------------------------------------')
		for r in weighted_runs:
			print r
		print ('--------------------------- selected by weighted formula -------------------------------' )
		print weighted_runs[0]
		print
else:
#		if args.select=='weights':
		if args.get_dir:
			print weighted_runs[0].rundir+'/'+weighted_runs[0].id
		else:
			print weighted_runs[0]


if args.ensemble:
	final_pdb=args.ensemble
	print 'extract 10 lowest energy models from optimal run to %s...'%final_pdb
	low10_file=weighted_runs[0].id+'/fullatom_pool/_auto_analysis/low_10.out'
	subprocess.check_call('pack_pdbs -silent %(low10_file)s > %(final_pdb)s'%locals(), shell=True )
