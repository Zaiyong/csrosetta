#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from os.path import basename
import os
import string
import argparse
import sys
import library
import StringIO
from assignment import noesy
from math import *
from assignment.noesy import SpinSystem, get_strips,spinsystem_evaluation,ResonanceBase, ResonanceList, Atom
import fasta
import copy
from sets import Set
import BmrbAtomNames
from math import fabs
from chemical import AllAminoAcidLib
from assignment.noesy.util import deviant_res_atoms, miss_res_atoms
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-min_vc", default=0.2, type=float )
parser.add_argument("-min_sep", default=0, type=int )
parser.add_argument("-keep_eliminated", action='store_true', default=False );
parser.add_argument("-print",dest='prin',type=int,default=1 );
parser.add_argument("-ref",help='the file involves true faults',default=None );
parser.add_argument("-scramble",help='scramble type',  choices={'carbon','residue','methyl','coupled'}, default=None, required=True );
library.add_standard_args( parser )
args = parser.parse_args()

class CollectResonance(ResonanceBase):
	def __init__( self, id, atom=None ):
		ResonanceBase.__init__( self, id, atom, None )
		self._freqs=[]
		self._indirect_freqs=[]
	def __str__(self):
		return ResonanceBase.__str__(self)
	def len_freqs(self):
		return len(self._freqs)
	def len_indirect_freqs(self):
		return len(self._indirect_freqs)
	def _freq_str( self ):
		self._freqs.sort(key=lambda tup: tup[0])
		median=self._freqs[ int( floor( 0.5*len( self._freqs) )) ][0];
		std=0
		for f in self._freqs:
			std+=(median-f[0])*(median-f[0])
		std/=len(self._freqs)
#		figure_out_partners
		partners=list(Set([f[2] for f in self._freqs]))
		#return "%d"%len(partners)
		return "%8.3f %8.3f %5d [ %s ]"%(median,sqrt(std),len(partners),"; ".join([ "%s"%p for p in partners ]))
#		return "%8.3f %8.3f [ %s ]"%(median,sqrt(std)," ".join([ "%7.3f(%1d)"%(f[0],f[1]) for f in self._freqs]))
	def freq(self):
		return (self._freq_lb, self._freq_ub)
	def add_freq( self,freq,id, other_guy ):
		self._freqs.append( (freq,id,other_guy) )

	def add_indirect_freq( self,freq,id, other_guy ):
		self._indirect_freqs.append( (freq,id, other_guy) )
	def iter_freqs(self):
		for freq in self._freqs:
			yield freq

def detect_wired_residue(res_list,ref_list):
	sequence=ref_list.sequence()
	wired_residues=[]
	for resid in range(1,len(sequence)+1):
		inter_len=0
		try:
			H_len=res_list.by_atom(Atom('H',resid)).len_freqs()
			for freq in res_list.by_atom(Atom('H',resid)).iter_freqs():
				if freq[2].resid() !=resid:
					inter_len+=1
		except KeyError:
			H_len=0
		if sequence[resid-1] !='G':
			HA_len=res_list.by_atom(Atom('HA',resid)).len_freqs()
			for freq in res_list.by_atom(Atom('HA',resid)).iter_freqs():
				if freq[2].resid() !=resid:
					inter_len+=1
		else:
			HA_len=res_list.by_atom(Atom('HA2',resid)).len_freqs()+res_list.by_atom(Atom('HA2',resid)).len_freqs()
			for freq in res_list.by_atom(Atom('HA2',resid)).iter_freqs():
				if freq[2].resid() !=resid:
					inter_len+=1
			for freq in res_list.by_atom(Atom('HA3',resid)).iter_freqs():
				if freq[2].resid() !=resid:
					inter_len+=1
		if 1.0*inter_len/(HA_len+H_len)<0.35:
			wired_residues.append('%d %s'%(resid,sequence[resid-1]))
	return wired_residues
	#print 'thre following residue are wrong   ',wrong_res_residues

def detect_wired_methyl(res_list,ref_list):
	wired_methyls=[]
	sequence=ref_list.sequence()
	for r in ref_list.itervalues():
		if r.name()[0]=='Q':
			if sequence[r.resid()-1] in 'AITM':
				if res_list.by_atom(r.atom()).len_freqs()<6:
					wired_methyls.append(r.atom())
			elif sequence[r.resid()-1]=='L' and r.name()=='QD1':
				if res_list.by_atom(r.atom()).len_freqs() < min(6,res_list.by_atom(Atom('QD2',r.resid())).len_freqs()):
					wired_methyls.append(r.atom())
			elif sequence[r.resid()-1]=='L' and r.name()=='QD2':
				if res_list.by_atom(r.atom()).len_freqs() < min(6,res_list.by_atom(Atom('QD1',r.resid())).len_freqs()):
					wired_methyls.append(r.atom())
			elif sequence[r.resid()-1]=='V' and r.name()=='QG1':
				if res_list.by_atom(r.atom()).len_freqs() < min(6,res_list.by_atom(Atom('QG2',r.resid())).len_freqs()):
					wired_methyls.append(r.atom())
			elif sequence[r.resid()-1]=='V' and r.name()=='QG2':
				if res_list.by_atom(r.atom()).len_freqs() < min(6,res_list.by_atom(Atom('QG1',r.resid())).len_freqs()):
					wired_methyls.append(r.atom())
	return wired_methyls

def detect_wired_coupled(wired_res_atoms,ref_list,res_list):
	sequence=ref_list.sequence()
	wired_couples=[]
	for atom in wired_res_atoms:
		if atom.elem()=='H' and (atom.name() not in ['H','HA']):
			wired_couples.append(atom)
	for r in ref_list.itervalues():
		if r.atom().elem()=='H':
			aa=ref_list.sequence()[r.resid()-1]
			if aa=='-': continue
			stereo_name=all_aa_lib.aa_lib(aa).stereo(r.name())
			if stereo_name:
				try:
					self_assign_len=res_list.by_atom(r.atom()).len_freqs()
					stereo_assign_len=res_list.by_atom(Atom(stereo_name,r.resid())).len_freqs()
				except KeyError:
					self_assign_len=0
					try:
						stereo_assign_len=res_list.by_atom(Atom(stereo_name,r.resid())).len_freqs()
					except KeyError:
						stereo_assign_len=0
				if (max(stereo_assign_len,self_assign_len) < 5 ) and (stereo_assign_len+self_assign_len >0 ):
					if stereo_assign_len > self_assign_len:
						wired_couples.append(r.atom())
	return wired_couples

def detect_wired_carbon(wired_res_atoms,res_list):
	wired_carbons=[]
	for atom in wired_res_atoms:
		if (atom.elem() !='C') or (atom.name() =='C'): continue
		wired_carbons.append(atom)
# 	for r in res_list.itervalues():
# 		try:
# 			if r.atom().element() != 'H': continue
# 			if (r.len_freqs()>0) and (r.len_indirect_freqs()==0):
# 				aa=ref_list.sequence()[r.resid()-1]
# 				wired_carbons.append(Atom(BmrbAtomNames.get_source_atom(aa,r.name()),r.resid()))
# 		except KeyError:
# 			continue
	return wired_carbons

def map_assign_to_res_list(crosspeaks,ref_list,res_list):
	expids={}
	last_id=0
	for peak in crosspeaks:
		if not args.keep_eliminated and peak.eliminated(): #or peak.comment.find("elim")):
			continue
		freqs=peak.freq()
		atoms=Set()
		dP1=peak.spin_info(1).atom_col
		dP2=peak.spin_info(2).atom_col
		for assignment in peak:
			if assignment.volume_contribution()<args.min_vc: continue
			if args.min_sep and ( abs( assignment.atom( dP1 ).resid() - assignment.atom( dP2 ).resid() ) < args.min_sep ): continue
			for d in range(1,peak.dim()+1):
				other_d = None
				label_col=None
				if d==dP1:
					other_d=dP2
					label_col=peak.spin_info(1).label_col
				if d==dP2:
					other_d=dP1
					label_col=peak.spin_info(2).label_col
				atom=assignment.atom( d )
				if atom not in atoms:
					atoms.add( atom )
					try:
						resonance = res_list.by_atom( atom )
					except KeyError:
						resonance = res_list.add_resonance( CollectResonance( -1, atom ) )
					s,a=peak.col2spin(d)
					folder=peak.spin_info(s).folder(a)
					unfolded=freqs[d-1]
					if folder.is_folded():
						ref_res = ref_list.by_atom( atom )
						unfolded = folder.inverse( freqs[d-1], ref_res.freq() )
					if peak.info().experiment_id() in expids:
						id=expids[peak.info().experiment_id()]
					else:
						last_id+=1
						id=last_id
						expids[peak.info().experiment_id()]=last_id
					other_guy=None
					if other_d:
						other_guy = assignment.atom( other_d )
					resonance.add_freq( unfolded, id, other_guy )
					if label_col > -1:
						resonance.add_indirect_freq( unfolded, id, other_guy )
	return res_list

def plot_venn(wrongs,type):
	true_faults_NO=0
	pred_faults_NO=len(wrongs)
	true_pred_faults_NO=0
	ref_file_lines=open(args.ref,'r').readlines()
	for line in ref_file_lines:
		tags=line.split()
		if len(tags)>0:
			true_faults_NO+=1
			try:
				for atom in wrongs:
					if atom.name()==tags[0] and atom.resid()==int(tags[1]):
						true_pred_faults_NO+=1
						break
			except:
				for wrong in wrongs:
					residue=wrong.split()[0]
					aa=wrong.split()[1]
					if residue==tags[0] and aa==tags[1]:
						true_pred_faults_NO+=1
						break
	print 'actual wrong number is         %d'%true_faults_NO
	print 'predicted wrong number is      %d'%pred_faults_NO
	print 'predicted true wrong number is %d'%true_pred_faults_NO
	venn2(subsets = (true_faults_NO-true_pred_faults_NO,pred_faults_NO-true_pred_faults_NO, true_pred_faults_NO),set_labels = ('actural wrong number', 'predicted wrong number'))
	plt.title("Venn diagram of actrual and predicted wrong %s"%type)
	plt.show()

crosspeaks=noesy.read_peak_files(args.peaks)
raw_res_list=ResonanceList()
ref_list = noesy.ResonanceList.read_from_stream( open(args.prot,'r') )
for r in ref_list.itervalues():
	atom=r.atom()
	raw_res_list.add_resonance( CollectResonance( -1, atom ) )
out=[]
all_aa_lib=AllAminoAcidLib()
#sequence=ref_list.sequence()
wired_res_atoms=deviant_res_atoms(ref_list,0.01)
missed_atoms=miss_res_atoms(ref_list)
missed_atoms_by_residue={}
wired_res_atoms_by_residue={}

for atom in missed_atoms:
	try:
		missed_atoms_by_residue[atom.resid()].append(atom.name())
	except KeyError:
		missed_atoms_by_residue[atom.resid()]=[]
		missed_atoms_by_residue[atom.resid()].append(atom.name())
for i,atoms in missed_atoms_by_residue.iteritems():
	print 'the lost atoms in residues %d are: %s'%(i,"  ".join(atoms))

res_list=map_assign_to_res_list(crosspeaks,ref_list,raw_res_list)

if args.scramble=='coupled':
	wrong_couples=detect_wired_coupled(wired_res_atoms,ref_list,res_list)
	print 'the wrong couples are [ %s ]'%";".join([ "%s"%p for p in wrong_couples])
	if args.ref:
		plot_venn(wrong_couples,'coupled')
elif args.scramble=='residue':
	wrong_residues=detect_wired_residue(res_list,ref_list)
	print 'the wrong residues are [ %s ]'%";".join([ "%s"%p for p in wrong_residues])
	if args.ref:
		plot_venn(wrong_residues,'residues')
elif args.scramble=='methyl':
	wrong_methyls=detect_wired_methyl(res_list,ref_list)
	print 'the wrong methyls are [ %s ]'%";".join([ "%s"%p for p in wrong_methyls])
	if args.ref:
		plot_venn(wrong_methyls,'methyls')
elif args.scramble=='carbon':
	wrong_carbons=detect_wired_carbon(wired_res_atoms,res_list)
	print 'the wrong carbons are [ %s ]'%";".join([ "%s"%p for p in wrong_carbons])
	if args.ref:
		plot_venn(wrong_carbons,'carbons')


if args.prin==0:
	for r in ref_list.itervalues():
		if r.atom().element() != 'H': continue
		aa=ref_list.sequence()[r.resid()-1]
		try:
			match = res_list.by_atom( r.atom() )
			out.append("%20s %8.3f %5s %10d    %s"%(r.atom(),r.freq(), aa,len(match._freqs),BmrbAtomNames.get_source_atom(aa,r.name())))
		except KeyError:
			out.append('%20s %8.3f %5s %10d    %s'%(r.atom(),r.freq(), aa, 0,' '))

	for i,r in enumerate(out):
		if i<len(out)-1:
			tags1=r.split()
			tags2=out[i+1].split()
			if tags1[1]==tags2[1] and tags1[-1]==tags2[-1] and tags1[-1] != '.':
				print tags1[0],tags2[0],tags1[1],tags1[4],tags2[4],int(fabs(int(tags1[4])-int(tags2[4])))

if args.prin==1:
	for r in ref_list.itervalues():
		if r.atom().element() != 'H': continue
		aa=ref_list.sequence()[r.resid()-1]
		try:
			match = res_list.by_atom( r.atom() )
			if match.len_freqs()==0:
				print r,aa,'nomatch'
			else:
				print r, aa, match
		except KeyError:
			print r, aa, 'nomatch'
