#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
import math
from PDB.Polypeptide import one_to_three
from assignment.noesy import SpinSystem, Resonance


def compare_predict_to_true(ss_proton_name_list,true_protons,ss,aa):
	pd_ss=1
	for r in ss_proton_name_list:
		pd_st=0
		if r=='H':
			h_atom_name='N'
		elif r.find('A')>0:
			h_atom_name='CA'
		elif r.find('B')>0:
			h_atom_name='CB'
		elif r.find('G')>0:
			if aa in 'RQEKMPL':
				h_atom_name='CG'
			else:
				h_atom_name='C'+r[1:3]
		for g in true_protons[h_atom_name]:
			try:
				error=math.fabs(g.freq()-ss.strip(r).proton().freq())
			except:
				true_protons[r]=Resonance(freq=ss.strip(r).proton().freq())
				error=0.0
			if error<=0.05:
				pd_st=1
				break
		pd_st=pd_st*pd_ss
	if pd_st==1:
		determine='correct'
	else:
		determine='incorrect'
	return determine
def outprint(residue_validation,ss_proton_name_list,target,resid,seq,systems):
	aa=seq[resid-1]
	if residue_validation.amount_ss()<1:
		print "%(target)20s %(AA)5s none-spinsystem: residue: residue %(resid)3d has no spinsystems \n" \
				%{'target':target,'AA':one_to_three(seq[resid-1]),'resid':residue_validation.resid()}
	if residue_validation.amount_ss()==1:
 		iter_ss=residue_validation.iter_spinsystems()
 		spinsystems=[]
 		for r in iter_ss:
 			spinsystems.append(r)
 		true_protons=residue_validation.true_protons()
		print true_protons
		ss_determine=compare_predict_to_true(ss_proton_name_list,true_protons,spinsystem[0],aa)
# 		for r in ss_proton_name_list:
# 			pd=1
# 			try:
# 				error=math.fabs(true_protons[r].freq()-spinsystems[0].strip(r).proton().freq())
# 			except:
# 				true_protons[r]=Resonance(freq=spinsystems[0].strip(r).proton().freq())
# 				error=0.0
# 			if error<=0.05:
# 				determine='correct'
# 				pd=pd*1
# 			else:
# 				determine='incorrect'
# 				pd=pd*0
# 			print "%(target)20s %(AA)5s unique-spinsystem: strip: compare to native: true CS of %(name)3s is %(true)8.3f and valid version is %(valid)8.3f , the determine is %(dtm)10s\n" \
# 					%{'target':target,'AA':one_to_three(seq[resid-1]),'name':r,'true':true_protons[r].freq(),'valid':spinsystems[0].strip(r).proton().freq(),'dtm':determine}
# 		if pd==1:
# 			ss_determine='correct'
# 		else:
# 			ss_determine='incorrect'
		print "%(target)20s %(AA)5s unique-spinsystem: spinsystem: residue %(resid)3d , the validation are %(ss_dtm)10s\n" \
				%{'target':target,'AA':one_to_three(seq[resid-1]),'resid':residue_validation.resid(),'ss_dtm':ss_determine}

	elif residue_validation.amount_ss()>1:
		iter_ss=residue_validation.iter_spinsystems()
		spinsystems=[]
		for r in iter_ss:
			spinsystems.append(r)
		true_protons=residue_validation.true_protons()
		print true_protons
		resid_determine=[]
		for s in spinsystems:
			ss_determine=compare_predict_to_true(ss_proton_name_list,true_protons,s,aa)
# 		for s in spinsystems:
# 			pd=1
# 			for r in ss_proton_name_list:
# 				try:
# 					error=math.fabs(true_protons[r].freq()-s.strip(r).proton().freq())
# 				except:
# 					true_protons[r]=Resonance(freq=s.strip(r).proton().freq())
# 					error=0.0
# 				if error<=0.05:
# 					determine='correct'
# 					pd=pd*1
# 				else:
# 					determine='incorrect'
# 					pd=pd*0
# 				print "%(target)20s %(AA)5s multi-spinsystem: strip: compare to native: true CS of %(name)3s is %(true)8.3f and valid version is %(valid)8.3f , the determine is %(dtm)10s\n" \
# 						%{'target':target,'AA':one_to_three(seq[resid-1]),'name':r,'true':true_protons[r].freq(),'valid':s.strip(r).proton().freq(),'dtm':determine}
# 			if pd==1:
# 				ss_determine='correct'
# 			else:
# 				ss_determine='incorrect'
			print "%(target)20s %(AA)5s multi-spinsystem: spinsystem: residue %(resid)3d , the validation are %(ss_dtm)10s\n" \
					%{'target':target,'AA':one_to_three(seq[resid-1]),'resid':residue_validation.resid(),'ss_dtm':ss_determine}
			resid_determine.append(ss_determine)
		print "%(target)20s %(AA)5s multi-spinsystem: residue: residue %(resid)3d has %(total_ss)5d spinsystems , %(correct_ss)3d spinsytems are correct and %(incorrect_ss)3d spinsystems are incorrect\n" \
				%{'target':target,'AA':one_to_three(seq[resid-1]),'resid':residue_validation.resid(),'total_ss':len(resid_determine),'correct_ss':resid_determine.count('correct'),'incorrect_ss':resid_determine.count('incorrect')}

# 	systems=sorted(systems, key=SpinSystem, reverse=True )
# 	for sys in systems:
# 		print 'spinsystem list: %5.3f %s'%(sys.score(1), sys)
