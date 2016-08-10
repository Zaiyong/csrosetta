#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
import string
import amino_acids
import library
import sys
def cut_sequence(sequence,start,end,verbose):
	if ( end==0 ):
		end=len(sequence)
	if ( end > len(sequence) ):
		sys.stderr.write("WARNING: your end value is out of range\n");
	new_seq=sequence[(start-1):end];
	if verbose:
		print "\nold fasta:\n %s"%sequence
		print "\nnew fasta:\n %s%s%s"%('-'*(start-1),new_seq,'-'*(len(sequence)-end))
	return new_seq, end

def read_fasta(file):
	return read_fasta_stream(open(file,'r'))

#compare to sequences, consider '-' as matches anything, and don't fail if length mismatches
def compare_fasta(fasta1,fasta2, strict_length=False ):
	max_length=max(len(fasta1),len(fasta2))
	if not strict_length:
		_fasta1=fasta1+'-'*(max_length-len(fasta1))
		_fasta2=fasta2+'-'*(max_length-len(fasta2))
	else:
		if len(fasta1) != len(fasta2): return False
		_fasta1=fasta1
		_fasta2=fasta2
	for a,b in zip(_fasta1,_fasta2):
		if a=='-' or b=='-': continue
		if a!=b: return False
	return True

def fill_gaps( fasta1, fasta2 ):
#	print 'FASTA_FILL_GAPS (IN1): ',fasta1
#	print 'FASTA_FILL_GAPS (IN2): ',fasta2
	for a,i in enumerate(fasta1):
		if a=='-' and len(fasta2)>i:
			fasta1[i]=fasta2[i]
#	print 'FASTA_FILL_GAPS: ',fasta1
	return fasta1

def read_fasta_stream(fd):
	for line in fd:
		if line[0]==">":
			continue
		return line[:-1];

def write_fasta(file, fasta, tag=None):
	fd=open(file,'w')
	if not tag:
		tag='t000_'
	fd.write(">%s\n"%tag)
	fd.write("%s\n"%fasta)
	fd.close()

def molecular_weight(fasta):
	import amino_acids
	tot=0
	for aa in fasta:
		tot+=amino_acids.molecular_weight[aa]
	amide_bonds_correction=(len(fasta)-1)*(15.9994+1.00794)
	return tot-amide_bonds_correction

def find_fasta_offset( fasta1, fasta2, verbose=1 ):
	try:
		return _find_fasta_offset( fasta1, fasta2, verbose )
	except library.InconsistentInput:
		return -_find_fasta_offset( fasta2, fasta1, verbose )

def _find_fasta_offset( fasta1, fasta2, verbose=1 ):
	if verbose>0:
		print "trying to match the following sequences"
		print fasta1
		print fasta2

	#find first usable pieces of sequenc
	last_gap=0
	first_aa=0
	matchable_pieces={}
	fasta1=fasta1+'-'
	for res_upl in range( 0, len(fasta1) ):
		#if this is a gap --> then is the end of the next matchable_piece
		if fasta1[res_upl]=='-':
			if last_gap<first_aa and res_upl-first_aa>4:
				matchable_pieces[fasta1[first_aa:res_upl]]=first_aa
			last_gap=res_upl
		elif last_gap == res_upl-1:
			first_aa=res_upl

	#start with longer pieces first to find match for offset
	sorted_pieces=sorted(matchable_pieces,key=len,reverse=True)

	#go through matchable pieces and see if it can be found in fasta
	inconsistent=True
	sequence_offset=0
	for pi in sorted_pieces:
		for l in range( 5, len(pi)):
#			print ' try to find '+pi[0:l]+' in ',fasta2
			if pi[0:l] in fasta2:
				sequence_offset=matchable_pieces[pi]-fasta2.index( pi[0:l] )
				inconsistent=False
				break

	if len(sorted_pieces):
		if verbose>0: print pi, " is in ", fasta2
		if verbose>0: print "with offset ", sequence_offset
		#now go from back to fron along the original fasta2 sequence to look for a matching piece
		#if not inconsistent:
		#	print "found offset: %d\n"%sequence_offset
		if inconsistent:
			print "cannot determine offset\n"
			raise library.InconsistentInput("source fasta is inconsistent with the target fasta sequence")
	else:
		raise library.MissingInput("not sufficient information to determine offset")
	return sequence_offset

def upl2fasta( lines ):
  #work out offset first
	dict={}
	for line in lines:
		if line[0]=="#": continue
		cols=line.split()
		if len(cols)==0: continue
		dict[int(cols[0])]=cols[1][0:3]
		dict[int(cols[3])]=cols[4][0:3]

	#residue numbers as "keys"
	keys=dict.keys()
	keys.sort()
#print "highest residue number", keys[-1] #last residue
	max_upl_res=keys[-1]

	#create fasta sequence where AA is known -- otherwise
	upl_fasta=list("-"*keys[-1])
	for key in keys:
		upl_fasta[key-1]=amino_acids.longer_names[ dict[key] ]

	upl_fasta="".join(upl_fasta)
	return upl_fasta

def pdb2fasta( file ):
	lines=open(file,'r').readlines()
	fasta=[];
	for line in lines:
		ss=string.split( line );
		if ss[0]!='ATOM': continue
		if len(ss)>5 and ss[2]=="CA":
			fasta.append(amino_acids.longer_names[ss[3]]);
	return "".join(fasta)

def dssp2fasta( file ):
	lines=open(file,'r').readlines()
	fasta=[];
	for line in lines:
		ss=string.split( line );
		if ss[-1]=='.' or ss[0]=='#': continue
		#if len(ss)>5 and ss[2]=="CA":
		fasta.append(ss[3]);
	return "".join(fasta)
