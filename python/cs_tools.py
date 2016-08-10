#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
## make mammoth structure alignments

def get_flexible_residues( pred, threshold ):
	filelist=open(pred).readlines()
	flexible_residues=[]
	for r in filelist:
		tags=r.split()
		if len(tags)==10:
			if float(tags[7])<threshold:
				flexible_residues.append( int(tags[0]) )
	return flexible_residues
