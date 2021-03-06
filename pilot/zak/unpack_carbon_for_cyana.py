#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

from sys import argv
import BmrbAtomNames
from cs import ProtCSFile
import cs
from assignment.noesy import Atom,Resonance,ResonanceList

assert( len(argv)>2)
infile=argv[1]
outfile=argv[2]

#print infile

tab=ProtCSFile()
tab.read_stream(  open(infile,'r') )
sequence=tab.sequence
res_in=ResonanceList.read_from_prot( tab )

def unpack_resonances(res_in):
	new_res=[]
	old_res=[]
	for resid,resonances in res_in.iter_residues():
		aa=res_in.sequence()[resid-1]
		for res in resonances:
			name=res.name()
			if res.atom().elem() not in 'HC':
				old_res.append(res)
			else:
				origin_atoms = BmrbAtomNames.anti_degenerate(aa,name)
				for atom in origin_atoms:
#					print resid,aa,name,origin_atoms
					newr = Resonance(-1,Atom(atom,resid),res.freq(),res.error())
					newr.ambiguity=res.ambiguity
					new_res.append(newr)
#					res.atom().set_name(origin_atoms[1])
	res_out=ResonanceList()
	res_out.set_sequence(res_in.sequence())
	for res in new_res:
		res_out.add_resonance(res)
	for res in old_res:
		res_out.add_resonance(res)
	return res_out
res_in = unpack_resonances(res_in)
res_out = unpack_resonances(res_in)

prot_data=res_out.generate_dict()
ambiguity=[]
for r in res_out.itervalues():
	if r.ambiguity:
		ambiguity.append( r.ambiguity )
if len(ambiguity):
	prot_data['AMBIGUITY']=ambiguity

prot_file = ProtCSFile().from_table( cs.NIH_table().from_dict( prot_data ) )
prot_file.write( open(outfile,'w'), True )
