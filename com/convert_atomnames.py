#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-

#application specific headers

#default headers
import argparse
from basic.options import ExampleArgumentParser
from os.path import basename
import traceback, sys
import amino_acids
import fasta
#toolbox headers
import library
import BmrbAtomNames

parser = ExampleArgumentParser(prog=basename(__file__),
															 description="Convert atomnames between BMRB/NMR Style and PDB/Xray/Rosetta Style.",
															 examples=[('%(prog)s decoy.pdb converted.pdb','convert atom names in decoy.pdb and store in convert.pdb'),
																				 ('%(prog)s in.prot converted.prot','convert atom names in in.prot and store in convert.prot')]
															 )

parser.add_argument("infile", help="an input file");
parser.add_argument("outfile",  help="an output file with converted atoms",nargs="?",default="stdout");
parser.add_argument("-format", choices=['pdb','prot'], default='pdb' );
parser.add_argument("-fasta", help="fasta file to fill in sequence information if required" );
parser.add_argument("-seq", help="sequence file to fill in sequence information if required" );
parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )
parser.add_argument("-verbose", default=1, type=int )


args = parser.parse_args()

#output:
verbose=args.verbose
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	if verbose: library.hello( __file__ )

sequence=None
if args.fasta:
	sequence=fasta.read_fasta( args.fasta )
if args.seq:
	sequence=library.read_aa3_sequence( args.seq)


class AtomConverter:
	def __init__( self, atom_col=None, aa_col=None, resid_col=None, sequence=None, verbose=1 ):
		self.atom = atom_col
		self.aa_col = aa_col
		self.resid_col = resid_col
		self.sequence = sequence
		self.reverse = None
		self.verbose = verbose

	def convert( self, line ):
		tags=line.split()

		if self.resid_col and sequence:
			aa=sequence[int(tags[self.resid_col])]
		elif self.aa_col:
			aa=tags[self.aa_col]
			if len(aa)==3:
				aa = amino_acids.long2short(aa)
		else:
			raise 'required aa_col or resid_col and sequence to convert...'

		name = tags[self.atom]
		try:
			new_name = BmrbAtomNames.translate( aa, name, not self.reverse is None )
		except KeyError:
			if self.reverse is None:
				if self.verbose>2: print 'cannot translate %s, try the other direction'%name
				try:
					new_name = BmrbAtomNames.translate( aa, name, reverse=True )
					self.reverse = True
				except KeyError:
					new_name=name
					pass
			else:
				#raise
				new_name=name
				pass



		if new_name != name:
			max_len=max(len(name),len(new_name))
			name=name.ljust(max_len)
			new_name=new_name.ljust(max_len)
			import re
			if max_len<4 and 'X' in re.sub('[0-9]','X',new_name) and re.sub('[0-9]','',new_name)[0]=='H':
				if re.sub('[0-9]','X',new_name)[0]=='X':
					new_name=new_name+' '
					name=' '+name
				else:
					new_name=' '+new_name
					name=name+' '
			if verbose>9: print 'replace >%s< with >%s< in line %s'%(name,new_name,line[:-1])
			line=line.replace(name,new_name)
		return line

def convert_generic( infile,
										 outfile,
										 verbose,
										 allowed_tags=None,
										 translated_tags=None,
										 converters=[]
):
	reverse = None
	for line in infile:
		tags=line.split()
		if not allowed_tags or tags[0] in allowed_tags:
			if not translated_tags or tags[0] in translated_tags:
				for converter in converters:
					line = converter.convert(line)
			outfile.write(line)
		elif verbose>3:
			print 'ignoring line: %s'%line[:-1]

def convert_pdbs( infile, outfile, verbose ):
	reverse = None
	for line in infile:
		tags=line.split()
		if tags[0] in ['REMARK','ATOM','HETATM','MODEL','TER','ENDMDL']:
			if tags[0] in ['ATOM']:
				aa = amino_acids.long2short(tags[3])
				name = tags[2]
				try:
					new_name = BmrbAtomNames.translate( aa, name, not reverse is None )
					if new_name == "H1" and name == "1H": new_name = "H"
				except KeyError:
					if reverse is None:
						if verbose>2: print 'cannot translate %s, try the other direction'%name
						new_name = BmrbAtomNames.translate( aa, name, reverse=True )
						if new_name == "H1" and name == "1H": new_name = " H"
						reverse = True
					else:
						raise
				if new_name != name:
					max_len=max(len(name),len(new_name))
					name=name.ljust(max_len)
					new_name=new_name.ljust(max_len)
					import re
					if max_len<4 and 'X' in re.sub('[0-9]','X',new_name) and re.sub('[0-9]','',new_name)[0]=='H':
						if re.sub('[0-9]','X',new_name)[0]=='X':
							new_name=new_name+' '
							name=' '+name
						else:
							new_name=' '+new_name
							name=name+' '
					if verbose>9: print 'replace >%s< with >%s< in line %s'%(name,new_name,line[:-1])
					line=line.replace(name,new_name)
			outfile.write(line)
		else:
			if verbose>3:
				print 'ignoring line: %s'%line[:-1]

try:
	infile=open(args.infile,'r')
	if args.format == "pdb":
		convert_pdbs( open(args.infile,'r'), outfile, verbose )
	if args.format == 'prot':
		convert_generic( infile,outfile,verbose,converters=[AtomConverter(atom_col=3,aa_col=5,verbose=verbose)] )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
