## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

import subprocess
from os import dup2,path
from os.path import exists
from operator import add
from math import sqrt
from os.path import *
import argparse
import sys
import copy
import shutil
import amino_acids
### toolbox library
#import library
from os import mkdir , makedirs

from warnings import *
import warnings

import traceback

from PDB import PDBParser
from PDB import PDBIO

from library import square

from library import mkdirp
import library

class NIH_table:
	def __init__(self):
		self.vars=None

		#format string
		self.format=""

		#line entries, indexed by (resi, atom)
		self.table={}

		#additional data entries:
		#for instance DATA SEQUENCE VAL would be data['SEQUENCE']='VAL'
		self.data={}

		self.has_atom=False

		self.header_prefix=''

	def read_file( self, file):
		fd=open(file,'r')
		self.read( fd )
		return self

	def read( self, fd ):
		self.__init__()
		self.read_header( fd )
		if not self.vars or not self.format:
			raise MissingInput("cannot find VARS and FORMAT entries in file %s"%file )
		self.read_body( fd )
		return self

	def iteritems(self):
		return self.table.iteritems()

	def get_header_prefix( self, line ):
		header_types=['DATA','VARS','REMARK','FORMAT']
		for t in header_types:
			if t in line:
				start=line.index(t)
				self.header_prefix=line[0:start].replace(' ','')
#				print 'found header prefix: ->%s<-'%self.header_prefix
				return

	def read_header( self, fd ):
		if self.header_prefix=='':
			self.header_prefix=None
		for line in fd:
#			print 'header line: ',line[:-1]
			if not self.header_prefix:
				self.get_header_prefix(line)
			if self.header_prefix:
				line=line.replace(self.header_prefix,'')
			tags=line.split()
			if not tags: continue
			if "DATA"==tags[0]:
				key=tags[1];
				tags[1]=""
				if key in self.data.keys():
					tags[1]=self.data[key]
				self.data[key]="".join(tags[1:])
			if "VARS"==tags[0]:
				self.vars=tags[1:]
			if "FORMAT"==tags[0]:
				self.format=line.strip()[6:]
			if self.vars and self.format:
			 	break #finish reading-loop, do some more setup and resume reading file
		if not self.vars:
			raise library.MissingInput("cannot find VARS entry in header")
		if not self.format:
			raise library.MissingInput("cannot find FORMAT entry in header")
		if len(self.format.split())!=len(self.vars):
			raise library.InconsistentInput("VARS and FORMAT entry have different length: \n%s \n%s"%(self.vars, self.format))

	def read_body( self, fd ):
		floats=[]
		ints=[]
		lists=[]
		optionals=[]
		f_tags=self.format.split()
		min_tags = len(self.vars)
		for f in f_tags:
			floats.append( "f" in f or "g" in f or "e" in f );
			ints.append( "d" in f );
			lists.append( "l" in f );
			optionals.append( "o" in f );
			if "o" in f:
				min_tags -= 1

		col_resid=self.vars.index("RESID")
		try:
			col_atom=self.vars.index("ATOMNAME")
			self.has_atom=True
		except:
			self.has_atom=False

			self.table={}
		for line in fd:
			tags=line.split()
			if len(tags)>=min_tags and len(tags)<=len(self.vars):
				resi=int(tags[col_resid])
				try:
					atom=tags[col_atom]
					if atom=="HN":
						atom="H"
					key=(resi,atom)
				except:
					key=resi
				value=tags
#				print tags, floats
				try:
					for ct in range(0, len(tags) ):
						if floats[ct]: value[ct]=float(tags[ct])
						elif ints[ct]: value[ct]=int(tags[ct])
						elif lists[ct]: value[ct]=[ int(t) for t in tags[ct].replace("[","").replace("]","").split(",")]
					self.table[key]=value
				except:
#					print "Problem parsing line ", line
#					print "Expected columns: ",self.vars
					raise

			elif len(tags)>0:
				raise library.InconsistentInput("expected %d entries (from VARS) but found only %d in line %s.\n VARS: %s"%(len(self.vars),len(tags),line," ".join(self.vars)))

	def from_dict( self, dict, exclude=[] ):
		import sets
		self.vars=dict.keys()
		for bad in exclude:
			try:
				pos = self.vars.index( bad )
				del self.vars[ pos ]
			except IndexError:
				pass
		self.table={}
		#got through data line by line (i - stores line-number)
		for i,resi in enumerate(dict['RESID']):
			try:
				atom=dict['ATOMNAME'][i]
				self.has_atom=True
				key=(resi,atom)
			except IndexError:
				key=resi
#			print key,[ tab[i] for column,tab in dict.iteritems() if not column in exclude ]
			self.table[key]=[ tab[i] for column,tab in dict.iteritems() if not column in exclude ]

		if 'RESNAME' in dict:
			fasta=get_sequence_from_column( self.get_slice('RESNAME') )
			self.data['SEQUENCE']=fasta
		return self

	def set_sequence( self, sequence ):
		sys.stderr.write('set sequence... \n')
		self.data['SEQUENCE']=sequence
		try:
			col_res1 = self.get_col('RESNAME')
		except ValueError:
			self.vars.append('RESNAME')
			self.format+=" %3s"
			col_res1 = None

		try:
			col_res3 = self.get_col('RESNAME3')
		except ValueError:
			self.vars.append('RESNAME3')
			self.format+=" %5s"
			col_res3 = None

		for key,val in self.table.iteritems():
#				val=self.table[key]
				try:
					resi=key[0]
				except:
					resi=key
				try:
					aa=sequence[resi-1]
				except IndexError:
					raise library.InconsistentInput('trying to find residue %d in fasta-string %s but string is too short'%(resi,sequence))
				aa3=amino_acids.short_to_long[aa]
				if col_res1: val[col_res1]=aa
				else: val.append(aa)
				if col_res3: val[col_res3]=aa3
				else: val.append(aa3)

		#renumber resid indices according to start...end
	def renumber( self, start, end ):
			col_resid=self.vars.index("RESID")
			new_table={}
			for key in self.table.keys():
				val=self.table[key]
				try:
					k=key[0]
					new_key=(key[0]-start+1, key[1])
					val[col_resid]=new_key[0]
				except:
					k=key
					new_key=key-start+1
					val[col_resid]=new_key
				if k<start or k>end:
					continue
				new_table[new_key]=val
			self.table=new_table

	def write_file( self, file ):
		fd=open(file,'w')
		self.write( fd )

	#write to file descriptor
	def write( self, fd ):
		self.write_header( fd )
		self.write_body( fd )

	def write_header( self, fd ):
		if not self.vars or not self.format: return
		for d in self.data.keys():
			out_str = self.data[d]
			for i in range(0,len(out_str),10):
				if i%50==0: fd.write("\n%sDATA %s"%(self.header_prefix,d));
				fd.write(" "+out_str[i:(i+10)]);
			fd.write("\n");

		if self.data:
			fd.write("\n")
		fd.write("%sVARS %s\n"%(self.header_prefix," ".join(self.vars)))
		fd.write("%sFORMAT %s\n\n"%(self.header_prefix,self.format))
#		print self.vars
#		print self.format

	def write_body( self, fd ):
		keys=self.table.keys()
		keys.sort()

		for key in keys:
			entry=self.table[key]
			try:
				fd.write((self.format%tuple(entry))+"\n")
			except TypeError:
				print 'cannot write ', entry, ' with the given format ', self.format
				raise

	def get( self, resi, atom, key ):
		col=self.vars.index( key )
		return self.table[ (resi, atom) ][ col ]

	def put( self, resi, atom, key, value ):
		col=self.vars.index( key )
		self.table[ (resi, atom)][col]= value

	#get column index
	def get_col( self, key ):
		return self.vars.index( key )

	def has_col( self, key ):
		return key in self.vars

	def __contains__(self, key):
		return key in self.vars

	def copy( self, table_in ):
#		print self.vars
		col_map = [ table_in.get_col(v) for v in self.vars ]
		self.table={}
		for key in table_in.table.keys():
			self.table[key]=[ table_in.table[ key ][ col ] for col in col_map ]

	#extract data from specific column
	def get_slice( self, column ):
		col=self.get_col( column )
		slice={}
		for key in self.table.keys():
			slice[key]=self.table[key][col]
		return slice

	def _add_column_field_and_return_index( self, column, in_front_of, format ):
		col=len(self.vars)
		if in_front_of:
			col=self.get_col( in_front_of )
		self.vars.insert( col, column)
		if format:
			tags=self.format.split('%')
			tags.insert(col+1, format[1:]+' ' )
			self.format='%'.join( tags )
		return col

	#insert data to specific column
	#column: column key (string)
	#values: (resid,atom)->value pairs
	def add_slice( self, column, values, in_front_of=None, format='' ):
		if not column in self.vars:
			col = self._add_column_field_and_return_index( column, in_front_of, format )
			insert = True
		else:
			col=self.get_col( column )
			insert=False

		for key in values.keys():
			if key in self.table:
				if insert:
					self.table[key].insert( col, values[key] )
				else:
					self.table[key][col]=values[key]
			else:
				fields=[0]*len(self.vars)
				colr=self.vars.index("RESID")
				colrn=None
				try:
					fields[colr]=key[0]
					cola=self.vars.index("ATOMNAME")
					fields[cola]=key[1]
				except:
					fields[colr]=key
				if "RESNAME" in self.vars:
					colrn=self.vars.index("RESNAME")
				if colrn and "SEQUENCE" in self.data:
					fields[colrn]=self.data["SEQUENCE"][fields[colr]-1]
				self.table[key]=fields
		return self

	#insert data to specific column
	#column: column key (string)
	#data: data (! must be exactly same length as existing table and ordered already in correct way)
	def add_column( self, column, data, in_front_of=None, format='' ):
		if len(self.table) != len(data):
			raise library.InsonsistentInput("data inserted as new column must be same length as existing table")
		if not column in self.vars:
			col = self._add_column_field_and_return_index( column, in_front_of, format )
			insert = True
		else:
			col=self.get_col( column )
			insert=False

		for key,value in zip(self.table.keys(),data):
			if insert:
				self.table[key].insert( col, value )
			else:
				self.table[key][col]=value
		return self

class TalosCSFile(NIH_table):
	def __init__(self):
		NIH_table.__init__(self)
		self.sequence=None
#		self=None

	def read_file(self, file ):
#		self=NIH_table()
		NIH_table.__init__(self)
		NIH_table.read_file( self, file )
		if 'SEQUENCE' in self.data:
			self.sequence=self.data['SEQUENCE']

	def write( self, fd ):
#lf:
#			return
		if self.sequence:
			self.data['SEQUENCE']=self.sequence
		NIH_table.write( self, fd )

	def write_file( self, file ):
		fd=open( file, 'w' )
		self.write(fd)

	def renumber( self, start, end ):
		self.sequence= self.sequence[(start-1):end]
		NIH_table.renumber( self, start, end )
		self.data['SEQUENCE']=self.sequence

	#make TalosCSFile from a general instance of NIH_Table
	def from_table( self, table_in, sequence=None ):
#		self.table = NIH_table()
		NIH_table.__init__(self)
		if 'RESNAME' in table_in.vars:
			self.vars = ['RESID','RESNAME','ATOMNAME','SHIFT']
		else:
			raise library.MissingInput("require sequence information to make talos-table out of prot-table")
		self.format='%4d %1s %4s %8.3f'
		if 'AMBIGUITY' in table_in.vars:
			self.vars += ['AMBIGUITY']
			self.format+=' %2d'

		self.copy( table_in )
		self.sequence = sequence
		if not self.sequence and 'SEQUENCE' in table_in.data:
			self.sequence = table_in.data['SEQUENCE']
			self.data['SEQUENCE']=self.sequence
		#figure out if RESNAME needs translating from aa3-->aa1
		resn=self.get_slice( 'RESNAME' )
		if len(resn[resn.keys()[0]])==3:
			resn_translated=table_op( resn, resn, lambda x,y: amino_acids.longer_names[x.upper()] )
			self.add_slice( 'RESNAME', resn_translated )
		elif not len(resn[resn.keys()[0]])==1:
			raise library.InconsistentInput("RESNAME column in input table should have either aa3 or aa1 format, i.e., ALA or A for residue names")
		return self

def get_sequence_from_column( data ):
	min_resi = 1000000
	max_resi = -min_resi
	for key,a in data.iteritems():
		if min_resi > key[0]: min_resi = key[0]
		if max_resi < key[0]: max_resi = key[0]
	sequence=list("-"*(max_resi))
	for key,a in data.iteritems():
		if len(a)==1:
			sequence[key[0]-1]=a
		else:
			sequence[key[0]-1]=amino_acids.long2short( a )
	return ''.join(sequence)

class ProtCSFile(NIH_table):
	def __init__(self):
		NIH_table.__init__(self)
		self.sequence=None
		self.header=True
#		self.table=None
	def auto_header( self, fd, has_aa3, has_aa1 ):
		headers = [ ( ['INDEX','SHIFT','SIGMA','ATOMNAME','RESID'], '%8d %10.3f %10.3f %7s %8d' ),
								( ['INDEX','SHIFT','SIGMA','ATOMNAME','AMBIGUITY','RESID'], '%8d %10.3f %10.3f %7s %5d %8d' ) ]

		for header in headers:
			self.vars = header[0]
			self.format = header[1]
#			print self.vars, self.format
			if has_aa3:
				self.vars.append('RESNAME3')
				self.format=self.format+' %6s'
			if has_aa1:
				self.vars.append('RESNAME')
				self.format=self.format+' %4s'
		  #read body
			auto_add = ['INTENSITY','STEREO','DONE']
			auto_fmt = [' %4d',' %ol','%z']
			done = False
			last_exc = None
			for auto_header,fmt in zip(auto_add,auto_fmt):
				try:
#					print 'Try: ', self.vars
					self.read_body( fd )
					done = True
					return
				except (ValueError, library.InconsistentInput) as exc:
					self.vars.append(auto_header)
					self.format+=fmt
					last_exc = exc
					fd.seek(0)

		if not done:
			raise last_exc

	def read(self, fd, has_aa3, has_aa1, sequence=None, header=False ):
		#self.table=NIH_table()
		NIH_table.__init__(self)
		if not header:
			self.header=False
			self.auto_header( fd, has_aa3, has_aa1 )

		else: #has header
			self.header=True
			self.read_header( fd )
			self.read_body( fd )
			has_aa3 = 'RESNAME3' in self.vars
			has_aa1 = 'RESNAME' in self.vars

		if not has_aa1:
			if has_aa3:
				aa3=self.get_slice('RESNAME3')
				aa={}
				for key,a in aa3.iteritems():
					aa[key]=amino_acids.longer_names[ a.upper() ]
				self.add_slice('RESNAME', aa)
				self.format=self.format+" %4s"
			elif sequence:
				aa3=self.get_slice('INDEX')
				aa={}
				for key,a in aa3.iteritems():
					aa[key]=sequence[key[0]-1]
				self.add_slice('RESNAME', aa)
				self.format=self.format+" %4s"
		if not sequence and (has_aa1 or has_aa3):
			aa1=self.get_slice('RESNAME')
			sequence = get_sequence_from_column( aa1 )

		if 'SEQUENCE' in self.data:
			sequence = self.data['SEQUENCE']
		return sequence

	def renumber( self, start, end ):
		if self.sequence:
			self.sequence= self.sequence[(start-1):end]
			self.data['SEQUENCE']=self.sequence
		NIH_table.renumber( self, start, end )



  def read_file( self, file, sequence=None ):
		fd=open(file,'r')
		self.read_stream( fd, sequence )

	def read_stream( self, fd, sequence=None ):
		has_header = False
		has_aa1 = False
		has_aa3 = False
		for l in fd:
			tags=l.split()
			if len(tags)<1: continue
			if 'VARS' in l and 'SHIFT' in l:
				has_header = True
				has_aa3 = 'RESNAME3' in tags
				has_aa = 'RESNAME' in tags
				break

#			print len(tags), has_aa3, has_aa1, len(tags[5]), len(tags[len(tags)-1])
		if not has_header:
			fd.seek(0)
			for l in fd:
				tags=l.split()
				num_string = 0
				index_string = []
				for ind,t in enumerate(tags):
					try:
						d = float(t)
					except:
						num_string += 1
						index_string.append(ind)
					if num_string >= 2 and len(tags[index_string[1]])==3:
						has_aa3 = True
				if len(tags)>=6 and len(tags[-1])==1:
					has_aa1 = True
				if len(tags)>6 and len(tags[-2])==3:
					has_aa3 = True

	#	assert( fd.seekable() )
		fd.seek(0)
		self.sequence = self.read(fd, has_aa3, has_aa1, sequence=sequence, header=has_header )


	#make ProtCSFile from a general instance of NIH_Table
	def from_table( self, table_in, sequence=None ):
		def insert_var( self, pos, var, format_str ):
			self.vars.insert(pos,var)
			fs = self.format.split()
			fs.insert( pos, ' '+format_str )
			self.format = " ".join(fs)

		NIH_table.__init__(self)

		self.vars = ['SHIFT','ATOMNAME','RESID']
		self.format='%10.3f %7s %8d'
		if 'RESNAME' in table_in.vars and 'RESNAME3' in table_in.vars:
			self.vars += ['RESNAME3','RESNAME']
			self.format += ' %6s %4s'
		elif 'RESNAME' in table_in.vars:
			self.vars += ['RESNAME']
			self.format += ' %4s'
		elif 'RESNAME3' in table_in.vars:
			self.vars += ['RESNAME3']
			self.format += ' %6s'
		if 'SIGMA' in table_in.vars:
			insert_var( self, 1, 'SIGMA', '%4.2f')
		if 'AMBIGUITY' in table_in.vars:
			pos = self.vars.index('RESID')
			insert_var( self, pos, 'AMBIGUITY', '%2d' )
		if 'STEREO' in table_in.vars:
			self.format += ' %s'
			self.vars.append('STEREO')
		if 'INDEX' in table_in.vars:
			insert_var( self, 0, 'INDEX', '%8d' )

#		print self.vars
#		print table_in.vars

		self.copy( table_in )
		self.sequence = sequence

		lines=self.get_slice( 'RESID' )
		indices={}
		sigmas={}
		keys=lines.keys()
		keys.sort()
		for id,key in enumerate(keys):
			indices[key]=id+1
			sigmas[key]=0.00
		if not 'INDEX' in self.vars:
			self.add_slice( 'INDEX', indices, in_front_of='SHIFT', format='%8d' )
		if not 'SIGMA' in self.vars:
			self.add_slice( 'SIGMA', sigmas, in_front_of='ATOMNAME', format='%10.3f' )
		if not self.sequence and 'SEQUENCE' in table_in.data:
			self.sequence = table_in.data['SEQUENCE']
			self.data['SEQUENCE']=self.sequence
		return self

	def write( self, fd, header=None ):
		if header is None:
			header=self.header
		if header:
			if self.sequence:
				self.data['SEQUENCE']=self.sequence
			self.write_header( fd )
		self.write_body( fd )

	def write_file( self, filename, header=None ):
		self.write( open( filename, 'w'), header )

#SpartaCalculator: allows to spawn Sparta calculations:
#	  refCS as reference chemical shifts
#   pdb the input pdb structure
#   pred_dir a directory name to store all Sparta data files, these are automatically reused
#   for subsequent runs with same pdb structure, unless overwrite=True.
#   if the Sparta class is instantiated for same pdb multiple times it will be copied from
#   previous instance using the cached 'SpartaStore'

class SpartaCalculator:
	def __init__( self, refCS, pred_dir, overwrite=False ):
		self.refCS=refCS
		self.pred_dir=pred_dir
		self.overwrite=overwrite
		self.SpartaStore={}
	def get_sparta( self, pdb ):
		return self.Sparta( pdb, self.refCS, self.pred_dir, self.overwrite, self.SpartaStore )

	class Sparta:
	  def __init__( self, pdb, refCS, pred_dir, overwrite, SpartaStore ):
		  if pdb in SpartaStore:
				old=SpartaStore[pdb]
#				print "retrieve: ", pdb, id(old)
				self.__copy_constr__( old )
				self.refCS = refCS
			else:
				library.mkdirp( pred_dir )
				self.refCS = refCS
				self.pred_dir = pred_dir
				self.pdb = pdb
				self.name= splitext( basename(pdb) )[0]
				self.overwrite = overwrite
				self.pred_fn = join( self.pred_dir, "pred_"+self.name+".tab" )
				self.struct = join( self.pred_dir, "struct_"+self.name+".tab" )
				self.csCalc_fn = join( self.pred_dir, "csCalc_"+self.name+".tab" )

				self.pred = NIH_table()
				self.csCalc = NIH_table()

				if not self.results_exist() or self.overwrite:
					self.runSparta()
				self.read_results()
				SpartaStore[pdb]=self;
			#print "generate: ", pdb, id(self)


		def __copy_constr__( self, old ):
			self.refCS = old.refCS
			self.pred_dir = old.pred_dir
			self.pdb = old.pdb
			self.name = old.name
			self.overwrite = old.overwrite
			self.pred_fn = old.pred_fn
			self.struct = old.struct
			self.csCalc_fn = old.csCalc_fn
			self.pred = old.pred
			self.csCalc = old.csCalc

		def __str__( self ):
			return "SPARTA: refCS: %s, pred_dir: %s, pdb: %s\n pred: %s, struct: %s, csCalc: %s\n"%(self.refCS,\
                                                         self.pred_dir, self.pdb, self.pred_fn, self.struct, self.csCalc_fn )

		def runSparta( self ):
			commands = ["sparta+",
									"-in %s -out %s -outS %s -ref %s -offset -outCS %s"%( self.pdb, self.pred_fn, self.struct, self.refCS, self.csCalc_fn ) ]
			print 'Execute: ',' '.join(commands)

			subprocess.call(commands, stderr=subprocess.PIPE)

		def read_results( self ):
			if exists( self.pred_fn ):
				self.pred.read_file( self.pred_fn )
			if exists( self.csCalc_fn ):
				self.csCalc.read_file( self.csCalc_fn )

		def results_exist( self ):
			return exists( self.pred_fn ) and ( exists( self.csCalc_fn ) or not self.refCS )

		def get_slice( self, column ):
			try:
				return self.pred.get_slice( column )
			except ValueError:
				try:
					return self.csCalc.get_slice( column )
				except ValueError:
					raise library.MissingInput("field %s is neither found in pred_ nor in csCalc_ files of Sparta+ output"%(column) )

#End Sparta class definition


def table_op( table1, table2, func, selection=None, default=0 ):
    table={}
    for key in table1.keys():
        try:
					if not selection or ( key[0] in selection.resi ):
            table[key]=func( table1[key], table2[key] )
					else:
						table[key]=default
        except KeyError:
            pass
    return table


def average_shifts( column ):
    global args
    pdblist=args.pdbs;
    pdb = pdblist[ 0 ]
    sparta = Sparta( pdb )
    tab=sparta.get_slice( column );
    #setup empty tables
    av_shifts=table_op( tab, tab, lambda x,y: 0.0 );
    av2_shifts=av_shifts;
    inv_N = 1.0 / len(pdblist)
    for pdb in pdblist:
        sparta = Sparta( pdb )
#        print "got : ", id( sparta )
        new_tab=sparta.get_slice( column )
        av_shifts=table_op( av_shifts, new_tab, lambda x,y: x+y*inv_N )
        av2_shifts=table_op( av2_shifts, new_tab, lambda x,y: x+y*y*inv_N )
    #subtract square of average : std = <x^2>-<x>^2
    std_shifts=table_op( av2_shifts, av_shifts, lambda x,y: x-y*y )
    return av_shifts,std_shifts

#by default no extra weighting per atom, since Sparta+ results, for instance, are already weighted according to sigma
def write_to_bfactor( pdbin, pdbout, result, sum_residue=True, \
                      func=(lambda x: x), atoms={"N":1.0,"H":1.0,"C":1.0,"CB":1.0,"CA":1.0,"HA":1.0}, \
                      response=(lambda x: x)
                      ):

	with warnings.catch_warnings(record=True) as w:
		parser=PDBParser()
		s=parser.get_structure( "t000_", pdbin )

	final_values={}
	for chain in s[0]:
		for residue in chain:
			resid=residue.get_id()[1]
			sum=0;
			if sum_residue:
				for atom in atoms:
					try:
						sum=sum+func(result[(resid,atom)])*atoms[atom]
					except KeyError:
						pass
				for atom in residue:
					atom.set_bfactor( response( sum ) )
				final_values[resid]=response( sum )
			else:
				for atom in residue:
					try:
						v=response( result[(resid,atom.get_name())] )
						atom.set_bfactor( v )
						final_values[atom]=v
					except KeyError:
						atom.set_bfactor( 0.0 )



	io=PDBIO()
	io.set_structure( s )
	io.save( pdbout )
	return final_values

def sum_table( table, selection ):
	sum=0
	for key in table.keys():
		if key[0] in selection.resi:
			sum = sum + table[key]
	return sum

class AtomSelection:
	def __init__( self, selection ):
		self.selection=selection
		self.resi=[]
		if selection:
			tagsp=selection.split('+')
			for t in tagsp:
				if '-' in t:
					tags=t.split('-')
					for i in range(int(tags[0]),int(tags[1])+1):
						self.resi.append(i)
				else:
					self.resi.append(int(t))
		else:
			self.resi=range(1, 1000)


