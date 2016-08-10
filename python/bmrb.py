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
import string
### toolbox library
import library
from os import mkdir , makedirs

from warnings import *

import traceback

#from Bio.PDB.PDBParser import *
#from Bio.PDB import PDBIO
#from numpy import square

from library import mkdirp
import library

class BMRB_Type :
	def __init__( self, tags ):
		self.tags=tags
		self.type=int

class BMRB_Int( BMRB_Type ):
	def __init__( self, tags ):
		BMRB_Type.__init__( self, tags )
		self.type=int

class BMRB_Str( BMRB_Type ):
	def __init__( self, tags ):
		BMRB_Type.__init__( self, tags )
		self.type=str

class BMRB_Float( BMRB_Type ):
	def __init__( self, tags ):
		BMRB_Type.__init__( self, tags )
		self.type=float


_tag_dictionary = {'ID': BMRB_Int(['_Atom_chem_shift.ID','_Atom_shift_assign_ID']),
         'RESID': BMRB_Int(['_Atom_chem_shift.Comp_index_ID','_Residue_seq_code','_Atom_one_residue_seq_code','Seq_ID_1','Seq_ID','Comp_index_ID']),
         'RESNAME': BMRB_Str(['_Atom_chem_shift.Comp_ID','_Residue_label','_Atom_one_residue_label','Comp_ID_1','Comp_ID']),
         'ATOMNAME': BMRB_Str(['_Atom_name','_Atom_chem_shift.Atom_ID','_Atom_one_atom_name','Atom_ID_1','Atom_ID']),
         'SHIFT': BMRB_Float(['_Chem_shift_value','_Atom_chem_shift.Val','_Chem_shift_value','Val','Chem_shift_val']),
         'SIGMA': BMRB_Float(['_Chem_shift_value_error','_Atom_chem_shift.Val_err','Val_err','_Residual_dipolar_coupling_value_error','Chem_shift_val_err']),
				 'RDC_TYPE': BMRB_Str(['_Residual_dipolar_coupling_ID','RDC_code']),
				 'RESID2': BMRB_Int(['_Atom_two_residue_seq_code','Seq_ID_2']),
				 'RESNAME2': BMRB_Str(['_Atom_two_residue_label','Comp_ID_2']),
				 'ATOMNAME2': BMRB_Str(['_Atom_two_atom_name','Atom_ID_2']),
				 'AMBIGUITY': BMRB_Int(['_Atom_chem_shift.Ambiguity_code','Ambiguity_code','_Chem_shift_ambiguity_code','_Chem_shift_ambiguity_type']),
				 'RDC': BMRB_Float(['_Residual_dipolar_coupling_value','Val'])
}

cs_loop_cols=['ID','RESID','RESNAME','ATOMNAME','SHIFT','SIGMA','AMBIGUITY']
rdc_loop_cols=['RDC_TYPE', 'RESID','RESNAME','ATOMNAME','RESID2','RESNAME2','ATOMNAME2','RDC','SIGMA']

#after reading the BmrbFile content is organized in frames of name <frame> which are organized in self._frames
#according to their category <frame-category>. A frame goes from 'save_<frame>' to 'save_'
#according to this entry:
#save_<frame>
#   _Saveframe_category      <frame-category>
# .... content
#save_
#From each frame we currently store only the 'loop', that is lines between 'loop_' and 'stop_'. The beginning of a
#loop has entries starting with '_XXX' which give the column names.
#we store each loop of a frame as object of class 'Loop' which has the data-members cols (the column headers) and data
#which is just a list with an entry per line. Each line-entry is a list of tags for this line.
#to extract data from a frame we use 'process_frame( <frame_category>, columns )'
#the columns are translated using the _tag_dictionary which takes care of ambiguous nameing found in different BMRB styles.
#the output is a dictionary {'colname1': [data1, data2, data3, ... , dataN ], 'colname2':[data1, data2, ..., dataN] }
#if multiple frames of the same category have fitting columns these will be appended to the dictionary...
from basic.Tracer import Tracer
tr = Tracer( "bmrb" )
#reads a BmrbFile into _frames
class BmrbFile:

	###################################################################
	# class Loop
	class Loop:
		def __init__(self, cols, data):
			if len(data):
				self.n_cols_in_data = len(data[0])
				self.n_cols = len(cols)
				self.n_repeats_per_line = self.n_cols_in_data/self.n_cols
				if self.n_cols_in_data % self.n_cols:
					msg='Loop-Header is inconsistent with data-entries, check number of columns\n'
					msg+='-'*80+'\n'
					msg+='Header: \n  '
					msg+='\n  '.join(cols)
					msg+='\n\nFirst-data entry: '+' | '.join(data[0])+'\n'
					msg+='-'*80
					raise library.InconsistentInput(msg)
			else:
				msg='Loop has no data-entries\n'
				msg+='-'*80+'\n'
				msg+='Header: \n'
				msg+='\n  '.join(cols)
				msg+='-'*80
				raise library.InconsistentInput(msg)

			self.cols=cols
			self.data=data

		def size( self ):
			return len(self.cols)

		def __str__( self ):
			return "[ "+", ".join(self.cols)+" ]"

		def __repr__( self ):
			return self.__str__()

		def process( self, cols ):
			trloop = Tracer( 'loop', tr )
			ind={}
			types={}

			for col in cols:
				for tag in _tag_dictionary[col].tags:
					if tag in self.cols:
						ind[col]=self.cols.index(tag);
						types[col]=_tag_dictionary[col].type
						break

			#should have now the indices of the requested columns and their type in ind and types
			#if no fitting columns return
			if len(ind)==0: return
			trloop.Debug('labels (bmrb):  ', self.cols)
			trloop.Debug('column indices: ', ind)
			#extract output dictionary
			output={}
			for col in cols:
				output[col]=[]

			#lines are already split into columns
#figure out which bmrb-type columns fit to the requested columns (in cols) in this loop

			for line in self.data:
				if line[0][0]=='#': continue
				for i,col in enumerate(cols):
#					print 'F', i, col, line, self.n_repeats_per_line
					if not col in ind.keys() or line[ind[col]]=='.':
						trloop.Debug('found empty columns: ',cols[i], ' and will remove it')
						del cols[i]
						return self.process( cols )
					for repeat in range(0,self.n_repeats_per_line):
#						print output[col]
#						print types[col]
#						print ind[col]
#						print 'colsindata',self.n_cols_in_data
#						print ind[col]+self.n_cols_in_data*repeat, repeat
						output[col].append(types[col](line[ind[col]+self.n_cols_in_data/self.n_repeats_per_line*repeat]))
			return output


  ####################################################################
	# class Frame:
	# a frame starts with save_XXXX and ends with save_
	# frames contain 'fields' and 'loops'. A loop is basically a table, first define what the columns are, then the data
	class Frame:
		def __init__( self, name, fields, loops ):
			self.name=name
			self.fields=fields
			self.loops=loops
		def __repr__( self ):
			str="Frame %s: \n"%self.name
			for k,f in self.fields.iteritems():
				str=str+"%s: %s\n"%(k,f)
			str=str+"and %d loops\n"%len(self.loops)
			return str

	###########################################################
	#BmrbFile
	def __init__(self, file, errors=None ):
		self._frames={}
		self.parse_file( open(file,'r'), errors )
		self.star3=False

	def get_frame_category(self,line):
		tags=line.split()
		if len(tags)<1: return None
		try:
			s3tag=tags[0].split('.')[1]
			if s3tag!='Sf_category': return None
			self.star3=True
		except IndexError:
			if tags[0]!='_Saveframe_category': return None
			self.star3=False
		return tags[1]

	#-------------------------------------------------------
	# main capture routine= read a frame from save_XXXX to save_
	# find next save_<frame>
	# return the name, i.e., <frame> and the category <frame-category>
	def Next_Save_Frame( self, file ):
		name=''
    for line in file:
			tags=string.split(line)
			tr.Debug('SKIP BETWEEN FRAMES: ',line[:-1])
			if len(tags)>0 and len(tags[0])>=5 and tags[0][:5]=='save_':
				tr.Debug('READ FRAME: ',line[:-1])
				name=tags[0][5:]
				for line in file:
					category=self.get_frame_category( line)
					if category: break
				return category, name
		return 'NO_CATEGORY', 'NO_NAME'


	#-----------------------------------------------------
	# read fields and loops of current frame
	# return as Loops and Fields (classes see above )
	def capture_loops( self, file, errors ):
		loops=[]
		fields={}
		multi_line_field=None
		col_nr=-1;
		for line in file:
			tr.Debug('READ FRAME :',line[:-1])
			if "'" in line or '"' in line:
				tags=[]
				within_field=False
				current_field=""
				for c in line:
					if c == "'" or c=='"':
						if within_field and len(current_field):
							tags.append(current_field)
							current_field=''
							within_field=False
						else:
							within_field=True
						continue
					if within_field:
						current_field+=c
					elif c==' ' and len(current_field):
						tags.append(current_field)
						current_field=''
					elif c in string.whitespace:
						continue
					else:
						current_field+=c
#				print 'T', tags, line
			else:
				tags=string.split(line)
			if len(tags) and tags[0][0]=='#':
				continue
			for i,t in enumerate(tags):
				if t[0]=='#':
					tags=tags[0:i]
					break
			if len(tags)>0 and tags[0]=='loop_':
				#this is the sart of a loop
        col_nr=0;
        col_id=[];
        data=[];
        continue
			if col_nr<0 and len(tags)>0 and tags[0][0]=='_':
				#this is a field for the frame
				if self.star3:
					fkey=tags[0].split('.')[1]
				else:
					fkey=tags[0]
				if len(tags)>1:
					fval=tags[1]
					fields[fkey]=fval
				else: multi_line_field='START'
				continue
			if col_nr>=0 and len(tags)>0 and tags[0][0]=='_':
				#this is a field for the current loop
				if self.star3:
					name=tags[0].split('.')[1]
				else:
					name=tags[0]
        col_id.append(name);
        col_nr+=1;
        continue
			if col_nr>=0 and len(tags)>0 and tags[0]=='stop_':
				#end of a loop
				try:
					loops.append( self.Loop( col_id, data ))
				except library.InconsistentInput as exc:
					exc.add("This loop will be ignored. Hopefully nothing important !!!")
					if not errors is None: errors.append(exc)
					else: raise exc

				col_nr=-1
        continue
			if col_nr>=0 and len(tags)>0:
				#data entry of current loop
        data.append(tags);
        continue
			if len(tags)>0 and tags[0]=='save_':
				return loops, fields
			if len(tags)>0 and tags[0][0]==';' and multi_line_field:
				if multi_line_field=='START':
					multi_line_field='CURRENT'
					mlf_data=[]
				elif multi_line_field=='CURRENT':
					multi_line_field=None
					fields[fkey]=mlf_data
				continue
			if len(tags)>0 and multi_line_field=='CURRENT':
				mlf_data.append(line[:-1])

		if col_nr>=0: #col_nr >= 0 means we are still in a loop --- EOF without getting 'stop_'
			try:
				loops.append( self.Loop( col_id, data ) )
			except library.InconsistentInput as exc:
				exc.add("This loop will be ignored. Hopefully nothing important !!!")
				if not errors is None: errors.append(exc)
				else: raise exc

		return loops, fields

	#go through all frames and store the respective loops
	def parse_file( self, file, errors=None ):
		MAX_FRAMES=1000  #to avoid hanging in corrupted files
		while MAX_FRAMES>0:
			MAX_FRAMES-=1
			SAVE_FRAME_CATEGORY, name=self.Next_Save_Frame( file );
			tr.Info( 'reading frame ( %s ): %s ...'%(SAVE_FRAME_CATEGORY, name))
			if SAVE_FRAME_CATEGORY=='NO_CATEGORY':
				break

			loops, fields = self.capture_loops( file, errors )
			self.add_frame(name, loops, fields, SAVE_FRAME_CATEGORY )

	def add_frame(self, name, loops, fields, CATEGORY='GENERIC'):
#		print 'ADD: ', name, loops, fields, CATEGORY
		if not CATEGORY in self._frames:
			self._frames[CATEGORY]={}
		self._frames[CATEGORY][name]=self.Frame(name, fields, loops )


	#how many different frame categories have been found ?
	def nframes( self ):
		return len( self._frames )

	#extract columns according to _tag_dictionary from loop

	def get_frame( self, categories ):
		for category in categories:
			try:
				frames=self._frames[category]
				return category, frames
			except KeyError:
				pass
		raise KeyError

#process frames of a certain category and return output that fits the requested columns
  def process_frame( self, categories, cols ):
		outputs=[]
		frames=None
		for category in categories:
			try:
				frames=self._frames[category]
				break
			except KeyError:
				pass
		if not frames:
			raise library.MissingInput("Cannot find category %s. Maybe not a proper BMRB file?"%category)

		for frame in frames.itervalues():
			for loop in frame.loops:
				#print 'L: ', loop
				tr.Debug("process loop with %3d columns in frame %s"%(loop.size(),frame.name))
				tr.Debug("loop: %s"%loop)
				output=loop.process( copy.copy(cols) )
				if output:
					outputs.append((len(output),frame.name,output))

		outputs=sorted(outputs,key=lambda x: -len(x[2].keys()))
		return outputs


_SEQUENCE_FIELD_KEYS=['_Mol_residue_sequence','Polymer_seq_one_letter_code']


def get_sequence( bmrb_file, verbose=0 ):
	categories=['monomeric_polymer','entity']
#	try:
	category,seq_frames=bmrb_file.get_frame(categories)
#	except KeyError:
#		raise library.InconsistentInput('Cannot find frame describing the molecule: %s'%' '.join(categories))

	sequences={}
	for frame in seq_frames.itervalues():
		for key in _SEQUENCE_FIELD_KEYS:
			try:
				sequences[frame.name]=''.join([ x.strip() for x in frame.fields[key]])
			except KeyError:
				pass
	if len(sequences.keys())<1:
		raise  library.InconsistentInput("Cannot find any of sequence fields '%s'. Maybe not a proper BMRB file?"%"','".join(_SEQUENCE_FIELD_KEYS))
	return sequences[sequences.keys()[0]], sequences, len(sequences)


def read_fullbmrb_or_trunkfile( filename, types, verbose=0, errors=None ):
	fasta=None
	if verbose==0:
		tr.set_priority( 1 )
	bmrb_file = BmrbFile( filename, errors )
	try:
		seq_frame=bmrb_file.get_frame(['monomeric_polymer','entity'])
	except KeyError:
		seq_frame=None

	cs_frame=None
	for type in types:
		try:
			cs_frame=bmrb_file._frames[type]
			break
		except KeyError:
			pass

	if not cs_frame:
		type=types[0]
		#reading under assumption of full STAR structure didn't work, reread from start this time assume that loops are at top-level
		loops, fields = bmrb_file.capture_loops( open(filename, 'r' ), errors )
		tr.Debug('L', loops)
		tr.Debug('F', fields)
		bmrb_file.add_frame('generic',loops,fields,type);
		seq_key=None
		seq_key=set(fields).intersection( set( _SEQUENCE_FIELD_KEYS ))

		if not seq_frame and len(seq_key)>0:
			fasta=''.join([ x.strip() for x in fields[seq_key.pop()]])
			nmol=1
	if verbose and seq_frame:
		tr.Info("*"*30)
		tr.Info("Molecule(s) found in BMRB: %s"%(', '.join(seq_frame[1].keys())))
		tr.Info("*"*30)

	if not fasta:
		try:
			(fasta, fastas, nmol)=get_sequence( bmrb_file, verbose )
		except KeyError:
			fasta=None
			nmol=1

	if nmol>1:
		print "WARNING: found multiple molecules in BMRB"
		print fasta
		exit()

	return bmrb_file, fasta
