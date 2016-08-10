#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'

from os.path import basename
from os import mkdir , makedirs
import string
import amino_acids
import textwrap
import hashlib

def square(x):
	return x*x

def hashfile(afile, hasher, blocksize=65536):
	buf = afile.read(blocksize)
	while len(buf) > 0:
		hasher.update(buf)
		buf = afile.read(blocksize)
	return hasher.digest()

def diff_files_md5(file1,file2):
	return hashfile(open(file1,'r'),hashlib.md5() ) == hashfile(open(file2,'r'),hashlib.md5())


def hello( name ):
	print '-'*75, """
---------             CS-Rosetta 3.0   (toolbox)                -----------
-----                                                                 -----
--     website: www.csrosetta.org                                        --
--     copyright: Oliver Lange                                           --
--     Reference: Lange et al. PNAS 2012,                                --
--         www.pnas.org/cgi/doi/10.1073/pnas.1203013109                  --
--                                                                       --
--     program name: %30s                      --
--                                                                       --
-----                                                                 -----
--------                                                            -------"""%basename(name)
	print '-'*75

def add_standard_args( parser ):
	parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )
	parser.add_argument("-tracer", nargs='*', default=None )

def init( args ):
	from basic.Tracer import init_tracer_from_cmdline
	init_tracer_from_cmdline(args.tracer)

def print_exception( print_traceback = False ):
	import traceback, sys
	if print_traceback:
		print traceback.print_tb( sys.exc_info()[2] )
#		print sys.exc_info()[1]
		print 'Exception: \n','\n'.join(traceback.format_exception_only(sys.exc_info()[0],sys.exc_info()[1]))
	else:
		#		print traceback.print_exception(sys.exc_type,None, None)
		print 'Excception dected: \n'
		print '\n'.join(traceback.format_exception_only(sys.exc_info()[0],sys.exc_info()[1]))

def extract_exception_msg():
	import traceback, sys
	exc=sys.exc_info()[1]
	exc_type=sys.exc_info()[0]
	msg=''
	if exc.message: msg+='%s\n'%str(exc.message)
	else:	msg+='%s\n'%' '.join(traceback.format_exception_only(exc_type,exc))
	return msg

def augment_exception_msg(more_msg):
	import traceback, sys
	msg=extract_exception_msg()
	msg+=more_msg
	sys.exc_info()[1].message=msg
	raise sys.exc_info()[0], sys.exc_info()[0](msg), sys.exc_info()[2]

def exception2str():
	import traceback, sys
	exc=sys.exc_info()[1]
	exc_type=sys.exc_info()[0]
	return (' '.join(traceback.format_exception_only(exc_type,exc)))[:-1]


class LibException(Exception):

	def __init__(self,msg):
		super(LibException,self).__init__(msg)
		self.message=msg

	def __str__(self):
		#        return "\n------------------      %s      --------------------------------------------\n"%(self.msg)
#		s='%s'%str(type(self)).replace("<class '","").replace("'>","").replace('library.','')
		s=''
		for line in self.message.split('\n'):
			s += "\n"+"\n".join( textwrap.wrap( line, width=100 ) )
		return s

	def add( self, str ):
		self.message=self.message+'\n'+str

class StubbedOut(LibException):
	pass

class MissingInput(LibException):
	pass

class InconsistentInput(LibException):
	pass

class RunException(LibException):
	pass

class MethodException(LibException):
	def __init__(self,method,msg):
		LibException.__init__(self,msg)
		self.method=method

	def __str__(self):
		return "\n"+"-"*80+"\n"\
			 +"in "+self.method.__str__()+"\n-----------------------------------\n"\
			 +self.method.module_description()+"\n"\
			 +"-"*80+"\n"\
			 +LibException.__str__(self)+"\n\n"\
			 +"-"*80+"\n"\

class ProgramError(LibException):
	def __str__(self):
		return """
\n\n------------------      %s      --------------------------------------------\n
\n please send a description of what you were doing, the cmd-line and the stack-trace (run with -traceback) to
 to : oliver.lange@tum.de \n
"""%(self.msg)

#returns the full name of prog prog.default.linuxgccrelease
# takes the fist present, limit using the option extras if you only like mpi for instance
def rosetta_executable( prog, extras=['default','static','mpi'] ):
	import subprocess
	platforms=['linux','macox','windows']
	cc=['gcc','icc']
	debug=['release','debug']
	testcmd='type %s >/dev/null 2>&1'
	for di in debug:
		for ci in cc:
			for ei in extras:
				for pi in platforms:
					name='%(prog)s.%(ei)s.%(pi)s%(ci)s%(di)s'%locals()
					try:
						subprocess.check_call(testcmd%name,shell=True)
						return name
					except subprocess.CalledProcessError as exc:
						pass

def read_rigid_file(rigid_file):
	lines=open(rigid_file,'r').readlines()
	start=None
	end=None
	for line in lines:
		tags=line.split()
		if not len(tags): continue
		if tags[0]!='RIGID': raise InconsistentInput('expected RIGID at start of line %s in file %s'%(line,rigid_file))
		s=int(tags[1])
		e=int(tags[2])
		print tags, s, e
		if not start or start>s: start=s
		if not end or end<e: end=e
	return start,end

def obj_is_list( val ):
	is_list=not isinstance(val,basestring)
	try:
		for file in val:
			pass
	except:
		is_list=False
	return is_list



def cst_is_centroid(file):
	lines = open(file).readlines()
	for line in lines:
		if (line[0]=="#"): continue;
		l=string.split(line);
		try:
			resi1 = int( l[2] );
			resi2 = int( l[4] );
		except:
			print '\n\n expected integer in col 3 and 5 of:\n'
			print l[2],l[4],'\n',l
			exit()

		atom1=l[1];
		atom2=l[3];
		centroid=["H","HA","CA","CB","C","O","N","ZN"]
		if not atom1 in centroid: return False
		if not atom2 in centroid: return False
	return True

def read_aa3_sequence(file):
	lines=open(file,'r').readlines()
	sequence=""
	for line in lines:
		try:
			tags=line.split()
			index = int(tags[1])
			if index<len(sequence):
				raise BadInput('require ordered sequence indices in .seq file. Problem at line: %s'%line)
			sequence=sequence.ljust(index-1,'-')+amino_acids.longer_names[ tags[0] ]
		except: #cannot obtain int(tags[1])
			sequence=sequence+amino_acids.longer_names[ tags[0] ]
	return sequence

def write_aa3_sequence(file,sequence):
	fd=open(file,'w')
	write_aa3_sequence_stream(fd,sequence)

def write_aa3_sequence_stream(fd,sequence):
	for ct,c in enumerate(sequence):
		if c=='-': continue
		fd.write('%5s %5d\n'%(amino_acids.short2long(c),ct+1))

def mkdirp( target_dir ):
	try:
#		if ".." in target_dir:
#			mkdir( target_dir );
#		else:
		makedirs( target_dir );
	except OSError as err:
		if not err[1]=="File exists":
			print err;
			raise

def _flatten_list(val, list ):
	if isinstance(val,basestring):
		list.append(val)
		return
	try:
		for i in val:
			_flatten_list(i, list)
	except:
		list.append(val)

def flatten_list( val ):
	list=[]
	_flatten_list( val, list)
	return list

def is_gzipped( file ):
	return file[-3:]==".gz"

def cst_has_HA_atoms( file ):
	lines = open(file).readlines()
	for line in lines:
      if (line[0]=="#"): continue
      l=string.split(line);
      try:
			resi1 = int( l[2] );
			resi2 = int( l[4] );
      except:
			print '\n\n expected integer in col 3 and 5 of:\n'
			print l[2],l[4],'\n',l
			exit()

		atom1=l[1];
      atom2=l[3];
      ha_atom=["HA"]
      if atom1 in ha_atom: return True
      if atom2 in ha_atom: return True
	return False

def _upl2mini_line( line, sequence_offset, max_res, sep=4, pad=0.15 ):
	if line[0]=="#": return ""
	cols=line.split()
	if (len(cols)<5): return ""
	resi=int(cols[0])
	resj=int(cols[3])
	atomi=cols[2]
	atomj=cols[5]
	atomi=atomi.replace("HN","H").replace("M","Q")
	atomj=atomj.replace("HN","H").replace("M","Q")
	if  ( (resi-resj)>=sep or (resj-resi)>=sep ) and (resi-sequence_offset)>0 and (resj-sequence_offset)>0 and (resi-sequence_offset<=max_res) and (resj-sequence_offset <= max_res):
		return "AmbiguousNMRDistance %5s %5d %5s %5d BOUNDED 1.5 %5.3f 0.3 NOE; rawdata %5.3f\n"%(atomi,int(cols[0])-sequence_offset,atomj,int(cols[3])-sequence_offset,float(cols[6])+pad,float(cols[6]))
	else: return ""


# def upl2mini_find_offset( fasta, lines ):
# 	#work out offset first

# 	#build dictionary res-number --> aa3
# 	dict={}
# 	for line in lines:
# 		if line[0]=="#": continue
# 		cols=line.split()
# 		if len(cols)<5: continue
# 		try:
# 			cols[1]
# 			cols[1][0:3]
# 		except:
# 			print 'Exception in line: %s'%line, cols

# 		if int(cols[0]) not in dict:
# 			dict[int(cols[0])]=cols[1][0:3]
# 		else:
# 			if dict[int(cols[0])]!=cols[1][0:3]:
# 				raise InconsistentInput("residue names within the upl file are inconsistent -- maybe to files with different sequence offset have been merged ? ")

# 		if int(cols[3]) not in dict:
# 			dict[int(cols[3])]=cols[4][0:3]
# 		else:
# 			if dict[int(cols[3])]!=cols[4][0:3]:
# 				raise InconsistentInput("residue names within the upl file are inconsistent -- maybe to files with different sequence offset have been merged ? ")



# 	#residue numbers as "keys"
# 	keys=dict.keys()
# 	keys.sort()
# 	print "highest residue number", keys[-1] #last residue
# 	max_upl_res=keys[-1]

# 	#create fasta sequence: use - for positions where aa3 is not known
# 	upl_fasta=list("-"*keys[-1])
# 	for key in keys:
# 		try:
# 			upl_fasta[key-1]=amino_acids.longer_names[ dict[key] ]
# 		except:
# 			print "WARNING: could not understand residue name %s"%dict[key]
# 	upl_fasta="".join(upl_fasta)
# 	print "TARGET FASTA: %s\n   UPL FASTA: %s\n"%(fasta, upl_fasta)

# 	#find offset
# 	start_upl=keys[0]-1
# 	inconsistent=False

# 	#find first usable pieces of sequence
# 	last_gap=0
# 	first_aa=0
# 	matchable_pieces={}
# 	for res_upl in range( 0, len(upl_fasta) ):
# 		#if this is a gap --> then is the end of the next matchable_piece
# 		if upl_fasta[res_upl]=='-':
# 			if last_gap<first_aa and res_upl-first_aa>4:
# 				matchable_pieces[upl_fasta[first_aa:res_upl]]=first_aa
# 			last_gap=res_upl
# 		elif last_gap == res_upl-1:
# 			first_aa=res_upl

# 	#start with longer pieces first to find match for offset
# 	sorted_pieces=sorted(matchable_pieces,key=len,reverse=True)

# 	#go through matchable pieces and see if it can be found in fasta
# 	inconsistent=True
# 	sequence_offset=0
# 	for pi in sorted_pieces:
# 		if pi in fasta:
# 			sequence_offset=matchable_pieces[pi]-fasta.index( pi )
# 			inconsistent=False
# 			break

# 	if len(sorted_pieces):
# 		print pi, " is in ", fasta
# 		print "with offset ", sequence_offset
# 		#now go from back to fron along the original fasta sequence to look for a matching piece
# 		if not inconsistent:
# 			print "found offset: %d\n"%sequence_offset
# 		else:
# 			print "cannot determine offset\n"
# 			raise InconsistentInput("upl fasta is inconsistent with the target fasta sequence")
# 	else:
# 		print "not sufficient sequence information in upl file -- assuming no offset"

# 	start_upl=sequence_offset-1
# 	#final-check:
# 	#all residue positions known from upl-fasta against the given fasta sequence using our offset
# 	inconsistent=False
# 	for ct in range( 0, min(len(fasta), max_upl_res-sequence_offset) ):
# 		start_upl=start_upl+1
# 		if upl_fasta[start_upl]=='-':
# 			continue
# 		if not fasta[ct]==upl_fasta[start_upl]:
# 			inconsistent=True
# 			print "found inconsistent aa %c-%c at position %d"%(fasta[ct],upl_fasta[start_upl],ct+1)


# 	if inconsistent:
# 		raise InconsistentInput("upl fasta is inconsistent with the target fasta sequence")

# 	return sequence_offset

def upl2mini( file, output_fn, fasta=None, sequence_offset=0, QFall_file=None, sep=4, pad=0.15, verbose=1 ):
	lines=open(file).readlines()
	length=200000
	from fasta import upl2fasta, find_fasta_offset
	if fasta:
#		sequence_offset=upl2mini_find_offset( fasta, lines )
		upl_fasta=upl2fasta( lines )
		sequence_offset=find_fasta_offset( upl_fasta, fasta, verbose=verbose )
		length=len(fasta)

	QFall=[]
	noQF=[]
	for line in lines:
		if "#QF" in line:
			QFall.append(line)
		else:
			noQF.append(line)

	if len(noQF) > 0 or not QFall_file:
		if isinstance(output_fn, str):
			output=open(output_fn,'w')
		else:
			output=output_fn
		for line in noQF:
			output.write(_upl2mini_line( line, sequence_offset, length, sep, pad  ))


   if len(QFall) > 0 and QFall_file:
		output=open(QFall_file,'w')
		for line in QFall:
			output.write(_upl2mini_line( line, sequence_offset, length, sep, pad ))

	return len(noQF),len(QFall)

