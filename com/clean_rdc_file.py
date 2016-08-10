#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-
from os.path import basename
import sys
import argparse
from basic.options import ExampleArgumentParser
import bmrb
### toolbox library
import library
from os.path import exists
import cs
import traceback
#############################

parser = ExampleArgumentParser(prog=basename(__file__), description="make Rosetta readable RDC list from any column based format",
examples=['%(prog)s input.txt -cols 1 2 3 4 5 6 7 -names resi1 resn1 atom1 resi2 resn2 atom2 rdc -out proper.rdc'])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("input", help="A bmrb file or the chemical shift section of an bmrb file");
parser.add_argument("-out", metavar="tab", help="chemical shift file",required=True);
parser.add_argument("-cols", help="which columns should be used from the input data ?", type=int,nargs='*', default=[1,2,3,4,7] );
parser.add_argument("-names", help="what do the columns mean?", choices={'expid','resi1','resn1','atom1','resi2','resn2','atom2','rdc'}, nargs='*', required=True );
parser.add_argument("-skip", help='skip x lines as header', default=0 , type=int );
parser.add_argument("-offset", help='an offset for the residue number', default=0, type=int );
#parser.add_argument("-header", help="write a header into the file", action='store_true', default=False )
library.add_standard_args( parser )

args = parser.parse_args()

#output:
verbose=1
library.hello(__file__)

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

class Interpreter:
	def __init__(self, order ):
		self._order=order
	def __call__(self,tags):
		return None
	def header(self ):
		return ''
	def expid(self):
		return 0
	def set_offset(self,offset):
		pass

	def __lt__(self,other):
		return self._order < other._order

class ExpID(Interpreter):
	def __init__(self, id_col=None, name='id' ):
		Interpreter.__init__(self, 1 )
		self._col=id_col
		self._nr=0
		self._last_exp_id = -100
	def __call__(self,tags):
		self._last_exp_id = int(tags[self._col])
		return ''
	def expid(self):
		return self._last_exp_id

class Constant(Interpreter):
	def __init__(self, astr, target_pos ):
		Interpreter.__init__(self,target_pos)
		self._str = astr
	def __call__(self, tags ):
		return self._str

class BadLine():
	pass

class ResID(Interpreter):
	def __init__(self, col, name ):
		if '1' in name:
			pos=1
		else:
			pos=10
		self._col = col
		Interpreter.__init__(self, 10 )

	def set_offset(self, offset):
		self._offset=offset;

	def __call__(self, tags ):
		id=int(tags[self._col])-self._offset
		if id<2:
			raise BadLine()
		return '%8d'%id

class ResName(Interpreter):
	def __init__(self, col, name ):
		self._col = col
		Interpreter.__init__(self, 10 )

	def __call__(self, tags ):
		return ''

class AtomName(Interpreter):
	def __init__(self, col, name ):
		if '1' in name:
			pos=2
		else:
			pos=11
		self._col = col
		Interpreter.__init__(self, 10 )
	def __call__(self, tags ):
		return tags[self._col]

class RDC(Interpreter):
	def __init__(self, col, name ):
		self._col = col
		Interpreter.__init__(self, 30 )

	def __call__(self, tags ):
		return '%8.3f'%float(tags[self._col])

interpreter_factory={'expid':ExpID,'resi1':ResID,'resi2':ResID,
										 'resn1':ResName,'resn2':ResName,
										 'atom1':AtomName,'atom2':AtomName,
										 'rdc':RDC }

try:
	interpreters=[]
	nr_res=0

#setup the interpreter
	for name,col in zip(args.names,args.cols):
		interpreters.append(interpreter_factory[name](col-1,name))

	interpreters.sort()
	for i in interpreters:
		i.set_offset( args.offset)
#generate the header
#main loop -- reading input file and interpreting each line
	skip=args.skip

	current_expid=0
	all_content={}

	for line in open(args.input,'r'):
		if skip>0:
			if verbose: print 'ignore header line due to -skip settings: %s'%line[:-1]
			skip-=1
			continue
		tags=line.split()
		if len(tags)==0: continue
		if tags[0][0]=='#':
			if verbose: print 'ignore comment line: %s'%line[:-1]
			continue

		out_line=''
		#run interpreters
		expid=0
		try:
			for i in interpreters:
				out_line+=' '+i(tags)
				expid+=i.expid()
			all_content.setdefault(expid,'')
			all_content[expid]+=out_line+'\n'
		except BadLine:
			pass

	for expid,content in all_content.iteritems():
		outfile=args.out
		if expid>0:
			import os
			base,ext=os.path.splitext(outfile)
			outfile='%s_%d%s'%(base,expid,ext)
		open(outfile,'w').write(content)




except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)



