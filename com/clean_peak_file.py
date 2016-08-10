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

parser = ExampleArgumentParser(prog=basename(__file__), description="make autoNOE-Rosetta readable peak list from any column based format",
examples=['%(prog)s input.txt -cols 1 3 4 6 7 -names id N h H V -tol 0.3 0.04 0.03 > proper.peaks',
					'%(prog)s input.txt proper.peaks'])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("input", help="A bmrb file or the chemical shift section of an bmrb file");
parser.add_argument("outfile", metavar="tab", help="chemical shift file",nargs="?",default="stdout");
parser.add_argument("-cols", help="which columns should be used from the input data [default: 1 2 3 4 7]?", type=int,nargs='*', default=[1,2,3,4,7] );
parser.add_argument("-names", help="what do the columns mean?", choices={'id','N','n','h','H','C','c','V','I'}, nargs='*', required=True );
#negative peaks are made positive when reading a peak-file.... ,'negV','negI','absI'}, nargs='*', required=True );
parser.add_argument("-tol", help="tolerances in sequence of appearance of resonances", nargs="*", default=None, type=float);
parser.add_argument("-skip", help='skip x lines as header', default=0 , type=int );
#parser.add_argument("-header", help="write a header into the file", action='store_true', default=False )
library.add_standard_args( parser )

args = parser.parse_args()

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
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
	def tolerances( self ):
		return ''
	def __lt__(self,other):
		return self._order < other._order

class PeakID(Interpreter):
	def __init__(self, id_col=None, name='id', reslist=[], tol_list=None ):
		Interpreter.__init__(self, 1 )
		self._col=id_col
		self._nr=0
	def __call__(self,tags):
		if self._col!=None:
			return '%8d'%int(tags[self._col])
		self._nr+=1
		return '%8d'%self._nr

class Constant(Interpreter):
	def __init__(self, astr, target_pos ):
		Interpreter.__init__(self,target_pos)
		self._str = astr
	def __call__(self, tags ):
		return self._str

class Resonance(Interpreter):

	def __init__(self, col, elem_letter, reslist=[], tol_list=None ):
		self._col=col
		self._elem_letter=elem_letter
		reslist.append(self)
		self._reso_id=len(reslist)
		Interpreter.__init__(self, 1+self._reso_id )
		if tol_list:
			self._tol=tol_list[self._reso_id-1]
		else:
			if self._elem_letter in 'NnCc':
				self._tol=0.3
			elif self._elem_letter=='H':
				self._tol=0.03
			else:
				self._tol=0.04

	def __call__(self, tags ):
		return '%8.3f'%float(tags[self._col])

	def header(self ):
		return '#INAME %d %s\n'%(self._reso_id,self._elem_letter)

	def tolerance( self ):
		return self._tol

	def name(self):
		return self._elem_letter

class Intensity(Interpreter):
	def __init__(self, col, name, reslist=[], tol_list=None ):
		self._col = col
		self._f=lambda x: x
		if 'neg' in name:
			self._f=lambda x: -x
		elif 'abs' in name:
			import math
			self._f=math.fabs
		Interpreter.__init__(self, 10 )

	def __call__(self, tags ):
		return '%8.3e'%(self._f(float(tags[self._col])))


interpreter_factory={'id':PeakID,'N':Resonance, 'C':Resonance, 'h':Resonance, 'H':Resonance,'c':Resonance, 'n':Resonance,
										 'I':Intensity, 'V':Intensity }
#,'negI':Intensity, 'negV':Intensity,'absI':Intensity }

try:
	interpreters=[]
	nr_res=0
	resonance_list=[]

#setup the interpreter
	for name,col in zip(args.names,args.cols):
		interpreters.append(interpreter_factory[name](col-1,name,resonance_list,args.tol))

	if not 'id' in args.names:
		interpreters.append(PeakID())
	if len(set(['I','V','negI','negV']).intersection(set(args.names)))==0:
		if verbose: print '\n','*'*80,'\n   WARNING: no intensity column specified... adding constant intensity of 1\n','*'*80,'\n'
		interpreters.append(Constant('1.00e+00',10))
	interpreters.append(Constant('0 U     ',8))
	interpreters.append(Constant('0.00e+00',11))

	interpreters.sort()

#generate the header
	header='# Number of dimensions %d\n'%len(resonance_list)
	for interpreter in interpreters:
		header+=interpreter.header()
	header+='#CYANAFORMAT '
	for res in resonance_list:
		header+='%c'%res.name()
	header+='\n'
	header+='#TOLERANCE '
	for res in resonance_list:
		header+=' %4.2f'%res.tolerance()
	header+='\n'
	if verbose: print header[:-1]
	outfile.write(header)


#main loop -- reading input file and interpreting each line
	skip=args.skip
	for line in open(args.input,'r'):
		if skip>0:
			if verbose: print 'ignore header line due to -skip settings: %s'%line[:-1]
			skip-=1
			continue
		tags=line.split()
		if len(tags)==0:
			continue
		if tags[0][0]=='#':
			if verbose: print 'ignore comment line: %s'%line[:-1]
			continue
		out_line=''
		for i in interpreters:
			out_line+=' '+i(tags)
		outfile.write(out_line+'\n')


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)



