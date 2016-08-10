#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments
from __future__ import division
from __future__ import generators
import sys
sys.path.append('/Applications/PyMOLX11Hybrid.app/pymol/modules/')
from Tkinter import *
#from pymol import *
from pymol import cmd
from tkFileDialog import *
import os,math,re
import tkMessageBox
import Pmw
#from pymol import cmd,selector
#from pymol.cmd import _feedback,fb_module,fb_mask,is_list,_cmd
#from pymol.cgo import *
#from pymol import stored
import tkColorChooser
#from pymol.vfont import plain
sys.path.append('/Users/zak/csrosetta3/pymol')
from Noesy_Data_Collection import Noesy_Data_Collection
from arrow import *
pdbfile=['']
def __init__(self):
        self.menuBar.addcascademenu('Plugin', 'ZakPlugin',
                                 'Launch CrossPeakViewer',
                                 label='CrossPeakViewer')
				self.menuBar.addmenuitem('ZakPlugin','command', 'Load NOESY file',
																 label='load NOESY file',
																 command=lambda s=self : Load_NOE(s) )
				self.menuBar.addmenuitem('ZakPlugin','command','Load pdb file',
																 label='load pdb file',
																 command=lambda s=self : Load_pdb(s) )
				self.menuBar.addmenuitem('ZakPlugin','command','line between atom for fun',
																 label='line between atom for fun',
																 command=lambda s=self : Select_atoms(s) )
class Load_NOE:
	def __init__(self,app):
		self._app=app
		self._output=self.input_file()
	def return_app(self):
		return self._app
	def input_file(self):
		import tkMessageBox
		import tkFileDialog
		import os
		import string
		NOE_file=tkFileDialog.askopenfilename(parent=self.return_app().root, title='Open the NOE file')
		print NOE_file
		noe_class=Noesy_Data_Collection.read_from_file(NOE_file)
		print noe_class
		print noe_class._peaklists[2]
		print noe_class._peaklists[1]
		print noe_class._peaklists[2]._crosspeaks[2]
		print noe_class._peaklists[2]._crosspeaks[0]
		print noe_class._peaklists[2]._crosspeaks[0]._assignments[0]
		print noe_class._peaklists[2]._crosspeaks[0]._assignments[3]
		return noe_class

class Load_pdb:
	def __init__(self,app):
		self._app=app
		self.input_file()
	def return_app(self):
		return self._app
	def input_file(self):
		import tkMessageBox
		import tkFileDialog
		import os
		import string
		global pdbfile
		pdbfile=tkFileDialog.askopenfilename(parent=self.return_app().root, title='Open the pdb file')
		#pdbfile='/Users/zak/1cwa.pdb'
		print pdbfile
		#pdb_input=open(pdb_file)
# 		wd=os.path.dirname(pdbfile)
# 		stuff = pdbfile.split('/')
# 		id = stuff[len(stuff)-1].split('.')
# 		pdbfile = wd + os.sep + id[0] + '.pdb'
# 		print pdbfile
# 		pocfile = wd + os.sep + id[0] + '.poc'
# 		pocInfofile = wd + os.sep + id[0] + '.pocInfo'
		cmd.load(pdbfile)
# 		pocNums={}
# 		pocDict={}
# 		pocin=open(pocInfofile, "r")
# 		print pocin

#Load_pdb(1)

class Select_atoms:
	def __init__(self,app):
		self._app=app
		self.draw_line()
	def return_app(self):
		return self._app
	def draw_line(self):
		import tkSimpleDialog
		import tkMessageBox
		import urllib
		import os
		import string
		from arrow import cgo_arrow
		atoms = tkSimpleDialog.askstring('select two atoms','please input two atoms for example "1,CA,2,CB" means CA in residue 1 and CB in residue 2 ',parent=self.return_app().root)
		terminals=atoms.split(',')
		print terminals
		pdb_input=open(pdbfile)
		pdblist=pdb_input.readlines()
		coord1=[]
		coord2=[]
		for r in pdblist:
			if len(r)>5:
				if r.split()[5]==terminals[0] and r.split()[2]==terminals[1]:
					coord1.append(float(r.split()[6]))
					coord1.append(float(r.split()[7]))
					coord1.append(float(r.split()[8]))
				if r.split()[5]==terminals[2] and r.split()[2]==terminals[3]:
					coord2.append(float(r.split()[6]))
					coord2.append(float(r.split()[7]))
					coord2.append(float(r.split()[8]))
		print coord1
		print coord2
		cgo_arrow(coord1,coord2)
		#cmd.do("select atom1, resi "+terminals[0]+" and name "+terminals[1])
		#cmd.do("select atom2, resi "+terminals[2]+" and name "+terminals[3])
		#coord1=cmd.iterate_state(1, atom1, "stored.sel2.append([x,y,z])")
		#print coord1
		#coord2=cmd.iterate_state(1, atom2, "stored.sel2.append([x,y,z])")
		#print coord2
