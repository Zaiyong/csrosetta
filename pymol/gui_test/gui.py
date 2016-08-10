#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments
import sys
sys.path.append('/Applications/PyMOLX11Hybrid.app/pymol/modules/')
from Tkinter import *
from tkFileDialog import askopenfilename
import Pmw
import os,math,re
filename=''

class FileDialogButtonClassFactory:
    def get(fn,filter='*'):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class FileDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.

            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                apply(Tkinter.Button.__init__, (self, master, cnf), kw)
                self.configure(command=self.set)
            def set(self):
                fd = PmwFileDialog(self.master,filter=filter)
                fd.title('Please choose a file')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)
class button:
	def __init__(self,parent):
		self.parent=parent
		self.press_button()
	def press_button(self):
		f=Frame(self.parent)
		f.pack(padx=155,pady=55)
		#self.entry = Entry(f,text="enter your choice")
		#self.entry.pack(side= TOP,padx=10,pady=12)
		self.peak_location=Pmw.EntryField(self.parent,
																			labelpos='w',
																			value="NOE_out.dat",
																			label_pyclass = FileDialogButtonClassFactory.get(self.set_peak_filename(filename),filter="*.dat"),
																			validate = {'validator':self.quickFileValidation,},
																			label_text="Browse:")

		entries = (self.peak_location)
		self.peak_location.pack(side=TOP, padx=15, pady=15)
		#Pmw.alignlabels(entries)
		#self._any.component('entry').focus_set()
		self.button=Button(f,text="browse",command=self.browse)
		self.button.pack(side=BOTTOM,padx=10,pady=10)
	def browse(self):
		#global filenmame
		#self.entry=askopenfilename()
		self.peak_location.value=askopenfilename()
		print self.peak_location.value
		#print filename
	def quickFileValidation(self,s):
		if s == '': return Pmw.PARTIAL
		elif os.path.isfile(s): return Pmw.OK
		elif os.path.exists(s): return Pmw.PARTIAL
		else: return Pmw.PARTIAL
    def set_peak_filename(self,filename):
        self.peak_location.setvalue(filename)
root=Tk()
root.title('blablablabla')
button=button(root)
root.mainloop()
