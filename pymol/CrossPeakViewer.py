#!/usr/bin/env python
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
#================================================================
from __future__ import division
from __future__ import generators
import string
import os,math,re
import Tkinter
from Tkinter import *
import tkMessageBox
import Pmw
from pymol import cmd,selector
import sys
from pymol.cmd import _feedback,fb_module,fb_mask,is_list,_cmd
from pymol.cgo import *
from pymol import stored
import tkColorChooser
from pymol.vfont import plain
sys.path.append('/Users/zak/csrosetta3/pymol')
sys.path.append('/Users/zak/csrosetta3/python')
from Noesy_Data_Collection import Noesy_Data_Collection
from Model import Model
from noe_tools import parse_NMR_name
from amino_acids import longer_names

NOE_class=[]

fastaid=[]

#================================================================

# try:
#     import pymol
#     REAL_PYMOL = True
# except ImportError:
#     REAL_PYMOL = False
#     class pymol:
#         class cmd:
#             def load(self,name,sel=''):
#                 pass
#             def get_names(self):
#                 return ['mol1','mol2','map1','map2']
#             def get_type(self,thing):
#                 if thing.startswith('mol'):
#                     return 'object:molecule'
#                 else:
#                     return 'object:map'
#                 f.close()
#         cmd = cmd()
#     pymol = pymol()

#================================================================

def __init__(self):
        self.menuBar.addmenuitem('Plugin', 'command',
                                 'Launch CrossPeakViewer',
                                 label='CrossPeakViewer',
                                 command = lambda s=self: CrossPeakViewer(s))


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
                fd.title('Please IHATE YOU BAD USER choose a file')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)

        return FileDialogButton
    get = staticmethod(get)

#================================================================
#================================================================


class CrossPeakViewer:

    def __init__(self,app):

        parent = app.root
        self.parent = parent

        self.main_dialog = Pmw.Dialog(parent,
                                 buttons = ('Exit',),
                                 title = 'CrossPeakViewer Plugin',
                                 command = self.buttonPressed)
        self.main_dialog.withdraw()
        Pmw.setbusycursorattributes(self.main_dialog.component('hull'))

        # the title
        main_label = Tkinter.Label(self.main_dialog.interior(),
                                text = 'CrossPeakViewer Plugin\nOliver Lange & Zaiyong Zhang\n<http://langelab.ch.tum.de>',
                                background = 'navy',
                                foreground = 'white',
                                #pady = 20,
                                )
        main_label.pack(expand = 0, fill = 'both', padx = 4, pady = 4)

        # the basic notebook

        self.main_notebook = Pmw.NoteBook(self.main_dialog.interior())
        self.main_notebook.pack(fill='both',expand=1,padx=3,pady=3)

        # Files Card

        main_page = self.main_notebook.add('Files')


        # path and file names

        main_group1 = Pmw.Group(main_page,tag_text='rosetta HALLO ZAK autoNOE output files')
        main_group1.pack(fill = 'both', expand = 0, padx = 10, pady = 5)

        self.peak_file = StringVar()
        self.peak_file.set("NOE_out.dat");
        self.peak_location = Pmw.EntryField(main_group1.interior(),
                                            labelpos='w',
                                            label_pyclass = FileDialogButtonClassFactory.get(self.set_peak_filename,filter="*.dat"),
                                            validate = {'validator':self.quickFileValidation,},
                                            value = "/Users/zak/NOE_out.dat",
                                            label_text = 'Browse:')

        self.peak_location.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.load_buttonbox = Pmw.ButtonBox(main_group1.interior(), padx=0)
        self.load_buttonbox.pack(side=LEFT,expand = 1, padx = 10, pady = 5)
        self.load_buttonbox.add('Load NOE',command=self.load_peakfile)

        #=============================================================

        # path and file names

        main_group2 = Pmw.Group(main_page,tag_text='pdb files')
        main_group2.pack(fill = 'both', expand = 0, padx = 10, pady = 5)

        self.pdb_file = StringVar()
        self.pdb_file.set("native.pdb");
        self.pdb_location = Pmw.EntryField(main_group2.interior(),
                                           labelpos='w',
                                           label_pyclass = FileDialogButtonClassFactory.get(self.set_pdb_filename,filter="*.pdb"),
                                           validate = {'validator':self.quickFileValidation,},
                                           value = "/Users/zak/native.pdb",
                                           label_text = 'Browse:')

        self.pdb_location.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.load_buttonbox2 = Pmw.ButtonBox(main_group2.interior(), padx=0)
        self.load_buttonbox2.pack(side=LEFT,expand = 1, padx = 10, pady = 5)
        self.load_buttonbox2.add('Load pdb',command=self.load_pdbfile)



        #=============================================================



        self.pages = {}
        self.on_display = {}
        self.on_zoom = False
        self.myview = False
        self.colorDic = {}
        self.solv_vals = {}


        self.main_notebook.setnaturalsize()
        self.showAppModal()
    def set_peak_filename(self,filename):
	    self.peak_location.setvalue(filename)

    def set_pdb_filename(self,filename):
	    self.pdb_location.setvalue(filename)
    def set_atom(self,atomname):
	    self.atom_location.setvalue(atomname)
    def buttonPressed(self,result):

	    if hasattr(result,'keycode'):
		    if result.keycode == 36:
			    if self.main_notebook.getcurselection()=='Files':
				    self.load_files()
	    elif result == 'Exit' or result == None:
		    self.main_dialog.withdraw()

    def quickFileValidation(self,s):
	    if s == '': return Pmw.PARTIAL
	    elif os.path.isfile(s): return Pmw.OK
	    elif os.path.exists(s): return Pmw.PARTIAL
	    else: return Pmw.PARTIAL

    def load_peakfile(self):
	    peak_file = self.peak_location.get()
	    if ( peak_file ):

		    print "HALLO would have loaded file... "+peak_file
		    Load_NOE(peak_file)
		    peaklist_num=len(NOE_class._peaklists)
		    num=0
		    global peak_dialog

		    for r in NOE_class._peaklists:
			    show_peak_dialog(r,num)
			    num=num+1

    def load_pdbfile(self):
        pdb_file = self.pdb_location.get()
        if ( pdb_file ):

            print "HALLO would have loaded file... "+pdb_file
            Load_pdb(pdb_file)
				lines=open(pdb_file).readlines()
				count=0
				global fastaid
				removechain=0
				oldresnum = '   '
				for line in lines:
					if (len(line)>20): # and (chainid == line[21]):
						line_edit = line
						if line[0:3] == 'TER':
							break
						elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
							line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
							if (line_edit[12:14] == 'SE'):
								line_edit = line_edit[0:12]+' S'+line_edit[14:]
							if len(line_edit)>75:
								if (line_edit[76:78] == 'SE'):
									line_edit = line_edit[0:76]+' S'+line_edit[78:]
						if line_edit[0:4] == 'ATOM':
							resnum = line_edit[23:26]
							if not resnum == oldresnum:
								count = count + 1
								longname = line_edit[17:20]
								if longer_names.has_key(longname):
									fastaid.append( longer_names[longname] );
								else:
									fastaid.append( 'X')
							oldresnum = resnum
							newnum = '%3d' % count
							line_edit = line_edit[0:23] + newnum + line_edit[26:]
							if removechain:
								line_edit = line_edit[0:21]+' '+line_edit[22:]
            #print fastaid
				pass
    def atom_arrow(self):
	    atoms = self.atom_location.get()
	    if (atoms ):
		    print "give a arrow"
		    Select_atoms(atoms)

    def showAppModal(self):
	    self.main_dialog.show()


class show_peak_dialog:
	def __init__(self,peaklist,num):
		self.crosspeaks_list_select=[]
		self.assignments_list_select=[]
		self.crosspeaks_item=[]
		self.assignments_item=[]
		self._num=num
		self._peaklist=peaklist
		self.parent=[]

		for r in self._peaklist._crosspeaks:
			self.crosspeaks_item.append(str(r._id)+' '+r._volume)
		self.peak_dialog=Pmw.Dialog(self.parent,
															#buttons=('list assign','show assign',),
															#defaultbutton='list assign',
															title=self._peaklist._nfn
															#command=assignments(),
															)
		self.peak_dialog.withdraw()
		Pmw.setbusycursorattributes(self.peak_dialog.component('hull'))
		self.notebook_peaks=[]
		self.notebook_peaks = Pmw.NoteBook(self.peak_dialog.interior())
		self.notebook_peaks.pack(fill='both',expand=1,padx=3,pady=3)
		page_peaks=self.notebook_peaks.add(str(self._num))
		group1=Pmw.Group(page_peaks,tag_text="crosspeaks")
		group1.pack(fill='both',expand=1,padx=11,pady=6)
		self.crosspeaks_list=Pmw.ScrolledListBox(group1.interior(),
															  items=self.crosspeaks_item,
															  labelpos='w',
															  listbox_height=7,
															  #dblclickcommand=list_assignments(peaklist, assignments_item),
															  usehullsize = 2,
															  hull_width=101,
															  hull_height=201
															  )
		self.crosspeaks_list.pack(fill='both',expand=2,padx=6, pady=6)
		self.crosspeaks_list_select=(self.crosspeaks_list.getcurselection())
		group2=Pmw.Group(page_peaks,tag_text="assignments")
		group2.pack(fill='both',expand=1,padx=11,pady=6)
		self.assignments_list=Pmw.ScrolledListBox(group2.interior(),
																items=[],
																labelpos='w',
																listbox_height=7,
																usehullsize=2,
																hull_width=101,
																hull_height=201
																)
		self.assignments_list.pack(fill='both',expand=2,padx=6,pady=6)
		self.assign_buttonbox = Pmw.ButtonBox(group1.interior(), padx=0)
		self.assign_buttonbox.pack(side=LEFT,expand = 1, padx = 10, pady = 5)
		self.assign_buttonbox.add('list assign',command=self.list_assign)
		self.assign_buttonbox2 = Pmw.ButtonBox(group2.interior(), padx=0)
		self.assign_buttonbox2.pack(side=LEFT,expand = 1, padx = 10, pady = 5)
		self.assign_buttonbox2.add('show assign',command=self.show_assign)
		self.pages_peaks = {}
		self.notebook_peaks.setnaturalsize()
		self.peak_dialog.show()

	def list_assign(self):
		print "blabla list assign"
		assignments=[]
		print self.crosspeaks_list.getcurselection()[0]
		for r in self._peaklist._crosspeaks:
			if str(r._id)==(self.crosspeaks_list.getcurselection()[0]).split()[0]:
				for i in r._assignments:
					#print i._eliminated
					#if i._eliminated==0:
					assignments.append(i._atom_string)
		if assignments==[]:
			print "this crosspeak is not assigned"
		else:
			print assignments
			self.assignments_list.setlist(assignments)
					 #break

	def show_assign(self):
		print "blabla show assign"
		select=self.assignments_list.getcurselection()[0]
		print select
		tag1=select[1:-1]
		tag2=tag1.split()
		print tag2
		#print len(select)
		dist_name=select
		name1=parse_NMR_name(tag2[6],fastaid[int(tag2[5])-1],int(tag2[5]))[0]
		selection1='resid '+tag2[5]+' and '+'name '+name1
		print selection1
		name2=parse_NMR_name(tag2[10],fastaid[int(tag2[9])-1],int(tag2[9]))[0]
		selection2='resid '+tag2[9]+' and '+'name '+name2
		print selection2
		cmd.distance("test",selection1,selection2)


class Load_NOE:
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
        global NOE_class
        NOE_class=Noesy_Data_Collection.read_from_file(self.return_app())
        print NOE_class._nafn
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
        pdbfile=self.return_app()
        print pdbfile
        cmd.load(pdbfile)

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

        atoms=self.return_app()
        terminals=atoms.split(',')
        print "hello terminals"
        print terminals
        dist_name=terminals[0]+'_'+terminals[1]+'_'+terminals[2]+'_'+terminals[3]
        selection1='resid '+terminals[0]+' and '+'name '+terminals[1]
        selection2='resid '+terminals[2]+' and '+'name '+terminals[3]
        cmd.distance(dist_name,selection1,selection2)

################################################################################
import os,fnmatch,time
import Tkinter,Pmw
#Pmw.setversion("0.8.5")

def _errorpop(master,text):
    d=Pmw.MessageDialog(master,
                        title="Error",
                        message_text=text,
                        buttons=("OK",))
    d.component('message').pack(ipadx=15,ipady=15)
    d.activate()
    d.destroy()

class PmwFileDialog(Pmw.Dialog):
	"""File Dialog using Pmw"""
	def __init__(self, parent = None, **kw):
		# Define the megawidget options.
		optiondefs = (
			('filter',    '*',              self.newfilter),
			('directory', os.getcwd(),      self.newdir),
			('filename',  '',               self.newfilename),
			('historylen',10,               None),
			('command',   None,             None),
			('info',      None,             None),
			)
		self.defineoptions(kw, optiondefs)
        # Initialise base class (after defining options).
		Pmw.Dialog.__init__(self, parent)

		self.withdraw()

	# Create the components.
		interior = self.interior()

		if self['info'] is not None:
			rowoffset=1
			dn = self.infotxt()
			dn.grid(row=0,column=0,columnspan=2,padx=3,pady=3)
		else:
			rowoffset=0

		dn = self.mkdn()
		dn.grid(row=0+rowoffset,column=0,columnspan=2,padx=3,pady=3)
		del dn

		# Create the directory list component.
		dnb = self.mkdnb()
		dnb.grid(row=1+rowoffset,column=0,sticky='news',padx=3,pady=3)
		del dnb

	# Create the filename list component.
		fnb = self.mkfnb()
		fnb.grid(row=1+rowoffset,column=1,sticky='news',padx=3,pady=3)
		del fnb

		# Create the filter entry
		ft = self.mkft()
		ft.grid(row=2+rowoffset,column=0,columnspan=2,padx=3,pady=3)
		del ft

		# Create the filename entry
		fn = self.mkfn()
		fn.grid(row=3+rowoffset,column=0,columnspan=2,padx=3,pady=3)
		fn.bind('<Return>',self.okbutton)
		del fn

		# Buttonbox already exists
		bb=self.component('buttonbox')
		bb.add('OK',command=self.okbutton)
		bb.add('Cancel',command=self.cancelbutton)
		del bb

		Pmw.alignlabels([self.component('filename'),
							  self.component('filter'),
							  self.component('dirname')])

	def infotxt(self):
		""" Make information block component at the top """
		return self.createcomponent(
			'infobox',
			(), None,
			Tkinter.Label, (self.interior(),),
			width=51,
			relief='groove',
			foreground='darkblue',
			justify='left',
			text=self['info']
			)

	def mkdn(self):
		"""Make directory name component"""
		return self.createcomponent(
			'dirname',
			(), None,
			Pmw.ComboBox, (self.interior(),),
			entryfield_value=self['directory'],
			entryfield_entry_width=40,
			entryfield_validate=self.dirvalidate,
			selectioncommand=self.setdir,
			labelpos='w',
			label_text='Directory:')

	def mkdnb(self):
		 """Make directory name box"""
		 return self.createcomponent(
			 'dirnamebox',
			 (), None,
			 Pmw.ScrolledListBox, (self.interior(),),
			 label_text='directories',
			 labelpos='n',
			 hscrollmode='none',
			 dblclickcommand=self.selectdir)

	def mkft(self):
		"""Make filter"""
		return self.createcomponent(
			'filter',
			(), None,
			Pmw.ComboBox, (self.interior(),),
			entryfield_value=self['filter'],
			entryfield_entry_width=40,
			selectioncommand=self.setfilter,
			labelpos='w',
			label_text='Filter:')

	def mkfnb(self):
		"""Make filename list box"""
		return self.createcomponent(
			'filenamebox',
			(), None,
			Pmw.ScrolledListBox, (self.interior(),),
			label_text='files',
			labelpos='n',
			hscrollmode='none',
			selectioncommand=self.singleselectfile,
			dblclickcommand=self.selectfile)

	def mkfn(self):
		"""Make file name entry"""
		return self.createcomponent(
			'filename',
			(), None,
			Pmw.ComboBox, (self.interior(),),
			entryfield_value=self['filename'],
			entryfield_entry_width=40,
			entryfield_validate=self.filevalidate,
			selectioncommand=self.setfilename,
			labelpos='w',
			label_text='Filename:')

	def dirvalidate(self,string):
		if os.path.isdir(string):
			return Pmw.OK
		else:
			return Pmw.PARTIAL

	def filevalidate(self,string):
		if string=='':
			return Pmw.PARTIAL
		elif os.path.isfile(string):
			return Pmw.OK
		elif os.path.exists(string):
			return Pmw.PARTIAL
		else:
			return Pmw.OK

	def okbutton(self):
			"""OK action: user thinks he has input valid data and wants to
			proceed. This is also called by <Return> in the filename entry"""
			fn=self.component('filename').get()
			self.setfilename(fn)
			if self.validate(fn):
				self.canceled=0
				self.deactivate()

	def cancelbutton(self):
		"""Cancel the operation"""
		self.canceled=1
		self.deactivate()

	def tidy(self,w,v):
		"""Insert text v into the entry and at the top of the list of
		the combobox w, remove duplicates"""
		if not v:
			return
		entry=w.component('entry')
		entry.delete(0,'end')
		entry.insert(0,v)
		list=w.component('scrolledlist')
		list.insert(0,v)
		index=1
		while index<list.index('end'):
			k=list.get(index)
			if k==v or index>self['historylen']:
				list.delete(index)
			else:
				index=index+1
				w.checkentry()

	def setfilename(self,value):
		if not value:
			return
		value=os.path.join(self['directory'],value)
		dir,fil=os.path.split(value)
		self.configure(directory=dir,filename=value)

		c=self['command']
		if callable(c):
			c()

	def newfilename(self):
		"""Make sure a newly set filename makes it into the combobox list"""
		self.tidy(self.component('filename'),self['filename'])

	def setfilter(self,value):
		self.configure(filter=value)

	def newfilter(self):
		"""Make sure a newly set filter makes it into the combobox list"""
		self.tidy(self.component('filter'),self['filter'])
		self.fillit()

	def setdir(self,value):
		self.configure(directory=value)

	def newdir(self):
		"""Make sure a newly set dirname makes it into the combobox list"""
		self.tidy(self.component('dirname'),self['directory'])
		self.fillit()

	def singleselectfile(self):
		"""Single click in file listbox. Move file to "filename" combobox"""
		cs=self.component('filenamebox').curselection()
		if cs!=():
			value=self.component('filenamebox').get(cs)
			self.setfilename(value)

	def selectfile(self):
		"""Take the selected file from the filename, normalize it, and OK"""
		self.singleselectfile()
		value=self.component('filename').get()
		self.setfilename(value)
		if value:
	    self.okbutton()

	def selectdir(self):
		"""Take selected directory from the dirnamebox into the dirname"""
		cs=self.component('dirnamebox').curselection()
		if cs!=():
			value=self.component('dirnamebox').get(cs)
			dir=self['directory']
		if not dir:
				dir=os.getcwd()
		if value:
			if value=='..':
				dir=os.path.split(dir)[0]
			else:
				dir=os.path.join(dir,value)
		self.configure(directory=dir)
		self.fillit()

	def askfilename(self,directory=None,filter=None):
		"""The actual client function. Activates the dialog, and
	   returns only after a valid filename has been entered
		(return value is that filename) or when canceled (return
		value is None)"""
		if directory!=None:
			self.configure(directory=directory)
		if filter!=None:
			self.configure(filter=filter)
		self.fillit()
		self.canceled=1 # Needed for when user kills dialog window
		self.activate()
		if self.canceled:
			return None
		else:
			return self.component('filename').get()

	lastdir=""
	lastfilter=None
	lasttime=0
	def fillit(self):
		"""Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
		if self.lastdir==self['directory'] and self.lastfilter==self['filter'] and self.lasttime>os.stat(self.lastdir)[8]:
			return
		self.lastdir=self['directory']
		self.lastfilter=self['filter']
		self.lasttime=time.time()
		dir=self['directory']
		if not dir:
			dir=os.getcwd()
		dirs=['..']
		files=[]
		try:
			fl=os.listdir(dir)
			fl.sort()
		except os.error,arg:
			if arg[0] in (2,20):
				return
			raise
		for f in fl:
			if os.path.isdir(os.path.join(dir,f)):
				dirs.append(f)
			else:
				filter=self['filter']
				if not filter:
					filter='*'
				if fnmatch.fnmatch(f,filter):
					files.append(f)
		self.component('filenamebox').setlist(files)
		self.component('dirnamebox').setlist(dirs)

	def validate(self,filename):
		"""Validation function. Should return 1 if the filename is valid,
		0 if invalid. May pop up dialogs to tell user why. Especially
		suited to subclasses: i.e. only return 1 if the file does/doesn't
		exist"""
		return 1

class PmwExistingFileDialog(PmwFileDialog):
    def filevalidate(self,string):
        if os.path.isfile(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL

    def validate(self,filename):
        if os.path.isfile(filename):
            return 1
        elif os.path.exists(filename):
            _errorpop(self.interior(),"This is not a plain file")
            return 0
        else:
            _errorpop(self.interior(),"Please select an existing file")
            return 0


#===============================================================
#
# stuff to deal with dlg files

def sortByEnergy(mlist):
    en = []
    idx = []
    for m in mlist:
        en.append(m.energy)
        idx.append(mlist.index(m))
    changed = True
    while changed:
        changed = False
        for i in range(len(en)-1):
            if en[i] > en[i+1]:
                dum = en[i+1]
                en[i+1] = en[i]
                en[i] = dum
                dum = idx[i+1]
                idx[i+1] = idx[i]
                idx[i] = dum
                changed = True
    new = []
    for i in range(len(idx)):
        new.append(mlist[idx[i]])
        new[i].num=i+1
    return new


# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    widget = Autodock(app)
    exitButton = Tkinter.Button(app.root, text = 'Exit', command = app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
