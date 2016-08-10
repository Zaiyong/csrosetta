#!/usr/bin/env python
#================================================================
from __future__ import division
from __future__ import generators

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
from Numeric import *
import tkColorChooser
from pymol.vfont import plain
from cnc.cnc_parser import *
from cnc.interaction import *
from cnc.HBond import *
from cnc.CstParser import *
from cnc.ViolParser import *
#================================================================

# Python backward-compatibility...
try:
    True
except:
    True = 1
try:
    False
except:
    False = 0
#
# Cheap hack for testing purposes
#

try:
    import pymol
    REAL_PYMOL = True
except ImportError:
    REAL_PYMOL = False
    class pymol:
        class cmd:
            def load(self,name,sel=''):
                pass
            def get_names(self):
                return ['mol1','mol2','map1','map2']
            def get_type(self,thing):
                if thing.startswith('mol'):
                    return 'object:molecule'
                else:
                    return 'object:map'
                f.close()
        cmd = cmd()
    pymol = pymol()

#================================================================

def __init__(self):
        self.menuBar.addmenuitem('Plugin', 'command',
                                 'Launch tCONCOORD',
                                 label='tCONCOORD',
                                 command = lambda s=self: Concoord(s))




# set the defaults


file_defaults = {
    "hb":"hbonds.dat",
    "ct":"contab.dat",
    "ctp":"tdist.ctp",
    "cst":"tdist.dat",
    'dist':'~/gromacs/src/concoord/tdist',
    'disco':'~/gromacs/src/concoord/dtisco',
#    'wdir':'~/projects/plugin/new'
    'wdir':os.path.abspath('.')
    }



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


#================================================================
#================================================================


class Concoord:

    def __init__(self,app):

        parent = app.root
        self.parent = parent

        self.dialog = Pmw.Dialog(parent,
                                 buttons = ('Exit',),
                                 title = 'tConcoord Plugin',
                                 command = self.buttonPressed)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        # the title

        w = Tkinter.Label(self.dialog.interior(),
                                text = 'tCONCOORD Plugin\nDaniel Seeliger\n<http://www.mpibpc.mpg.de/groups/grubmueller/start/people/dseelig/index.html>',
                                background = 'navy',
                                foreground = 'white',
                                #pady = 20,
                                )
        w.pack(expand = 0, fill = 'both', padx = 4, pady = 4)

        # the basic notebook

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',expand=1,padx=3,pady=3)

        # Files Card

        page = self.notebook.add('Files')


        # path and file names

        group1 = Pmw.Group(page,tag_text='cdist output files')
        group1.pack(fill = 'both', expand = 0, padx = 10, pady = 5)


        self.hbfile = StringVar()
        self.ctfile = StringVar()
        self.ctpfile = StringVar()
        self.dist_exec = StringVar()
        self.disco_exec = StringVar()
        self.cstfile = StringVar()
        self.violfile = StringVar()
        self.wdir = StringVar()
        self.hbfile.set(file_defaults['hb'])
        self.ctfile.set(file_defaults['ct'])
        self.ctpfile.set(file_defaults['ctp'])
        self.cstfile.set(file_defaults['cst'])
        self.dist_exec.set(file_defaults['dist'])
        self.disco_exec.set(file_defaults['disco'])
        self.wdir.set(file_defaults['wdir'])
        self.violfile.set("S.viol");
        self.hblocation = Pmw.EntryField(group1.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_dlgfilename,filter="*.dat"),
                                         validate = {'validator':self.quickFileValidation,},
                                         value = file_defaults['hb'],
                                         label_text = 'Browse:')

        self.hblocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.cstlocation = Pmw.EntryField(group1.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_dlgfilename,filter="*.dat"),
                                         validate = {'validator':self.quickFileValidation,},
                                         value = file_defaults['cst'],
                                         label_text = 'Browse:')

        self.cstlocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.ctlocation = Pmw.EntryField(group1.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_dlgfilename,filter="*.dat"),
                                         validate = {'validator':self.quickFileValidation,},
                                         value = file_defaults['ct'],
                                         label_text = 'Browse:')

        self.ctlocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.ctplocation = Pmw.EntryField(group1.interior(),
                                          labelpos='w',
                                          label_pyclass = FileDialogButtonClassFactory.get(self.set_dlgfilename,filter="*.pdb"),
                                          validate = {'validator':self.quickFileValidation,},
                                          value = file_defaults['ctp'],
                                          label_text = 'Browse:')

        self.ctplocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.viollocation = Pmw.EntryField(group1.interior(),
                                          labelpos='w',
                                          label_pyclass = FileDialogButtonClassFactory.get(self.set_violfilename,filter="*.viol"),
                                          validate = {'validator':self.quickFileValidation,},
                                          value = "S.viol",
                                          label_text = 'Browse:')

        self.viollocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.load_buttonbox = Pmw.ButtonBox(group1.interior(), padx=0)
        self.load_buttonbox.pack(side=LEFT,expand = 1, padx = 10, pady = 5)
        self.load_buttonbox.add('Load Files',command=self.load_files)


        group2 = Pmw.Group(page,tag_text='work directory')
        group2.pack(fill = 'both', expand = 0, padx = 10, pady = 5)

        self.wdirlocation = Pmw.EntryField(group2.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_dlgfilename,filter="*.dat"),
                                         validate = {'validator':self.quickFileValidation,},
                                         value = file_defaults['wdir'],
                                         label_text = 'Browse:')

        self.wdirlocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)


        group3 = Pmw.Group(page,tag_text='tCONCOORD executables')
        group3.pack(fill = 'both', expand = 0, padx = 10, pady = 5)

        self.distlocation = Pmw.EntryField(group3.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_dlgfilename,filter="*.dat"),
                                         validate = {'validator':self.quickFileValidation,},
                                         value = file_defaults['dist'],
                                         label_text = 'Browse:')

        self.distlocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.discolocation = Pmw.EntryField(group3.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_dlgfilename,filter="*.dat"),
                                         validate = {'validator':self.quickFileValidation,},
                                         value = file_defaults['disco'],
                                         label_text = 'Browse:')

        self.discolocation.pack(fill = 'both', expand = 1, padx = 10, pady = 5)



        #=============================================================

        self.pages = {}
        self.on_display = {}
        self.on_zoom = False
        self.myview = False
        self.colorDic = {}
        self.solv_vals = {}

        page = self.notebook.add('Interactions')
        rgroup = Pmw.Group(page,tag_text='Interaction Table')


        rgroup.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.rgroup = Pmw.NoteBook(rgroup.interior())
        self.rgroup.pack(fill='both',expand=1,padx=3,pady=3)

        self.status_line = Label(rgroup.interior(), #relief='groove',
                                 relief='sunken',
                                 font='helvetica 12', anchor='w',fg='yellow',bg='black')
#        bg='#ddddff')
        self.status_line.pack(side='left', fill='x', expand=True)



        #=============================================================

        # the exclusion card

        page = self.notebook.add('Exclusions')
        self.exclusions = []
        group = Pmw.Group(page, tag_text='Exclusions')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        text = ''


        lfre=Frame(group.interior())
        bar=Scrollbar(lfre,)
        self.extext=Text(lfre,yscrollcommand=bar.set,background="white")
        bar.config(command=self.extext.yview)

        self.extext.insert(END,text)
        self.extext.pack(side=LEFT,expand="yes",fill="both")
        bar.pack(side=LEFT,expand="yes",fill="y")
        lfre.pack(expand="yes",fill="both")

        #=============================================================
        # the about card

        page = self.notebook.add('About')
        group = Pmw.Group(page, tag_text='About tCONCOORD Plugin')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        text = """
        This plugin help setting up tCONCOORD runs
        """

        lfre=Frame(group.interior())
        bar=Scrollbar(lfre,)
        ll=Text(lfre,yscrollcommand=bar.set,background="#ddddff",font="Times 14")
        bar.config(command=ll.yview)

        ll.insert(END,text)
        ll.pack(side=LEFT,expand="yes",fill="both")
        bar.pack(side=LEFT,expand="yes",fill="y")
        lfre.pack(expand="yes",fill="both")

        self.notebook.setnaturalsize()

        self.showAppModal()


    #===========================================================

    # functions

    def tk_color_dialog(self):
        color = tkColorChooser.Chooser(
            initialcolor='white',title='Choose color').show()
        if color[0] is not None:
            return [color[0][0]/100.,
                    color[0][1]/100.,
                    color[0][2]/100.]
        else: return None


    def show_meta_info(self, value):

        name = self.rgroup.getcurselection()
        objl = self.nameDic[name][value]
        self.pages[name]['text'].clear()
        s = ''
        for obj in objl:
            if obj.type in HBOND_TYPES:
                s += obj.hb_meta_str()
            elif obj.type=='PHO':
                s+=obj.hphob_meta_str()
            elif obj.type=='NET':
                s+=obj.net_meta_str()

        self.pages[name]['text'].insert('end',s)

    def set_outfilename(self,filename):
        self.outlocation.setvalue(filename)

    def set_dlgfilename(self,filename):
        self.dlglocation.setvalue(filename)

    def set_violfilename(self,filename):
        self.viollocation.setvalue(filename)

    def buttonPressed(self,result):
        print 'button', result.keycode
        if hasattr(result,'keycode'):
            if result.keycode == 36:
                if self.notebook.getcurselection()=='Files':
                    self.load_files()
        elif result == 'Exit' or result == None:
            self.dialog.withdraw()


    def quickFileValidation(self,s):
        if s == '': return Pmw.PARTIAL
        elif os.path.isfile(s): return Pmw.OK
        elif os.path.exists(s): return Pmw.PARTIAL
        else: return Pmw.PARTIAL


    def load_files(self):

        hbf = self.hblocation.get()
        ctpf = self.ctplocation.get()
        ctf = self.ctlocation.get()
        cstf = self.cstlocation.get()
        violf = self.viollocation.get()
        ctp = CTPParser().parse(ctpf)
        hb = HBondParser().parse(hbf)
        ct = ContabParser().parse(ctf)
        cst = CstParser().parse(cstf,ctp)
        if ( violf != "S.viol" ):
            viols = ViolParser().parse(violf,ctp)
        self.model = ctp.get_model()
        self.ObjDic = {}
        self.nameDic = {}
        self.pdbstr = IObj().build_ca_str(self.model)
        cmd.read_pdbstr(self.pdbstr,'CA-Trace')
        self.ObjDic['UNC'] = []
        self.ObjDic['PHO'] = []
        self.nameDic['UNC'] = {}
        self.nameDic['PHO'] = {}
        self.ObjDic['NET'] = []
        self.nameDic['NET'] = {}

        for hbond in hb.hbonds:
            obj = IObj().build_hbond_obj(hbond,ct,self.model)
            if not obj.constr:
                self.ObjDic['UNC'].append(obj)
                if self.nameDic['UNC'].has_key(obj.hash):
                   self.nameDic['UNC'][obj.hash].append(obj)
                else:
                    self.nameDic['UNC'].update({obj.hash:[obj]})
            else:
                if not self.ObjDic.has_key(obj.type):
                    self.ObjDic[obj.type] = [obj]
                else:
                    self.ObjDic[obj.type].append(obj)
                if not self.nameDic.has_key(obj.type):
                    self.nameDic[obj.type] = {obj.hash:[obj]}
                else:
                    if self.nameDic[obj.type].has_key(obj.hash):
                        self.nameDic[obj.type][obj.hash].append(obj)
                    else:
                        self.nameDic[obj.type].update({obj.hash:[obj]})


#        for r1, dic in ct.con.items():
#            for r2, tpl in dic.items():
#                if 'PHO' in tpl:
#                    obj = IObj().build_hphob_obj(r1,r2,self.model,'PHO')
#                    if not self.ObjDic.has_key(obj.type):
#                        self.ObjDic[obj.type] = [obj]
#                    else:
#                        self.ObjDic[obj.type].append(obj)
#                    if not self.nameDic.has_key(obj.type):
#                        self.nameDic[obj.type] = {obj.hash:[obj]}
#                    else:
#                        if self.nameDic[obj.type].has_key(obj.hash):
#                            self.nameDic[obj.type][obj.hash].append(obj)
#                        else:
#                            self.nameDic[obj.type].update({obj.hash:[obj]})


# read NETWORK items
        net = cst.data[ CST_NET ]
        for bond in net:
            obj = IObj().build_hphob_obj(bond.res1,bond.res2,self.model,'NET')
            if not self.ObjDic.has_key(obj.type):
                self.ObjDic[obj.type] = [obj]
            else:
                self.ObjDic[obj.type].append(obj)
                if not self.nameDic.has_key(obj.type):
                    self.nameDic[obj.type] = {obj.hash:[obj]}
                else:
                    if self.nameDic[obj.type].has_key(obj.hash):
                        self.nameDic[obj.type][obj.hash].append(obj)
                    else:
                        self.nameDic[obj.type].update({obj.hash:[obj]})
# read VIOLATIONS
        if ( violf != "S.viol" ):
            for bond in viols.data:
                print bond
                obj = IObj().build_hphob_obj(bond.res1,bond.res2,self.model,VIOL_TYPES[bond.type])
                if not self.ObjDic.has_key(obj.type):
                    self.ObjDic[obj.type] = [obj]
                else:
                    self.ObjDic[obj.type].append(obj)
                    if not self.nameDic.has_key(obj.type):
                        self.nameDic[obj.type] = {obj.hash:[obj]}
                    else:
                        if self.nameDic[obj.type].has_key(obj.hash):
                            self.nameDic[obj.type][obj.hash].append(obj)
                        else:
                            self.nameDic[obj.type].update({obj.hash:[obj]})

#read LR items
#         net = cst.data[ CST_LR ]
#         for bond in net:
#             obj = IObj().build_hphob_obj(bond.res1,bond.res2,self.model,'LR')
#             if not self.ObjDic.has_key(obj.type):
#                 self.ObjDic[obj.type] = [obj]
#             else:
#                 self.ObjDic[obj.type].append(obj)
#                 if not self.nameDic.has_key(obj.type):
#                     self.nameDic[obj.type] = {obj.hash:[obj]}
#                 else:
#                     if self.nameDic[obj.type].has_key(obj.hash):
#                         self.nameDic[obj.type][obj.hash].append(obj)
#                     else:
#                         self.nameDic[obj.type].update({obj.hash:[obj]})


        for key in INT_TYPES:
            if self.ObjDic.has_key(key):
                try:
                    self.colorDic[key]=self.ObjDic[key][0].AO.color
                except:
                    self.colorDic[key]=[1,0,0]
                self.update_combo(key)


    def Objs_as_list(self,name):
        l = []
        for obj in self.ObjDic[name]:
            if not obj.hash in l:
                l.append(obj.hash)
        return l

    def update_combo(self,name):
        try:
            self.rgroup.delete(name)
        except:
            pass
        list = self.Objs_as_list(name)
        self.on_screen=name
        self.on_display[name] = False

        self.pages[name] = {'name':self.rgroup.add(name)}
        self.pages[name].update({'structs':list})
        self.buttonbox = Pmw.ButtonBox(self.pages[name]['name'], padx=3)
        self.buttonbox.pack(fill='x',side=TOP)
        self.buttonbox.add('Show/Hide',command=self.arrow_display)
#        self.buttonbox.add('Hide Inter.',command=self.hide_arrows)
        self.buttonbox.add('Change Color',command=self.change_color)
        self.buttonbox.add('Remove Inter.',command=self.remove)
        self.buttonbox.add('Zoom in/out',command=self.zoom_inter)
        self.buttonbox.add('Apply Solv.',command=self.apply_solv)

        self.pages[name]['group1'] = Pmw.Group(self.pages[name]['name'],tag_text = "Interactions")
        self.pages[name]['group1'].pack(side=LEFT,expand=0)
        self.pages[name]['group2'] = Pmw.Group(self.pages[name]['name'],tag_text = "Information")
        self.pages[name]['group2'].pack(side=LEFT,expand=0)


        if name in ['BBH','SSH','SBH']:
            self.pages[name]['group3'] = Pmw.Group(self.pages[name]['name'],tag_text = "Hydrophobic Prot.")
            self.pages[name]['group3'].pack(side=LEFT,expand=0)

        self.pages[name]['combo'] = Pmw.ComboBox(self.pages[name]['group1'].interior(),
                                                 label_text=name,
                                                 labelpos='nw',
                                                 scrolledlist_items= list,
                                                 selectioncommand=self.show_meta_info,
                                                 listbox_height=17,
                                                 listbox_width=1,

                                                 dropdown=False)

        self.pages[name]['combo'].pack(side='left')


        self.pages[name]['text'] = Pmw.ScrolledText(self.pages[name]['group2'].interior(),
                                                    borderframe=5,
                                                    vscrollmode='dynamic',
                                                    hscrollmode='dynamic',
                                                    labelpos='n',
                                                    label_text=name,
                                                    text_width=50,
                                                    text_height=12,
                                                    text_wrap='none',
                                                    text_background='#000000',
                                                    text_foreground='green'
                                                    )
        self.pages[name]['text'].pack(expand=0,side=TOP)

        self.pages[name]['stat'] = Pmw.ScrolledText(self.pages[name]['group2'].interior(),
                                                    borderframe=5,
                                                    vscrollmode='dynamic',
                                                    hscrollmode='dynamic',
                                                    labelpos='n',
                                                    label_text='STATISTICS',
                                                    text_width=50,
                                                    text_height=10,
                                                    text_wrap='none',
                                                    text_background='white',
                                                    text_foreground='black'
                                                    )
        self.pages[name]['stat'].pack(expand=0,side=TOP)

        if name in ['BBH','SBH','SSH']:
            self.pages[name]['pval'] = DoubleVar()
            try:
                self.pages[name]['pval'].set(self.solv_vals[name])
            except:
                self.pages[name]['pval'].set(2.2)
                self.solv_vals[name]=2.2
            self.slider(self.pages[name]['group3'].interior(), self.pages[name]['pval'], 1.8, 2.6, 'Threshold')


        self.rgroup.selectpage(name)
#        if self.pages[name].has_key('group3'):
#            Pmw.alignlabels((self.pages[name]['group1'],
#                             self.pages[name]['group2'],
#                             self.pages[name]['group3'],
#                             self.buttonbox))
#        else:
#            Pmw.alignlabels((self.pages[name]['group1'],
#                             self.pages[name]['group2'],
#                             self.buttonbox))

        self.pages[name]['stat'].clear()
        string = self.get_statistics()
        self.pages[name]['stat'].insert('end',string)
        self.status_line.configure(text ='Loading %s' % name)
       # self.arrow_display()


    def get_statistics(self):

        text='INTERACTION STATISTICS:\n'
        text+='=======================\n'
        for key in INT_TYPES:
            try:
                text+='# %-7s = %d\n' % (key,len(self.ObjDic[key]))
            except:
                pass
        return text


    def zoom_inter(self):

        name = self.rgroup.getcurselection()
        if self.on_zoom:
            sel =  self.pages[name]['combo'].get()
            if sel:
                obj = self.nameDic[name][sel][0]
                if not obj.AO.on_screen:
                    cmd.zoom(sel)
                    obj.AO.on_screen=True
                    self.on_zoom = True
                else:
                    if self.myview:
                        cmd.set_view(self.myview)
                        self.on_zoom = False
                    else:
                        cmd.zoom()
            else:
                cmd.zoom()
        else:
            for key in self.nameDic.values():
                for l in key.values():
                    for obj in l:
                        obj.AO.on_screen = False

            sel =  self.pages[name]['combo'].get()
            if sel:
                obj = self.nameDic[name][sel][0]
                obj.AO.on_screen = True
                self.myview = cmd.get_view()
                cmd.zoom(sel)
                self.on_zoom = True


    def apply_solv(self):

        cur = self.rgroup.getcurselection()
        for name in ['BBH','SBH','SSH']:
            self.solv_vals[name]=self.pages[name]['pval'].get()

        text='Resorting with BBH=%2.2f, SBH=%2.2f, SSH=%2.2f' %\
              (self.solv_vals['BBH'],self.solv_vals['SBH'],self.solv_vals['SSH'])

        self.reset_hbonds()
        self.sort_hbonds(self.solv_vals)

        for key in INT_TYPES:
            if self.ObjDic.has_key(key):
                self.update_combo(key)

        for name in ['BBH','SBH','SSH']:
            self.pages[name]['pval'].set(self.solv_vals[name])


        self.status_line.configure(text = text)
        self.rgroup.selectpage(cur)

    def reset_hbonds(self):

        if self.ObjDic.has_key('UNC'):
            for obj in self.ObjDic['UNC']:
                self.ObjDic[obj.type].append(obj)
                if self.nameDic[obj.type].has_key(obj.hash):
                    self.nameDic[obj.type][obj.hash].append(obj)
                else:
                    self.nameDic[obj.type].update({obj.hash:[obj]})
            self.ObjDic['UNC'] = []
            self.nameDic['UNC'] = {}
        for key, l in self.ObjDic.items():
            if key!='EXCLUDED':
                for obj in l:
                    obj.AO.color = self.colorDic[key]
                    obj.constr=True
                    obj.status = 'CONSTRAINT'
#        for key in self.ObjDic.keys():
#            self.update_combo(key)

    def sort_hbonds(self,vals):
        uncl = []
        for key, val in vals.items():
            conl = []
            for obj in self.ObjDic[key]:
                if obj.type in HBOND_TYPES:
                    if obj.I.prot > val:
                        obj.constr=False
                        obj.status='UNCONSTRAINT'
                        obj.AO.color=self.colorDic['UNC']
                        obj.AO.build_arrow()
                        uncl.append(obj)
                    else:
                        conl.append(obj)
                        obj.constr=True
                        obj.AO.color = self.colorDic[key]
                        obj.AO.build_arrow()
            self.ObjDic[key] = conl
            self.nameDic[key] = {}
            for obj in conl:
                if self.nameDic[key].has_key(obj.hash):
                    self.nameDic[key][obj.hash].append(obj)
                else:
                    self.nameDic[key].update({obj.hash:[obj]})

        self.ObjDic['UNC'] = uncl
        self.nameDic['UNC'] = {}
        for obj in uncl:
            if self.nameDic['UNC'].has_key(obj.hash):
                self.nameDic['UNC'][obj.hash].append(obj)
            else:
                self.nameDic['UNC'].update({obj.hash:[obj]})



    def arrow_display(self):
        view = cmd.get_view()
        name = self.rgroup.getcurselection()
        if not self.on_display[name]:
            for obj in self.ObjDic[name]:
                cmd.load_cgo(obj.AO.cgo,obj.hash,state=1)
                obj.AO.on_screen = True

            self.status_line.configure(text = 'Showing all %s' % name)
            self.on_display[name] = True
        else:
            for obj in self.ObjDic[name]:
                cmd.delete(obj.hash)
                obj.AO.on_screen = False
            self.status_line.configure(text = 'Hiding all %s' % name)
            self.on_display[name] = False
        cmd.set_view(view)



    def change_color(self):
        name = self.rgroup.getcurselection()
        color = self.tk_color_dialog()
        if color:
            self.colorDic = color
            for obj in self.ObjDic[name]:
                obj.AO.color = color
                obj.AO.build_arrow()
            self.arrow_display()
            self.arrow_display()


    def remove(self):
        name = self.rgroup.getcurselection()
        sel =  self.pages[name]['combo'].get()
        if not self.colorDic.has_key('EXCLUDED'):
            self.colorDic['EXCLUDED']=[1,0,0]

        if sel:
            objl = self.nameDic[name][sel]
            for obj in objl:
                r1 = obj.r1
                r2 = obj.r2
                if (r1,r2) not in self.exclusions:
                    self.exclusions.append((r1,r2))
                    self.update_exclusion_card(r1,r2)
                obj.status='EXCLUDED'

                if self.ObjDic.has_key('EXCLUDED'):
                    self.ObjDic['EXCLUDED'].append(obj)
                else:
                    self.ObjDic['EXCLUDED'] = [obj]
                if self.nameDic.has_key('EXCLUDED'):
                    if self.nameDic['EXCLUDED'].has_key(obj.hash):
                        self.nameDic['EXCLUDED'][obj.hash].append(obj)
                    else:
                        self.nameDic['EXCLUDED'].update({obj.hash:[obj]})
                else:
                    self.nameDic['EXCLUDED']={obj.hash:[obj]}

                self.nameDic[name][sel].remove(obj)

                if not self.nameDic[name][sel]:
                    del self.nameDic[name][sel]
                self.ObjDic[name].remove(obj)

            for key in INT_TYPES:
                if self.ObjDic.has_key(key):
                    if key in ['BBH','SBH','SSH']:
                        self.solv_vals[key] = self.pages[key]['pval'].get()
                    self.update_combo(key)
        self.rgroup.selectpage(name)


    def update_exclusion_card(self,r1,r2):
        text='%5d  %5d\n' % (r1,r2)
        self.extext.insert(END,text)

    def slider(self, parent, variable, low, high, label):
        """Make a slider [low,high] tied to variable."""
        widget = Scale(parent, orient='vertical',
          from_=high, to=low,  # range of slider
          # tickmarks on the slider "axis":
          tickinterval=(high-low)/5.0,
          # the steps of the counter above the slider:
          resolution=(high-low)/100.0,
          label=label,    # label printed above the slider
          length=300,     # length of slider in pixels
          variable=variable)  # slider value is tied to variable
        widget.pack(side='right')
        return widget




    def fileopen(self, filename, mode):
        try:
            fp = open(filename,mode)
            return fp
        except:
            tkMessageBox.showerror('Error','Could not open file %s' % filename)
            return None

    def showAppModal(self):
        #self.dialog.activate(geometry = 'centerscreenalways', globalMode = 'nograb')
        self.dialog.show()
        #self.dialog.activate(geometry = 'centerscreenalways')


#
# The classes PmwFileDialog and PmwExistingFileDialog and the _errorpop function
# are taken from the Pmw contrib directory.  The attribution given in that file
# is:
################################################################################
# Filename dialogs using Pmw
#
# (C) Rob W.W. Hooft, Nonius BV, 1998
#
# Modifications:
#
# J. Willem M. Nissink, Cambridge Crystallographic Data Centre, 8/2002
#    Added optional information pane at top of dialog; if option
#    'info' is specified, the text given will be shown (in blue).
#    Modified example to show both file and directory-type dialog
#
# No Guarantees. Distribute Freely.
# Please send bug-fixes/patches/features to <r.hooft@euromail.com>
#
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

class Model:
    def __init__(self):
        self.atoms = []
        self.energy = 0.
        self.info = []
        self.num = 0
        self.as_string = ''
        self.info_as_string = ''
    def str2mod(self,string):
        list = string.split('\n')
        for line in list:
            if 'ATOM' in line:
                self.atoms.append(line.split(':')[1].strip())
                self.as_string+=line.split(':')[1].strip()+'\n'
            elif 'USER' in line:
                self.info.append(line.split(':')[1].strip())
                self.info_as_string+=line.split(':')[1].strip()+'\n'
            elif 'MODEL' in line:
                self.num = int(line.split()[2])
        for line in self.info:
            if 'Docked Energy' in line:
                x = line.split('=')[1]
                self.energy = float(x.split()[0])

    def writeMod(self,fp):
        print >>fp,'MODEL%8d' % self.num
        for at in self.atoms:
            print >>fp, at
        print >>fp,'ENDMDL'


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
