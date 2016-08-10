#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

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
