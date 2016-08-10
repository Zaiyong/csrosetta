#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
from sys import argv

assert( len(argv)>2)
rdc_file = argv[1]
max_resid=int(argv[2])
list = open(rdc_file,'r').readlines()

