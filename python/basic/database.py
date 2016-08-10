#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

DATABASE_ENV='csrosettaDir'
bltin_open=open

import os

def open(file):
    if not DATABASE_ENV in os.environ:
        raise library.RunException('Environment variable "csrosettaDir" not found. Please source csrosetta3/com/init')
    filename=os.environ[DATABASE_ENV]+"/database/"+file
    return bltin_open(filename,'r')
