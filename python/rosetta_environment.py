#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

#check environment to set some defaults
#ENV_ROSETTA_PATH='ROSETTA3_PATH'
ENV_ROSETTA_BIN='ROSETTA3_BIN'
ENV_ROSETTA_DB='ROSETTA3_DB'
ENV_ROSETTA_PLATFORM='ROSETTA3_PLATFORM'
ENV_ROSETTA_BUILD='ROSETTA3_BUILD'

bin=None
db=None
platform=None
build=None

import os
try:
    #	rosetta_path=os.environ[ENV_ROSETTA_PATH]
	bin=os.environ[ENV_ROSETTA_BIN]
	db=os.environ[ENV_ROSETTA_DB]
except KeyError:
	print '\n'+'*'*60
	print 'WARNING: Cannot find ROSETTA_XXX Environment variables'
	print '*'*60+'\n'
	db=None
	bin=None

try:
#	rosetta_path=os.environ[ENV_ROSETTA_PATH]
	platform=os.environ[ENV_ROSETTA_PLATFORM]
except KeyError:
	platform='linux'

try:
#	rosetta_path=os.environ[ENV_ROSETTA_PATH]
	build=os.environ[ENV_ROSETTA_BUILD]
except KeyError:
	build='mpi'


def rosetta_exec(cmd,mpi=None):
	if mpi:
		return cmd+'.mpi.'+platform+'gccrelease'
	else:
		return cmd+'.'+build+'.'+platform+'gccrelease'

