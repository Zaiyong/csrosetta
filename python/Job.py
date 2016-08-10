#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
## make mammoth structure alignments

import string
from os.path import basename,dirname

class Job:
	def __init__(self, fn_job ):
		self.job=splitext(basename(fn_job))[0]
		self.jobpath=dirname(fn_job)
		jobbase=self.jobpath+"/"+self.job
            #platform_file = $jobpath/$( echo $jobname | awk -v FS="_" '{print $1}' ).generic
		tags=string.split( self.job, "_" )
		self.platform_file= self.jobpath+"/"+tags[0]+".generic"

	def __str__(self ):
		return "job %s  -- for platform %s\n"%(self.job, self.platform_file )
