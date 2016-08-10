#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


import glob
from os import path
import os
import fnmatch
import subprocess
import library
import shutil

ignore_global=['.svn','.DS_Store','_IGNORE','_RELEASE','*.pyc','*~','*.pybk']

release_tag='ver3.0'
patch_tag='ver1.5'

def write_copyrighted_file( file_name, root ):
	copyright='''###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2013 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
### version: %s
### 
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###
'''%( release_tag.replace('ver','') )

	biopython_extra_notice='''
#                 Biopython License Agreement
#
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.
#
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE.
'''

	in_file=open(file_name,'r')
	in_mode=os.stat(file_name).st_mode

	target_name=root+file_name
	print 'copyrighting %s --> %s'%(file_name, target_name)
	library.mkdirp(path.dirname(target_name))
	out_file=open(target_name,'w')

	copyrighted=False
	biopython='python/PDB/' in file_name
	seen_bio_header=False
	count_bio_header=5
	for l in in_file:
		head1='#!'==l[0:2] and ( 'python' in l or 'bash' in l )
		head2='##-*-' in l and 'python' in l
		if not head1 and not head2 and not copyrighted:
			out_file.write(copyright)
			copyrighted=True
		if 'Thomas Hamelryck' in l:
			seen_bio_header=True
		if seen_bio_header: count_bio_header-=1
		if count_bio_header==0:
			out_file.write(biopython_extra_notice)
		if 'make mammoth structure alignments' in l: continue
		out_file.write(l)
	out_file.close()
	os.chmod(target_name, in_mode)

def remove( d, accept, ignore ):
#	print d, ignore
	if accept and not d in accept:
		return True
	if d in ignore:
		return True
	for i in ignore:
		if fnmatch.fnmatch(d,i):
			return True
	return False

os.chdir('/home/olange/')
archived_files=[]
for root,dirs,files in os.walk('csrosetta3'):
#	print "R",root
#	print "D", dirs
	accept=None
	if '_RELEASE' in files:
        #read all lines and remove whitespace include /n
		accept=map( lambda s: s.strip(), open(root+'/_RELEASE','r').readlines() )

	ignore=[]
	if '_IGNORE' in files:
		ignore=map( lambda s: s.strip(), open(root+'/_IGNORE','r').readlines() )
#		print ignore
	#the little [:] makes it an inplace modification of the list
	dirs[:]=[x for x in dirs if not remove( x, accept, ignore+ignore_global ) ]
	files[:]=[root+"/"+x for x in files if not remove( x, accept, ignore+ignore_global ) ]
	archived_files += files

fn_release_list='csrosetta3/maintainance/release_files_%s.txt'%release_tag
fd=open(fn_release_list,'w')
for f in archived_files:
	fd.write('%s\n'%f)
	print f
fd.close()


print 'make copyrighted versions of files'

release_tree='csrosetta3/release_files/'
library.mkdirp(release_tree+'csrosetta3')

for f in archived_files:
	if not 'tutorials/inputs' in f and not 'csrosetta3/database/' in f:
		write_copyrighted_file(f,release_tree)
	else:
		library.mkdirp(path.dirname(release_tree+'/'+f))
		shutil.copy(f,release_tree+'/'+f)

archive='csrosetta3/maintainance/csrosetta_toolbox_%s.tgz'%release_tag
print '\nstore files in archive %s...'%archive
subprocess.call("cd %s; tar -cvzf %s -T %s"%(release_tree, '../../'+archive, '../../'+fn_release_list), shell=True, stdout=subprocess.PIPE);

import sys
# upload to website must go through web interface because drupal remembers the file-size...
#
sys.stdout.write('\nupload to website...')
sys.stdout.flush()
site='di34ron@webdev1.lrz.de'
feeds='/nfs/web_tum/www/n/di34ron/webserver/drupal/sites/default/files/private/downloads/'
cmd='scp %s %s:%s/%s'%(archive,site,feeds,os.path.basename(archive))
print cmd
print 'run command and update file-size via phpmyadmin or upload manually'
subprocess.call(cmd, shell=True, stdout=subprocess.PIPE);

os.chdir('/home/olange/csrosetta3/maintainance')
site='di34ron@webdev1.lrz.de'
feeds='/nfs/web_tum/www/n/di34ron/webserver/drupal/sites/default/files/private/downloads/'
cmd='scp patch_rosetta3.4_to_CSROSETTA3_%s.txt %s:%s/'%(patch_tag,site,feeds)
print cmd
print 'run command and update file-size via phpmyadmin or upload manually'
subprocess.call(cmd, shell=True, stdout=subprocess.PIPE);

#
# if files are uploaded via csp their size needs to be fixed.
# create appropriate content by uploading a 0 length file
# overwrite file using scp
# go to table drupal_file_managed and run SQL Query like this
# UPDATE  `drupal_file_managed` SET filesize =221210759 WHERE fid =62
# where fid is the correct file-id for your file

print '...done'
