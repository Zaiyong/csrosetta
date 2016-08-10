#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

import argparse
import sys
import string


parser=argparse.ArgumentParser(description="transform argparse help output into formatted html ready for putting on website",
                                add_help=True)
parser.add_argument("-method", help="only output information that is specific for the method", default=None)

args = parser.parse_args()

input=sys.stdin
usage=[]
desc=[]
help={}
help_key_list=[]
last_help=[]
last_help_key=None
help_offset=0
progname=None
usage_block_active=False
desc_block_active=False
positional_block_active=False
help_block_active=False
spec_method_help_active=False
example_block_active=False
examples=[]
#parse input line by line,
#collect information in fields like usage, desc, help and examples
def scan_usage(start_usage, usage):
	for l in sys.stdin:
		tags=l.split()
		if len(tags)==0:
			return usage
		usage.append(l[start_usage:-1])

def scan_description():
	desc=[]
	for l in sys.stdin:
		tags=l.split()
		if len(tags)==0:
			return desc
		desc.append(l[:-1])
	return desc

def scan_examples():
	examples=[]
	expos=9999
	last_ex=None
	first=True
	for l in sys.stdin:
		tags=l.split()
		if len(tags)==0:
			if first: continue
			if last_ex:
				examples.append((last_ex, exdesc))
				last_ex=None
				continue
			return examples

		if expos<len(l) and len(l[:expos+2].strip())==0:
			exdesc.append((string.strip(l)))
			continue

		if last_ex:
			examples.append((last_ex, exdesc))

		expos=l.index(tags[0])
		last_ex=l.strip()
		exdesc=[]
		first=False
	return examples

def scan_positional( help ):
	help_pos=99999
	last_key=None
	for l in sys.stdin:
		tags=l.split()
		if len(tags)==0:
			if last_key:
				help[last_key]=last_help
			return help
		if help_pos<len(l) and len(l[:help_pos].strip())==0:
			last_help.append(l[help_pos:-1])
			continue
		if last_key:
			help[last_key]=last_help
		help_pos=l.index(tags[1])
		last_key=tags[0]
		last_help=[l[help_pos:-1]]
	return help

for l in sys.stdin:
	tags=l.split()
	if len(tags)==0: continue
	if "usage:" in tags[0]:
		progname=tags[1]
		start_usage=len("usage: ")+len(progname)+1
		usage=[l[start_usage:-1]]
		usage=scan_usage( start_usage, usage )
		desc=scan_description()
		continue
	if "positional arguments:" in l:
		help = scan_positional( help )
		help_key_list=help.keys()
#		print 'after positionals', help
		help_block_active
		continue

	if help_block_active and tags[0]=='examples:':
		help_block_active=False
		examples=scan_examples()
		continue

	if "optional arguments:" in l:
		help_block_active=True
		continue
	if not help_block_active:
		continue
	if tags[0]=="-h,":
		help_offset=l.index('show')
		last_help_key=string.strip(l[0:help_offset])
		last_help=[string.strip(l[help_offset:-1])]
		continue
	if "extra options for method" in l or " extra options selected when" in l:
		spec_method_help_active=True
		continue
	if args.method and not spec_method_help_active:
		continue
	if string.strip(l)[0]=='-':
		help[last_help_key]=last_help
		help_key_list.append(last_help_key)
		#have to figure out if help starts in same line or next one...
		#the 2nd token is the argument, unless it is a boolean flag (-overwrite)
		last_help=[]
		last_help_key=string.strip(l[0:help_offset])

		tag_end=l.index(tags[0])+len(tags[0])
#		if tag_end>help_offset:
#			last_help_key=string.strip(l[0:last_tag_index])
		if len(tags)<2:
			continue
		def find_next_tag( start, line):
			def find_end_brace( start, cstart, cend, line):
				level=1
				i=start+1
				while i<len(line) and level>0:
					if line[i]==cstart: level+=1
					if line[i]==cend: level-=1
					i+=1
				return i
			begin=start
			while line[begin]==' ': begin+=1
			c=line[begin]
			if c=='[': return begin,find_end_brace(begin,'[',']',line)
			elif c=='{': return begin,find_end_brace(begin,'{','}',line)
			return begin, find_end_brace(begin,None,' ',line)
		last_tag_end=tag_end
		tag_begin, tag_end = find_next_tag( tag_end, l )
		#all tags that end before help_offset are definitely part of key.
		#if the first tag after help_offset is starting exactly at help_offset it is likely start of the description
		while tag_end<help_offset and tag_end<len(l):
			last_tag_end=tag_end
			tag_begin, tag_end = find_next_tag( tag_end, l )
		#this is a line that also has description starting here
		if tag_begin==help_offset:
			last_help_key=l[0:last_tag_end]
			st=string.strip(l[help_offset:-1])
			if len(st):
				last_help.append(st)
		else:
			#a line where the cmd-section goes all the way to the end
			last_help_key=l[:-1]
		continue

	if len(l)>help_offset and string.strip(l[help_offset:-1])>0:
		last_help.append(string.strip(l[help_offset:-1]))


if last_help_key:
	help[last_help_key]=last_help
	help_key_list.append(last_help_key)

#print examples

if not args.method:
 print "<h2>Description:</h2>"
 for i in desc:
	print i
 print "<h2>Usage:</h2>"
 print "<table>"
 print "<tr><td>%10s</td><td>%s</td></tr>"%(progname,usage[0])
 for i in usage[1:]:
	print "<tr><td>%10s</td><td>%s</td></tr>"%(" "*8,i)
 print "</table>"

extra_msg=""
if args.method:
	extra_msg=" (only %s specific options)"%args.method
print "<table>"
print "<h2>Detailed Help:%s</h2>"%extra_msg
for i in help_key_list:
	print "<tr><td>%s</td><td>%s</td></tr>"%(i, " ".join(help[i]))
print "</table>"

if len(examples):
	print "<h2>Examples:</h2>"
#	print examples
	start="<pre>"
	str=""
	for ex in examples:
		str=str+start+ex[0]
		if len(ex[1]) and len(ex[1][0]):
			str=str+"</pre>"+" ".join(ex[1])
			start="<pre>"
			stop=""
		else:
			start="\n"
			stop="</pre>"
	print str+stop

#print help


