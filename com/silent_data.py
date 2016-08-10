#!/usr/bin/env python2.7
## make mammoth structure alignments

import string
from silent_lib import ReadSilentData
import sys

from os.path import basename
import argparse
from basic.options import ExampleArgumentParser
import traceback
import library


#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

parser = ExampleArgumentParser(prog=basename(__file__), description="obtain fasta sequences from provided pdb-files",
examples=['%(prog)s silent_file.out rms score > rms_score.txt',
          '%(prog)s scores_silent.fsc rms atom_pair_constraint rdc > rms_score.txt'],
aliases=['extract_scores','silent_data'])

parser.add_argument("silentfile", help="input silent file", default=None)
parser.add_argument('columns', nargs='*', help='columns for output')
parser.add_argument('-pipe', help="read from stdin", default=False, action='store_true')
parser.add_argument('-header', help="print header in first line", action="store_true", default = False )

library.add_standard_args( parser )

args = parser.parse_args()
names = args.columns

try:
  if args.pipe:
    infile=sys.stdin
    # if pipe, then first positional argument is actually a column... so args.silentfile is added to names
    names=[args.silentfile]+args.columns
  else:
    infile = open(args.silentfile,'r')
  read_rt = 'RT' in names
  sfd = ReadSilentData( names )
  if args.header:
      print " ".join( list ( x[0] for x in names ) )
  for l in infile:
    data=sfd.read_line( l )
    if data:
      if read_rt:
        rt_data = None
        while ( not rt_data ):
          try:
            next_line = infile.next()
            rt_data = sfd.read_line( next_line, data )
          except StopIteration:
            print 'not found RT until line ', next_line
            rt_data = data
        data=rt_data
      if isinstance( data, str ):
        print data
      else:
        print " ".join(data)


except library.LibException as inst:
  if args.traceback:
    print traceback.print_exc(inst )
  else:
    print traceback.print_exception(sys.exc_type, sys.exc_value, None)
