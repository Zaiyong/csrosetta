#!/usr/bin/env python2.7
## make mammoth structure alignments

from basic.options import ExampleArgumentParser
import sys
from os.path import basename
import argparse
import string
import fasta


parser = ExampleArgumentParser(prog=basename(__file__), description="obtain molecular weight for fasta sequences",
examples=['%(prog)s target1.fasta > mol.weight',
          'cat target1.fasta | %(prog)s > mol.weight'])

parser.add_argument("fasta", nargs='?', help="pdb files", default=None );
args = parser.parse_args()

if not args.fasta or len(args.fasta)==0:
    fasta_file=sys.stdin
else:
    fasta_file=open(args.fasta,'r')


seq=fasta.read_fasta_stream(fasta_file)
print "%5.1f"%(fasta.molecular_weight(seq)/1000)

