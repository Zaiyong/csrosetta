## make mammoth structure alignments

import string
from glob import glob
#from sys import argv,stderr,exit
#import sys
from os import popen,system,fdopen,mkdir,makedirs
from os import dup2,path
from os.path import exists
from operator import add
from math import sqrt
from os.path import basename
import argparse
import sys
import shutil
### toolbox library


import traceback

# definition of options for the method RASREC
sub_method_code = flag_lib+"/methods/_rasrec_base/rasrec_options.py"
if exists( sub_method_code ):
    exec open(sub_method_code, 'r' )
else:
    print "CANNOT FIND METHOD CODE %s"%sub_method_code
    exit()

class RasrecMethod(RasrecBaseMethod):
    def __init__(self,name,path):
        RasrecBaseMethod.__init__(self,name,path)

method = RasrecMethod(method_name, method_path)
