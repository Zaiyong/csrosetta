##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'


from os import path
import argparse
import numpy

### toolbox library
from library import MethodException

import traceback
from silent_lib import ReadSilentData


# definition of options for the method RASREC
sub_method_code = flag_lib+"/methods/_rasrec_base/rasrec_options.py"
if path.exists( sub_method_code ):
	exec open(sub_method_code, 'r' )
else:
	print "CANNOT FIND METHOD CODE %s"%sub_method_code
	exit()

if 'group' in locals():
	group.add_argument('-contacts', nargs=1, metavar='<target>.cmp', help='contact map prediction')

if 'run_group' in locals():
	run_group.add_argument('-cutoff', help='prediction threshold for contacts', type=float, required=True )
	run_group.add_argument('-combine', help='use combination of restraints', action='store_true', default=False)
	run_group.add_argument('-upper_dist', help='set upper bound', type=float, default=10)
	run_group.add_argument('-lower_dist', help='set lower bound', type=float, default=1.5)
	run_group.add_argument('-min_sep', help='minimum sequence separation for generated restraints', type=int, default=4)

class ContactMethod(RasrecBaseMethod):
	def __init__(self,name,path):
		RasrecBaseMethod.__init__(self,name,path)
        #self.option2dir['contacts']='bioinformatics'

    # builds a rosetta constraint file given a contactmap
	def generate_contact_restraints( self, contactmap_file, cutoff, min_sep, upper_bound, lower_bound , run_dir):
      #output cst file will be located in actual run directory
		outfile = run_dir + "/inputs/contactmap_restraints.cst"

		#parse textfile and store contact map in a numpy matrix
		count = 0
		for line in open( contactmap_file, 'r'):
			split = line.split()
			if (count == 0):
				ma = numpy.zeros((len(split), len(split)), numpy.float)
			for i in range (0, len(split)):
				ma[count, i] = numpy.float(split[i])
			count += 1
		print "SEQUENCE: %5d %s"%(count,self.fasta)
        #check whether contactmap length matches the given fasta file
		assert ( len(self.fasta) == ma.shape[0] ), \
		"Length of fasta sequence (" + str(len(self.fasta)) + ") " + \
		"is different from length of contactmap (" + str(ma.shape[0]) + ")."

        #check if matrix is symmetric
		assert ( ( ma == ma.T ).all() ), \
		"Contactmap is not symmetric."

		o = open( outfile, 'w')
		for i in range(0,ma.shape[0]):
			#Matrix is symmetrical -> no need to run through whole matrix
			for j in range( i+1, ma.shape[1]):
				# only create restraints if the value in the matrix is higher than the cutoff and if the two residues are not closer than the min. distance.
				if ( ma[i,j] > cutoff ) and abs(i-j) > min_sep:
					o.write("AtomPair CA %4d CA %4d BOUNDED %3.1f %3.1f 1 ContactPreds %5.1f\n"%(i+1,j+1,lower_bound,upper_bound,ma[i][j]))

	def make_target_flags(self, run, setup, filename, flags, subs ):
		RasrecBaseMethod.make_target_flags( self, run, setup, filename, flags, subs )
		args=self.get_args()

		fl = self.file_library
		flags.write("-broker:setup @@setup_contacts.tpb")

		if not args.combine:
			fl.override("others", "setup_contacts.tpb","COMBINE_RATIO 1")

		self.generate_contact_restraints( setup.abspath(args.contacts), args.cutoff, args.min_sep, args.upper_dist, args.lower_dist, run.rundir() )
		self.add_restraint_scores()#add restraint scores

	def setup_file_library( self ):
		RasrecBaseMethod.setup_file_library( self )
		args=self.get_args()
		fl = self.file_library
		path = flag_lib+"/methods/contacts/"
#        fl.provide_file( "flags", path, "flags_contactsbroker" )
		fl.provide_file( "others" , path, "setup_contacts.tpb")
        #fl.show()

	def motd( self, rundir ):
		RasrecBaseMethod.motd( self, rundir )


#echo "Starting from beta.top stage3 .... Make sure it is not a helical protein, eventually remove this line from flags_noe_assign"
#echo -iterative:initial_beta_topology init_phaseII/beta.top >> flags_noe_assign
#fi

method = ContactMethod(method_name, method_path)
