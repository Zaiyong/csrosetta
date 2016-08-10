#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'


import string
from glob import glob
#from sys import argv,stderr,exit
#import sys
from os import popen,system,fdopen
from os import dup2
from os.path import exists
from operator import add
from math import sqrt
from os.path import basename
import argparse
import sys
import traceback
import cs
from assignment import noesy
### toolbox library
import library
import fasta
from sets import Set
import BmrbAtomNames
import math


class FloatingResonance( noesy.Resonance ):

	def __init__(self):
		noesy.Resonance.__init__(self)
		self._float_partners=Set()

	def add_partner_id(self, id):
		self._float_partners.add( id )

	def add_partner_ids(self, ids):
		for i in ids:
			self.add_partner_id( i )

	@classmethod
	def from_resonance(obj,resonance):
		obj=FloatingResonance()
		obj.__dict__ = resonance.__dict__
		obj._float_partners=Set()
		return obj

	def float_partners(self):
		return self._float_partners

	def float_partners_str(self):
		str='%s'%list(self._float_partners)
		return str.replace(' ','')


	def __str__(self):
		return "FLOAT: "+noesy.Resonance.__str__(self)+' %s'%list(self._float_partners)


class ProcessProtForAutoNOE:

	@staticmethod
	def get_parser():
		from basic.options import ExampleArgumentParser
		parser = ExampleArgumentParser(prog=basename(__file__), description="renumber residues in chemical shift file",
																	 examples=['%(prog)s original.prot prepared.prot -combine_threshold 0.01 -noheader'])

		parser.add_argument("infile", help="chemical shift file");
		parser.add_argument("outfile", help="chemical shift file",nargs="?",default="stdout");
		parser.add_argument("-fasta",help="figure out trimming from given sequence");
		parser.add_argument("-combine_threshold", type=float,
												help="if resonances of proR/proS protons are less different than threshold, they are combined",default=0.01)
		parser.add_argument("-combine_heavy_threshold", type=float,
												help="if resonances of proR/proS protons are less different than threshold, they are combined",default=0.0)
		parser.add_argument("-combine_stereo_threshold", type=float,
												help="stereo protons with frequency differences less than threshold are combined", default=0.01)
		parser.add_argument("-v", help='verbosity 0,1,2,3', type=int, default=1)
		parser.add_argument("-noheader", dest='header', help='output header or not', action='store_false', default=True)
		parser.add_argument("-nostereo", dest='stereo', help='with stereo column or not', action='store_false', default=True)
		parser.add_argument("-cyana_ssa", help='produce cyana_ssa constraint file to keep stereospecific assignments', default=None)
		parser.add_argument("-ambiguity", help='add ambiguity column in output file', action='store_true', default=False)
		parser.add_argument("-missing", action='store_true', help='if atoms are missing the respective ambiguity group is not formed', default=False)

		library.add_standard_args( parser )
		return parser

	def __init__( self, args_in ):
		self._args=args_in
	# take a resonance and get the pool-name of combined protons:
	# HB2 --> QB
	# HB --> QB ( already a pool )...
	# do we need to be aa-type dependent ?
	@classmethod
	def from_options( obj, **kwargs ):
		import argparse
		args=obj.get_parser().parse_args(['infile','outfile'])
		print args
		for k,v in kwargs.iteritems():
			setattr( args, k, v)

		obj=ProcessProtForAutoNOE( args )
		return obj

	def pool_name( self, name, aa ):
		try:
			pn, s = BmrbAtomNames.get_combine( aa, name )
			if not pn:
				pn = name
			rejects = None
			if 'QQ' in pn:
				rejects = BmrbAtomNames.get_rejects( aa, name )
			return pn, s, rejects
		except KeyError:
			return name, None, None

	# obtain the combinable pools of atoms in the list of resonances
	# the pools will be lists of tuples (r, s) with r= resonance and s=source atom -- e.g., the heavy atom (s) pool atoms are bound to
	def get_pools( self, resonances, aa ):
		pools={}
		for r in resonances:
			name=r.name()
			pn, s, rejects = self.pool_name(name,aa )
			if rejects:
				#check if any of the reject atoms is present in the original resonances
				if len([ True for rs in resonances if rs.name() in rejects ])>0:
					# global args
					if self._args.v>=2: print 'pool %s is not formed because also atoms with name %s are present'%(pn," ".join(rejects))
					pools.setdefault( name, [] ).append( (r,None) )
					continue
			pools.setdefault( pn, [] ).append( (r,s) )
			#fix the problem in files with missing protons:
			#if only an HE1 (for instance) is present and it is bound to a non-degenerate CE1 we cannot combine these
		for key, alist in pools.items():
			if 'Q' in key and len(alist)==1 and aa in 'FYW':
#				print key, alist
				if '1' in alist[0][0].name() or '2' in alist[0][0].name() or '3' in alist[0][0].name():
					del pools[key]
					pools[ alist[0][0].name() ] = alist

		return pools

	# calculate the mean and std of resonances in the pool
	# should maybe also take the error propagation into account
	def freq_stats_of_pool( self, resonances):
		std=0.0
		if len(resonances)>1:
			mean=0.0
			std=0.0
			for r,s in resonances:
				f=r.freq()
				mean+=f
				std+=f*f
	#		mean/=len(resonances)
	#		std/=len(resonances)
			#print mean, std, len(resonances)
			std-=mean*mean/len(resonances)
			std/=len(resonances)-1
			mean/=len(resonances)
			if std<1e-12: std=0.0
			return mean,math.sqrt(std)
		else:
			return resonances[0][0].freq(), 0.0

	def check_sources( self, pools, sources, proton_group, aa, resid, combined_pools, cannot_combine ):
		from sets import Set

		#if no sources are given, things are fine
		if len( sources ) == 0: return True

		# should sources also be combined into a pool ? PHE CD1, CD2 --> CD
		meta_pools = Set([ self.pool_name( s, aa)[0] for s in sources ])
		assert (len( meta_pools ) < 2)
		m_key = list(meta_pools)[0]  # name of the combined pool: e.g., CD
		if self._args.v >=3: print 'check heavy-atom pool %s for proton-group %s'%(m_key,proton_group)
	#	if self._args.missing:
	#		if self.atoms_missing_from_pool( aa, m_key, sources ):
	#			if self._args.v >= 2: print 'not combined because missing atoms'
	#			return False
		try:
			meta_pool = pools[ m_key ]
			m_mean,m_std = self.freq_stats_of_pool( meta_pool )
			if self._args.v >=3: print '             has std=%5.3f length=%5d sources[0]=%5s'%(m_std,len([ True for r,s in meta_pool if r.name() not in sources ]),sources[0])
			#check that all members of pool are also sources, otherwise no combination of meta_pool necessary
			# --> e.g., HD11, HD12, HD13 -->source: CD1 --> meta_pool CD ( CD1, CD2 ) but CD2 is not part of the sources, so we shouldn't combine to CD
			# --> but PHE QD (already combined in input) --> source CD --> pool CD ( CD1, CD2 ) --> source is pool-name --> combination required
			if len([ True for r,s in meta_pool if r.name() not in sources ])==0 or m_key==sources[0]:

				# okay, we should combine the heavy-atom pool because we need to integrate multiple sources: QQG -> CG1, CG2 -> CG, aromatics: HD1/CD1, HD2/CD2->QD/CD
				if m_std <= self._args.combine_heavy_threshold:
					if self._args.cyana_ssa or ( self._args.missing and not 'QQ' in proton_group ):
						return True
					cum_error=math.sqrt(sum([r.error()*r.error() for r,s in meta_pool]))
					newrcp = noesy.Resonance( id=meta_pool[0][0].id(), atom=noesy.Atom( m_key.replace( 'P','' ), resid ), freq=m_mean, error=m_std+cum_error )
					newrcp.ambiguity = meta_pool[0][0].ambiguity
					combined_pools[ m_key ] = newrcp
					if self._args.v >= 2:	print m_key.replace('P',''), meta_pool[0][0].id(), m_std
				else:
					if self._args.v >= 2: print '%5s %7.3f'%(m_key,m_std), " | ".join([ '%s'%v[0] for v in meta_pool ]), len(sources), sources
					cannot_combine.add( m_key )

					# couldn't combine this pool --> next
					return False
		except KeyError: #couldn't find heavy-atom for sources -- this might just be missing data.... report and move on
			if self._args.v >= 1: print 'Resonance %4s %3d is missing as heavy-atom for group %4s'%(m_key, resid, proton_group )
		return True

	def atoms_missing_from_pool( self, aa, pool_name, pool ):
		import BmrbAtomNames
		correct_size = len( BmrbAtomNames.anti_degenerate( aa, pool_name.replace('PH','Q').replace('P','') ) )
		return len(pool) < correct_size

	# figure out whether atom orgainzed in pools should be combined into a single resonance HD11, HD12, HD13 --> QD1
	def combine_pools( self, pools, aa): #, resid, reslist ):
		from sets import Set
		# global args
		combined_pools = {}
		cannot_combine = Set()
		cyana_stereo_constraints = []
		for key, value in pools.iteritems():
			resid=value[0][0].resid()
			combine_threshold=self._args.combine_threshold
			name=key
			if "PH" in key:
				name=key.replace("PH","Q") # make ProR/ProS - H also  into QH if frequencies are close
				combine_threshold=self._args.combine_stereo_threshold

			# compute mean and stddev of pool-resonances
			mean,std = self.freq_stats_of_pool( value )
			if self._args.v >= 2:	print '%5s %7.3f'%(key,std), " | ".join([ '%s'%v[0] for v in value ]), " @ "+" | ".join([ '%s'%v[1] for v in value ])
			# if pool has already been determined skip
			if key in combined_pools or key in cannot_combine:
				continue

			# if pool is haevy atom pool, don't combine directly, but wait for proton pools to trigger combination
			if key[0:2] == 'PN' or key[0:2] == 'PC': #shouldn't combine heavy-atoms if not driven by proton-combination
				continue
			if key[0] in 'CN':
				continue

			if self._args.missing:
				if self.atoms_missing_from_pool( aa, key, value ):
					if self._args.v >= 2: print 'not combined because missing atoms'
					continue
			# do not combine pool if stddev is high
			cyana_atoms = []
			if std > combine_threshold:
#				print 'should I combine stuff ---> ', key, value[0][0], value[0][0].ambiguity
				if ( "PH" in key ) or ( not value[0][0].ambiguity is None and ( "QQ" in key or "QB"==key or "QD"==key or "QE"==key ) ):
					ids=[]
					for r,s in value:
						if r.ambiguity>=2 or ( r.ambiguity is None ):
							ids.append(r.id())
						else:
							cyana_atoms.append( r.atom() )
					for k,(r,s) in enumerate(value):
						if len(ids):
							newr=FloatingResonance.from_resonance(r)
						#print '-->%s'%newr
							newr.add_partner_ids( ids )
							value[k]=(newr,s)
				if len(cyana_atoms):
					cya_string = ""
					for atom in cyana_atoms:
						cya_string += ' %s'%atom.name()
					cyana_stereo_constraints.append( 'atom stereo "'+cya_string+' %5d'%cyana_atoms[0].resid()+'"   # %s'%aa )
				continue

			# which heavy atoms are bound to pool atoms, are there different sources --> e.g., PHE HD1->CD1, HD2->CD2 --> combine to QD, CD ?
			sources = Set()
			for r,s in value:
				sources.add( s )
			# remove None
			sources = [ s for s in sources if s ]

			# can we combine protons without violating sources: e.g., if different source-heavy-atoms exist, they should have compatible frequencies
			if not self.check_sources( pools, sources, name, aa, resid, combined_pools, cannot_combine ):
				cannot_combine.add( key )
				if self._args.v >= 1:	print '--> cannot combine pool %4s %5d with frequency stdev = %5.3f because (heavy) source-atoms are different'%(key,resid,std)
				if len(value)==1: 		# is this a forced combine due to QX in input data ? --> raise exception...
					meta_pool=pools[ self.pool_name( list(sources)[0], aa)[0] ]
					raise library.InconsistentInput('combined proton-group %(key)s %(resid)d in input file is inconsistent with different frequencies for heavy atoms'%locals()+
																 "\n"+"\n".join([ '%s'%v[0] for v in meta_pool]))
				continue

			# this pool is okay for combination
			# new error should be std+sqrt(sum(std^2))
			cum_error=math.sqrt(sum([r.error()*r.error() for r,s in value]))
			newrc = noesy.Resonance( id=value[0][0].id(), atom=noesy.Atom( name, resid ), freq=mean, error=std+cum_error )
			newrc.ambiguity = value[0][0].ambiguity
			combined_pools[ key ] = newrc

		return combined_pools,cyana_stereo_constraints

	def generate_combined_resonances( self, pools, combined_pools ):
		new_resonances = []
		for key, value in pools.iteritems():
			try:
				new_resonances.append( combined_pools[ key ] )
	#			print 'gcr1: %s'%combined_pools[key]
			except KeyError:
				for v in value:
	#				print 'gcr2: %s'%v[0]
					new_resonances.append( v[0] )
		return new_resonances


	def clean_up_names( self, res_in ):
		for resid,resonances in res_in.iter_residues():
			aa = res_in.sequence()[resid-1]
			for r in resonances:
				name = r.name()
				# replace some names...
				if name == "HN": name ="H";
				# can we do this better using the information in BmrbAtomNames ?
				if aa == "L" and name == "HD1": name = "QD1";
				if aa == "L" and name == "HD2": name = "QD2";
				if aa == "V" and name == "HG1": name = "QG1";
				if aa == "V" and name == "HG2": name = "QG2";
				if aa == "T" and name == "HG2": name = "QG2";
				if aa == "A" and name == "HB" : name = "QB";
				if aa == "G" and name == "HA" : name = "QA";
				if aa == "I" and name == "HG1": name = "QG1";
				if aa == "I" and name == "HG2": name = "QG2";
				if aa == "I" and name == "HD1": name = "QD1";
				if aa == "M" and name == "HE" : name = "QE";
				if r.name() != name:
					if self._args.v >= 1: print 'replace name %s --> %s on aa: %s'%(r.name(),name,aa)
					r.atom()._name=name

		####
		### start of main program
		###
	def __call__(self, infile, outfile, fasta=None ):
			target_fasta=0

			from cs import ProtCSFile
			tab=ProtCSFile()
			tab.read_stream( infile )
			sequence=tab.sequence

			if not sequence and fasta:
				sequence=fasta
				tab.set_sequence(sequence)
			if not sequence and self._args.fasta:
				sequence=fasta.read_fasta(self._args.fasta)
				tab.set_sequence(sequence)

			#combine atoms into QX if possible/necessary
			res_in=noesy.ResonanceList.read_from_prot( tab )
			self.clean_up_names( res_in )

			res_out=noesy.ResonanceList()
			res_out.set_sequence( sequence )
			cyana_ss_constraints = []
			for resid,resonances in res_in.iter_residues():
		#		print 'reso: ',"\n".join(["%s"%r for r in resonances])
				aa=res_in.sequence()[resid-1]

				#copy heavy atoms
				if self._args.v >= 2: print 'residue %d %s round 1...'%(resid,aa)
				pools = self.get_pools( resonances, aa )
				combined_pools, cya_ss = self.combine_pools( pools, aa )
				new_resonances = self.generate_combined_resonances( pools, combined_pools )

				if self._args.v >=2: print 'residue %d %s round 2...'%(resid,aa)
				pools = self.get_pools( new_resonances, aa )
				combined_pools, cya_ss = self.combine_pools( pools, aa )
				new_resonances = self.generate_combined_resonances( pools, combined_pools )
				cyana_ss_constraints += cya_ss

				for r in new_resonances:
					res_out.add_resonance( r )


			prot_data=res_out.generate_dict()
			floats=[]
			ambiguity=[]
			for r in res_out.itervalues():
				ambiguity.append( r.ambiguity )
				try:
					floats.append( r.float_partners_str() )
				except AttributeError as exc:
		#			print exc
					floats.append( None )
			if self._args.stereo:
				prot_data['STEREO']=floats
			if self._args.ambiguity:
				prot_data['AMBIGUITY']=ambiguity

		#	print floats
			nih_table = cs.NIH_table().from_dict( prot_data )
#			print nih_table.vars
#			print nih_table.table
#			print 'convert to ProtCS-File'
			prot_file = cs.ProtCSFile().from_table( nih_table )
			prot_file.write( outfile, header=self._args.header )

			if self._args.cyana_ssa:
				fd = open( self._args.cyana_ssa, 'w')
				for line in cyana_ss_constraints:
					fd.write('%s\n'%line)

