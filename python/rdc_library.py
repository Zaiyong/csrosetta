#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
## make mammoth structure alignments
from cs import NIH_table

class RDCFile(NIH_table):
	def __init__(self):
		NIH_table.__init__(self)
		self.sequence=None
#		self=None

	def read_file(self, file ):
#		self=NIH_table()
		NIH_table.__init__(self)
		NIH_table.read_file( self, file )
		if 'SEQUENCE' in self.data:
			self.sequence=self.data['SEQUENCE']

	def write( self, fd ):
#lf:
#			return
		if self.sequence:
			self.data['SEQUENCE']=self.sequence
		NIH_table.write( self, fd )

	def write_file( self, file ):
		fd=open( file, 'w' )
		self.write(fd)

	def renumber( self, start, end ):
		self.sequence= self.sequence[(start-1):end]
		NIH_table.renumber( self, start, end )
		self.data['SEQUENCE']=self.sequence
		resi2=self.get_slice( 'RESID2' )
		resi_shifted=table_op( resi2, resi2, lambda x,y: x-start+1 )
		self.add_slice( 'RESID2', resi_shifted )


	def from_table( self, table_in, sequence=None ):
#		self.table = NIH_table()
		NIH_table.__init__(self)
		if 'RESNAME' in table_in.vars:
			self.vars = ['RESID','RESNAME','ATOMNAME','RESID2','RESNAME2','ATOMNAME2','RDC']
		else:
			raise library.MissingInput("require sequence information to make talos-table out of prot-table")
		self.format='%4d %1s %2s %4d %1s %2s %8.3f'
		self.copy( table_in )
		self.sequence = sequence
		if not self.sequence and 'SEQUENCE' in table_in.data:
			self.sequence = table_in.data['SEQUENCE']
			self.data['SEQUENCE']=self.sequence
		#figure out if RESNAME needs translating from aa3-->aa1
		resn=self.get_slice( 'RESNAME' )
		if len(resn[resn.keys()[0]])==3:
			resn_translated=table_op( resn, resn, lambda x,y: amino_acids.longer_names[x.upper()] )
			self.add_slice( 'RESNAME', resn_translated )
		elif not len(resn[resn.keys()[0]])==1:
			raise library.InconsistentInput("RESNAME column in input table should have either aa3 or aa1 format, i.e., ALA or A for residue names")

		resn=self.get_slice( 'RESNAME2' )
		if len(resn[resn.keys()[0]])==3:
			resn_translated=table_op( resn, resn, lambda x,y: amino_acids.longer_names[x.upper()] )
			self.add_slice( 'RESNAME2', resn_translated )
		elif not len(resn[resn.keys()[0]])==1:
			raise library.InconsistentInput("RESNAME2 column in input table should have either aa3 or aa1 format, i.e., ALA or A for residue names")

#a Rosetta RDC-file is modelled by RDC_Data
class RDC_Line:
    def __init__( self, res1, atom1, res2, atom2, rdc ):
        self.atom1=atom1
        self.atom2=atom2
        self.res1=res1
        self.res2=res2
        self.rdc=rdc

class RDC_Data:
    def __init__( self ):
		 self.data=[]
		 pass

    def read_file( self, file ):
        self.data=[]
        lines = open( file  ).readlines()
        for line in lines:
            if line[0]=='#': continue
            if len(line) < 5: continue
            cols = line.split()
            self.data.append( RDC_Line( int( cols[0] ), cols[1], int( cols[2] ), cols[3], float( cols[4] ) ) );
		  return self

    #return list of the RDC values
	 def size(self):
		 return len(self.data)

    def rdcs( self ):
        r =[]
        for l in self.data:
            r.append( l.rdc )
        return r

    def estimate_Da_and_R_hist( self, rdcs=None, binwidth=None ):
        if not rdcs:
            rdcs=self.rdcs()

        mindata = min(rdcs)
        maxdata = max(rdcs)

        if not binwidth:
            binwidth = (maxdata - mindata) / 15.0;

        histogram = []
        bincenter = []
        currentbincenter = mindata + binwidth/2.0;
        numbins =  int((maxdata-mindata)/binwidth)
        for bin in range(numbins):
            histogram.append( 0.0 )
            bincenter.append( currentbincenter )
            currentbincenter += binwidth

        total = 0
        for rdc in rdcs:
            bin = int( (rdc - mindata)/binwidth )
            if bin<0       : bin = 0
            if bin>=numbins: bin = numbins - 1
            histogram[bin] += 1
            total = total + 1


        for bin in range(numbins):
            histogram[bin] = histogram[bin]/total * 100.0

		  max_pop=0
		  max_pos=-1
		  for bin,h in enumerate(histogram):
			  if max_pop<h:
				  max_pop=h
				  max_pos=bin

#		  print histogram
#		  print bincenter
#		  print max_pos

#compute Da and R
		  Dxx=bincenter[max_pos]
		  Dzz=maxdata
		  Dyy=mindata
#		  print Dxx,Dyy,Dzz

#what kind of weird math is that?
# only the else case seems to be visited
# and that is basically taking the number as is ...
# okay, this is no using the histogram data at all,
# what shoudl be done is a least-squares fit using Dxx, Dzz, Dyy from above, so that
# all 3 data points influence the final choice of R,D and not just the maximum and minimum of the distribution
        if max(abs(bincenter[0]),abs(bincenter[numbins-1]))==abs(bincenter[0]):
			  Dzz=abs(bincenter[0])*abs(bincenter[0])/bincenter[0]
			  Dyy=abs(bincenter[numbins-1])*abs(bincenter[numbins-1])/bincenter[numbins-1]
        else:
			  Dzz=abs(bincenter[numbins-1])*abs(bincenter[numbins-1])/bincenter[numbins-1]
			  Dyy=abs(bincenter[0])*abs(bincenter[0])/bincenter[0]
        Dxx=-Dzz-Dyy
        R=(1+2*Dxx/Dzz)*2/3
        rdc_range=maxdata-mindata
        inv_rdc_range2=1/(rdc_range*rdc_range)
        return Dzz/2,R,rdc_range,inv_rdc_range2

