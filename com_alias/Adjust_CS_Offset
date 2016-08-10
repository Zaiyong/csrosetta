awk '   
	NF~4 && $1!="DATA"{
		if( $3==atomName) printf(" %4d %1s %4s %8.3f\n", $1,$2,$3,$4+offset);
		else printf(" %4d %1s %4s %8.3f\n", $1,$2,$3,$4);
	}
	NF!=4 || $1=="DATA"{ print };
' atomName=$2 offset=$3 $1 > temp.cs_offset_adj

 mv temp.cs_offset_adj $1


# CS-ROSETTA: System for Chemical Shifts based protein structure prediction using ROSETTA
# (C) Shen and Bax 2007-2008, Lab of Chemical Physics, NIDDK, NIH
#
#
# Add_CS_offset.com: add offset to the checmical shifts of give atom type
#	syntax: Adjust_CS_offset.com cs.tab CA -2.7
#		add -2.7ppm offset to all CA chemical shift in cs.tab
