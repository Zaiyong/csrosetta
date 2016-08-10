#!/bin/bash

unpack_pdbs final.pdb

#have to get rid of extra residues for RDC calculations
for i in $( ls final_*.pdb ); do renumber_pdb $i $i.trim -fasta $CM_FASTA ; done

rm -f rescored.out
$CM_BINPATH/score_jd2.$CM_EXEC_EXT -in:file:s final_*.pdb.trim -out:file:silent rescored.out @flags_rescore -out:level 200 -out:levels core:100 >> log_rescore.txt

