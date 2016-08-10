#!/bin/bash

mkdir -p pick_fragments
cd pick_fragments
#1. prepare input files
wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15339 -O 2jrm.bmrb
bmrb2talos 2jrm.bmrb -ignore_errors > 2jrm.tab

#3.1 figure out flexible regions in original protein sequence
echo running talos+ to figure out flexible regions
mkdir -p untrimmed
cd untrimmed
ln -sf ../2jrm.tab
talos+ -in 2jrm.tab > talos.log
cat pred.tab
cd ..

#3.2 from the pred.tab file we found that we want to trim off residue 1-6 at N-terminal and 53-60 at C-terminal
renumber_talos -s 7 -e 52 2jrm.tab 2jrm_trim.tab
talos2fasta 2jrm_trim.tab > 2jrm_trim.fasta

#3.3 now pick the fragments
echo picking fragments...
mkdir -p fragments
cd fragments
ln -sf ../2jrm_trim* .
pick_fragments -cs 2jrm_trim.tab -hom
cd ..


