#!/bin/bash

#overwrite to keep tutorial targets out of users target-lib
CS3_BENCH_TARGETLIB=../tutorial_targets

mkdir -p abrelax_notrim
cd abrelax_notrim

#run without trimming
setup_target -target 2jrm -method abrelax -frags ../inputs/2jrm.frags*dat.gz -fasta ../inputs/2jrm.fasta -cs ../inputs/2jrm.tab -flexible_residues $(seq 1 6) $( seq 53 60 )

setup_run -target 2jrm -method abrelax -dir . -cycle_factor 10 -nstruct 10 -job interactive

cd 2jrm/run

# single processor run
source production.interactive.job

# or using multiple cores
source production.interactive.job -n 10

#
echo TUTOR: Extract scores and low-energy decoys

echo CMD: extract_scores score rms chem_shifts
extract_scores score rms
extract_decoys decoys.out -score 2 > low2.out

echo TUTOR
extract_pdbs -score 2
cd ../..
