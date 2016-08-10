#!/bin/bash
rm -Rf init_phaseII
mkdir -p init_phaseII
#cat $( basename $CM_FLAGFILE ) | grep -v iterative | grep -v \@ > init_phaseII/flags_autoNOE_init
for i in $( cat flag_list.txt  | grep flag ); do cat $i | grep -v \@flag |  grep -v iterative >> init_phaseII/flags_autoNOE_init; done

# --- run one round of auto-assignment
cd init_phaseII

$CM_BINPATH/r_pdb2top.$CM_EXEC_EXT -in:file:silent $CM_AUTONOE_DECOYS -out:top beta.top

$CM_BINPATH/r_noe_assign.$CM_EXEC_EXT \@flags_autoNOE_init -noesy:in:decoys $CM_AUTONOE_DECOYS -noesy_weights:dcut $CM_AUTONOE_DCUT  -cycle 7 -noesy:out:cst assigned.cst.filter -noesy:no_network > LOG_2

$CM_BINPATH/r_noe_assign.$CM_EXEC_EXT \@flags_autoNOE_init -noesy:in:decoys $CM_AUTONOE_DECOYS -noesy_weights:dcut $CM_AUTONOE_DCUT  -cycle 7 -noesy:out:cst assigned.cst -out:min_seq_sep 5 -noesy:no_network > LOG_5

cat LOG_2 | grep assign.parameters > NOESY_PARAM
cd ..
# --- finished preparing the first auto-assignment


