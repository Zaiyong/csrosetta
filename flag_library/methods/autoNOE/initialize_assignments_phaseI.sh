rm -Rf initial_assignment
mkdir -p initial_assignment
#cat $CM_FLAGFILE | grep -v iterative | grep -v \@ > initial_assignment/flags_autoNOE_init
for i in $( cat flag_list.txt  | grep '^flags_' | grep -v flag_list.txt | grep -v flags_phase | grep -v initnoe | grep -v nmr_relax_patch ); do cat $i | grep -v \@flag |  grep -v iterative | sed s@\ ../@\ ../../@g >> initial_assignment/flags_autoNOE_init; done

cd initial_assignment
cp ../flags_initnoe* .
$CM_BINPATH/r_noe_assign.$CM_EXEC_EXT \@flags_autoNOE_init @flags_initnoe_sampling -out:levels core.io.pdb:warning core.io.database:warning core.pack:warning -noesy:out:cst noe_auto_assign.cst -noesy:out:min_seq_sep 5  > cyana.log

$CM_BINPATH/r_noe_assign.$CM_EXEC_EXT \@flags_autoNOE_init @flags_initnoe_filter -out:levels core.io.pdb:warning core.io.database:warning core.pack:warning -noesy:out:cst noe_auto_assign.cst.filter -noesy:out:min_seq_sep 2 -noesy:out:worst_prob_class 5 >> cyana.log

cd ..
