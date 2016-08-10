cat *.flags | grep -v iterative | grep -v @ > flags_cyana_test
for i in $( cat *.flags | grep @ | awk -v FS="/" '{print $2}'); do cat $i | grep -v iterative >> flags_cyana_test; done
cp -Rf $CM_METHODPATH/tools .
cp -Rf $CM_METHODPATH/rdc_relax_patch .
cp -Rf $CM_METHODPATH/*_patch .
cp -Rf $CM_METHODPATH/rdc_pool_patch .


$HOME/rosetta_src/mini/bin/r_noe_assign.default.linuxgccrelease @../flags_cyana_test  -mute_warning core.io.pdb -mute_warning core.io.database -mute_warning core.pack -noesy:out:cst noe_auto_assi
gn.cst -noesy:out:min_seq_sep 5 >& cyana.log
$HOME/rosetta_src/mini/bin/r_noe_assign.default.linuxgccrelease @../flags_cyana_test  -mute_warning core.io.pdb -mute_warning core.io.database -mute_warning core.pack -noesy:out:cst noe_auto_assi
gn.cst.filter -noesy:out:min_seq_sep 2 >>& cyana.log

