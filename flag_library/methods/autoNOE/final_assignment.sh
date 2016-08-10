# input the final structures you want to run assignment on... [default: fullatom_pool/low_10.out]
decoys=$1

rm -Rf final_assignment
mkdir -p final_assignment
cd final_assignment

cp ../initial_assignment/flags_autoNOE_init .
cp ../flags_phaseIII .

if [ "$decoys" = "" ]; then
   decoys=../fullatom_pool/low_10.out
   if [ ! -e $decoys ]; then
       if [ ! -e ../fullatom_pool/decoys.out ]; then
	   echo 'Cannot find any decoys, is the run completed ?'
	   exit 1
       fi
       extract_decoys -formula score-atom_pair_constraint-rdc -N 10 -verbose 0 > $decoys
   fi
fi

r_noe_assign.mpi.linuxgccrelease \@flags_autoNOE_init @flags_phaseIII -out:levels core:error -noesy:out:cst noe_auto_assign.cst -noesy:in:decoys $decoys -noesy:out:min_seq_sep 5 -noesy:out:peaks NOE_final.out