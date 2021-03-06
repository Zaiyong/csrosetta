#!/bin/bash

mode=$1
if [ "$mode" != "verbose" ]; then
echo
echo "tutorial started in silent mode"
echo 'call as "tutorial_autoNOE.sh verbose" to get output of commands'
fi



# ----------------------------------------------------------------------------------------------------------------
# PREAMBLE
# the code in the PREAMBLE is specific to make these scripts work as tutorials.
# this code is not usually carried out by a user who performs the steps explained in the TUTORIAL part of this script
#
#
function tail_log() {
#parameter file, N
echo tail of output: $1
echo + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
tail -n $2 $LOGS/$1
echo + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
echo
}

# overwrite to keep tutorial targets out of users target-lib,
# note that a relative path here would mean that setup_target and setup_run
# have to be started from the same directory. With absolute paths this restriction does not apply
CS3_BENCH_TARGETLIB=$(pwd)/tutorial_targets


if [ ! -e inputs ]; then
    echo
    echo ERROR: the inputs directory of the CSROSETTA toolbox is required
    echo Solution1: either start the tutorials in the folder they are distributed in
    echo Solution2: issue the command 'ln -s <patch_to_csrosetta3>/tutorials/inputs' in your local folder
    echo Trying second solution automatically: ln -s $( dirname $( dirname `which setup_target` ) )/tutorials/inputs
    echo
    ln -fs $( dirname $( dirname `which setup_target` ))/tutorials/inputs
    if [ -e inputs ]; then
	if [ -e inputs/2jrm_trim.fasta ]; then
	    echo 'Succeeded in creating a link. Starting the tutorial'
	    echo
	else
	    echo 'Failed in creating the link'
	    echo
	    exit 1
	fi
    else
	echo 'Failed in creating the link'
	echo
	exit 1
    fi
fi

function cmd() {
#parameters cmd, log, N, crucial, exit_cmd
log=$2
N=$3
if [ $# -lt 4 ]; then
  exit_cmd='exit 1'
else
  exit_cmd="$4"
fi

if [ "$log" != "/dev/null" ]; then
    _log=$LOGS/$log.log
else
    _log=/dev/null
fi

if [ "$mode" == "verbose" ]; then
  echo CMD: $1
  if [ "$log" != "/dev/null" ]; then
    echo Redirecting output to $( echo $LOGS | sed s@$(pwd)/@@ )/$log.log
    $1 > $_log \
	&& tail_log $log.log $N || echo $exit_cmd
  else
    $1 > $_log || echo $exit_cmd
  fi
else
  echo CMD: $1
  $1 > $_log
fi

}


#overwrite to keep tutorial targets out of users target-lib
CS3_BENCH_TARGETLIB=$(pwd)/tutorial_targets

mkdir -p run_autoNOE_logs
LOGS=$(pwd)/run_autoNOE_logs

# END_PREAMBLE
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# START TUTORIAL
#
# NOTE 1: The redirecting of output with > log file && tail_log || exit
#    is also tutorial specific as a user you probably run these scripts interactively
#
echo
echo TUTOR: Tutorial \'autoNOE\' is starting...
echo

IN=inputs/NeR103A/
RUNDIR=run_autoNOE

#3. setup a ROSETTA run
if [ -e $RUNDIR ]; then
    echo
    echo RUN directory $RUNDIR exists already... Just checking for completed runs
    echo
    CONTINUE=yes
else
    CONTINUE=no
fi

if [ "$CONTINUE" == "no" ]; then
#1 prepare target from scratch
echo TUTOR: Setup Target with fragments, fasta and chemical shifts...
echo TUTOR: Assuming fragments have already been generated, see fragment picking tutorial...
cmd "setup_target -target ner103a -method autoNOE -frags $IN/ner*dat.gz \
    -fasta $IN/ner103.fasta -cs $IN/ner103.autotrim.tab" setup_target.1 15

#2a prepare prot file
echo TUTOR: We also need the chemical shift assignments in the .prot file format
cmd "bmrb2prot $IN/ner103a.bmrb ner103a.prot" make_prot 15
echo TUTOR: and the resonance assignment file needs to be trimmed to match the trimming of the fragments
cmd "renumber_prot ner103a.prot ner103a.trim.prot -fasta $IN/ner103.fasta" trim_prot 15

#2b prepare peak files
echo TUTOR: Prepare the peak files for reading by autoNOE-Rosetta
echo TUTOR: ali peaks...
cmd "clean_peak_file $IN/NeR103A_CASD_ali_noesy.list ali.peaks -cols 1 2 3 4 -skip 1 -names h C H I " clean_peaks.1 15

echo TUTOR: Prepare the peak files for reading by autoNOE-Rosetta
echo TUTOR: ali peaks...
cmd "clean_peak_file $IN/NeR103A_CASD_aro_noesy.list aro.peaks -cols 1 2 3 4 -skip 1 -names h C H I " clean_peaks.2 15

echo TUTOR: Prepare the peak files for reading by autoNOE-Rosetta
echo TUTOR: 15n peaks...
cmd "clean_peak_file $IN/NeR103A_CASD_n_noesy.list n.peaks -cols 1 2 3 4 -skip 1 -names h N H I " clean_peaks.3 15

#3 add peaks to target setup
echo TUTOR: Add the peaks to the setup of the target...
cmd "setup_target -target ner103a -method autoNOE -peaks *peaks -shifts ner103a.trim.prot" setup_target.2 15

rm -f ali.peaks n.peaks aro.peaks ner103a.prot ner103a.trim.prot
fi

if [ "$CONTINUE" == "no" ]; then
 echo TUTOR: Generate a CS-ROSETTA RUN directory
 echo
 cmd "setup_autoNOE -target ner103a -method autoNOE -dir $RUNDIR -job slurm -cycle_factor 0.05" setup_run 15
fi

#select here to run test or production (test is much faster, but doesn't produce meaningful results)
#RDIR=test
#RUN=test
 RUN=production
 RDIR=run
echo INFO: Selected to run $RUN version of protocol -- uncomment/comment in script to select between test and production

#change to run directory
cmd "cd $RUNDIR/ner103a/"
if [ "$CONTINUE" == "no" ]; then
for cst in $( ls -d cst* ); do
 # run the pre-assignment step
 cmd "cd $cst/$RDIR" /dev/null 0
 cmd ". initialize_assignments_phaseI.sh" preassign 15


 # start a 48-core job (200-500 cores are the recommended range)
queue=none
 echo
 echo TUTOR: 'Start a multiprocessor (48 core) job using the SLURM queuing system ...'
 echo CMD: sbatch -n 48 $RUN.slurm.job
 type sbatch >/dev/null 2>&1  && \
 sbatch -n 48 $RUN.slurm.job && queue=slurm || \
     echo 'sbatch not found, you probably have a different queuing system'
 echo
 if [ "$queue" == "none" ]; then
   echo TUTOR: 'Start a multiprocessor (64 core) job using the MOAB (as e.g., on JUROPA/Juelich) queuing system ...'
   echo CMD: msub -l nodes=4:ppn=16 $RUN.slurm.job
   type msub >/dev/null 2>&1  && \
   msub -l nodes=4:ppn=16 $RUN.slurm.job && queue=moab || \
   echo 'msub not found, you probably have a different queuing system'
 fi

 cmd "cd ../.." /dev/null 0
done
fi

echo
echo TUTOR: check log-files in directory logs_'<jobid>/'
echo TUTOR: the most important log is logs_'<jobid>/log_2', which contains the Pool-Manager process
echo TUTOR: log_0 is from the dedicated IO process and usually contains no errors
echo TUTOR: log_1 is from the job-distributor process and also usually without errors
echo TUTOR: log_3 and higher is from the worker processes. These are all identical.
echo TUTOR: look for errors in log_2 and log_3. If you get a specific message that job crashed with rank X,
echo TUTOR: look in log_X.
echo

echo TUTOR: check for completion - there should be a directory fullatom_pool_stage8

all_complete=1
pwd
for cst in $( ls -d cst* ); do
cmd "cd $cst/$RDIR"

if [ -e fullatom_pool_stage8 ] && [ -e fullatom_pool_stage8/decoys.out ]; then
    echo $cst $RUN 'run is completed'
else
    all_complete=0
    echo $RUN 'run is not completed'
    echo 'Wait (several hours/days)'
    if [ -e batch_000001 ]; then
	echo Batches generated so far:
	ls -d batch*
    fi
    if [ -e centroid_pool_stage1 ]; then
	echo Stages finished
	ls -d *pool_stag*
    fi
fi
cmd "cd ../../" /dev/null 0
done

if [ $all_complete -eq 0 ]; then
    echo some runs are not completed... wait more
    exit 0
fi

echo TUTOR: all runs are completed
echo TUTOR: run analysis script to check which restraint weight is optimal, and extract final structures
cmd "autoNOE_select_final_run -ensemble final.pdb" analysis 30


