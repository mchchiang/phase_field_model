#!/bin/bash
#
# SGE (Sun Grid Engine) submission script for an array of runs
#
# usage: qsub -t job_range -N job_name ./qsub_script.sh job_index_file
#
#$ -cwd              # run from current directory
#$ -V                # use all shell environment variables
#
# Choose a queue:
#$ -q sopa.1.day      # cdt.7.day, cm.7.day, sopa.1.day 
#
# Set job runtime
#$ -l h_rt=10:00:00  # time limit = 7days  (in sopa.1.day the time limit = 24h)
#
# Specify the standard output and error log file 
#$ -e stderr_$JOB_NAME.log
#$ -o stdout_$JOB_NAME.log
#
# Choose a parallel environment:
#$ -pe omp 4

job_index_file=$1
sim_name=$(sed "$SGE_TASK_ID"'q;d' $job_index_file)

exe="/Disk/feltz_staging/s1309877/phase_field_model/c/benjamin_model/bin/exe/run_phase_field_model"

local_dir=$PWD
job_in_dir=$local_dir/input/$sim_name
job_out_dir=$local_dir/output/$sim_name

node_job_dir=/scratch/s1309877_sim/$sim_name/

# Find the remote scratch directory when the job is running
remote_scratch=/Disk/$(hostname -s)_staging
# Replace 'scratch' with the correct remote scratch directory
remote_job_dir=${node_job_dir/\/scratch/$remote_scratch}

# Stagein data
stagein() {
    echo "Staging data from $job_in_dir to $remote_job_dir"
    mkdir -p $remote_job_dir \
	&& rsync -ravg $job_in_dir/* $remote_job_dir
    echo
}

# Stageout data
stageout() {
    echo "Staging data from $remote_job_dir to $job_out_dir"
    mkdir -p $job_out_dir \
	&& rsync -ravg $remote_job_dir/* $job_out_dir \
	&& rm -rf $remote_job_dir
    echo
}

# Actual computation
compute() {
    echo "Running computation"
    cd $node_job_dir
    export OMP_NUM_THREADS=${NSLOTS}
    $exe params_${sim_name}.txt > ${sim_name}.log
    cd $local_dir
    echo
}

###############################################################################

# Run the stagein, compute, and stageout steps

stagein || {
    echo "ERROR! Initial data stagein failed!"
    echo "Aborting job."
    exit 1
}
compute || {
    echo "ERROR! Compute job failed!"
    echo "Job's staging area has been left alone."
    echo "You can access this at: $remote_job_dir"
    exit 1
}
stageout || {
    echo "ERROR! Final data stageout failed!"
    echo "You may need to manually rescue data from job's staging area."
    echo "You can access this at: $remote_job_dir"
    exit 1
}
