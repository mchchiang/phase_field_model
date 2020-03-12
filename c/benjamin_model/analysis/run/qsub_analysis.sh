#!/bin/bash
#
# SGE (Sun Grid Engine) submission script for an array of runs
#
# Usage: qsub -t job_range -N job_name ./qsub_analysis.sh analysis_script.sh
#        [script_args ...]
#
#$ -cwd              # run from current directory
#$ -V                # use all shell environment variables
#
# Choose a queue:
#$ -q cdt.7.day      # cdt.7.day, cm.7.day, sopa.1.day softcm.7.day
#
# Set job runtime
#$ -l h_rt=72:00:00  # time limit = 7days  (in sopa.1.day the time limit = 24h)
#
# -l h=phcomputecm0[12]
# -l h=statix
#
# Specify the standard output and error log file 
#$ -e stderr_$JOB_NAME.log
#$ -o stdout_$JOB_NAME.log
#
# Choose a parallel environment:
#$ -pe omp 64

analysis_sh=$1

export OMP_NUM_THREADS=${NSLOTS}

bash $analysis_sh ${@:2}
