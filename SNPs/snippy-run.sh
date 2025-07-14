#!/bin/bash
#===============================================================================
# File Name    : snippy_02_run.sh
# Description  : This script runs each command in the file (snippy_commands.txt) generated from snippy_01_generate.sh in parallel and outputs them into their respective cluster folders
# Author       : Lindsey Hall, hall.l.r@wustl.edu
# Created      : 230629
#===============================================================================
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=snippy_02_run
#SBATCH --array=1-102%50
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/snippy/z_snippy_%a_%A.out
#SBATCH --error=slurm_out/snippy/z_snippy_%a_%A.out

eval $( spack load --sh miniconda3 )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate /ref/gdlab/software/envs/snippy

basedir="$PWD"
snippyclusters_output="${basedir}/d14_snippy/snippyclusters/snippy_commands.txt"

set -x

# Replace instances of ${SLURM_CPUS_PER_TASK} with 8 in snippy_commands.txt
sed -i 's/${SLURM_CPUS_PER_TASK}/8/g' "$snippyclusters_output"

# Read the commands from the modified text file into an array
mapfile -t commands < "$snippyclusters_output"

# Submit each command as a separate job in the array
command="${commands[$SLURM_ARRAY_TASK_ID - 1]}"
$command
RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in $command!"
  exit $RC
fi
