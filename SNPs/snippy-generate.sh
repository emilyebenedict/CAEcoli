#!/bin/bash
#===============================================================================
# File Name    : snippy_01_generate.sh
# Description  : This script runs snippy_00_structure.py
# Usage        : python3 snippy_clusters.py <full path to input csv> <full path to output directory> <full path to assemblies directory> <full path to reads directory>
# Author       : Lindsey Hall, hall.l.r@wustl.edu
# Created      : 230629
#===============================================================================
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=snippy_01_generate
#SBATCH --array=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/snippy/z_snippy_%a_%A.out
#SBATCH --error=slurm_out/snippy/z_snippy_%a_%A.out

basedir="$PWD"
assemblies="${basedir}/tmp_checkM_eeb"
reads="${basedir}/d01_trimmedreads"
outdir="${basedir}/d14_snippy/snippyclusters"
snippyclusters="${basedir}/snippy_00_structure.py"
snippyclusters_input="${basedir}/clusters_for_snippy.csv"

set -x
#python3 snippy_00_structure.py <input csv file path> <output directory> <assemblies directory> <reads directory>
python3 ${snippyclusters} ${snippyclusters_input} ${outdir} ${assemblies} ${reads}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi
