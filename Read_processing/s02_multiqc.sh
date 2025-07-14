#!/bin/bash

#===============================================================================
# Name         : s02_multiqc.sh
# Description  : Consolidates output from s01_fastqc.sh
# Usage        : sbatch s02_multiqc.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.4
# Created On   : 2019_01_08
# Modified On  : 2024_12_12
#===============================================================================

#Submission script for HTCF
#SBATCH --job-name=multiqc
#SBATCH --mem=4G
#SBATCH --output=slurm_out/multiqc/z_multiqc_%A.out
#SBATCH --error=slurm_out/multiqc/z_multiqc_%A.out

eval $( spack load --sh py-multiqc )

basedir="$PWD"
indir="${basedir}/d02_read_qc"
rawin="${indir}/raw_reads/fastqc"
cleanin="${indir}/clean_reads/fastqc"
rawout="${indir}/raw_reads/multiqc"
cleanout="${indir}/clean_reads/multiqc"

set -x

time multiqc ${rawin} -o ${rawout} -n raw_multiqc
time multiqc ${cleanin} -o ${cleanout} -n clean_multiqc

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi
