#!/bin/bash
#===============================================================================
# File Name    : s08_mlst.sh
# Usage        : sbatch s08_mlst.sh
# Author       : Erin Newcomer, erin.newcomer@wustl.edu
# Version      : 1.0
# Created On   : 2020-06-03
# Last Modified: 2024-12-14
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=MLST
#SBATCH --mem=16G
#SBATCH --array=1-102%50
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/mlst/x_new_mlst_%a.out
#SBATCH --error=slurm_out/mlst/y_new_mlst_%a.err
#SBATCH --mail-user=ebenedict@wustl.edu
#SBATCH --mail-type=END

#load the module
eval $( spack load --sh mlst@2.22.1)
 
#store the base directory
basedir="$PWD"
indir="${basedir}/d04_spades"
outdir="${basedir}/d07_mlst"

#make the output directory
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/all_seq_list.txt`

mlst --minscore 50 ${indir}/${sample}/scaffolds.fasta >>  ${outdir}/EC_mlst_eeb.txt
