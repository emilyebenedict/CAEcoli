#!/bin/bash
#===============================================================================
#
# File Name    : s04_checkm.sh
# Description  : This script will assess the quality of genome assemblies.
# Usage        : sbatch s04_checkm.sh
# Author       : Jian Ryou, j.ryou@wustl.edu
# Version      : 1.0
# Created On   : Wed Aug 17 17:37:20 CST 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=checkm
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G
#SBATCH --output=slurm_out/checkm/x_checkm_%a.out

eval $( spack load --sh py-checkm-genome )
#module load checkm

basedir="$PWD"
indir="${basedir}/d04_spades"
tmpdir="${basedir}/tmp_checkM_eeb"
outdir="${basedir}/d05_checkM"

#make output directory
mkdir -p ${outdir}

#read in the slurm array task
#sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/all_seq_list.txt`

# copy scaffold files from unicycler to a temporary directory (checkM requires all your files being in a single folder)
# IMPORTANT: the find command is meant to find the output dirs from flye.
mkdir -p ${tmpdir}

#cp ${indir}/${sample}/scaffolds.fasta ${tmpdir}/${sample}.fasta

set -x


time checkm lineage_wf -f ${outdir}/EC_checkm_output_eeb.txt\
        -t ${SLURM_CPUS_PER_TASK} \
        -x fasta \
        --tab_table \
        ${tmpdir} ${outdir}

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
