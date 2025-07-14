#!/bin/bash

#===============================================================================
# File Name    : s09_amrfinder.sh
# Description  : Annotations ARGs and VFs in assembled genomes or contigs
# Usage        : sbatch s09_amrfinder.sh
# Author       : Luke Diorio-Toth
# Version      : 1.3
# Created On   : Tue Oct 11 14:49:35 CDT 2022
#===============================================================================

#SBATCH --job-name=amrfinder
#SBATCH --array=1-102%50
#SBATCH --cpus-per-task=4
#SBATCH --mem=500M
#SBATCH --output=slurm_out/amrfinder/z_amrfinder_%a_%A.out
#SBATCH --error=slurm_out/amrfinder/z_amrfinder_%a_%A.out

eval $( spack load --sh amrfinder@3.10.42 )

basedir="$PWD"
indir="${basedir}/d09_prokka_hq"
outdir="${basedir}/d14_amrfinder_escherichia2"

# because there is inconsistent formatting to GFF files, specify the format
# options include prokka, bakta, pgap, etc. (see docs for the full list)
annotation_format="prokka"

#make output directory and read in the slurm array task
mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/all_seq_list.txt`

set -x
time amrfinder --plus \
  -n ${indir}/${sample}/${sample}.fna \
  -p ${indir}/${sample}/${sample}.faa \
  -g ${indir}/${sample}/${sample}.gff* \
  -a ${annotation_format} \
  --name ${sample} \
  -o ${outdir}/${sample}_out.tsv \
  --threads ${SLURM_CPUS_PER_TASK} \
  --organism Escherichia
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred in ${sample}!"
  exit $RC
fi
