#!/bin/bash
#===============================================================================
# File Name    : spades.sh
# Usage        : sbatch spades.sh
# Author       : Lindsey Hall, hall.l.r@wustl.edu
# Version      : 1.0
# Created On   : 23-05-08
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=spades
#SBATCH --array=1-102%50
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --output=slurm_out/spades/x_spades_%a.out
#SBATCH --error=slurm_out/spades/y_spades_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

eval $(spack load --sh spades@3.15.3)

#module load spades/3.14.0-python-2.7.15

basedir="$PWD"
indir="${basedir}/d01_trimmedreads"
outdir="${basedir}/d04_spades"

mkdir -p ${outdir}

#sample list (no extensions)
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/all_seq_list.txt`

set -x

#run spades for each sample
spades.py \
  --isolate  \
  -1 ${indir}/${sample}_FW_clean.fastq \
  -2 ${indir}/${sample}_RV_clean.fastq  \
  -t ${SLURM_CPUS_PER_TASK}  \
  --tmp-dir /tmp/${sample}_tmp \
  -o ${outdir}/${sample}

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully!"
else
  echo "Error occured for Sample ${sample}"
  exit $RC
fi
