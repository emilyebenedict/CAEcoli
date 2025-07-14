#!/bin/bash
#===============================================================================
# File Name    : snippy_03_core.sh
# Description  : This script runs snippy-core on all clusters
# Author       : Lindsey Hall, hall.l.r@wustl.edu
# Created      : 230629
#===============================================================================
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=snippy_03_core
#SBATCH --array=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/snippy/z_snippy_%a_%A.out
#SBATCH --error=slurm_out/snippy/z_snippy_%a_%A.out

eval $( spack load --sh snp-dists )
eval $( spack load --sh miniconda3 )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate /ref/gdlab/software/envs/snippy

basedir="$PWD"
assemblies="${basedir}/tmp_eeb_checkM"
outdir="${basedir}/d14_snippy/snippyclusters"

set -x

# Loop through each cluster folder in the output directory
for cluster_dir in "$outdir"/*/; do
    # Get the cluster name from the cluster folder path
    cluster_name=$(basename "$cluster_dir")

    # Create the "core_snps" folder in the cluster folder
    core_snps_dir="$cluster_dir/core_snps"
    mkdir -p "$core_snps_dir"

    # Change to the "core_snps" folder
    cd "$core_snps_dir" || exit

    # Get the reference file path from the "reference.txt" file
    reference_file="../reference.txt"
    reference="${assemblies}/$(cat "$reference_file").fasta"

    # Get the samples file path from the "samples.txt" file
    samples_file="$cluster_dir/samples.txt"

    # Format the sample paths for snippy-core command
    sample_paths=$(awk '{print "../"$0}' "$samples_file" | paste -sd " ")

    # Run the snippy-core command with the reference and samples
    snippy-core --ref "$reference" --prefix snps $sample_paths

    snp-dists *.aln > core_snps_matrix.tsv
    
    # Move back to the parent directory (cluster folder)
    cd ../../

done

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in $command!"
  exit $RC
fi
