#!/usr/bin/env python3
#===============================================================================
# File Name    : snippy_04_pairwise.py
# Description  : this script compiles all cluster core snp matrices into a pairwise csv
# Author       : Lindsey Hall, hall.l.r@wustl.edu
# Usage        : python3 snippy_04_pairwise.py <full path to folder containing cluster subfolders> <full path to output directory> 
# Created      : 230629
#===============================================================================
import os
import csv
import sys

# Get the input and output folder paths from command-line arguments
input_folder = sys.argv[1]
output_folder = sys.argv[2]

# Iterate over each folder in the input folder
for folder_name in os.listdir(input_folder):
    folder_path = os.path.join(input_folder, folder_name)
    if os.path.isdir(folder_path):  # Check if it's a directory
        snp_matrix_path = os.path.join(folder_path, "core_snps", "core_snps_matrix.tsv")
        output_file = os.path.join(output_folder, f"{folder_name}_pairwise_distances.csv")

        # Read the SNP distance matrix
        with open(snp_matrix_path, "r") as matrix_file:
            reader = csv.reader(matrix_file, delimiter="\t")
            snp_matrix = list(reader)

        # Modify the header row
        if snp_matrix and len(snp_matrix[0]) > 0:
            snp_matrix[0][0] = "Sample"

        # Extract sample names from the first row and column
        sample_names = snp_matrix[0][1:] if snp_matrix and len(snp_matrix) > 0 else []

        # Open the output file for writing
        with open(output_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["sample1", "sample2", "distance"])  # Write the header row

            # Iterate over each sample in the matrix
            for i in range(1, len(sample_names) + 1):
                for j in range(i + 1, len(sample_names) + 1):
                    sample1 = snp_matrix[i][0] if len(snp_matrix) > i else ""
                    sample2 = snp_matrix[j][0] if len(snp_matrix) > j else ""
                    distance = snp_matrix[i][j] if len(snp_matrix) > i and len(snp_matrix[i]) > j else ""
                    writer.writerow([sample1, sample2, distance])

print("CSV files successfully generated.")

