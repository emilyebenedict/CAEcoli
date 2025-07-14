
#!/usr/bin/env python3
#===============================================================================
# File Name    : snippy_00_structure.py
# Description  : This script reads in a csv with 3 columns: 1. sample_name 2.  cluster 3. reference. This script makes a folder named after each unique cluster value. Inside each folder, it should contain 2 text files: 1. samples.txt, containing a list of all of the sample names that have that cluster name; and 2. reference.txt, a text file with the unique reference name for that cluster. Then it goes through each folder and makes a final file in the output dir called snippy_commands.txt which contains snippy commands for every sample with the cluster reference, which will be outputted into the individual cluster folders
# Usage        : python3 snippy_00_structure.py <full path to input csv> <full path to output directory> <full path to assemblies directory> <full path to reads directory>
# Author       : Lindsey Hall, hall.l.r@wustl.edu
# Created      : 230629
#===============================================================================
import csv
import os
import sys

def create_folder_structure(csv_file, output_dir, assemblies_path, reads_path):
    # Read the CSV file
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        data = list(reader)

    # Create a dictionary to store cluster information
    clusters = {}

    # Iterate over each row in the data
    for row in data:
        sample_name, cluster, reference = row

        # Create a folder for each unique cluster if it doesn't exist
        if cluster not in clusters:
            folder_path = os.path.join(output_dir, cluster)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            clusters[cluster] = {
                'samples': [],
                'reference': set()
            }

        # Append the sample name to the samples list for the cluster
        clusters[cluster]['samples'].append(sample_name)

        # Add the reference to the reference set for the cluster
        clusters[cluster]['reference'].add(reference)

    print('Folder structure and text files created successfully.')

    # Create text files in each cluster folder
    for cluster, data in clusters.items():
        samples_file = os.path.join(output_dir, cluster, 'samples.txt')
        reference_file = os.path.join(output_dir, cluster, 'reference.txt')

        with open(samples_file, 'w') as file:
            file.write('\n'.join(data['samples']))

        with open(reference_file, 'w') as file:
            file.write('\n'.join(data['reference']))

    print('Text files created successfully.')

    # Create snippy_commands.txt file in the output directory
    snippy_commands_file = os.path.join(output_dir, 'snippy_commands.txt')
    with open(snippy_commands_file, 'w') as file:
        for cluster, data in clusters.items():
            folder_path = os.path.join(output_dir, cluster)
            samples_file = os.path.join(output_dir, cluster, 'samples.txt')
            reference_file = os.path.join(output_dir, cluster, 'reference.txt')

            # Read the reference name from the reference.txt file
            with open(reference_file, 'r') as ref_file:
                reference_name = ref_file.read().strip()

            with open(samples_file, 'r') as samples_file:
                samples = samples_file.read().splitlines()

            for sample in samples:
                assembly_file = os.path.join(assemblies_path, f"{sample}.fasta")
                fwd_read_file = os.path.join(reads_path, f"{sample}_FW_clean.fastq")
                rev_read_file = os.path.join(reads_path, f"{sample}_RV_clean.fastq")

                command = f"snippy --cpus ${{SLURM_CPUS_PER_TASK}} --outdir {folder_path}/{sample} --force --ref {assemblies_path}/{reference_name}.fasta --R1 {fwd_read_file} --R2 {rev_read_file}\n"
                file.write(command)


    print('snippy_commands.txt file created successfully.')


if __name__ == '__main__':
    # Check if the CSV file path, output directory, assemblies directory, and reads directory are provided as command-line arguments
    if len(sys.argv) != 5:
        print('Please provide the input CSV file path, output directory, assemblies directory, and reads directory as command-line arguments.')
        print('Usage: python3 snippy_clusters.py <input csv file path> <output directory> <assemblies directory> <reads directory>')
        sys.exit(1)

    # Get the command-line arguments
    csv_file_path = sys.argv[1]
    output_directory = sys.argv[2]
    assemblies_directory = sys.argv[3]
    reads_directory = sys.argv[4]

    # Check if the CSV file exists
    if not os.path.isfile(csv_file_path):
        print('The specified input CSV file does not exist.')
        sys.exit(1)

    # Check if the output directory exists or create it if it doesn't
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Call the function to create the folder structure, text files, and snippy_commands.txt file
    create_folder_structure(csv_file_path, output_directory, assemblies_directory, reads_directory)
