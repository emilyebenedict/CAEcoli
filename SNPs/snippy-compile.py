#!/usr/bin/env python3
#===============================================================================
# File Name    : snippy_05_compile.py
# Description  : this script combines all pairwise csvs into one, excluding rows that contain "Reference"
# Author       : Lindsey Hall, hall.l.r@wustl.edu
# Usage        : python3 snippy_05_compile.py <full path to folder containing cluster subfolders> 
# Created      : 230629
#===============================================================================
import os
import csv
import sys

# Get the input folder path from command-line argument
input_folder = sys.argv[1]

# Get the list of CSV files in the input folder
csv_files = [file for file in os.listdir(input_folder) if file.endswith(".csv")]

# Initialize the concatenated data list
concatenated_data = []

# Iterate over each CSV file
for csv_file in csv_files:
    file_path = os.path.join(input_folder, csv_file)

    # Read the data from the current CSV file
    with open(file_path, "r", newline="") as csvfile:
        reader = csv.reader(csvfile)
        data = list(reader)

    # Exclude rows that contain "Reference"
    data = [row for row in data if "Reference" not in row]

    # Append the data to the concatenated data list
    concatenated_data.extend(data[1:])  # Exclude the header row for all files except the first

# Prepend the header row to the concatenated data
concatenated_data.insert(0, data[0])

# Define the output file path
output_file = os.path.join(input_folder, "concatenated_data.csv")

# Write the concatenated data to the output file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(concatenated_data)

print("Concatenated CSV file successfully generated.")
