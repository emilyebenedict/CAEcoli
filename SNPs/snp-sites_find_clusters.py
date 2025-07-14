# Author: Lindsey Hall, hall.l.r@wustl.edu
# This script generates clusters for snippy based on snp-dists output and a user-determined threshold

import pandas as pd
import os
import sys

def assign_core_snp_groups(edgelist_df):
    core_snp_groups = {}
    group_counter = 0
    
    for index, row in edgelist_df.iterrows():
        sample1 = row['sample1']
        sample2 = row['sample2']
        
        if sample1 in core_snp_groups:
            group = core_snp_groups[sample1]
            core_snp_groups[sample2] = group
        elif sample2 in core_snp_groups:
            group = core_snp_groups[sample2]
            core_snp_groups[sample1] = group
        else:
            group_counter += 1
            group = f"core_snp_group_{group_counter}"
            core_snp_groups[sample1] = group
            core_snp_groups[sample2] = group
    
    return core_snp_groups


def process_pairwise_core_snps(input_file, core_snp_cutoff):
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    # Remove rows with 'pseudo_reference_sequence' in any column
    df = df[~df.apply(lambda row: 'pseudo_reference_sequence' in row.values, axis=1)]
    
    df.columns = ['sample1', 'sample2', 'core_snp_dist']
    df_swapped = df[['sample2', 'sample1', 'core_snp_dist']].copy()
    df_combined = pd.concat([df, df_swapped])
    df_combined.drop_duplicates(inplace=True)
    df_combined['core_snp_dist'] = pd.to_numeric(df_combined['core_snp_dist'])
    
    df_filtered = df_combined[df_combined['core_snp_dist'] <= core_snp_cutoff]
    
    core_snp_groups = assign_core_snp_groups(df_filtered[['sample1', 'sample2']])
    
    df_groups = pd.DataFrame(core_snp_groups.items(), columns=['sample', 'core_snp_group_lte'])
    
    df_result = pd.merge(df_filtered, df_groups, left_on='sample1', right_on='sample', how='left')
    df_result.drop('sample', axis=1, inplace=True)
    
    output_dir = os.path.dirname(input_file)
    intermediate_output_file = os.path.join(output_dir, f"pairwise_core_snps_lte_{core_snp_cutoff}.tsv")
    df_result.to_csv(intermediate_output_file, sep='\t', index=False)
    
    unique_samples = df_result[['sample1', 'core_snp_group_lte']].drop_duplicates()
    reference_samples = unique_samples.groupby('core_snp_group_lte')['sample1'].apply(lambda x: sorted(x)[0]).reset_index()
    reference_samples.columns = ['core_snp_group_lte', 'reference']
    
    final_output = pd.merge(unique_samples, reference_samples, on='core_snp_group_lte')
    
    final_output.rename(columns={'sample1': 'sample_name'}, inplace=True)
    final_output_file = os.path.join(output_dir, "clusters_for_snippy.csv")
    
    # Remove groups with only one sample
    final_output = final_output[final_output.groupby('core_snp_group_lte')['sample_name'].transform('size') > 1]
    
    final_output[['sample_name', 'core_snp_group_lte', 'reference']].to_csv(final_output_file, sep=',', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <core_snp_cutoff>")
    else:
        input_file = sys.argv[1]
        core_snp_cutoff = int(sys.argv[2])
        process_pairwise_core_snps(input_file, core_snp_cutoff)
