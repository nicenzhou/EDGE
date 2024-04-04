##################################################################################################################################
## 																																
## 	Script Name: step9.Merge_results.py													
## 	Description: This script contains the codes to merge alpha and GWAS results for multiple phenotypes and ancestries.  
## 	Requirement: No.		
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      Codes allow merging the alpha and GWAS summary statistics for each SNP. It allows for multiple conditions at once.
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
##	      EDGE GWAS: 
##	        Zhou, J.  et al. Flexibly encoded GWAS identifies novel nonadditive SNPs in individuals of African and European ancestry.
##	        medRxiv 2023.06.01.23290857; doi: https://doi.org/10.1101/2023.06.01.23290857 
## 				
##################################################################################################################################

# ------------------------------------- #
## Merge EDGE alpha calculation
## AND
## EDGE GWAS (alpha application)
## For all chromosomes
## For all combinations of ancestries and phenotypes separately using one chunk of code
# ------------------------------------- #
import os
import pandas as pd

# Define the directories
directories = ['/mnt/project/Results/EUR_BMI/', '/mnt/project/Results/EUR_CAD/', '/mnt/project/Results/AFR_BMI/', '/mnt/project/Results/AFR_CAD/']

# Mapping of directory names to corresponding suffixes
directory_suffixes = {
    '/mnt/project/Results/EUR_BMI/': 'EUR_bmi',
    '/mnt/project/Results/EUR_CAD/': 'EUR_cad',
    '/mnt/project/Results/AFR_BMI/': 'AFR_bmi',
    '/mnt/project/Results/AFR_CAD/': 'AFR_cad'
}

# Loop through each directory
for directory in directories:
    # Get a list of all CSV files in the directory
    csv_files = [file for file in os.listdir(directory) if file.endswith('.csv')]

    # Initialize lists 
    alpha_files = []
    result_files = []

    # Iterate through each CSV file and append its data to the corresponding list
    for csv_file in csv_files:
        file_path = os.path.join(directory, csv_file)
        
        if "accumulated_alpha" in csv_file:
            alpha_df = pd.read_csv(file_path)
            alpha_files.append(alpha_df)
        elif "accumulated_result" in csv_file:
            result_df = pd.read_csv(file_path)
            result_files.append(result_df)

    # Merge alpha and GWAS summary first
    alpha_files = pd.concat(alpha_files, ignore_index=True)
    result_files = pd.concat(result_files, ignore_index=True)

    # Remove rows with NA
    alpha_files = alpha_files.dropna()
    result_files = result_files.dropna()

    # Save the merged files
    alpha_files.to_csv('chr_all_merged_alpha.csv', index=False)
    result_files.to_csv('chr_all_merged_result.csv', index=False)

    # Reload the merged files
    alpha_data = pd.read_csv('chr_all_merged_alpha.csv')
    result_data = pd.read_csv('chr_all_merged_result.csv')

    # Merge the two tables based on the Variant ID/SNPs columns
    final_merged_data = pd.merge(alpha_data, result_data, left_on='Variant.ID', right_on='snp', how='inner')

    # Remove rows with NA
    final_merged_data = final_merged_data.dropna()

    # Extract chromosome, position, reference, and alternate alleles from Variant.ID
    final_merged_data[['chr', 'pos', 'ref', 'alt']] = final_merged_data['Variant.ID'].str.split(':', expand=True)

    # Add 'chr' column to integers 
    final_merged_data['chr'] = final_merged_data['chr'].str.replace('chr', '').astype(int)

    # Sort the DataFrame first by chromosome and then by position within each chromosome
    # Allows to run LocusZoom if necessary
    final_merged_data.sort_values(by=['chr', 'pos'], inplace=True)

    # Drop intermediate columns
    final_merged_data.drop(columns=['chr', 'pos', 'ref', 'alt'], inplace=True)

    # Get the suffix for the output filename
    suffix = directory_suffixes.get(directory, 'unknown')

    # Construct the output filename
    output_filename = f'chr_all_final_merged_{suffix}.csv'
    
    # Save the final merged data 
    final_merged_data.to_csv(output_filename, index=False)

    # Print the first two rows as a preview
    print(f"Preview of merged data for {directory}:")
    print(final_merged_data.head(2))
    print()

    print(f"Merged data saved to {output_filename}")

print("All operations completed.")
