##################################################################################################################################
## 																																
## 	Script Name: step8.2.EDGE_GWAS_Linear_BMI.py													
## 	Description: EDGE alpha was calculated with the training set and EDGE GWAS was performed with the test set for a continuous outcome (BMI).
## 	Requirement: This script needs a specific Python package called pandas_plink to read the genotypes. 				
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script showed the EDGE alpha calculation and EDGE GWAS for applying the calculated alpha for a continuous outcome.
##	      log-transformed BMI as an example. 
##	      The following files will be generated as outputs:
##	        1) A file contains alpha values
##	        2) A file contains GWAS results
##	        3) A log file contains SNPs that are either not converged from the model or failed in the optimization
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
##	      EDGE GWAS: 
##	        Zhou, J. et al. Flexibly encoded GWAS identifies novel nonadditive SNPs in individuals of African and European ancestry.
##	        medRxiv 2023.06.01.23290857; doi: https://doi.org/10.1101/2023.06.01.23290857 
##	      EDGE GxG: 
##	        Hall, M. A. et al. Novel EDGE encoding method enhances ability to identify genetic interactions. PLoS Genetics 17, e1009534 (2021).
##	      Pandas-plink in python:
##	          https://pandas-plink.readthedocs.io/en/latest/index.html
## 				
##################################################################################################################################

# ------------------------------------- #
## EDGE alpha calculation
## AND
## EDGE GWAS (alpha application)
## An illustration using chromosome 1
# ------------------------------------- #
import pandas as pd
import numpy as np
import time
from pandas_plink import read_plink1_bin

G = read_plink1_bin("ukb21007_c1_bmi_{eur/afr}_hardcall8_qced_reordered.bed", 
                    "ukb21007_c1_bmi_{eur/afr}_hardcall8_qced_reordered.bim", 
                    "ukb21007_c1_bmi_{eur/afr}_hardcall8_qced_reordered.fam", verbose=True)

train_df = pd.read_csv("train_df_bmi_{eur/afr}.csv", header=0, index_col=0)
train_df.set_index('IID', inplace=True)
test_df = pd.read_csv("test_df_bmi_{eur/afr}.csv", header=0, index_col=0)
test_df.set_index('IID', inplace=True)
covariates_shared = pd.read_csv("covariates_shared_bmi_{eur/afr}.txt", header=0, index_col=0, sep="\t")

import statsmodels.api as sm
import joblib
from joblib import Parallel, delayed

outcome = "log10bmi"
covariates = ['age_baseline', 'sex', 'genotype_batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']

# Initialize empty DataFrames to accumulate results
accumulated_alpha_df = pd.DataFrame(columns=['Variant.ID', 'Alpha.Value', 'Ref.Allele', 'Alt.Allele', 'EAF'])
accumulated_result_df = pd.DataFrame(columns=['snp', 'coef', 'std_err', 't', 'P>|t|', 'conf_int_low', 'conf_int_high'])

samples = G.coords['sample'].values

def fit_linear_regression_codom(het, hom, outcome_column, covariates):
    # Merge the data with covariates_shared
    data = pd.merge(het, hom, left_index=True, right_index=True, suffixes=('_het', '_hom'))
    merged_df = pd.merge(data, covariates_shared, left_index=True, right_index=True)
    merged_df = merged_df.dropna()
    
    # Use the provided outcome column
    outcome_column = merged_df[outcome]

    # Function to fit linear regression on a subset of data
    def fit_linear_regression_codom_internal(x_column):
        # Adjust the X_chunk to include both _het and _hom columns
        X_chunk = merged_df[[x_column + '_het', x_column + '_hom'] + covariates]
        y_chunk = outcome_column

        non_missing_mask = ~X_chunk.isnull().any(axis=1)
        X_chunk = X_chunk[non_missing_mask]
        y_chunk = y_chunk[non_missing_mask]

        X_chunk = sm.add_constant(X_chunk)

        model = sm.OLS(y_chunk, X_chunk)
        result = model.fit()

        # Extracting the coefficients, p-values, and confidence intervals only for the identifier
        result_df = pd.DataFrame({
            'snp': x_column,
            'coef_het': result.params[x_column + '_het'],
            'coef_hom': result.params[x_column + '_hom'],
            'std_err_het': result.bse[x_column + '_het'],
            'std_err_hom': result.bse[x_column + '_hom'],
            't_het': result.tvalues[x_column + '_het'],
            't_hom': result.tvalues[x_column + '_hom'],
            'P>|t|_het': result.pvalues[x_column + '_het'],
            'P>|t|_hom': result.pvalues[x_column + '_hom'],
            'conf_int_low_het': result.conf_int().loc[x_column + '_het', 0],
            'conf_int_high_het': result.conf_int().loc[x_column + '_het', 1],
            'conf_int_low_hom': result.conf_int().loc[x_column + '_hom', 0],
            'conf_int_high_hom': result.conf_int().loc[x_column + '_hom', 1],
        }, index=[0])

        return result_df

    # Get the list of X variables
    x_columns = merged_df.filter(regex='_het$').columns.str.replace('_het', '').tolist()

    # Parallelize linear regression fitting for each X variable
    results = Parallel(n_jobs=64)(delayed(fit_linear_regression_codom_internal)(x_column) for x_column in x_columns)

    # Concatenate the results into a single DataFrame
    result_df = pd.concat(results, ignore_index=True)

    return result_df


# Define a function for linear regression
def fit_linear_regression(data, outcome, covariates):
    # Merge the data with covariates_shared
    merged_df = pd.merge(data, covariates_shared, left_index=True, right_index=True)
    merged_df = merged_df.dropna()

    # Use the provided outcome column
    outcome_column = merged_df[outcome]

    # Function to fit linear regression on a subset of data
    def fit_linear_regression_internal(x_column):
        X_chunk = merged_df[[x_column] + covariates]
        y_chunk = outcome_column

        non_missing_mask = ~X_chunk.isnull().any(axis=1)
        X_chunk = X_chunk[non_missing_mask]
        y_chunk = y_chunk[non_missing_mask]

        X_chunk = sm.add_constant(X_chunk)

        model = sm.OLS(y_chunk, X_chunk)
        result = model.fit()

        # Extracting the coefficients, p-values, and confidence intervals only for the identifier
        result_df = pd.DataFrame({
            'snp': x_column,
            'coef': result.params[x_column],
            'std_err': result.bse[x_column],
            't': result.tvalues[x_column],
            'P>|t|': result.pvalues[x_column],
            'conf_int_low': result.conf_int().loc[x_column, 0],
            'conf_int_high': result.conf_int().loc[x_column, 1]
        }, index=[0])

        return result_df

    # Get the list of X variables
    x_columns = merged_df.drop([outcome] + covariates, axis=1).columns

    # Parallelize linear regression fitting for each X variable
    results = Parallel(n_jobs=64)(delayed(fit_linear_regression_internal)(x_column) for x_column in x_columns)

    # Concatenate the results into a single DataFrame
    result_df = pd.concat(results, ignore_index=True)

    return result_df

# Iterate over each column chunk in df
# Specify the path pattern for the CSV files (adjust as needed)
#csv_path_pattern = 'chr9_{}.csv'

# Specify the number of chunks (adjust based on the total number of chunks)
#num_chunks = len(snp) // chunk_size + 1

# Loop through the CSV files and apply "Later code" to each DataFrame
#for i in range(num_chunks):
#    # Read the CSV file
#    df = pd.read_csv(csv_path_pattern.format(i), index_col='IID')
#start_number = 20001  # replace with your desired starting number
#end_number = len(G.coords['snp'].values)  # replace with your desired ending number

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

skipped_snps = []

# Loop over SNPs
#for i in range(start_number, end_number + 1):
for i in range(len(G.coords['snp'].values)):
    variant_name = "variant" + str(i)

    try:
        test = G.sel(variant=variant_name).values
    except KeyError:
        print(f"Variant {variant_name} not found. Stopping the loop.")
        break
        
    snp = G.coords['snp'].values[i]

    df = pd.DataFrame(test, index=samples, columns=[snp])
    df.index.name = 'IID'
    df.index = df.index.astype('int64')
    column_chunks = [df.iloc[:, i:i+1] for i in range(0, df.shape[1])]

    # Loop over chunks
    for chunk in column_chunks:
        for column in chunk.columns:
            het = chunk[column].replace(2, 0)
            hom = chunk[column].replace({0:1, 1: 0, 2: 0})

            train_het_data = het.loc[het.index.isin(train_df.index)]
            train_hom_data = hom.loc[hom.index.isin(train_df.index)]
            test_geno_data = df.loc[hom.index.isin(test_df.index)]
            test_geno_data_chunk = test_geno_data[chunk.columns].loc[test_geno_data.index.isin(test_df.index)]

            try:
                train_result_df = fit_linear_regression_codom(train_het_data, train_hom_data, outcome, covariates)
            except (np.linalg.LinAlgError, ConvergenceWarning):
                print(f"Skipping SNP {snp} due to Singular matrix or ConvergenceWarning error")
                skipped_snps.append(snp)
                continue
            
            valid_values = test[~np.isnan(test)]
            num_zeros = np.sum(valid_values == 0)

            eaf = 2 * (num_zeros) / (2 * len(valid_values)) if len(valid_values) > 0 else np.nan
            
            # Varinat name
            variant_name = f'variant{i}'
            a0_values = G.a0.sel(variant=variant_name).values
            a1_values = G.a1.sel(variant=variant_name).values
            
            # Calculate alpha_df
            alpha_df = pd.DataFrame({
                'Variant.ID': train_result_df['snp'],
                'Alpha.Value': train_result_df['coef_het'] / train_result_df['coef_hom'] if not train_result_df.empty else np.nan,
                'Ref.Allele': a1_values,
                'Alt.Allele': a0_values,
                'EAF': eaf
            })

            accumulated_alpha_df = pd.concat([accumulated_alpha_df, alpha_df], ignore_index=True)

            identifier = alpha_df['Variant.ID'].iloc[0]
            alpha_value = alpha_df['Alpha.Value'].iloc[0] if not alpha_df.empty else np.nan

            test_geno_data_chunk_r = test_geno_data_chunk.replace({identifier: {0: 1, 1: alpha_value, 2: 0}})

            test_result_df = fit_linear_regression(test_geno_data_chunk_r, outcome, covariates)

            accumulated_result_df = pd.concat([accumulated_result_df, test_result_df[test_result_df['snp'] == identifier]], ignore_index=True)

            accumulated_alpha_df.to_csv('chr1_accumulated_alpha.csv', index=False)
            accumulated_result_df.to_csv('chr1_accumulated_result.csv', index=False)
            
# Write skipped SNPs into a log file
with open('chr1_skipped_snps.log', 'w') as f:
    f.write('\n'.join(skipped_snps))
