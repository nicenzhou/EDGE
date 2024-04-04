##################################################################################################################################
## 																																
## 	Script Name: step2.Apply_alpha_to_MVP_for_CAD_and_BMI.py													
## 	Description: It performs an EDGE GWAS of MVP cohort using a pre-calculated alpha from the UK Biobank dataset. 
## 	Requirement: This script needs a specific Python package called pandas_plink to read the genotypes.		
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script shows the way to apply calculated alphas to a new cohort, as running the EDGE GWAS. 
##	      Both phenotypes (CAD and log-transformed BMI) were considered and codes were set up as below.
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
##	        https://pandas-plink.readthedocs.io/en/latest/index.html
## 				
##################################################################################################################################

# ------------------------------------- #
## Binary outcome
## Logistic regression
## Example: CAD 
## An illustration using chromosome 1
# ------------------------------------- #
import pandas as pd
import numpy as np
import time
from pandas_plink import read_plink1_bin
import statsmodels.api as sm
from joblib import Parallel, delayed
from scipy.special import expit

G = read_plink1_bin("chr1.bed", 
                    "chr1.bim", 
                    "chr1.fam", verbose=True)

test_df = pd.read_csv("test_df.csv", header=0, index_col=0)
test_df.set_index('IID', inplace=True)
covariates_shared = pd.read_csv("covariates_shared_cad.txt", header=0, index_col=0, sep="\t")

outcome = "CAD"
covariates = ['age_baseline', 'sex', 'genotype_batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']

alpha_df = pd.read_csv('chr_all_final_merged_african_cad.csv', usecols=['Variant.ID', 'Alpha.Value'])

accumulated_info_df = pd.DataFrame(columns=['Variant.ID', 'Alpha.Value', 'Ref.Allele', 'Alt.Allele', 'EAF'])
accumulated_result_df = pd.DataFrame(columns=['snp', 'coef', 'std_err', 'z', 'P>|z|', 'conf_int_low', 'conf_int_high'])

samples = G.coords['sample'].values

def fit_logistic_regression(data, outcome, covariates):
    merged_df = pd.merge(data, covariates_shared, left_index=True, right_index=True)

    def fit_logistic_regression_internal(x_column):
        X_chunk = merged_df[[x_column] + covariates]
        y_chunk = merged_df[outcome]

        non_missing_mask = ~X_chunk.isnull().any(axis=1)
        X_chunk = X_chunk[non_missing_mask]
        y_chunk = y_chunk[non_missing_mask]

        X_chunk = sm.add_constant(X_chunk)

        model = sm.Logit(y_chunk, X_chunk)
        result = model.fit(method='bfgs', maxiter=1000)

        result_df = pd.DataFrame({
            'snp': x_column,
            'coef': result.params[x_column],
            'std_err': result.bse[x_column],
            'z': result.tvalues[x_column],
            'P>|z|': result.pvalues[x_column],
            'conf_int_low': result.conf_int().loc[x_column, 0],
            'conf_int_high': result.conf_int().loc[x_column, 1]
        }, index=[0])

        return result_df

    x_columns = merged_df.drop([outcome] + covariates, axis=1).columns

    results = Parallel(n_jobs=64)(delayed(fit_logistic_regression_internal)(x_column) for x_column in x_columns)

    result_df = pd.concat(results, ignore_index=True)
    return result_df

#start_number = 0  # replace with your desired starting number if you would like to perform the analysis within certain region(s)
#end_number = 20  # replace with your desired ending number if you would like to perform the analysis within certain region(s) 

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

skipped_snps = []

# Loop over SNPs
#for i in range(start_number, end_number + 1):
for i in range(len(G.coords['snp'].values)):
    variant_name = "variant" + str(i)
    test = G.sel(variant=variant_name).values
    snp = G.coords['snp'].values[i]

    df = pd.DataFrame(test, index=samples, columns=[snp])
    df.index.name = 'IID'
    df.index = df.index.astype('int64')
    column_chunks = [df.iloc[:, i:i+1] for i in range(0, df.shape[1])] 
        
    # Loop over chunks
    for chunk in column_chunks:
        for column in chunk.columns:
            identifier = snp
            alpha_value_row = alpha_df.loc[alpha_df['Variant.ID'] == identifier, 'Alpha.Value']
    
            if not alpha_value_row.empty and not np.isnan(alpha_value_row.iloc[0]):
                alpha_value = alpha_value_row.iloc[0]
            else:
                # Skip this SNP and move to the next one
                skipped_snps.append(snp)
                continue
            
            valid_values = test[~np.isnan(test)]
            num_zeros = np.sum(valid_values == 0)
            num_ones = np.sum(valid_values == 1)
            
            eaf = (2 * (num_zeros) + num_ones) / (2 * len(valid_values)) if len(valid_values) > 0 else np.nan
            
            # Varinat name
            variant_name = f'variant{i}'
            a0_values = G.a0.sel(variant=variant_name).values
            a1_values = G.a1.sel(variant=variant_name).values
            
            # Calculate alpha_df
            info_df = pd.DataFrame({
                    'Variant.ID': [identifier],
                    'Alpha.Value': [alpha_value],
                    'Ref.Allele': [a1_values],
                    'Alt.Allele': [a0_values],
                    'EAF': [eaf]
                    })

            accumulated_info_df = pd.concat([accumulated_info_df, info_df], ignore_index=True)
            
            test_geno_data = df.loc[df.index.isin(test_df.index)]
            test_geno_data_chunk = test_geno_data[chunk.columns].loc[test_geno_data.index.isin(test_df.index)]
            test_geno_data_chunk_r = test_geno_data_chunk.replace({identifier: {0: 1, 1: alpha_value, 2: 0}})

            test_result_df = fit_logistic_regression(test_geno_data_chunk_r, outcome, covariates)

            accumulated_result_df = pd.concat([accumulated_result_df, test_result_df[test_result_df['snp'] == identifier]], ignore_index=True)

            accumulated_info_df.to_csv('chr1_accumulated_info.csv', index=False)
            accumulated_result_df.to_csv('chr1_accumulated_result.csv', index=False)
            
# Write skipped SNPs into a log file
with open('chr1_skipped_snps.log', 'w') as f:
    f.write('\n'.join(skipped_snps))


# ------------------------------------- #
## Continuous outcome
## Linear regression
## Example: log-transformed BMI
## An illustration using chromosome 1
# ------------------------------------- #
import pandas as pd
import numpy as np
import time
from pandas_plink import read_plink1_bin
import statsmodels.api as sm
from joblib import Parallel, delayed
from scipy.special import expit

G = read_plink1_bin("chr1.bed", 
                    "chr1.bim", 
                    "chr1.fam", verbose=True)

test_df = pd.read_csv("test_df_bmi.csv", header=0, index_col=0, sep="\t")
test_df.set_index('IID', inplace=True)
covariates_shared = pd.read_csv("covariates_shared_bmi.csv", header=0, index_col=0)

outcome = "log10bmi"
covariates = ['age_baseline', 'sex', 'genotype_batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']

alpha_df = pd.read_csv('chr_all_final_merged_african_cad.csv', usecols=['Variant.ID', 'Alpha.Value'])

accumulated_info_df = pd.DataFrame(columns=['Variant.ID', 'Alpha.Value', 'Ref.Allele', 'Alt.Allele', 'EAF'])
accumulated_result_df = pd.DataFrame(columns=['snp', 'coef', 'std_err', 'z', 'P>|z|', 'conf_int_low', 'conf_int_high'])

samples = G.coords['sample'].values

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

#start_number = 0  # replace with your desired starting number if you would like to perform the analysis within certain region(s)
#end_number = 20  # replace with your desired ending number if you would like to perform the analysis within certain region(s) 

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

skipped_snps = []

# Loop over SNPs
#for i in range(start_number, end_number + 1):
for i in range(len(G.coords['snp'].values)):
    variant_name = "variant" + str(i)
    test = G.sel(variant=variant_name).values
    snp = G.coords['snp'].values[i]

    df = pd.DataFrame(test, index=samples, columns=[snp])
    df.index.name = 'IID'
    df.index = df.index.astype('int64')
    column_chunks = [df.iloc[:, i:i+1] for i in range(0, df.shape[1])] 
        
    # Loop over chunks
    for chunk in column_chunks:
        for column in chunk.columns:
            identifier = snp
            alpha_value_row = alpha_df.loc[alpha_df['Variant.ID'] == identifier, 'Alpha.Value']
    
            if not alpha_value_row.empty and not np.isnan(alpha_value_row.iloc[0]):
                alpha_value = alpha_value_row.iloc[0]
            else:
                # Skip this SNP and move to the next one
                skipped_snps.append(snp)
                continue
            
            valid_values = test[~np.isnan(test)]
            num_zeros = np.sum(valid_values == 0)
            num_ones = np.sum(valid_values == 1)

            eaf = (2 * (num_zeros) + num_ones) / (2 * len(valid_values)) if len(valid_values) > 0 else np.nan
            
            # Varinat name
            variant_name = f'variant{i}'
            a0_values = G.a0.sel(variant=variant_name).values
            a1_values = G.a1.sel(variant=variant_name).values
            
            # Calculate alpha_df
            info_df = pd.DataFrame({
                    'Variant.ID': [identifier],
                    'Alpha.Value': [alpha_value],
                    'Ref.Allele': [a1_values],
                    'Alt.Allele': [a0_values],
                    'EAF': [eaf]
                    })

            accumulated_info_df = pd.concat([accumulated_info_df, info_df], ignore_index=True)
            
            test_geno_data = df.loc[df.index.isin(test_df.index)]
            test_geno_data_chunk = test_geno_data[chunk.columns].loc[test_geno_data.index.isin(test_df.index)]
            test_geno_data_chunk_r = test_geno_data_chunk.replace({identifier: {0: 1, 1: alpha_value, 2: 0}})

            test_result_df = fit_linear_regression(test_geno_data_chunk_r, outcome, covariates)

            accumulated_result_df = pd.concat([accumulated_result_df, test_result_df[test_result_df['snp'] == identifier]], ignore_index=True)

            accumulated_info_df.to_csv('chr1_accumulated_info.csv', index=False)
            accumulated_result_df.to_csv('chr1_accumulated_result.csv', index=False)
            
# Write skipped SNPs into a log file
with open('chr1_skipped_snps.log', 'w') as f:
    f.write('\n'.join(skipped_snps))
