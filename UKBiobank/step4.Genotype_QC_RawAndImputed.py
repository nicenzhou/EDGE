##################################################################################################################################
## 																																
## 	Script Name: step4.Genotype_QC_RawAndImputed.py									
## 	Description: This script contains codes for QCing the raw genotypes (bim/bed/fam) and imputed genotypes (BGEN).
## 	Requirement: PLINK1.9 and/or PLINK2 needs for running the genotype QC.				
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      Genotype QC for both raw genotype calls and TOPMed imputed genotype data.
##	      We recommend only selecting individuals who have the TOPMed imputed genotype data, 
##	        since later analyses will be done with that TOPMed imputed data only.
##	      The raw genotype QC is required for PC estimation for the AFR population. 
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
##	      PLINK1.9 and PLINK2 installation and usage:
##	        https://www.cog-genomics.org/plink/1.9/
##	        https://www.cog-genomics.org/plink/2.0/
## 				
##################################################################################################################################

# ------------------------------------- #
## Package installations and loading 
# ------------------------------------- #
#IF RUN PLINK1.9/PLINK2 within Python Notebook
!wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240318.zip # replace with the version you preferred
!unzip plink2_linux_x86_64_20240318.zip

import os
import pandas as pd


# ------------------------------------- #
## Revise the BGEN sample files
# ------------------------------------- #
!cp /mnt/project/Bulk/Imputation/Imputation\ from\ genotype\ \(TOPmed\)/*.sample . ; 

folder_path = '.'  # Replace with the actual folder path

# Iterate over files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.sample'):  # Process only files with the ".sample" extension
        file_path = os.path.join(folder_path, filename)
        
        # Read sample file
        sample_file = pd.read_csv(file_path, sep=' ')
        
        # Fix the datatype for sex
        sample_file['sex'] = sample_file['sex'].astype(str)
        sample_file.at[0, 'sex'] = 'D'
        
        # Write out the modified data frame as a space-delimited file without including the row indices
        sample_file.to_csv(file_path, sep=' ', index=False, encoding='utf-8-sig', header=True)

!dx upload *.sample --path /{PATH}/ --brief


# ------------------------------------- #
## QC: BGEN files
# ------------------------------------- #
# Run through chr1 to chr22
plink2 --threads 32 
         --bgen ukb21007_*.bgen ref-first 
         --sample ukb21007_*.sample 
         --keep keep_full_afr_noheader.txt 
         --autosome --minimac3-r2-filter 0.3 
         --hwe 10e-6 midp
         --set-missing-var-ids @:#:\$r:\$a 
         --new-id-max-allele-len 100 
         --rm-dup 
         --export bgen-1.2 "bits=8"
         --out ukb21007_c*_qced; 

# IF you forget to input "bits=8" for UKB, you will need to regenerate the BGEN file with 8 bits using qctool which has been pre-installed in SAK.
# For example,
# qctool -g ukb21007_c1_qced.bgen -bgen-bits 8 -og ukb21007_c1_qced_8bits.bgen;

# BEGN to plink file 
plink2 --threads 32 
         --bgen ukb21007_*.bgen ref-first 
         --sample ukb21007_*.sample 
         --keep keep_full_afr_noheader.txt
         --autosome --minimac3-r2-filter 0.3 
         --hwe 10e-6 midp 
         --set-missing-var-ids @:#:\$r:\$a 
         --new-id-max-allele-len 100 
         --rm-dup  
         --make-bed 
         --out ukb21007_c1_qced;

# Obtain a list of individuals who remained after the QC
!cp /mnt/project/Geno/bgensample/bgen_to_bfile/ukb21007_c1_qced.fam .
import pandas as pd
fam_file = pd.read_csv("ukb21007_c1_qced.fam", sep='\t', header=None)

keep_file = fam_file.iloc[:,:2]
keep_file.to_csv('keep_afr_frombgen_noheader.txt', header=False, sep='\t', na_rep='NA', index=False, quoting=3)
!dx upload *.txt --path /{PATH}/ --brief


# ------------------------------------- #
## Genotype calls QC
# ------------------------------------- #
# Merge genotype calls
cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ; 
ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt;
plink --threads 32 
         --merge-list files_to_merge.txt 
         --keep keep_afr_frombgen_noheader.txt
         --autosome-xy 
         --make-bed 
         --out geno_merged_afr_full_r; 
         
rm ukb22418_c[1-9]*; 
rm *.txt;
