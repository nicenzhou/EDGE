##################################################################################################################################
## 																																
## 	Script Name: step1.Identify_shared_SNPs_with_UKB.R													
## 	Description: Identify the shared SNPs between UKB and MVP data. 
## 	Requirement: Packages in base R are required. 				
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script allows identifying the shared SNPs between UKB and MVP. 
##	      EDGE GWAS requires a pre-calculated alpha to recode heterozygous genotypes. 
##	      Only SNPs with a pre-calculated alpha need to be considered for running through the EDGE GWAS. 
##	      The shared SNPs between UKB and MVP were extracted to apply the alpha from the UKB to the MVP cohort.
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
##	      EDGE algorithm in GxG: 
##	        Hall, M. A. et al. Novel EDGE encoding method enhances ability to identify genetic interactions. PLoS Genetics 17, e1009534 (2021).
## 				
##################################################################################################################################

# ------------------------------------- #
## Extract the shared SNPs
# ------------------------------------- #
