##################################################################################################################################
## 																																
## 	Script Name: step4.3.Phenotype_covariates_AFR.py										
## 	Description: This script prepares the CAD phenotype with PCs from PC-AiR for the AFR population.
## 	Requirement: The phenotype file and PC file from previous steps. 			
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script contains steps for merging the AFR phenotype file with the PCs from PC-AiR.
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
##	      NONE
## 				
##################################################################################################################################

# ------------------------------------- #
## Create a pheno file for AFR
# ------------------------------------- #
import pandas as pd
Pheno = pd.read_table("/{PATH}/pheno_nopcs_full_afr.txt",sep="\t")

pc_pcair = pd.read_table("/{PATH}/afr_cad_pcair_pcs.txt",sep="\t",header=None)
pc_pcair = pc_pcair.rename(columns={0: "IID", 1: "PC1", 2: "PC2", 3: "PC3", 4: "PC4", 5: "PC5", 6: "PC6", 7: "PC7", 8: "PC8", 9: "PC9", 
                       10: "PC10", 11: "PC11", 12: "PC12", 13: "PC13", 14: "PC14", 15: "PC15", 16: "PC16", 17: "PC17", 18: "PC18", 19: "PC19", 
                       20: "PC20", 21: "PC21", 22: "PC22", 23: "PC23", 24: "PC24", 25: "PC25", 26: "PC26", 27: "PC27", 28: "PC28", 29: "PC29", 
                       30: "PC30", 31: "PC31", 32: "PC32" })

Pheno_pcs_afr = pd.merge(left=Pheno, right=pc_pcair, how="left", on=['IID'])[["IID","CAD","sex","age_baseline","genotype_batch",'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']]
Pheno_pcs_afr = Pheno_pcs_afr[Pheno_pcs_afr['PC1'].notna()]
Pheno_pcs_afr.to_csv('pheno_pcairpc_afr.txt', sep='\t', na_rep='NA', index=False, header=True, quoting=3)
!dx upload pheno_pcairpc_afr.txt --path /{PATH}/ --brief
