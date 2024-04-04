##################################################################################################################################
## 																																
## 	Script Name: step3.Determine_EUR_and_AFR_ancestry.py												
## 	Description: Ancestry assigned based on the UMAP results for UKB participants (for both cases and controls).
## 	Requirement: This script needs tables from previous steps as input. 				
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:
##	      This script helps to assign the ancestry groups based on the UMAP results.
##	      The EUR was selected based on the UMAP results using n_neighbors = 200, min_distance = 0.1.
##	        Individuals were selected as EUR using UMAP1 < 10 and UMAP2 < 15.
## 		    The AFR was selected based on the UMAP results using n_neighbors = 250, min_distance = 0 (more distinct groups).
##	        Individuals were selected as EUR using UMAP1 > 12.5 and UMAP2 > 5.5.
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:																														
## 		    NONE
## 				
##################################################################################################################################


# ------------------------------------- #
## Package installations and loading 
# ------------------------------------- #
# Update seaborn and bokeh packages (need to restart the Kernel after)
!pip install seaborn --upgrade
!pip install bokeh

# Import the libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# ------------------------------------- #
## Read info and UMAP results files
# ------------------------------------- #
case = pd.read_table("/{PATH}/case_output_full.csv",index_col=0)
phenotype = pd.read_table("/{PATH}/phenotype_full.csv",index_col=0)

umap_result_200_01 = pd.read_table("/{PATH}/proj_umap_pca_200_01.csv",sep="\t")
umap_result_250_00 = pd.read_table("/{PATH}/proj_umap_pca_250_00.csv",sep="\t")
ukbb_info = pd.read_table("/{PATH}/ukbb_info_full.csv",index_col=0)


# ------------------------------------- #
## Determine the EUR and AFR
# ------------------------------------- #
umap_result_eur = umap_result_200_01.query('UMAP1 < 10 & UMAP2 < 15')
umap_result_eur["IID"] = umap_result_eur["FID"]
umap_result_eur.rename(columns={'FID': 'FID'}, inplace=True)

umap_result_afr = umap_result_250_00.query('UMAP1 > 12.5 & UMAP2 > 5')
umap_result_afr["IID"] = umap_result_afr["FID"]
umap_result_afr.rename(columns={'FID': 'FID'}, inplace=True)

ukbb_full_eur = pd.merge(left=ukbb_info, right=umap_result_eur, how="left", on=['FID', 'IID'])
ukbb_full_eur = ukbb_full_eur[ukbb_full_eur['UMAP1'].notna()]
ukbb_full_afr = pd.merge(left=ukbb_info, right=umap_result_afr, how="left", on=['FID', 'IID'])
ukbb_full_afr = ukbb_full_afr[ukbb_full_afr['UMAP1'].notna()]


# ------------------------------------- #
## Plot UMAP results after selections
# ------------------------------------- #
plt.rcParams['figure.dpi']= 300
sns.set_style('white')

race="self_reported_ethnic_group"
reorder_palette = ['#8c564b','#ff7f0e','#9467bd','#e377c2','#7f7f7f','#1f77b4','#2ca02c']

sns.set_palette(palette=reorder_palette)

g = sns.jointplot(
    x = "UMAP1",
    y = "UMAP2",
    s = 3,
    hue = race,
    marginal_kws={ "common_norm":False },
    ratio = 4,
    linewidth = 0,
    data = ukbb_full_{eur/afr} # Input the specific ancestry group for visualization
    )

# Place legend outside of plot
sns.move_legend(g.ax_joint, "upper left", bbox_to_anchor=(1.25, 0.75), frameon=False)

# Export data as PNG
plt.savefig('umap_result_[eur/afr].png', bbox_inches='tight')

# Print counts by race and sex and first lines from data frame 
print( ukbb_full.groupby(["sex", "self_reported_ethnic_group"])["self_reported_ethnic_group"].count() )


# ------------------------------------- #
## Remove individuals who were dropped
# ------------------------------------- #
#  Remove the individuals who have been withdrawn from the study
list_removed = pd.read_table("/{PATH}/w52374_2023-04-25.csv",sep="\t",header = None)
list_removed = list_removed.rename(columns={0: "IID"})
 
iid_col = 'IID'
mask = ~umap_result_afr[iid_col].isin(list_removed[iid_col])
umap_result_afr = umap_result_afr[mask][["FID","IID"]]

# Remove the individuals who did not pass the phenotype QC
iid_col = 'IID'
mask = umap_result_{eur/afr}[iid_col].isin(phenotype[iid_col])
umap_result_{eur/afr} = umap_result_{eur/afr}[mask][["FID","IID"]]


# ------------------------------------- #
## Identify CAD cases and controls
# ------------------------------------- #
# Identify cases and controls for different ancestry groups
mask_case = umap_result_{eur/afr}[iid_col].isin(case[iid_col])
mask_control = ~umap_result_{eur/afr}[iid_col].isin(case[iid_col])

case_umap_{eur/afr} = umap_result_{eur/afr}[mask_case][["FID","IID"]]
control_umap_{eur/afr} = umap_result_{eur/afr}[mask_control][["FID","IID"]]

case_umap_{eur/afr}.loc[:, 'CAD'] = 1
control_umap_{eur/afr}.loc[:, 'CAD'] = 0

# Obtain the sex and genotype_batch for all selected individuals
umap_{eur/afr} = pd.concat([case_umap_{eur/afr}, control_umap_{eur/afr}])
Pheno = pd.merge(umap_{eur/afr},phenotype)[["FID","IID","CAD","sex","age_baseline","genotype_batch"]]

# write the Pheno dataframe to a new CSV file
Pheno.to_csv('pheno_nopcs_full_{eur/afr}.csv', sep='\t', na_rep='NA', index=False, quoting=3)

# Prepare a file with only FID and IID which will remind for later analyses
ind = Pheno[["FID","IID"]]
ind.to_csv('keep_full_{eur/afr}.txt', sep='\t', na_rep='NA', index=False, quoting=3)
ind.to_csv('keep_full_{eur/afr}_noheader.txt', sep='\t', na_rep='NA', index=False, header= False, quoting=3)

# Prepare a file with phenotype with PCs
Pheno_pcs = pd.merge(umap_{eur/afr},phenotype)[["IID","CAD","sex","age_baseline","genotype_batch"]]
ukbb_full_pcs_{eur/afr} = ukbb_full_{eur/afr}[["FID","IID",'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20']]
Pheno_pcs = pd.merge(Pheno_pcs,ukbb_full_pcs_{eur/afr})[["IID","CAD","sex","age_baseline","genotype_batch",'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20']]
Pheno_pcs.to_csv('pheno_pcs_full_{eur/afr}.txt', sep='\t', na_rep='NA', index=False, quoting=3)

# Upload all result files
!dx upload *.txt *.csv *.png --path /{PATH}/ --brief
