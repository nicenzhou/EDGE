##################################################################################################################################
## 																																
## 	Script Name: step2.UMAP_analysis.py														
## 	Description: UMAP analysis for ancestry re-classifications, especially to obtain European- and African-ancestry populations.
## 	Requirement: This script needs tables from step1. umap-learn Python package is required. No need for the spark Python. 				
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
## 		    This script used UMAP to re-classify individuals into different ancestry groups based on their PCs from genetics.
##	      UMAP is recommended to run on everyone. 
##	      Using the Data-Field 22009 as input for UMAP: https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=22009
##	      All individuals were included: 19,525,280 items of data are available, covering 488,132 participants.
##	      Array indices run from 1 to 40.
##	      Units of measurement are Units.
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
## 			UMAP-learn in python:
## 			    https://umap-learn.readthedocs.io/en/latest/index.html
## 			Supplementary Materials S3.3.2 Details of PCA for release:
## 			    https://www.biorxiv.org/content/10.1101/166298v1
## 				
##################################################################################################################################

# ------------------------------------- #
## Package installations and loading 
# ------------------------------------- #
# Update seaborn and bokeh packages (need to restart the Kernel after)
!pip install seaborn --upgrade
!pip install bokeh

# Install umap-learn package (need to restart the Kernel after)
!conda install -c conda-forge umap-learn -y

# Import libraries
%matplotlib inline

import matplotlib.pyplot as plt
import collections
from collections import defaultdict
import gzip
import itertools
import numpy as np
import pandas as pd
import os
import time
import seaborn as sns

from ipywidgets import interact
import bokeh
import bokeh.io
from bokeh.io import push_notebook
from bokeh.plotting import figure, show, save, output_notebook, output_file

# Import colour palettes for later on
from bokeh.palettes import Category20b
from bokeh.palettes import Purples
from bokeh.palettes import Greens
from bokeh.palettes import YlOrBr
from bokeh.palettes import YlOrRd
from bokeh.palettes import PuOr
from bokeh.palettes import RdGy


# ------------------------------------- #
## Run UMAP
# ------------------------------------- #
# Read PC table
import pandas as pd
ukbb_pcs = pd.read_table(
    "/{PATH}/ukbb_pcs_full.csv",
    index_col="IID")

# Run UMAP and tSNE on PCs
import umap.umap_ as umap

# Number of principal components to use
n_pc = 20

# PC table 
proj_pca = ukbb_pcs.drop(['FID'], axis=1).to_numpy()

# Project the principal components via UMAP to 2 dimensions
umap_name =['UMAP1','UMAP2']
proj_umap_pca = umap.UMAP(n_components=2, n_neighbors=250, min_dist=0).fit_transform(proj_pca[:,:n_pc])
proj_umap_pca_df = pd.DataFrame(proj_umap_pca_20_15_05, columns = umap_name, index=ukbb_pcs['FID'])
proj_umap_pca_df.to_csv('proj_umap_pca_20_500_00.csv', sep='\t', na_rep='NA', quoting=3)


# ------------------------------------- #
## Visualize the UMAP results
# ------------------------------------- #
# Obtain self-reported ethnic group
ukbb_info = pd.read_table("/{PATH}/ukbb_info_full.csv",index_col="FID")

ukbb_full = pd.concat( [ukbb_info, proj_umap_pca_df], axis=1 )

ukbb_full

plt.rcParams['figure.dpi']= 300
sns.set_style('white')

race="self_reported_ethnic_group"
reorder_palette = ['#1f77b4', '#ff7f0e', '#2ca02c','#e377c2','#8c564b','#d62728','#9467bd','#7f7f7f'] 
sns.set_palette(palette=reorder_palette)

g = sns.jointplot(
    x = "UMAP1",
    y = "UMAP2",
    s = 3,
    hue = race,
#     marginal_kws={ "multiple":"fill" },
    marginal_kws={ "common_norm":False },
    ratio = 4,
    linewidth = 0,
    data = ukbb_full
    )

# Place legend outside of plot
sns.move_legend(g.ax_joint, "upper left", bbox_to_anchor=(1.25, 0.75), frameon=False)

# Export data as PNG
plt.savefig('proj_umap_pca_20_500_00.png', bbox_inches='tight')

# Print counts by race and sex and first lines from data frame 
print( ukbb_full.groupby(["sex", "self_reported_ethnic_group"])["self_reported_ethnic_group"].count() )

# Upload all results
!dx upload *.csv *.png --path /{PATH}/ --brief
