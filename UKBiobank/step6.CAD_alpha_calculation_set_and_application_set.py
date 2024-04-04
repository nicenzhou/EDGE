##################################################################################################################################
## 																																
## 	Script Name: step6.CAD_alpha_calculation_set_and_application_set.py
## 	Description: This script helps to do a 50/50 split in a similar distribution for later alpha calculation and test.
## 	Requirement: This step needs the phenotype files from previous steps.
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script contains steps to generate a training set with 50% of individuals and a test set with the remaining individuals.
##	      The distribution of the phenotype is desired to be similar between the training set and the test set. 
##	      The training set is for EDGE alpha calculation and the test set is for EDGE alpha application. 
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
##	      EDGE algorithm in GxG: 
##	        Hall, M. A. et al. Novel EDGE encoding method enhances ability to identify genetic interactions. PLoS Genetics 17, e1009534 (2021).
## 				
##################################################################################################################################

# ------------------------------------- #
## Package installations and loading
# ------------------------------------- #
!pip install scikit-learn

import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import train_test_split


# ------------------------------------- #
## AFR CAD
# ------------------------------------- #
pheno = pd.read_csv("/{PATH}/pheno_pcairpc_afr.txt", delim_whitespace=True) 

# Separate features (X) and target variable (y)
X = pheno.drop(columns=['CAD'])  # Features
y = pheno['CAD']  # Target variable

# Split the data into training and test sets (50% training, 50% test) with stratification
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42, stratify=y)

# Create DataFrames for training and testing
train_df = pd.concat([X_train, y_train], axis=1)
test_df = pd.concat([X_test, y_test], axis=1)

train_df.to_csv('train_df_afr.csv', index=True, header=True)
test_df.to_csv('test_df_afr.csv', index=True, header=True)

!dx upload *.csv --path /{PATH}/

# ------------------------------------- #
## EUR CAD
# ------------------------------------- #
pheno = pd.read_csv("/{PATH}/pheno_pcs_full_eur.txt", delim_whitespace=True) 

# Separate features (X) and target variable (y)
X = pheno.drop(columns=['CAD'])  # Features
y = pheno['CAD']  # Target variable

# Split the data into training and test sets (50% training, 50% test) with stratification
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42, stratify=y)

# Create DataFrames for training and testing
train_df = pd.concat([X_train, y_train], axis=1)
test_df = pd.concat([X_test, y_test], axis=1)

train_df.to_csv('train_df_eur.csv', index=True, header=True)
test_df.to_csv('test_df_eur.csv', index=True, header=True)

!dx upload *.csv --path /{PATH}/
