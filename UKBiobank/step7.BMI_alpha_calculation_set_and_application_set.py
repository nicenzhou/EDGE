##################################################################################################################################
## 																																
## 	Script Name: step7.BMI_alpha_calculation_set_and_application_set.py
## 	Description: This script helps to extract the BMI and do a 50/50 split in a similar distribution for later alpha calculation and test.
## 	Requirement: This step needs the phenotype files from previous steps. Spark Python is required for the phenotype extraction.
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script contains steps to extract the BMI and to generate a training set with 50% of individuals 
##	        and a test set with the remaining individuals.
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
## Spark initialization 
# (Request more RAM for pyspark; Done only once; do not rerun this cell unless you select Kernel -> Restart kernel).
# ------------------------------------- #
from pyspark import SparkConf
conf = SparkConf().set("spark.kryoserializer.buffer.max", "2000m")
from pyspark.sql import SparkSession
spark = SparkSession.builder.config(conf=conf).getOrCreate()

# Import packages
import databricks.koalas as ks
import dxpy
import dxdata
import pandas as pd
import pyspark
import re

# Update seaborn and bokeh packages (need to restart the Kernel after)
!pip install seaborn --upgrade
!pip install bokeh
import matplotlib.pyplot as plt
import seaborn as sns


# ------------------------------------- #
## BMI phenotyping
# ------------------------------------- #
# Automatically discover dispensed database name and dataset id
dispensed_dataset = dxpy.find_one_data_object(
    typename='Dataset', 
    name='app*.dataset', 
    folder='/', 
    name_mode='glob')
dispensed_dataset_id = dispensed_dataset['id']
dataset = dxdata.load_dataset(id=dispensed_dataset_id)

# Select participants
participant = dataset['participant']

# Field ID
field_ids = ['21001_i0','21001_i1','21001_i2','21001_i3']

# Returns all field objects for a given UKB showcase field id
def fields_for_id(field_id):
    from distutils.version import LooseVersion
    field_id = str(field_id)
    fields = participant.find_fields(name_regex=r'^p{}(_i\d+)?(_a\d+)?$'.format(field_id))
    return sorted(fields, key=lambda f: LooseVersion(f.name))
  
fields = [participant.find_field(name='eid')] + [fields_for_id(f)[0] for f in field_ids]
field_description = pd.DataFrame({
    "Field": [f.name for f in fields],
    "Title": [f.title for f in fields],
    "Coding": [f.coding.codes if f.coding is not None else '' for f in fields]
})
field_description

participant_df = participant.retrieve_fields(fields = fields, engine=dxdata.connect()).to_koalas()
df = participant_df
df = df.rename(columns=
                         {'eid':'IID', 
                          'p21001_i0': 'BMI_Instance_0',
                          'p21001_i1': 'BMI_Instance_1',
                          'p21001_i2': 'BMI_Instance_2',
                          'p21001_i3': 'BMI_Instance_3'
                          })
df_phenotype = df.to_pandas()

df_phenotype.to_csv('BMI_forall.csv', sep='\t', na_rep='NA', quoting=3)

!dx upload BMI_forall.csv --path /{PATH}/ --brief


# ------------------------------------- #
## BMI training and test sets
# ------------------------------------- #
train_df = pd.read_csv("/{PATH}/train_df_{eur/afr}.csv", header=0, index_col=0)
test_df = pd.read_csv("/{PATH}//test_df_{eur/afr}.csv", header=0, index_col=0)
bmi_df = pd.read_csv("/{PATH}//BMI_forall.csv", header=0, index_col=0, sep="\t")

# Obtain the BMI at enrollment
bmi_df_at_enrollment = bmi_df[['IID', 'BMI_Instance_0']]

train_merged_df = pd.merge(train_df, bmi_df_at_enrollment, on='IID')
train_cleaned_df = train_merged_df.dropna()
print(train_cleaned_df)

test_merged_df = pd.merge(test_df, bmi_df_at_enrollment, on='IID')
test_cleaned_df = test_merged_df.dropna()
print(test_cleaned_df)


# ------------------------------------- #
## Distribution of BMI before and after
## log-transformation for skewness
# ------------------------------------- #
# Plot the distribution using Seaborn
sns.histplot(train_cleaned_df['BMI_Instance_0'], bins=30, kde=True, color='skyblue')
plt.title('Distribution of BMI_Instance_0 (train)')
plt.xlabel('BMI')
plt.ylabel('Frequency')
plt.show()

sns.histplot(test_cleaned_df['BMI_Instance_0'], bins=30, kde=True, color='skyblue')
plt.title('Distribution of BMI_Instance_0 (test)')
plt.xlabel('BMI')
plt.ylabel('Frequency')
plt.show()

sns.histplot(train_cleaned_df['log10bmi'], bins=30, kde=True, color='skyblue')
plt.title('Distribution of log10(BMI_Instance_0 + 1) (train)')
plt.xlabel('log10(BMI_Instance_0 + 1)')
plt.ylabel('Frequency')
plt.show()

sns.histplot(test_cleaned_df['log10bmi'], bins=30, kde=True, color='skyblue')
plt.title('Distribution of log10(BMI_Instance_0 + 1) (test)')
plt.xlabel('log10(BMI_Instance_0 + 1)')
plt.ylabel('Frequency')
plt.show()


# ------------------------------------- #
## Prepare the final tables
# ------------------------------------- #
train_cleaned_df = train_cleaned_df[['IID', 'log10bmi', 'age_baseline', 'sex', 'genotype_batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20']]
test_cleaned_df = test_cleaned_df[['IID', 'log10bmi', 'age_baseline', 'sex', 'genotype_batch', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20']]

train_cleaned_df.to_csv('train_df_bmi_{eur/afr}.csv', sep='\t', na_rep='NA', quoting=3)
test_cleaned_df.to_csv('test_df_bmi_{eur/afr}.csv', sep='\t', na_rep='NA', quoting=3)

!dx upload *.csv --path /{PATH}/ --brief
