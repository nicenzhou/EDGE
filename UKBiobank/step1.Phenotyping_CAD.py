##################################################################################################################################
## 																																
## 	Script Name: step1.Phenotyping_CAD.py														
## 	Description: Phenotype determination for UKB participants using Spark Python.
## 	Requirement: This has to be run on UK Biobank RAP using Spark Python. 				
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
## 		    This script was used to classify the individuals into CAD cases and controls using:
##	         ICD9, ICD10, OPCS4, Vascular_heart_problems_diagnosed_by_doctor, Non_cancer_illness_code, and Operation_code.
##        This script also contains codes that help to extract the first 20 PCs for step 2 to run a UMAP for further 
##           ancestry classification. A PCA plot was also generated for visualizing the PC1 and PC2. 
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:																														
## 			The first part of the codes for data dispensing is from UKB OpenBio:
## 			    https://github.com/dnanexus/OpenBio/blob/master/UKB_notebooks/ukb-rap-pheno-basic.ipynb
## 				
##################################################################################################################################


# ------------------------------------- #
## Python version check
# ------------------------------------- #
!python --version


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


# ------------------------------------- #
## Data dispensing at the individual level for selected entities
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
field_ids = ['31','22000','22001','22006','21022',
             '22009_a1','22009_a2','22009_a3','22009_a4','22009_a5','22009_a6','22009_a7','22009_a8','22009_a9','22009_a10',
             '22009_a11','22009_a12','22009_a13','22009_a14','22009_a15','22009_a16','22009_a17','22009_a18','22009_a19','22009_a20',
             '22019','41202','41204','41200','41272','41210','41203','41205','22021','21000_i0',
             '20002_i0','20004_i0','6150_i0',
             '20002_i1','20004_i1','6150_i1',
             '20002_i2','20004_i2','6150_i2']

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

# QC with sex, sex Chromosome aneuploidy, and kinship
df = participant_df
df_qced = df[
    (df['p31'] == df['p22001']) & # Filter in sex and genetic sex are the same
    (df['p22019'].isnull())     # Not Sex Chromosome aneuploidy
]
print(df_qced.shape) #(487288, 44)

# Annotation of the columns
df_qced = df_qced.rename(columns=
                         {'eid':'IID', 
                          'p31': 'sex',
                          'p22000': 'genotype_batch',
                          'p22001': 'genetic_sex',
                          'p22006': 'ethnic_group',
                          'p21022': 'age_baseline',
                          'p22009_a1': 'PC1',
                          'p22009_a2': 'PC2',
                          'p22009_a3': 'PC3',
                          'p22009_a4': 'PC4',
                          'p22009_a5': 'PC5',
                          'p22009_a6': 'PC6',
                          'p22009_a7': 'PC7',
                          'p22009_a8': 'PC8',
                          'p22009_a9': 'PC9',
                          'p22009_a10': 'PC10',
                          'p22009_a11': 'PC11',
                          'p22009_a12': 'PC12',
                          'p22009_a13': 'PC13',
                          'p22009_a14': 'PC14',
                          'p22009_a15': 'PC15',
                          'p22009_a16': 'PC16',
                          'p22009_a17': 'PC17',
                          'p22009_a18': 'PC18',
                          'p22009_a19': 'PC19',
                          'p22009_a20': 'PC20',
                          'p22019': 'sex_chromosome_aneuploidy',  
                          'p22021': 'kinship_to_other_participants',
                          'p21000_i0': 'self_reported_ethnic_group',
                          'p41202': 'ICD10_main',
                          'p41204': 'ICD10_secondary',
                          'p41200': 'OPCS4_main',
                          'p41272': 'OPCS4',
                          'p41210': 'OPCS4_secondary',
                          'p41203': 'ICD9_main',
                          'p41205': 'ICD9_secondary',
                          'p20002_i0': 'Non_cancer_illness_code_Instance_0',
                          'p20002_i1': 'Non_cancer_illness_code_Instance_1',
                          'p20002_i2': 'Non_cancer_illness_code_Instance_2',
                          'p20004_i0': 'Operation_code_Instance_0',
                          'p20004_i1': 'Operation_code_Instance_1',
                          'p20004_i2': 'Operation_code_Instance_2',
                          'p6150_i0': 'Vascular_heart_problems_diagnosed_by_doctor_Instance_0',
                          'p6150_i1': 'Vascular_heart_problems_diagnosed_by_doctor_Instance_1',
                          'p6150_i2': 'Vascular_heart_problems_diagnosed_by_doctor_Instance_2'})


# ------------------------------------- #
## Prepare the input table for the PCA plot
# ------------------------------------- #
# Add FID column 
df_qced['FID'] = df_qced['IID']
df_qced.dropna(subset=['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20'], how='any', inplace=True)

# Identify the individuals with null in self-reported ethnic group 
df_qced = df_qced.to_pandas()
import pandas as pd
import numpy as np
df_qced[["self_reported_ethnic_group"]] = df_qced[["self_reported_ethnic_group"]].fillna(0)
df_qced[["self_reported_ethnic_group"]] = df_qced[["self_reported_ethnic_group"]].applymap(str)

# Recode Race
df_qced = df_qced.replace( {
    "self_reported_ethnic_group": { 
        "1.0":"White",
        "1001.0":"White",
        "2001.0":"Mixed",
        "3001.0":"Asian or Asian British",
        "4001.0":"Black or Black British",
        "2.0":"Mixed",
        "1002.0":"White",
        "2002.0":"Mixed",
        "3002.0":"Asian or Asian British",
        "4002.0":"Black or Black British",
        "3.0":"Asian or Asian British",
        "1003.0":"White",
        "2003.0":"Mixed",
        "3003.0":"Asian or Asian British",
        "4003.0":"Black or Black British",
        "4.0":"Black or Black British",
        "2004.0":"Mixed",
        "3004.0":"Asian or Asian British",
        "5.0":"Chinese",
        "6.0":"Other ethnic group",
        "-1.0":"Do not know",
        "-3.0":"Prefer not to answer",
        "0.0":"Do not know"
        }
    } )

#Check ethnic groups
df_qced['self_reported_ethnic_group'].unique()

# write the merged data frame to a new CSV file
df_qced.to_csv('ukbb_info_full.csv', sep='\t', na_rep='NA', quoting=3)


# ------------------------------------- #
## Generate the PCA plot
# ------------------------------------- #
# Check the range of the first two PCs
print("mi: " , df_qced["PC1"].min() , " | max: ", df_qced["PC1"].max() )
print("mi: " , df_qced["PC2"].min() , " | max: ", df_qced["PC2"].max() )

# Plot PCA results
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi']= 300
sns.set_style('white')

race="self_reported_ethnic_group"
reorder_palette = ['#1f77b4', '#ff7f0e', '#2ca02c','#e377c2','#8c564b','#d62728','#9467bd','#7f7f7f'] 
sns.set_palette(palette=reorder_palette)

g = sns.jointplot(
    x = "PC1",
    y = "PC2",
    s = 3,
    hue = race,
#     marginal_kws={ "multiple":"fill" },
    marginal_kws={ "common_norm":False },
    ratio = 4,
    linewidth = 0,
    data = df_qced
    )

# Place legend outside of plot
sns.move_legend(g.ax_joint, "upper left", bbox_to_anchor=(1.25, 0.75), frameon=False)

# Export data as PNG
plt.savefig('ukbb_pca_full.png', bbox_inches='tight')

# Print counts by race and sex and first lines from dataframe 
print( df_qced.groupby(["sex", "self_reported_ethnic_group"])["self_reported_ethnic_group"].count() )

# Prepare PC tables for UMAP
ukbb_pcs = df_qced[['FID','IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20']]
ukbb_pcs = ukbb_pcs.set_index('IID')
ukbb_pcs.to_csv('ukbb_pcs_full.csv', sep='\t', na_rep='NA', quoting=3)


# ------------------------------------- #
# Identify the CAD cases
# ------------------------------------- #
# Create a phenotype table                
df_phenotype = df_qced[['FID','IID','sex','genotype_batch','age_baseline',
                        'ICD10_main','ICD10_secondary',
                        'OPCS4_main','OPCS4','OPCS4_secondary',
                        'ICD9_main','ICD9_secondary',
                        'Non_cancer_illness_code_Instance_0',
                        'Non_cancer_illness_code_Instance_1',
                        'Non_cancer_illness_code_Instance_2',
                        'Operation_code_Instance_0',
                        'Operation_code_Instance_1',
                        'Operation_code_Instance_2',
                        'Vascular_heart_problems_diagnosed_by_doctor_Instance_0',
                        'Vascular_heart_problems_diagnosed_by_doctor_Instance_1',
                        'Vascular_heart_problems_diagnosed_by_doctor_Instance_2']]
df_phenotype.to_csv('phenotype_full.csv', sep='\t', na_rep='NA', quoting=3)

# Using clinical and self-reported codes to identify the CAD cases
df = df_phenotype

# Define the codes dictionary to look up the CAD
icd10_codes = ['I210','I211','I212','I213','I214','I219','I220','I221','I228','I229','I231','I232','I233','I236','I214','I240','I241','I248','I249','I252']
icd9_codes = ['410','411','412']
opcs4_codes = ['K40','K41','K42','K43','K44','K45','K46','K49','K501','K75']
vhpdbd_codes = ['1','2']
#ncic_codes = ['1066','1074','1075','1076','1077','1471','1483','1484','1485','1486','1487','1078','1584','1585','1488','1489','1586','1587','1490','1588','1079','1080','1589','1590','1426']
ncic_codes = ['1074', '1075']
op_codes = ['1095','1523','1070']

def check_code(code, codes):
    code_str = str(code)
    for c in codes:
        if re.search(c, code_str):
            return True
    return False

# Extract the relevant columns
icd10_columns = [col for col in df.columns if 'ICD10' in col]
icd9_columns = [col for col in df.columns if 'ICD9' in col]
opcs4_columns = [col for col in df.columns if 'OPCS4' in col]
vhpdbd_columns = [col for col in df.columns if 'Vascular_heart_problems_diagnosed_by_doctor' in col]
ncic_columns = [col for col in df.columns if 'Non_cancer_illness_code' in col]
op_columns = [col for col in df.columns if 'Operation_code' in col]

relevant_columns = icd10_columns + icd9_columns + opcs4_columns + vhpdbd_columns + ncic_columns + op_columns

# Create a mask for the rows that contain at least one matching code
icd10_mask = df[icd10_columns].applymap(lambda x: check_code(x, icd10_codes)).any(axis=1)
icd9_mask = df[icd9_columns].applymap(lambda x: check_code(x, icd9_codes)).any(axis=1)
opcs4_mask = df[opcs4_columns].applymap(lambda x: check_code(x, opcs4_codes)).any(axis=1)
vhpdbd_mask = df[vhpdbd_columns].applymap(lambda x: check_code(x, vhpdbd_codes)).any(axis=1)
ncic_mask = df[ncic_columns].applymap(lambda x: check_code(x, ncic_codes)).any(axis=1)
op_mask = df[op_columns].applymap(lambda x: check_code(x, op_codes)).any(axis=1)

# Select the relevant rows from the data frame
icd10_df = df[icd10_mask]
icd9_df = df[icd9_mask]
opcs4_df = df[opcs4_mask]
vhpdbd_df = df[vhpdbd_mask]
ncic_df = df[ncic_mask]
op_df = df[op_mask]

# Add a column to indicate which codes have been identified
icd10_df['Source'] = 'ICD10'
icd9_df['Source'] = 'ICD9'
opcs4_df['Source'] = 'OPCS4'
vhpdbd_df['Source'] = 'VHPDBD'
ncic_df['Source'] = 'NCIC'
op_df['Source'] = 'OP'

# Merge the data frames
merged_df = pd.concat([icd10_df, icd9_df, opcs4_df, vhpdbd_df, ncic_df, op_df])
merged_df['Sources'] = merged_df.groupby(merged_df.index)['Source'].agg(lambda x: ', '.join(map(str, x)))
merged_df = merged_df.drop_duplicates(subset=['IID'])

# Select the final desired columns
output_df = merged_df[['FID', 'IID', 'sex', 'genotype_batch', 'age_baseline',
                       'ICD10_main', 'ICD10_secondary', 
                       'OPCS4_main', 'OPCS4', 'OPCS4_secondary', 
                       'ICD9_main', 'ICD9_secondary', 
                       'Non_cancer_illness_code_Instance_0',
                       'Non_cancer_illness_code_Instance_1',
                       'Non_cancer_illness_code_Instance_2',
                       'Operation_code_Instance_0',
                       'Operation_code_Instance_1',
                       'Operation_code_Instance_2',
                       'Vascular_heart_problems_diagnosed_by_doctor_Instance_0',
                       'Vascular_heart_problems_diagnosed_by_doctor_Instance_1',
                       'Vascular_heart_problems_diagnosed_by_doctor_Instance_2',
                       'Sources']]

def get_identified_codes(row, codes):
    identified_codes = []
    for col in row.index:
        if 'ICD10' in col:
            for code in codes['icd10']:
                if code in str(row[col]):
                    identified_codes.append(code)
        elif 'ICD9' in col:
            for code in codes['icd9']:
                if code in str(row[col]):
                    identified_codes.append(code)
        elif 'OPCS4' in col:
            for code in codes['opcs4']:
                if code in str(row[col]):
                    identified_codes.append(code)
        elif 'Non_cancer_illness_code' in col:
            for code in codes['ncic']:
                if code in str(row[col]):
                    identified_codes.append(code)
        elif 'Vascular_heart_problems_diagnosed_by_doctor' in col:
            for code in codes['vhpdbd']:
                if code in str(row[col]):
                    identified_codes.append(code)
        elif 'Operation_code' in col:
            for code in codes['op']:
                if code in str(row[col]):
                    identified_codes.append(code)
    return ','.join(identified_codes)

output_df['codes_identified'] = output_df.apply(lambda row: get_identified_codes(row, {'icd10': icd10_codes, 'icd9': icd9_codes, 'opcs4': opcs4_codes, 'vhpdbd': vhpdbd_codes, 'op':op_codes, 'ncic': ncic_codes}), axis=1)
output_df = output_df[['FID','IID','sex','age_baseline','genotype_batch','Sources','codes_identified']]

# Write the merged data frame 
output_df.to_csv('case_output_full.csv', sep='\t', na_rep='NA', quoting=3)

# Upload result files
!dx upload *.csv --path /{PATH}/ --brief
!dx upload *.png --path /{PATH}/ --brief
