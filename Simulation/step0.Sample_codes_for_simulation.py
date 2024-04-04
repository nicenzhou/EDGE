##################################################################################################################################
## 																																
## 	Script Name: step0.Sample_codes_for_simulation.py			
## 	Description: It allows performing the simulations by generating the desired SNPs with different combinations of parameters. 
## 	              Inheritance patterns (recessive, sub-additive, additive, super-additive, dominant, and DOMDEV) are also tested.
## 	              Also the simulation of covariates is implemented within the codes for both binary and continuous.
## 	Requirement: This script needs concurrent.futures, clarite, and pandas-genomics python packages.
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      The codes help to generate the SNPs with different combinations of parameters for alpha calculation and application.
##	      The EDGE algorithm was applied to calculate the alpha to the test set only.
##	      The traditional encodings (recessive, sub-additive, additive, super-additive, and dominant) were also tested.
##	      The most recent encoding (DOMDEV) was also tested.
##	      Binary phenotype was considered as below.
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
##	      Pandas-genomic in python:
##	        https://pandas-genomics.readthedocs.io/en/latest/
##	      CLARITE in python:
##	        https://clarite-python.readthedocs.io/en/latest/
## 				
##################################################################################################################################

# ------------------------------------- #
## Packages loading
# ------------------------------------- #
import concurrent.futures
import pandas as pd
import numpy as np
import clarite
import time
from pandas_genomics import sim, scalars


# ------------------------------------- #
## Define inheritance patterns
## Simulations
## GWAS with different encodings
# ------------------------------------- #
# Name conversation for parameters
def conv_u(i):
    switcher = {
        "REC": "RECESSIVE",
        "SUB": "SUB_ADDITIVE",
        "ADD": "ADDITIVE",
        "SUP": "SUPER_ADDITIVE",
        "DOM": "DOMINANT",
        "HET": "HET",
        "NUL": "ADDITIVE",
    }
    return switcher.get(i, "Invalid Group")


def conv_l(i):
    switcher = {
        "REC": "Recessive",
        "SUB": "Sub-Additive",
        "ADD": "Additive",
        "SUP": "Super_Additive",
        "DOM": "Dominant",
        "HET": "Heterozygous",
        "NUL": "NULL",
    }
    return switcher.get(i, "Invalid Group")


def simulations(seed):

    train_seed = seed
    test_seed = seed + 2000

    ALL_RESULTS_ENCODING = pd.DataFrame()
    ALL_RESULTS_EDGE_ALPHA = pd.DataFrame()
    #ALL_RESULTS_ENCODING_EDGE = pd.DataFrame()

    for ab1, ab2 in work_group:

        ab1u = conv_u(ab1)
        ab2u = conv_u(ab2)
        ab1l = conv_l(ab1)

        # BLOCK 1
        # Main Effect for SNP1 without interaction
        # Training data
        train_main = sim.BAMS.from_model(
            eff1=getattr(sim.SNPEffectEncodings, ab1u),
            eff2=getattr(sim.SNPEffectEncodings, ab2u),
            penetrance_base=PEN_BASE,
            penetrance_diff=PEN_DIFF,
            main1=1,
            main2=0,
            interaction=0,
            snp1=variant1,
            snp2=variant2,
            random_seed=train_seed)

        # BLOCK 2
        train_me = train_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB)
        train_me = train_me.sort_values(by="Outcome",ascending=False)
        
        np.random.seed(train_seed)
        Age_1 = pd.DataFrame(np.random.poisson(lam=31, size=(int(0.507*n_controls))).astype(int))
        Age_2 = pd.DataFrame(np.random.poisson(lam=54.5, size=(int(0.316*n_controls))).astype(int))
        Age_3 = pd.DataFrame(np.random.poisson(lam=82.5, size=(n_controls-int(0.507*n_controls)-int(0.316*n_controls))).astype(int))
        Age_4 = pd.DataFrame(np.random.poisson(lam=31, size=(int(0.144*n_cases))).astype(int))
        Age_5 = pd.DataFrame(np.random.poisson(lam=54.5, size=(int(0.435*n_cases))).astype(int))
        Age_6 = pd.DataFrame(np.random.poisson(lam=82.5, size=(n_cases-int(0.144*n_cases)-int(0.435*n_cases))).astype(int))
        Age = pd.concat([Age_1,Age_2,Age_3,Age_4,Age_5,Age_6], axis=0, ignore_index=True, sort=False)
        BMI_1 = pd.DataFrame(np.random.normal(loc=26.5, scale=1.0, size=(n_controls)).astype(int))
        BMI_2 = pd.DataFrame(np.random.normal(loc=35, scale=7.5, size=(n_cases)).astype(int))
        BMI = pd.concat([BMI_1,BMI_2], axis=0, ignore_index=True, sort=False)
        Sex_1 = pd.DataFrame(np.random.binomial(n=1, p=0.481, size=(n_controls)).astype(int))
        Sex_2 = pd.DataFrame(np.random.binomial(n=1, p=0.525, size=(n_cases)).astype(int))
        Sex = pd.concat([Sex_1,Sex_2], axis=0, ignore_index=True, sort=False)
        Smoking_1 = pd.DataFrame(np.random.binomial(n=1, p=0.125, size=(n_controls)).astype(int))
        Smoking_2 = pd.DataFrame(np.random.binomial(n=1, p=0.64, size=(n_cases)).astype(int))
        Smoking = pd.concat([Smoking_1,Smoking_2], axis=0, ignore_index=True, sort=False) 
        train_COV = pd.concat([Age,BMI,Sex,Smoking], axis=1, ignore_index=True, sort=False)
        train_cov = pd.concat([train_me, train_COV.reindex(train_me.index)], axis=1)
        train_cov.columns=['Outcome','SNP1','SNP2','Age', 'BMI', 'Sex','Smoking']
        train_cov = train_cov[['Outcome','Age', 'BMI', 'Sex','Smoking']]

        # BLOCK 3
        # Calculate weights from the training dataset
        edge_weights_me_t = train_me.genomics.calculate_edge_encoding_values(
            data=train_cov, outcome_variable="Outcome",covariates=['Age', 'BMI', 'Sex','Smoking'])
        edge_weights_me = edge_weights_me_t.copy()
        edge_weights_me.insert(loc=0, column="BioAct", value=ab1l)
        edge_weights_me.insert(loc=0, column="TrainSeed", value=train_seed)
        edge_weights_me.insert(loc=0, column="TestSeed", value=test_seed)

        # BLOCK 4
        # Test data
        test_main = sim.BAMS.from_model(
            eff1=getattr(sim.SNPEffectEncodings, ab1u),
            eff2=getattr(sim.SNPEffectEncodings, ab2u),
            penetrance_base=PEN_BASE,
            penetrance_diff=PEN_DIFF,
            main1=1,
            main2=0,
            interaction=0,
            snp1=variant1,
            snp2=variant2,
            random_seed=test_seed)
        test_me = test_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB)
        np.random.seed(test_seed)
        Age_1 = pd.DataFrame(np.random.poisson(lam=31, size=(int(0.507*n_controls))).astype(int))
        Age_2 = pd.DataFrame(np.random.poisson(lam=54.5, size=(int(0.316*n_controls))).astype(int))
        Age_3 = pd.DataFrame(np.random.poisson(lam=82.5, size=(n_controls-int(0.507*n_controls)-int(0.316*n_controls))).astype(int))
        Age_4 = pd.DataFrame(np.random.poisson(lam=31, size=(int(0.144*n_cases))).astype(int))
        Age_5 = pd.DataFrame(np.random.poisson(lam=54.5, size=(int(0.435*n_cases))).astype(int))
        Age_6 = pd.DataFrame(np.random.poisson(lam=82.5, size=(n_cases-int(0.144*n_cases)-int(0.435*n_cases))).astype(int))
        Age = pd.concat([Age_1,Age_2,Age_3,Age_4,Age_5,Age_6], axis=0, ignore_index=True, sort=False)
        BMI_1 = pd.DataFrame(np.random.normal(loc=26.5, scale=1.0, size=(n_controls)).astype(int))
        BMI_2 = pd.DataFrame(np.random.normal(loc=35, scale=7.5, size=(n_cases)).astype(int))
        BMI = pd.concat([BMI_1,BMI_2], axis=0, ignore_index=True, sort=False)
        Sex_1 = pd.DataFrame(np.random.binomial(n=1, p=0.481, size=(n_controls)).astype(int))
        Sex_2 = pd.DataFrame(np.random.binomial(n=1, p=0.525, size=(n_cases)).astype(int))
        Sex = pd.concat([Sex_1,Sex_2], axis=0, ignore_index=True, sort=False)
        Smoking_1 = pd.DataFrame(np.random.binomial(n=1, p=0.125, size=(n_controls)).astype(int))
        Smoking_2 = pd.DataFrame(np.random.binomial(n=1, p=0.64, size=(n_cases)).astype(int))
        Smoking = pd.concat([Smoking_1,Smoking_2], axis=0, ignore_index=True, sort=False) 
        test_COV = pd.concat([Age,BMI,Sex,Smoking], axis=1, ignore_index=True, sort=False)
        #test_cov = pd.concat([test_me, test_COV.reindex(test_me.index)], axis=1)
        #test_cov.columns=['Outcome','SNP1','SNP2','Age', 'BMI', 'Sex','Smoking']
        #test_cov = test_cov[['Outcome','Age', 'BMI', 'Sex','Smoking']]
        
        # Run Regression by using weightes from CLARITE
        for v_enc in encoding:
            #encode = str("encode_" + v_enc.lower())
            # test_me_enc = getattr(test_me.genomics, encode) # More Pythonics but didnt works
            if v_enc != "DOMDEV":
                if v_enc.lower() == "additive":
                    test_me_enc = test_me.genomics.encode_additive()
                elif v_enc.lower() == "recessive":
                    test_me_enc = test_me.genomics.encode_recessive()
                elif v_enc.lower() == "dominant":
                    test_me_enc = test_me.genomics.encode_dominant()
                elif v_enc.lower() == "codominant":
                    test_me_enc = test_me.genomics.encode_codominant()
                else:
                    test_me_enc = test_me.genomics.encode_edge(
                        encoding_info=edge_weights_me_t)
                
                test_me_enc = pd.concat([test_me_enc, test_COV.reindex(test_me_enc.index)], axis=1)
                test_me_enc.columns=['Outcome','SNP1','SNP2','Age', 'BMI', 'Sex','Smoking']
                
                results_me = clarite.analyze.association_study(
                    data=test_me_enc, outcomes="Outcome",covariates=['Age', 'BMI', 'Sex','Smoking'])
                results_me["odds ratio"] = np.exp(results_me["Beta"])
                results_me.insert(loc=0, column="Encoding", value=v_enc)
                results_me.insert(loc=0, column="BioAct", value=ab1l)
                results_me.insert(loc=0, column="TrainSeed", value=train_seed)
                results_me.insert(loc=0, column="TestSeed", value=test_seed)

                ALL_RESULTS_ENCODING = pd.concat(
                    [ALL_RESULTS_ENCODING, results_me], axis=0)

                """if v_enc.lower() == "edge":
                    ALL_RESULTS_ENCOGIND_EDGE = pd.concat(
                        [ALL_RESULTS_ENCODING_EDGE, results_me], axis=0)"""  
        
            else:
                # DOMDEV Encoding
                test_me_pb000_DOMDEV = test_me_enc
                test_me_pb000_DOMDEV["COV1"] = test_me_pb000_DOMDEV["SNP1"]
                test_me_pb000_DOMDEV["COV2"] = test_me_pb000_DOMDEV["SNP2"]

                test_me_pb000_DOMDEV["SNP1"] = test_me_pb000_DOMDEV["SNP1"].replace(
                    2, 0)
                test_me_pb000_DOMDEV["SNP2"] = test_me_pb000_DOMDEV["SNP2"].replace(
                    2, 0)

                test_me_pb000_DOMDEV_SNP1_t = test_me_pb000_DOMDEV[[
                    "Outcome", "SNP1", "COV1"]]
                test_me_pb000_DOMDEV_SNP1_t = pd.concat([test_me_pb000_DOMDEV_SNP1_t, test_COV.reindex(test_me_pb000_DOMDEV_SNP1_t.index)], axis=1)
                test_me_pb000_DOMDEV_SNP1_t.columns=['Outcome','SNP1','COV1','Age', 'BMI', 'Sex','Smoking']
                test_me_pb000_DOMDEV_SNP2_t = test_me_pb000_DOMDEV[[
                    "Outcome", "SNP2", "COV2"]]
                test_me_pb000_DOMDEV_SNP2_t = pd.concat([test_me_pb000_DOMDEV_SNP2_t, test_COV.reindex(test_me_pb000_DOMDEV_SNP2_t.index)], axis=1)
                test_me_pb000_DOMDEV_SNP2_t.columns=['Outcome','SNP2','COV2','Age', 'BMI', 'Sex','Smoking']

                DOMDEV_results_me_pb000_SNP1 = clarite.analyze.association_study(
                    data=test_me_pb000_DOMDEV_SNP1_t,
                    outcomes="Outcome",
                    covariates=["COV1",'Age', 'BMI', 'Sex','Smoking']
                )
                DOMDEV_results_me_pb000_SNP2 = clarite.analyze.association_study(
                    data=test_me_pb000_DOMDEV_SNP2_t,
                    outcomes="Outcome",
                    covariates=["COV2",'Age', 'BMI', 'Sex','Smoking']
                )

                DOMDEV_results_me = pd.concat(
                    [DOMDEV_results_me_pb000_SNP1, DOMDEV_results_me_pb000_SNP2])

                DOMDEV_results_me["odds ratio"] = np.exp(
                    DOMDEV_results_me["Beta"])
                DOMDEV_results_me.insert(loc=0, column="Encoding", value=v_enc)
                DOMDEV_results_me.insert(loc=0, column="BioAct", value=ab1l)
                DOMDEV_results_me.insert(
                    loc=0, column="TrainSeed", value=train_seed)
                DOMDEV_results_me.insert(
                    loc=0, column="TestSeed", value=test_seed)

                ALL_RESULTS_ENCODING = pd.concat(
                    [ALL_RESULTS_ENCODING, DOMDEV_results_me], axis=0)

        ALL_RESULTS_EDGE_ALPHA = pd.concat(
            [ALL_RESULTS_EDGE_ALPHA, edge_weights_me], axis=0)

    return ALL_RESULTS_EDGE_ALPHA, ALL_RESULTS_ENCODING


def mp_treads():
    FINAL_RESULTS_EDGE_ALPHA = pd.DataFrame()
    FINAL_RESULTS_ENCODING = pd.DataFrame()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        future_to = {executor.submit(simulations, n): n for n in n_loops}

        for future in concurrent.futures.as_completed(future_to):

            data = future.result()
            data_0 = pd.DataFrame(data[0])
            data_1 = pd.DataFrame(data[1])

            FINAL_RESULTS_EDGE_ALPHA = pd.concat(
                [FINAL_RESULTS_EDGE_ALPHA, data_0], axis=0)
            FINAL_RESULTS_ENCODING = pd.concat(
                [FINAL_RESULTS_ENCODING, data_1], axis=0)

            #FINAL_RESULTS_EDGE_ALPHA = pd.concat([FINAL_RESULTS_EDGE_ALPHA, data_0], axis=0)

    FINAL_RESULTS_EDGE_ALPHA.to_csv(
        f"{output_path}/EDGE_alpha_Results_{num_samples}_{case_control_ratio}_pb{PEN_BASE}_pd{PEN_DIFF}_maf{MAFA}.txt", sep=";")

    FINAL_RESULTS_ENCODING.to_csv(
        f"{output_path}/All_Results_{num_samples}_{case_control_ratio}_pb{PEN_BASE}_pd{PEN_DIFF}_maf{MAFA}.txt", sep=";")

    print(f"Finished in {round(time.perf_counter()-start,2)} second(s)")


if __name__ == "__main__":
    # ==> PARAMETERS:

    # Define the variant for two SNPs
    variant1 = scalars.Variant("1", 1, id="rs1", ref="A", alt=["C"])
    variant2 = scalars.Variant("1", 2, id="rs2", ref="G", alt=["T"])

    # Define Case-Control ratio
    num_samples = 2000
    case_control_ratio = '1:3'
    n_controls = int(num_samples*3/4)
    n_cases = num_samples - n_controls
    PEN_BASE = (1-0.1)/2
    PEN_DIFF = 0.1
    MAFA = 0.1
    MAFB = 0.1
    
    # Interations for Train and Test
    n_loops = range(1000)

    work_group = (
        ["REC", "ADD"],
        ["SUB", "ADD"],
        ["ADD", "ADD"],
        ["SUP", "ADD"],
        ["DOM", "ADD"],
        ["HET", "ADD"],
        ["NUL", "ADD"])

    encoding = ("Additive", "DOMDEV", "Recessive",
                "Dominant", "Codominant", "EDGE")

    output_path = "/storage/home/jpz5091/scratch/SimWithCovars/pen_diff01/maf01"
    #output_path = "files"

    start = time.perf_counter()

    mp_treads()
