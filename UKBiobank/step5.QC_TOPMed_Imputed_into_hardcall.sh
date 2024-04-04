##################################################################################################################################
## 																																
## 	Script Name: step5.QC_TOPMed_Imputed_into_hardcall.sh											
## 	Description: QC for the TOPMed Imputed data was performed and output as the hardcall genotypes. 
## 	Requirement: Both PLINKs are required. keep_afr_frombgen_noheader.txt from step 4.3 and keep_full_eur_noheader.txt from step 3.
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script contains steps to re-assign the missing variant IDs, genotype call rate QC, sample call rate QC, 
##	      HWE check, MAF check, LD pruning, and reorder with minor allele as effect allele.
##	      The reorder step is mandatory because the minor allele is required to be the effect allele for later analysis.  
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
## BEGN QC
# ------------------------------------- #
plink2 \
--threads 32 \
--bgen ukb21007_c*_b0_v1.bgen ref-first \ 
--sample ukb21007_c*_b0_v1.sample \
--missing-code 0 \
--set-missing-var-ids @:#:\$r:\$a \
--new-id-max-allele-len 100 \
--rm-dup \
--keep {keep_afr_frombgen_noheader.txt/keep_full_eur_noheader.txt} \
--geno 0.01 \
--mind 0.01 \
--hwe 10e-6 midp \
--maf 0.05 \
--snps-only just-acgt \
--make-bed --out ukb21007_qced; 

plink2 \ 
--threads 32 \
--bfile ukb21007_qced \
--indep-pairwise 500kb 0.8 \
--out ukb21007_qced_ld; 

plink2 --threads 32 \
--bfile ukb21007_qced \
--extract ukb21007_qced_ld.prune.in \
--keep {keep_afr_frombgen_noheader.txt/keep_full_eur_noheader.txt} \
--geno 0.01 \
--mind 0.01 \
--hwe 10e-6 midp \
--maf 0.05 \
--snps-only just-acgt \
--make-bed \
--out ukb21007_c{1..22}_cad_{afr/eur}_hardcall8_qced;


# ------------------------------------- #
## Reorder minor allele as effect allele
# ------------------------------------- #
cp /{PATH}/ukb21007_c[1-9]* . ; 
for file in *.bed *.bim *.fam; do 
  if [[ -f "$file" ]]; 
  then prefix="${file%.*}"; 
  plink --bfile "$prefix" --make-bed --out "${prefix}_reordered"; 
  fi; 
done
