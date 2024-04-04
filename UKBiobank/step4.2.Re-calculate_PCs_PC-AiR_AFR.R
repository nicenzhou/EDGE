##################################################################################################################################
## 																																
## 	Script Name: step4.2.Re-calculate_PCs_PC-AiR_AFR.R											
## 	Description: PCA was limited to unrelated individuals for AFR using PC-AiR. This script prepares the input for PC-AiR.
## 	Requirement: A QC through PLINK has to be done to the BEGN files and raw genotype calls. PC-AiR needs to be run in R. 				
## 	Authors: Jiayan Zhou <jyzhou@stanford.edu>																																												
## 																																
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Details:																														
##	      This script contains steps for PC calculation through PC-AiR.
## 																																	
## ---------------------------------------------------------------------------------------------------------------------------- 
## 
## 	Reference:
##	      EDGE GWAS: 
##	        Zhou, J. et al. Flexibly encoded GWAS identifies novel nonadditive SNPs in individuals of African and European ancestry.
##	        medRxiv 2023.06.01.23290857; doi: https://doi.org/10.1101/2023.06.01.23290857 
##	      EDGE GxG: 
##	        Hall, M. A. et al. Novel EDGE encoding method enhances ability to identify genetic interactions. PLoS Genetics 17, e1009534 (2021).
##	      PC-AiR/Population Structure and Relatedness Inference using the GENESIS Package:
##	        https://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html
## 				
##################################################################################################################################

# ------------------------------------- #
## PC calculation using PC-AiR
# ------------------------------------- #
## IF R IN RAP GIVES AN ERROR ABOUT IGRAPH, 
## YOU MIGHT RUN THIS SECTION IN THE LOCAL MACHINE.
# ------------------------------------- #
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install GENESIS R package
BiocManager::install(version='3.16')
BiocManager::install("GENESIS")

library(SNPRelate)
library(GENESIS)
library(GWASTools)

snpgdsBED2GDS(bed.fn = "geno_merged_afr_full_r.bed", 
              bim.fn = "geno_merged_afr_full_r.bim", 
              fam.fn = "geno_merged_afr_full_r.fam", 
              out.gdsfn = "geno_merged_afr_full_r.gds")

showfile.gds(closeall=TRUE)
#gdsfile <- system.file(filename = "geno_merged_afr_full.gds", package="GENESIS")
gds <- snpgdsOpen("~/geno_merged_afr_full_r.gds")
#snpgdsClose(gds)

# LD Pruning
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1), verbose=TRUE)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)

# KING input
king <- snpgdsIBDKING(gds,num.thread=32,verbose=TRUE)
KINGmat <- kingToMatrix(king)

library(GWASTools)
snpgdsClose(gds)
geno <- GdsGenotypeReader(filename = "geno_merged_afr_full_r.gds")
genoData <- GenotypeData(geno)

# run PC-AiR on pruned SNPs
mypcair <- pcair(genoData, kinobj = KINGmat, divobj = KINGmat, num.cores = 32, verbose = TRUE)
summary(mypcair)

# plot top 2 PCs
plot(mypcair)
# plot PCs 3 and 4
plot(mypcair, vx = 3, vy = 4)

pcair_results = as.data.frame(mypcair$vectors)
write.table(pcair_results, file="afr_cad_pcair_pcs.txt", sep="\t", col.names = F, row.names = T)
