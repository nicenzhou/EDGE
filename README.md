# EDGE: Elastic Data-Driven Encoding for GWAS

## Overview

EDGE (Elastic Data-Driven Encoding) is a novel method for genome-wide association studies (GWAS) that determines the inheritance model each SNP contributes to a given trait, allowing for unique and flexible SNP encoding. Unlike traditional GWAS that assume an additive inheritance model, EDGE can identify additive and also nonadditive effects, including dominant and recessive patterns.

## Background

Most GWAS assume an additive inheritance model, which assigns heterozygous genotypes half the risk of homozygous-alternate genotypes. This has led to a focus on additive genetic effects in complex disease research. Growing evidence indicates that many single-nucleotide polymorphisms (SNPs) have nonadditive effects, including dominant and recessive effects, which are missed by the additive model alone.

## Applications

EDGE GWAS analysis identified nonadditive inheritance patterns for more than 52% of genome-wide significant loci for:
- Coronary artery disease (CAD)
- Body mass index (BMI)

This research:
- Identifies novel disease-risk SNPs missed by traditional additive models
- May improve polygenic risk prediction in diverse populations
- Provides a framework for future applications to thousands of disease phenotypes

## Authors

**Jiayan Zhou**<sup>1,2</sup>, Andre Luis Garao Rico<sup>3,4</sup>, Lindsay Guare<sup>4</sup>, Million Veteran Program, Kyong-Mi Chang<sup>5,6</sup>, Philip S. Tsao<sup>1,2</sup>, Themistocles L. Assimes<sup>1,2*</sup>, Shefali Setia Verma<sup>7*</sup>, Molly Ann Hall<sup>3,4*</sup>

<sup>*</sup>Co-supervising senior authors

### Affiliations

1. VA Palo Alto Healthcare System, Palo Alto, CA 94304, USA
2. Department of Medicine, Division of Cardiovascular Medicine, Stanford University School of Medicine, Stanford, CA 94305, USA
3. Department of Genetics, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 19104, USA
4. Institute for Biomedical Informatics, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 19104, USA
5. Corporal Michael J. Crescenz VA Medical Center, Philadelphia, PA 19104, USA
6. Department of Medicine, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 19104, USA
7. Department of Pathology and Laboratory Medicine, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 19104, USA

## Contact

**Corresponding Author:**  
Molly Ann Hall - molly.hall@pennmedicine.upenn.edu

**For questions about the code:**  
Jiayan Zhou - jyzhou@stanford.edu

## Citation

If you use EDGE in your research, please cite: 
Zhou J, Rico ALG, Guare L, et al. Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits. [Journal TBD]

## Abstract
Abstract
Background: Most genome-wide association studies (GWAS) assume an additive inheritance model, which assigns heterozygous genotypes half the risk of homozygous-alternate genotypes. This has led to a focus on additive genetic effects in complex disease research. Growing evidence indicates that many single-nucleotide polymorphisms (SNPs) have nonadditive effects, including dominant and recessive effects, which are missed by the additive model alone. 
Results: To address this issue, we developed Elastic Data-Driven Encoding (EDGE) to determine the inheritance model each SNP contributes to a given trait, allowing for unique and flexible SNP encoding in GWAS. Simulation results demonstrate that EDGE provides higher power than additive and other genetic encoding models across a wide range of simulated inheritance patterns while maintaining a conserved false positive rate. EDGE GWAS on data from the UK BioBank and the Million Veteran Program, comprising more than 500,000 individuals, identified nonadditive inheritance patterns for more than 52% of the genome-wide significant loci for coronary artery disease and body mass index. 
Conclusions: This research lays the groundwork for integrating nonadditive genetic effects into GWAS workflows to identify novel disease-risk SNPs, which may ultimately improve polygenic risk prediction in diverse populations and provide a springboard for future applications to thousands of disease phenotypes.

