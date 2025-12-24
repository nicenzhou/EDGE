# EDGE: Elastic Data-Driven Encoding for GWAS

## Overview

Most GWAS assume an additive inheritance model, which assigns heterozygous genotypes half the risk of homozygous-alternate genotypes. This has led to a focus on additive genetic effects in complex disease research. Growing evidence indicates that many single-nucleotide polymorphisms (SNPs) have nonadditive effects, including dominant and recessive effects, which are missed by the additive model alone. EDGE (Elastic Data-Driven Encoding) is a novel method for genome-wide association studies (GWAS) that determines the inheritance model each SNP contributes to a given trait, allowing for unique and flexible SNP encoding. Unlike traditional GWAS that assume an additive inheritance model, EDGE can identify additive and also nonadditive effects, including dominant and recessive patterns.

## Statistical Model

EDGE employs a flexible encoding approach based on a regression model that separately estimates effects for heterozygous and homozygous alternate genotypes:

**Equation 1: Regression Model**

$$E(Y | SNP_{Het}, SNP_{HA}, COV_i) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i$$

Where:
- $Y$ = phenotype/outcome
- $SNP_{Het}$ = indicator for heterozygous genotype
- $SNP_{HA}$ = indicator for homozygous alternate genotype
- $COV_i$ = covariates
- $\beta_{Het}$ = effect size for heterozygous genotype
- $\beta_{HA}$ = effect size for homozygous alternate genotype
- $\beta_{cov_i}$ = effect sizes for covariates
- 
**Equation 2: Encoding Parameter**

$$\alpha = \frac{\beta_{Het}}{\beta_{HA}}$$

Where:
- $$\alpha$$ = encoding parameter representing the ratio of heterozygous to homozygous alternate effects
- $$\beta_{Het}$$ = effect size for heterozygous genotype (from Equation 1)
- $$\beta_{HA}$$ = effect size for homozygous alternate genotype (from Equation 1)

### Inheritance Model Classification

The encoding parameter `α` determines the inheritance pattern for each SNP:

| Inheritance Pattern | α Range | Description |
|-------------------|---------|-------------|
| **Under-recessive** | α < 0 | Heterozygous effect opposite direction of homozygous alternate |
| **Recessive** | 0 ≤ α < 0.125 | Effect primarily in homozygous alternate genotype |
| **Sub-additive** | 0.125 ≤ α < 0.375 | Heterozygous effect less than half of homozygous alternate |
| **Additive** | 0.375 ≤ α < 0.625 | Heterozygous effect approximately half of homozygous alternate |
| **Super-additive** | 0.625 ≤ α < 0.875 | Heterozygous effect more than half of homozygous alternate |
| **Dominant** | 0.875 ≤ α ≤ 1 | Heterozygous effect similar to homozygous alternate |
| **Over-dominant** | α > 1 | Heterozygous effect exceeds homozygous alternate effect |

This flexible approach allows EDGE to identify the optimal genetic encoding for each SNP based on the data, rather than assuming a fixed additive model.

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

