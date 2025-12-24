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

## Inheritance Model Classification

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

### Large-Scale GWAS Implementation

EDGE was applied to genome-wide association studies using:
1. **UK Biobank (UKB)**: Calculate encoding parameter $\alpha$ for each SNP
2. **Million Veteran Program (MVP)**: Independent validation
3. **Meta-analysis**: Combined UKB and MVP results

### Key Findings

EDGE GWAS meta-analysis identified **nonadditive inheritance patterns for more than 52% of genome-wide significant loci** for:
- **Coronary artery disease (CAD)**
- **Body mass index (BMI)**

### Impact

This research:
- Identifies novel disease-risk SNPs missed by traditional additive models
- May improve polygenic risk prediction in diverse populations
- Provides a framework for future applications to thousands of disease phenotypes

## Citation

If you use EDGE in your research, please cite: 
Zhou, J., Rico, A. L. G., Guare, L., Million Veteran Program, Chang, K. M., Tsao, P. S., Assimes, T. L., Verma, S. S., & Hall, M. A. (2023). Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits. *medRxiv*, 2023.06.01.23290857. https://doi.org/10.1101/2023.06.01.23290857

**BibTeX:**
```bibtex
@article{zhou2023flexibly,
  title={Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits},
  author={Zhou, Jiayan and Rico, Andre Luis Garao and Guare, Lindsay and Million Veteran Program and Chang, Kyong-Mi and Tsao, Philip S and Assimes, Themistocles L and Verma, Shefali Setia and Hall, Molly Ann},
  journal={medRxiv},
  pages={2023--06},
  year={2023},
  publisher={Cold Spring Harbor Laboratory Press},
  doi={10.1101/2023.06.01.23290857}
}```

MIT License

Copyright (c) 2023 Jiayan Zhou, Andre Luis Garao Rico, Lindsay Guare, Kyong-Mi Chang, Philip S. Tsao, Themistocles L. Assimes, Shefali Setia Verma, Molly Ann Hall

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
