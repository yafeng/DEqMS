# DEqMS
DEqMS is a statistical tool for testing differential protein expression in quantitative proteomic analysis, developed by Yafeng Zhu @ Karolinska Institutet. 

DEqMS is a published method, if you use it in your research, please cite:
Zhu et al. DEqMS: A Method for Accurate Variance Estimation in Differential Protein Expression Analysis. 
Mol Cell Proteomics 2020 Mar 23. [PMID: 32205417](https://pubmed.ncbi.nlm.nih.gov/32205417/)


## Installation
To install lastest DEqMS package (v1.6.0), start R (version "4.0") and enter:
```{r}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DEqMS")

```
## Introduction
DEqMS is developed on top of Limma. However, Limma assumes same prior variance for all genes. In proteomics, the accuracy of protein abundance estimates varies by the number of peptides/PSMs quantified in both label-free and labelled data. Proteins quantification by multiple peptides or PSMs are more accurate. DEqMS package is able to estimate different prior variances for proteins quantified by different number of PSMs/peptides, therefore achieving better accuracy. The package can be applied to analyze both label-free and labelled proteomics data.

## How to use it
Browse DEqMS online Vignettes [here](https://bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html)
