# Genomic Wide Association Study (GWAS)

## 01. Introduction
GWAS is commonly used to test the association between a phenotype and genotypes. For that several models are proposed. The GWAS protocol is based on [cloudfiel.github.io/GWASTutorial](cloudfield.github.io/GWASTutorial/06_Association_tests/#association-tests-using-plink). As we can see on the following picture, GWAS requires initialisation, filtering, QC and PCA steps. 

<img width="515" height="470" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GWAS_overview.png" />

For the study, 3 GWAS models are realized with a standard baseline GWAS (PLINK logistic), a mixed-model GWAS (GEMMA) corrected for kinship and small population effects, and a mixed-model GWAS (GEMMA) but only with males to compare between tuskless and tusked type males.

Similar as [Campbell et al. (2021)](https://www.science.org/doi/full/10.1126/science.abe7389), I focused the analyses on the whole genome and not only chromosome 1 and X, to identify the regions 

## 02. Standard baseline GWAS

The aim of this part is to build a first exploratory baseline of GWAS and to select the number of principal components to retain. PLINK logistic regression is used.
The model equation : 

> Tusk ~ SNP + sex + population + PCs

The script of this first GWAS model are available [SEE 08.GWAS_PLINK_logistic_PC5_PC10.slurm & 09:GWAS_QQ_Manhattan.R]. 



