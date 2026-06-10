# Genomic Wide Association Study (GWAS)

## 1. Introduction
GWAS is commonly used to test the association between a phenotype and genotypes. For that several models are proposed. The GWAS protocol is based on [cloudfiel.github.io/GWASTutorial](cloudfield.github.io/GWASTutorial/06_Association_tests/#association-tests-using-plink). As we can see on the _following picture_, GWAS requires initialisation, filtering, QC and PCA steps. 

<img width="515" height="470" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GWAS_overview.png" />

For the study, 3 GWAS models are realized with ***a standard baseline GWAS (PLINK logistic)***, ***a mixed-model GWAS (GEMMA)*** corrected for kinship and small population effects, and ***a mixed-model GWAS (GEMMA) for males only*** to compare between tuskless and tusked type males.

Similar as [Campbell et al. (2021)](https://www.science.org/doi/full/10.1126/science.abe7389), I focused the analyses on the whole genome and not only chromosome 1 and X, to identify the regions 

## 2. Standard baseline GWAS

The aim of this part is to build a first exploratory baseline of GWAS and to select the number of principal components to retain. PLINK logistic regression is used.
The **model equation** : 

> ***Tusk ~ SNP + sex + population + PCs***

The script of this first GWAS model are available [SEE 08.GWAS_PLINK_logistic_PC5_PC10.slurm & 09.GWAS_QQ_Manhattan.R]. 

### _2.1 Q-Q plots and Genomic Control Inflation Factor_

Visually, the observed line for Q-Q plot with 5PCs is under the expected line meaning that observed p-values are less extremes as expected. A Genomic Control Inflation Factor, ***lambda GC***,  equal to 1 is a perfect alignment between expected and observed. Over 1.1, it is inflation and under 1 the model is conservative, over-corrected or lack of power. 

The ***Genomic Control Inflation Factor*** is defined as : 
> lambda GC = Median(X²obs)/Median(X²1,nul)
> 
> Median(X²1,nul) = 0.455
> 
> X²obs = F-1|X²1 (1-p)

Here lambda GC = 0.167, meaning that I probably lacking of power and that GWAS model is too conservative (pvalues too high and test statistics too low). The second Q-Q plot is the opposite with an inflation (lambda GC = 4.286) > 1, a signal that covariable numbers is too high.

<img width="315" height="410" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/QQ_PC5.png" /> <img width="315" height="370" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/QQ_PC10.png" />



