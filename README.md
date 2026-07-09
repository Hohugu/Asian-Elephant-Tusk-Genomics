# Genomic project on tusk genes identification in Asian elephants (Elephas maximus)

## Introduction

In Asian elephants, it has been observed by Evans in the 1900s that Asian elephants display a strong sexual dimorphism by the only presence of tusk in males.
Unlike African elephants, only males display tusk, however the phenotype without tusk can also be found in males as well as females. 

This difference of phenotype is still unexplained, and consequently the aim of this project is to identify the genes which are potentially associated with the phenotype tusk in Asian elephants. 
For that, blood samples have been collected and DNA have been extracted and sequenced from the Myanmar Timber Enterprise (MTE) before 2021. Thanks to the India Institute Science (IISc) collaboration, the study covered 93 individuals. The analyses were based on analysis in [Campbell-Staton et al., 2021](https://www.science.org/doi/full/10.1126/science.abe7389) which have identified two genes involved in the expression and the growth of tusk in African elephants in an anthropogenetic selective pressure (i.e poaching and loss of tusk in African elephants). Theses genes are AMELX and MEP1A, which can also be find in mice. In this study, females were also taken into account for the analyses of gene identification.

This project is separated into four sections : **WGS Variant Calling**, **Pre-GWAS**, **GWAS** and **Genetic analyses**. Each section form one branch of the github. You can then navigate between branches by clicking on the dropdown menu. 

## Project Part

### WGS Variant Calling
The first part **WGS Variant Calling** is a workflow intitled "*Variant Calling Pipeline : GATK Best practise Germline short variant discovery (SNPs+Indels) workflow*" from these sources :
 - > https://www.youtube.com/watch?v=iHkiQvxyr5c&t=1803s
 - > https://www.cog-genomics.org/plink/1.9/strat
 - > https://vcftools.sourceforge.net/man_latest.html 
 - > https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

The aim of this part is to start with sequencing reads and perform a series of steps to determine a set of genetic variants. We start with data preprocesssing, then variant discovery and finally filtering and annotation, to end with a clean VCF file, Analysis-ready. 

 
  <img width="834" height="452" alt="image" src="https://github.com/user-attachments/assets/60897ab0-24cc-4cd1-8974-b2775ab4661b" />


For this part, I used ***.sra*** files from NCBI collaborator project to download all the individual files on **Puhti** server from CSC, in order to convert them into **.fastq** and **.fastq.gz** files and then into **.bam** files to combine them to those for Myanmar. Before the recalibrate base quality scores and the conversion into .bam files, ***boostrapping*** is usually used, especially when no known set of variants is available. However, no boostrapping is performed before this study. 

All the script are numeroted and explicitly intitled to make the protocol easier to follow. 

### Pre-GWAS
For this step, I based all the steps on the manual from Physalia, called "introduction_to_GWAS", but also on other format like GWAS tutorial available on https://cloufield.github.io/GWASTutorial/01_Dataset/
This step concerns data pre-processing with initial and explonatory data analysis (IDA & EDA). Then, a quality control is effectuated, followed by the imputation step. I require different software and packages to compute the pipeline for data pre-processing but also for GWAS later : RStudio (R) using the R package rrBLUP, PLINK 1.9 (PLINK is a very powerful open-source genome analysis toolset, with a lot of computationally efficient functions for data pre-processing, filtering, formatting, and analysis), [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html). A R script is also provided to prepare and explore the phenotype and genomic data. 

From the introduction GWAS manual : "
**Initial data analysis (IDA)** mainly focuses on data cleaning, a first screening, and transformation (if necessary) to ensure data quality and confirm that our data set meets the relevant distributional and model assumptions. Data cleaning may include steps such as elimination of duplicate records, handling of missing values, identification of systematic errors, or correction of coding inconsistencies.

**Exploratory data analysis (EDA)** is used to examine data sets and summarize their main characteristics. EDA helps to discover patterns in the data, spot anomalies and outliers, test a hypothesis and check our assumptions. EDA tells us what data can reveal beyond the formal modeling or hypothesis testing and provides a better understanding of data set variables and their interactions. It can also help determine if the statistical techniques you are considering for
data analysis are appropriate. IDA and EDA often employ data visualization methods to check distributional characteristics, identify relationships between variables, identify potential cofactors, and spot inconsistencies. 

We will run an IDA and EDA for both our phenotypic data and the genotypic data."

In [PreGWAS](https://github.com/Hohugu/Asian-Elephant-Tusk-Genomics/blob/02.Pre-GWAS/README.md#pre-genomic-wide-association-study---gwas) steps, I will also conduct the Genotype filtering. No imputation is required here. 

### GWAS, Selection analyses and Functional annotation

 [GWAS](https://github.com/Hohugu/Asian-Elephant-Tusk-Genomics/blob/03.GWAS/README.md#genome-wide-association-study-gwas), [Selection analyses](https://github.com/Hohugu/Asian-Elephant-Tusk-Genomics/blob/04.Selection-Differentiation-scans/README.md#selection-and-population-differentiation-scans), [Functional annotation](https://github.com/Hohugu/Asian-Elephant-Tusk-Genomics/blob/06.Functional-Regulatory-Annotation/README.md#functional-and-regulatory-proximal-annotation) and [Heritability](https://github.com/Hohugu/Asian-Elephant-Tusk-Genomics/blob/07.Pedigree-Heritability/README.md#pedigree-and-heritability) are declined by branches. Inside these branches, the detail of the analyses are provided. The aim is to identified SNPs that show strong association with the phenotype. For that different models are tested to find the fittest model. Selection analyses, including FST, nucleotide diversity, Tajima's D, LD, heterozygosity, are realized in a whole-genome scan approach, reproducting Campbell et al. protocol. Regions with significant and redondant results are considered as candidate genes and selected for functional annotation. The latter analysis allow to describe the content enrichment, location (codent or non-codent), and potential gene description. Rare and regulatory elements are also investigated but with Whole genome sequencing (WGS) limitations. Indeed, WGS can not provide informations as much as CHiP or ATAC-seq can do on functional verification and mechanisms on gene interactions. Heritability and pedigree analyses are also realized using MCMCglmm.
