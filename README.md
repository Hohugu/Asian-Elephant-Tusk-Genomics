# Genomic project on tusk genes identification in Asian elephants (Elephas maximus)

## Introduction

In Asian elephants, it has been observed by Evans in the 1900s that Asian elephants display a strong sexual dimorphism by the only presence of tusk in males.
Unlikely African elephants, only males display tusk, however the phenotype without tusk can also be found in males as well as females. 

This difference of phenotype still unexplained, and consequently the aim of this project is to identify the genes which are potentially associated with the phenotype tusk in Asian elephants. 
For that, blood samples have been collected and DNA have been extracted and sequenced from the Myanmar Timber Enterprise (MTE) before 2021. Thanks to the India Institude Science (IISc) collaboration, the study covered 93 individuals. The analysis were based on analysis in [Campbell-Staton et al., 2021](https://www.science.org/doi/full/10.1126/science.abe7389) which have identified two genes involved in the expression and the growth of tusk in African elephants in an anthropogenetic selective pressure (i.e poaching and loss of tusk in African elephants). Theses genes are AMELX and MEP1A, which can also be find in mices. In this study, females were also taken into account for the analysis. In this study on Asian elephants, we will also take females into analysis of gene identification given they is a dimorphism context. 

This project is separated into three sections : **WGS Variant Calling**, **Pre-GWAS** and **GWAS**. 

## Project Part

### WGS Variant Calling
The first part **WGS Variant Calling** is a workflow intitled "*Variant Calling Pipeline : GATK Best practise Germline short variant discovery (SNPs+Indels) workflow*" from these sources :
 - > https://www.youtube.com/watch?v=iHkiQvxyr5c&t=1803s
 - > https://www.cog-genomics.org/plink/1.9/strat
 - > https://vcftools.sourceforge.net/man_latest.html 
 - > https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

The aim of this part is to start with sequencing reads and perform a series of steps to determine a set of genetic variants. We start with data preprocesssing, then variant discovery and finally filtering and annotation, to end with a clean VCF file, Analysis-ready. 

 
  <img width="834" height="452" alt="image" src="https://github.com/user-attachments/assets/60897ab0-24cc-4cd1-8974-b2775ab4661b" />


For this part, I used ***.sra*** files from NCBI collaborator project to download all the individual files on **Puhti** server from CSC, in order to convert them into **.fastq** and **.fastq.gz** files and then into **.bam** files to combine them to those for Myanmar. Before the recalibrate base quality scores and the conversion into .bam files, I used ***boostrapping***, because no known set of variants were available before this study to respect the bio-informatic protocol of variant calling. 

All the script are numeroted and explicitly intitled to make the protocol easier to follow. 

### Pre-GWAS
For this step, I based all the steps on the manual from Physalia, called "introduction_to_GWAS", but also on other format ...
This step concerns data pre-processing with initial and explonatory data analysis (IDA & EDA). Then, a quality control is effectuated, followed by the imputation step. I require different software and packages to compute the pipeline for data pre-processing but also for GWAS later : RStudio (R) using the R package rrBLUP, PLINK 1.9 (PLINK is a very powerful open-source genome analysis toolset, with a lot of computationally efficient functions for data pre-processing, filtering, formatting, and analysis), [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html). A R script is also provided to prepare and explore the phenotype and genomic data. 

(From the introduction GWAS manual :) "
**Initial data analysis (IDA)** mainly focuses on data cleaning, a first screening, and transformation (if necessary) to ensure data quality and confirm that our data set meets the relevant distributional and model assumptions. Data cleaning may include steps such as elimination of duplicate records, handling of missing values, identification of systematic errors, or correction of coding inconsistencies.

**Exploratory data analysis (EDA)** is used to examine data sets and summarize their main characteristics. EDA helps to discover patterns in the data, spot anomalies and outliers, test a hypothesis and check our assumptions. EDA tells us what data can reveal beyond the formal modeling or hypothesis testing and provides a better understanding of data set variables and their interactions. It can also help determine if the statistical techniques you are considering for
data analysis are appropriate. IDA and EDA often employ data visualization methods to check distributional characteristics, identify relationships between variables, identify potential cofactors, and spot inconsistencies. 

We will run an IDA and EDA for both our phenotypic data and the genotypic data."

In PreGWAS step, I will also conduct the Genotype filtering and then the Imputation. For the Imputation, I will use Beagle.
