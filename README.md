# Pre-Genomic Wide Association Study - GWAS

## Introduction

Pre-GWAS step is a crucial step to generate the right file format but also to explore the genomic data. I have based this protocol mostly on the Physalia GWAS course but also on [tutorial](https://cloufield.github.io/GWASTutorial/01_Dataset/) and [plink webpage](www.cog-genomics.org/plink/1.9/filter). PLINK (v.1.90b6.24), R, and biokit were used during the Pre-GWAS steps. As mentionned in the main README file, the pre-gwas step concerns data pre-processing with initial and explonatory data analysis (IDA and EDA respectively). This two first steps generated and used the files for data cleaning, descriptive informations and transformation to ensure to the right conformity of the files for gwas models, but also to examine and summarise data, discovery patterns or structuration that have to be taken into account for the model configuration. 

## 01. IDA & EDA
### _1.1. Genomics data_
The raw genomic data called "all_contigs.raw.vcf.gz" [See 01.WGS-Variant-Calling] is used to generate the first genomic files as output of the first script [01_init_fullgenome_plink_qc.slurm]. PLINK conversion has been applied on the whole genome with a total genotyping rate of 0.989733. The genomics metrics, such as frequencies, allele counting, genomics frequencies, the missingness, Hardy-Weinberg and Heterozygosity files were computed. From these files, the corresponding explonatory figures and tables were generated. Contrary to the Physalia course manual, I do not preselected chromosomes to explore all the Asian elephant genomes. The EDA/quality control is done with R in HPC (i.e CSC/Puhti) [SEE 02.global_EDA_QC.R] to generate tables and visualize the genomics metrics. 

### _1.2. Phenotypic data_
The Phenotypic data were also explored in the same Rscript, but also verified in local machine [See 02.global_EDA_QC.R]. 


In PreGWAS step, I will also conduct the Genotype filtering and then the Imputation. For the Imputation, I will use Beagle.


Similar as [Campbell et al. (2021)](https://www.science.org/doi/full/10.1126/science.abe7389), I focused the analyses on the whole genome and not only chromosome 1 and X, to identify the regions 
