# Pre-Genomic Wide Association Study - GWAS

## Introduction

Pre-GWAS step is a crucial step to generate the right file format but also to explore the genomic data. I have based this protocol mostly on the Physalia GWAS course but also on [tutorial](https://cloufield.github.io/GWASTutorial/01_Dataset/) and [plink webpage](www.cog-genomics.org/plink/1.9/filter). PLINK (v.1.90b6.24), R, and biokit were used during the Pre-GWAS steps. As mentionned in the main README file, the pre-gwas step concerns data pre-processing with initial and explonatory data analysis (IDA and EDA respectively). This two first steps generated and used the files for data cleaning, descriptive informations and transformation to ensure to the right conformity of the files for gwas models, but also to examine and summarise data, discovery patterns or structuration that have to be taken into account for the model configuration. The Pre-GWAS also included the filtering step, and post filtering quality control on raw data (i.e including both sexes) and among males (i.e comparison between tusked and tuskless males). The LD pruning and PCA (population structure) are crucial steps to include control variables for GWAS models. 

## 01. IDA & EDA
_Genomics data_ : The raw genomic data called "all_contigs.raw.vcf.gz" [See 01.WGS-Variant-Calling] is used to generate the first genomic files as output of the first script [01_init_fullgenome_plink_qc.slurm]. PLINK conversion has been applied on the whole genome with a total genotyping rate of 0.989733. The genomics metrics, such as frequencies, allele counting, genomics frequencies, the missingness, Hardy-Weinberg and Heterozygosity files were computed. From these files, the corresponding explonatory figures and tables were generated. Contrary to the Physalia course manual, I do not preselected chromosomes to explore all the Asian elephant genomes. The EDA/quality control is done with R in HPC (i.e CSC/Puhti) [SEE 02.global_EDA_QC.R] to generate tables and visualize the genomics metrics. 
_Phenotypic data_ : The Phenotypic data were also explored in the same Rscript, but also verified in local machine [See 02.global_EDA_QC.R]. 

After EDA and QC, for 93 individuals, 28M variants are presents in the raw genotype data. The mean sample and variant missing are both equal to 1.03% with a median variant missing of 0, meaning that the majority of SNPs are complete and that sequencing and global genotyping is good. I detect that the phenotype is extremely correlated to sexe, as only males show a variability in tusk (i.e tusked or tuskless) while females are tuskless. This first observation allows to conlude on a strong sexual dimorphism on tusk phenotype. Statistically, this confounding factor may induced tusk signal to sexe signal. The three following graphs show that the heterozygosity rate between females and males is the same in average, despite a large variation for males. This figure means that both sexes dispose of an similare average level of heterozygosity. The same can be observed regarding the level of heterozygosity between tusk types males. 

<img width="330" height="470" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/02.Pre-GWAS/tusk_by_sex.png" /> <img width="330" height="470" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/02.Pre-GWAS/heterozygosity_by_sex.png" /> <img width="330" height="470" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/02.Pre-GWAS/heterozygosity_by_tusk.png" />

The measure of the level of heterozygosity between India and Myanmar shows a real difference between both populations.
Heterozygoty returns a mean F = 0.0557 and the median = 0.0357. These measures reveal an exces of heterozygoty which may reflect population structure, high diversity, repeated regions, copy number variation (CNV), or a biological variation. **Therefore, this reflections will have to be explore later**. 

## 02. Genotype filtering and Post EDQ/QC

### _2.1. Filtering_

The filtering step is realized at this moment of the pipeline based on filtering markers (minor allele frequency MAF, missing calls per marker and missing calls per individuals). MAF measures the frequence of the variant in the population (e.g. MAF = 0.5, both alleles are balanced; MAF = 0.01, rare variant). MAF filters out the rare variants because rare SNPs are often noised, few informatives, give false positives. F_miss measures the percentages of missing genotypes per individuals (e.g. F_MISS = 0.40, means that the individuals have 40% of missing genotype, so bad quality). The lmiss measures the missing rate per variant (e.g. F_MISS = 0.60, means that SNP absent on 60% of individuals). Based on the EDA results, I started the filtering step with : 
            
            > _Variants with a minor allele frequency below 5% : **geno = 0.05**_
            > _Samples with missing call rates exceeding 5% : **mind = 0.05**_
            > _Variants with missing call rates exceeding 1% : **maf = 0.01**_




In PreGWAS step, I will also conduct the Genotype filtering and then the Imputation. For the Imputation, I will use Beagle.


Similar as [Campbell et al. (2021)](https://www.science.org/doi/full/10.1126/science.abe7389), I focused the analyses on the whole genome and not only chromosome 1 and X, to identify the regions 
