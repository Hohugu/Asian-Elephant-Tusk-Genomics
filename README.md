# Campbell-style Analyses on windows level

## 1.Introduction

Campbell et al. have generates a VCF from aligned reads using GATK protocol as I did for this study. Most of the downstream analyses use filtering genotype/variants via vcftools, cyvcf2 or RAiSD. For FST and DXY, they use ANGSD and likelihood genotypes to generate uncertainty related to sequencing depth.

I did not reproduce the full read-mapping, genotype-likelihood, and filtering pipeline from Campbell et al. Instead, I performed a Campbell-style analysis from our already QC-filtered GWAS genotype dataset: GWAS_qc_named_pheno.

Here a resume of the protocol for this section : 
- FST     : windows of 10 kb
- DXY     : windows of 50 SNPs, step 10 SNPs
- muLD    : RAiSD μLD
- HetDev  : SNPs TX 15–85%, missing in TT, windows of 10 SNPs, step 2
- top 5%  : main threshold for Venn
- top 1%  : additional strict threshold

Campbell et al. used 13 high-depth samples with 7 tuskless and 6 tusked and then 18 samples for ANGSD FST / DXY.

For this study, I used 93 individuals : 53 TX, 40 TT. This dataset is therefore different from Campbell and already QC.
FST and DXY were calculated from called genotypes rather than ANGSD genotype likelihoods.

For RAiSD μLD, Campbell et al. used RAiSD on high-depth tuskless samples, while I used 53 tuskless. Campbell et al. explain that HetDev 10 SNPs windows corresponded to approximately 10kb. However, some HetDev windows in our dataset are several Mb long due to the rarety of filtering SNPs. I produced a version with span <= 10 kb.

To resume : We performed a Campbell-style analysis using FST, DXY, RAiSD μLD and tuskless heterozygosity deviation. Unlike Campbell et al., who used specific high-depth subsets and ANGSD genotype likelihoods for FST/DXY, our analysis was performed
from the QC-filtered GWAS genotype dataset. Therefore, this analysis is used as a comparative overlap framework rather than an exact reproduction.

## 2.


  
