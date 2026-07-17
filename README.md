# Campbell-style Analyses on windows level

## 1.Introduction

Campbell et al. have generates a VCF from aligned reads using GATK protocol as I did for this study. Most of the downstream analyses use filtering genotype/variants via vcftools, cyvcf2 or RAiSD. For FST and DXY, they use ANGSD and likelihood genotypes to generate uncertainty related to sequencing depth.

I did not reproduce the full read-mapping, genotype-likelihood, and filtering pipeline from Campbell et al. Instead, I performed a Campbell-style analysis from our already QC-filtered GWAS genotype dataset: GWAS_qc_named_pheno.

Here a resume of the protocol for this section : 
- FST     : fenêtres 10 kb
-
- DXY     : fenêtres 50 SNPs, step 10 SNPs
-
- muLD    : RAiSD μLD
-
- HetDev  : SNPs TX 15–85%, absents chez TT, fenêtres 10 SNPs, step 2
-
- top 5%  : seuil principal pour le Venn
-
- top 1%  : seuil strict additionnel
