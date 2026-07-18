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

To resume : We performed a Campbell-style analysis using FST, DXY, RAiSD μLD and tuskless heterozygosity deviation (i.e HetDev). Unlike Campbell et al., who used specific high-depth subsets and ANGSD genotype likelihoods for FST/DXY, our analysis was performed
from the QC-filtered GWAS genotype dataset. Therefore, this analysis is used as a comparative overlap framework rather than an exact reproduction.

2 differents approaches are done : one to realize a Venn diagram based on 10kb window/bin level and one Venn diagram on SNP-level.

## 2. Venn diagram 10kb window/bin level

This first Venn diagram overlap candidate windows/regions: **Campbell-style window/bin overlap analysis**

FST     : windows 10 kb top 5% ou top 1%
DXY     : windows 50 SNPs top 5% ou top 1%
muLD    : windows RAiSD top 5% ou top 1%
HetDev  : windows 10 SNPs top 5% ou top 1%

This diagram answer the question : Which genomics regions are common between analyses ? 

__Venn diagram 10kb on 5pct__ :
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/09.Campbell-style-Analyses/Campbell_corrected_Venn_10kb_bins_5pct_Campbell.png" />

__Venn diagram 10kb on 1pct__ : 
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/09.Campbell-style-Analyses/Campbell_corrected_Venn_10kb_bins_1pct_Campbell_threshold.png" />

## 3. Venn diagram SNP-level

A SNP would be “significant” for an analysis only if the SNP position is located inside a significant window : **SNP-level extraction from the significant Campbell-style windows**

FST SNP       = SNP inside windows FST top 5% / 1%
DXY SNP       = SNP inside windows DXY top 5% / 1%
muLD SNP      = SNP inside windows RAiSD μLD top 5% / 1%
HetDev SNP    = SNP inside windows HetDev top 5% / 1%

This new Venn diagram represents the SNPs overlapping retained for each analyses, and answer the question : Which SNPs are common between Fst, Dxy, muLD and HetDev ? 

- Windows-level : CM044048.1:14920001–14930000 is common to FST+DXY+muLD
- SNP-level : Which SNPs are common between Fst, Dxy, muLD and HetDev ?

  
  
