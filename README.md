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

To resume : We performed a Campbell-style analysis using FST, DXY, RAiSD μLD and tuskless heterozygosity deviation (i.e HetDev). Unlike Campbell et al., who used specific high-depth subsets and ANGSD genotype likelihoods for FST/DXY, our analysis was performed from the QC-filtered GWAS genotype dataset. Therefore, this analysis is used as a comparative overlap framework rather than an exact reproduction.

2 differents approaches are done : one to realize a Venn diagram based on 10kb window/bin level and one Venn diagram on SNP-level. For each metric, all the genomic windows are ordered from the strongest to the weakest signal. For example : I keep 5% windows with the highest Fst and same for the 1%, which makes this latter more restrictive.

## 2. Venn diagram 10kb window/bin level

This first Venn diagram overlap candidate windows/regions: **Campbell-style window/bin overlap analysis**

- FST     : windows 10 kb top 5% ou top 1%
- DXY     : windows 50 SNPs top 5% ou top 1%
- muLD    : windows RAiSD top 5% ou top 1%
- HetDev  : windows 10 SNPs top 5% ou top 1%

This diagram answer the question : Which genomics regions are common between analyses ? 

**<ins>Venn diagram 10kb top 5%<ins>** :

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/09.Campbell-style-Analyses/Campbell_corrected_Venn_10kb_bins_5pct_Campbell.png" />

556 bins common FST + DXY + muLD are find. But 0 common bins are returned between all metrics contrary to 305 found in Campbell et al. 2021.

**<ins>Venn diagram 10kb top 1%<ins>**: 

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/09.Campbell-style-Analyses/Campbell_corrected_Venn_10kb_bins_1pct_Campbell_threshold.png" />

For top 1% only 10 bin common FST + DXY + muLD are find.

## 3. Venn diagram SNP-level

A SNP would be “significant” for an analysis only if the SNP position is located inside a significant window : **SNP-level extraction from the significant Campbell-style windows**

- FST SNP       = SNP inside windows FST top 5% / 1%
- DXY SNP       = SNP inside windows DXY top 5% / 1%
- muLD SNP      = SNP inside windows RAiSD μLD top 5% / 1%
- HetDev SNP    = SNP inside windows HetDev top 5% / 1%

This SNP-level Venn diagram represents SNPs falling within the retained top-tail windows for each analysis. It answers the question: which SNPs are located in genomic windows shared by FST, DXY, μLD, and HetDev?

- Windows-level : CM044048.1:14920001–14930000 is common to FST+DXY+muLD
- SNP-level : Which SNPs are common between Fst, Dxy, muLD and HetDev ?

**For top 5%** : 
- FST&DXY&muLD              : 11,199 SNPs
- FST&DXY&HetDev_TX         : 20 SNPs
- FST&HetDev_TX             : 87 SNPs
- muLD&HetDev_TX            : 0 SNP

**For top 1%** : 
The total analysed SNPs is 21,295,223 with 1,794,892 SNPs in at least one analyse. 

- FST SNPs                  : 288,845
- DXY SNPs                  : 397,710
- muLD SNPs                 : 1,139,206
- HetDev_TX SNPs            : 68

- FST&DXY&muLD	            : 92 SNPs
- DXY&muLD	                : 3071 SNPs
- FST&muLD	                : 9922 SNPs
- FST&DXY	                  : 18036 SNPs

SNPs in 1 analyse       : 1,764,047
SNPs in 2 analyses      : 30,753
SNPs in 3 analyses      : 92
SNPs in 4 analyses      : 0

92 SNPs are returned within FST&DXY&muLD intersection, while none overlapping is detected with Fst, Dxy, muLD and HetDev for tuskless.

The main common signal is so Fst + Dxy + muLD, while HetDev is very rare and does not joined to muLD in the dataset.

## 4. LOC126069858 & AMELX

This two regions have been extrated and analysed to investigate the signals identified by Campbell et al. but also by the genome-wide scan (i.e. 03.GWAS & 04.Selection-Differentiation-scan).

**For LOC126069858**, the FST+DXY+muLD signal is nearby, but not inside the gene body. The bin 10 kb CM044048.1:14920001-14930000 overlap the edge of the gene, but can contain SNPs with the FST+DXY+muLD signal near to the gene.
So no signal is detected directly inside LOC126069858.

Inside the bin 10kb which overlap the beginning of the gene. The bin CM044048.1:14920001–14930000 contain 42 SNPs FST+DXY+muLD at 5%. The window is located between 14,920,064 and 14,922,466 so upstream of LOC126069858, which starts at 14,929,369.

**For AMELX**, No SNP have been found for signal FST+DXY+muLD in AMELX region. In the large region AMELX CM044047.1:167.5–171.0 Mb, no SNP with FST+DXY+muLD signal is detected. However, simple and double signals qre detected :

|**For top 5%** : | **For top 1%** :|
|-----------------------|----------------------|
|3083 SNPs muLD only  | 647 SNPs muLD seul|
|135 SNPs DXY only  | 92 SNPs FST seul|
|590 SNPs FST only  | 0 SNP DXY|
|307 SNPs FST + muLD  | 0 SNP FST + muLD|
|25 SNPs FST + DXY  | 0 SNP FST + DXY|
|0 SNP FST + DXY + muLD  | 0 SNP FST + DXY + muLD|

At the SNP level, no variant within the exact LOC126069858 gene body was shared by FST, DXY and μLD. However, 42 SNPs in the adjacent 10-kb bin immediately upstream of LOC126069858 were shared by FST, DXY and μLD at the 5% tail threshold. Within LOC126069858 itself, SNPs were mainly supported by FST and μLD, without DXY convergence. No FST+DXY+μLD SNP-level convergence was observed in the AMELX regional window.




