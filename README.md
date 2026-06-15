# Genomic Wide Association Study (GWAS)

## 1. Introduction
GWAS is commonly used to test the association between a phenotype and genotypes. For that several models are proposed. The GWAS protocol is based on **cloudfiel.github.io/GWASTutorial**. As we can see on the _following picture_, GWAS requires initialisation, filtering, QC and PCA steps. 

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

However, by trying with 2 and 3 PCs, results return different, with lambda GC of 0.936 and 0.904 respectively. 

<img width="315" height="410" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/QQ_PC2.png" /> <img width="315" height="370" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/QQ_PC3.png" />

Visually Q-Q plots inform that models are lightly conservative with correcte calibration for 99.99% of SNPs. Therefore, I will check the Manhattan plots:

### _2.2 Manhattan plots_

<img width="315" height="410" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/Manhattan_PC2.png" /> <img width="315" height="370" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/Manhattan_PC3.png" />

Manhattan plot pour 2 PCs shows a SNP very strongly differentiate between groups (P = 4.3x10^-198 and 3.2x10^-45), which is after verification a modelisation's artefact. The disappearance of these signals on the 3PCs plot supports the hypothesis that the PC2's model is distorted by artificial SNP. **GWAS model with 3 PCs remains the most likely model**.

## 3. GEMMA GWAS mixed-model

The mixed-model GWAS is implemented with GEMMA, which provides standard modern WGS mixed-model GWAS, efficient and well suited for structured populations. The aim is to correct for kinship and small population effects. Mixed models are important here because population-associated SNPs may appear falsely associated with  tusk phenotype. For this, I compute a genomic relationship matrix (GRM) and add covariables as following : 

> ***Tusk ~ SNP + sex + GRM + population + PCs***

### _3.1 Genomic Control Inflation Factor & Plots_

Le lambda GC is equal to 0.858, indicating a slightly conservative model and without extreme artefact of Plink PC2 as observing in Manhattan plot. The GEMMA QQ plot shows a strong deviation. This deviation would come from that some rare variants with strong effect are present or that artefacts remain. 

<img width="630" height="730" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GEMMA_PC3_Manhattan.png" /> <img width="315" height="370" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GEMMA_PC3_QQ.png" />

### _3.2 Bonferroni SNPs_

To identify significant SNP which may associated to phenotype, I used the ***Bonferroni threshold*** : 

> PBonf = 0.05 / Ntests  with Ntests = number of total variants tested, 21 283 317
> 
> PBonf = 2.35e-09
> 
> -log10(PBonf) = 8.63

So all SNPs with -log10(P) superior to 8.63 are genome-wide significants. **262 significants SNPs are superior to 8.63**. 
Each SNPs seem to be part of LD blocs by chromosomes. I am particularly interested by SNPs in chromosome 1 and X/Y. 

From 262 SNPs, only 121 SNPs are inside coding genes and therefore have gene description and functions available from the gff annotation file. 
Genes with highest significants SNPs are not directly related to tusk development as [Campbell et al. (2021)](https://www.science.org/doi/full/10.1126/science.abe7389) found for African elephants with AMELX and MEP1a genes. These genes are involved in enamel and dentin expression and development. Further analyses will be made following Campbell et al protocol. I detected real loci associated to tusk morphology. However, with 93 individuals and an unbalanced binary phenotype, some signals may reflected a family or geographic structure. Nethertheless, Lambda GC (0.858) suggests that GEMMA would over corrected than inflate statistics, which lead main signals very likely. 

### _3.3 Annotation & Candidate genes_

After identification of Bonferroni significant SNPs, I have searched for genes names and functions related to the significant SNPs. A total of 47 unique genes overlapped Bonferroni-significant GEMMA associations.

**Table of Bonferroni-significant genes identified by GEMMA mixed-model GWAS (PC1–PC3)**

| Gene name | Gene symbol | NC chr | CM chr | Start | End | Best SNP | Description | Best p-value |
|------------|--------|------------|------------|------------|------------|------------|------------------|------------|
| LOC126073692 | IL6R | NC_064821.1 | CM044022.1 | 192207359 | 192266911 | CM044022.1:192234282:C:G | interleukin 6 receptor | 2.115881e-19 |
| LOC126080233 | LOC126080233 | NC_064825.1 | CM044026.1 | 106830358 | 106864305 | CM044026.1:106837590:A:G | olfactory receptor 8K3-like | 4.617250e-18 |
| LOC126057365 | PDGFC | NC_064831.1 | CM044032.1 | 23468358 | 23734364 | CM044032.1:23712353:C:T | platelet derived growth factor C | 6.904120e-15 |
| LOC126058333 | GTF2F2 | NC_064832.1 | CM044033.1 | 34951228 | 35137192 | CM044033.1:35093113:A:C | general transcription factor IIF subunit 2 | 6.904120e-15 |
| C5H4orf19 | PGCKA1 | NC_064823.1 | CM044024.1 | 123449792 | 123553816 | CM044024.1:123473250:A:G | PDCD10 and GCKIII kinases associated 1 | 1.200328e-14 |
| LOC126087709 | MYZAP | NC_064831.1 | CM044032.1 | 47712552 | 47805874 | CM044032.1:47769846:A:G | myocardial zonula adherens protein | 8.652342e-13 |
| LOC126067331 | CNBD2 | NC_064843.1 | CM044044.1 | 33506136 | 33548234 | CM044044.1:33525949:A:C | cyclic nucleotide binding domain containing 2 | 8.542653e-13 |
| LOC126070942 | LOC126070942 | NC_064821.1 | CM044022.1 | 17357210 | 17815040 | CM044022.1:17533772:C:T | olfactory receptor 7G2-like | 8.542653e-13 |
| LOC126067332 | LOC126067332 | NC_064843.1 | CM044044.1 | 33549010 | 33578813 | CM044044.1:33549161:A:C | uncharacterized | 8.542653e-13 |
| LOC126064930 | CALB2 | NC_064839.1 | CM044040.1 | 39119360 | 39147775 | CM044040.1:39127397:C:T | calbindin 2 | 5.851082e-13 |
| LOC126066272 | EDIL3 | NC_064820.1 | CM044021.1 | 89441686 | 89973518 | CM044021.1:89746872:A:T | EGF like and discoidin domains 3 | 1.238987e-12 |
| LOC126071946 | GTPBP8 | NC_064819.1 | CM044020.1 | 1702397 | 1725352 | CM044020.1:1710266:A:G | GTP binding protein 8 (putative) | 1.296276e-12 |
| LOC126079052 | LOC126079052 | NC_064825.1 | CM044026.1 | 117237266 | 117238198 | CM044026.1:117237421:A:G | olfactory receptor 5M3-like | 1.314175e-12 |
| LOC126080471 | LOC126080471 | NC_064825.1 | CM044026.1 | 115144203 | 115163744 | CM044026.1:115160186:A:G | olfactory receptor 10AG1-like | 1.314175e-12 |
| LOC126073164 | LOC126073164 | NC_064819.1 | CM044020.1 | 90003377 | 90377483 | CM044020.1:90108700:A:C | patr class I histocompatibility antigen A126 | 1.292884e-13 |

You can find the rest of the table in **Table of Bonferroni-significant genes**. I noticed that some genes are uncharacterized while they are significant.

Then, I retrieved the candidates genes from Campbell et al., 2021 in the **GGF.file** reference genome of Asian elephants to find their positions and their associated p-values from the GEMMA mixed-model (i.e. GWAS GEMMA-PC3). None of the candidat genes do not return significants and so not directly associated to tusk phenotype from the GWAS GEMMA-PC3 (< 2.35e-09) [SEE 12.bis_Campbell_candidates_GEMMA_PC3.sh]. 

**Table of orthologue candidat genes in Asian elephants.** 

| Candidat Gene  | Gene symbol | CHR | Position | Description | Best p-value |
| ----- | ------ | ------ | ------| --------------- | ------- |
| ODAM  | ODAM_LOC126077139 | NC_064823.1 | 87832910-87841335 | odontogenic, ameloblast associated | 3.37e-03 |
| MEP1A | MEP1A_LOC126069657 | NC_064819.1 | 111679947-111709525 | meprin A subunit alpha | 1.52e-02 |
| AMTN  | AMTN_LOC126076721 | NC_064823.1 | 87416785-87432258 | amelotin | 5.02e-02 |
| ENAM  | ENAM_LOC126077098 | NC_064823.1 | 87316884-87332008 | enamelin | 5.10e-02 |
| AMELX | AMELX_LOC126069472 | NC_064846.1 | 168725165-168729418 | amelogenin X-linked | 8.00e-02 |
| AMBN  | AMBN_LOC126076719 | NC_064823.1 | 87354290-87366652 | ameloblastin | 9.55e-02 |

I checked the closest Bonferroni hits from the candidates genes in a window of 500k, 1M and 2M [SEE 12.ter_Annotate_nearby_Bonferroni_hits.sh]. 
For the 500kb window, despite the non-significant association with tusk morphology, a **suggestif signal** is observed near to ENAM/AMBN/AMTN cluster. This signal remains stronger than the majority of the genomic background below Bonferroni threshold and than genes located in the region containing enamel or teeth genes.

| Candidate Gene | Closest suggestive SNP | p-value | Distance from candidate gene | SNP chromosome | Gene containing SNP | Gene coordinates | Annotation |
|---|---:|---:|---:|---|---|---|---|
| ENAM | CM044024.1:87145701:G:T | 1.846930e-08 | 171,183 bp | NC_064823.1 / CM044024.1 | Intergenic | NA | Upstream of ENAM/AMBN/AMTN enamel gene cluster |
| AMBN | CM044024.1:87145701:G:T | 1.846930e-08 | 208,589 bp | NC_064823.1 / CM044024.1 | Intergenic | NA | Upstream of AMBN, within ENAM/AMBN/AMTN enamel gene cluster region |
| AMTN | CM044024.1:87145701:G:T | 1.846930e-08 | 271,084 bp | NC_064823.1 / CM044024.1 | Intergenic | NA | Upstream of AMTN, within ENAM/AMBN/AMTN enamel gene cluster region |
| AMELX | CM044047.1:168377782:C:T | 9.656607e-08 | 347,383 bp | NC_064846.1 / CM044047.1 | Intergenic | NA | Upstream of AMELX |
| ODAM | CM044024.1:88318180:A:G | 1.639486e-06 | 476,845 bp | NC_064823.1 / CM044024.1 | Intergenic | NA | Downstream of ODAM, near ENAM/AMBN/AMTN/ODAM enamel gene region |
| MEP1A | CM044020.1:111267040:C:T | 1.321988e-05 | 412,907 bp | NC_064819.1 / CM044020.1 | Intergenic | NA | Upstream of MEP1A |

Moreover, a SNP located on chromosome NC_064846.1 on gene lncRNA LOC126069593 (CM044047.1:170145312:G:T) is returned significant from GEMMA GWAS mixed-model (p = 1.56e-10). This region is at 1.42Mb from AMELX (amelogenin X-linked) candidate gene, on the same chromosome, previously involved in tusk morphology in African elephants (Campbell et al., 2021). No Bonferroni-significant SNP has been found directly in AMELX, but GEMMA GWAS signal is detected inside the same chromosomic bloc and might reflect a local association.

| Candidate gene | Closest Bonferroni SNP   |      p-value | Distance from candidate gene | SNP chromosome | Gene containing SNP | Gene coordinates    | Annotation      |
| -------------- | ------------------------ | -----------: | ---------------------------: | -------------- | ------------------- | ------------------- | --------------- |
| AMELX          | CM044047.1:170145312:G:T | 1.563143E-10 |                 1,415,894 bp | NC_064846.1    | LOC126069593        | 170107747–170332119 | uncharacterized |

A functional annotation on this SNP and a LD analysis are needed to asses its potential link with AMELX region.

## 4. GEMMA only-males GWAS mixed-model

The only-male GWAS is similar to the previous mixed-model GWAS. 

### _4.1 Genomic Control Inflation Factor & Plots_









