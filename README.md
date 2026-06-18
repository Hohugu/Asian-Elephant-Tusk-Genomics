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

> ***Tusk ~ SNP + sex + GRM + population + PC1 + PC2 + PC3***
> 
The script of this second GWAS model are available [SEE 10.prepare_GEMMA_input.R & 11.GEMMA_mixed_model_PC3.slurm]. 

### _3.1 Genomic Control Inflation Factor & Plots_

Le lambda GC is equal to 0.858, indicating a slightly conservative model and without extreme artefact of Plink PC2 as observing in Manhattan plot. The GEMMA QQ plot shows a strong deviation. This deviation would come from that some rare variants with strong effect are present or that artefacts remain. Script for the plots [SEE 12.GEMMA_Manhattan_QQ_PC3.R]

<img width="630" height="730" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GEMMA_PC3_Manhattan.png" /> <img width="315" height="370" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GEMMA_PC3_QQ.png" />

### _3.2 Bonferroni SNPs_

To identify significant SNP which may associated to phenotype, I used the ***Bonferroni threshold*** : 

> PBonf = 0.05 / Ntests  with Ntests = number of total variants tested, 21 283 317
> 
> PBonf = 2.35e-09
> 
> -log10(PBonf) = 8.63

So all SNPs with -log10(P) superior to 8.63 are genome-wide significants. **262 significants SNPs are superior to 8.63**. 
Each SNPs seem to be part of LD blocs by chromosomes. I am particularly interested by SNPs in chromosome 1, knowing that MEP1a is also present in chromosome 1. 

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
| LOC126067332 | LOC126067332 | NC_064843.1 | CM044044.1 | 33549010 | 33578813 | CM044044.1:33549161:A:C | uncharacterized | 8.542653e-13 |
| LOC126064930 | CALB2 | NC_064839.1 | CM044040.1 | 39119360 | 39147775 | CM044040.1:39127397:C:T | calbindin 2 | 5.851082e-13 |
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

Suggestive signals are identified close to the candidates genes. Suggestive signals mean that the Bonferroni threshold is not reached but the signals still strong. These signals were extracted by taking the id of SNPs to have the exact position and then used the distance from the candidate genes to verify the position and the related gene name. The distance is calculated compared to candidat gene coordonates. 
For the 500kb window, no Bonferroni-significant association overlapped the Campbell et al. candidate genes. ODAM (odontogenic ameloblast) associated showed a weaker local signal (best p = 1.64E-06 within 500 kb downtream, and best p = 1.85E-08 within 1Mb upstream). However, suggestive associations (p < 1E-05) were detected within 500 kb of all six candidate loci, with the strongest signal approximately 171–271 kb upstream of the ENAM/AMBN/AMTN enamel gene cluster (p = 1.85E-08). This latter signal remains stronger than the majority of the genomic background below Bonferroni threshold and than genes located in the region containing enamel or teeth genes. **Further analysis have to be made to detect LD between ENAM gene and the suggestive association**.

**Table of Closest suggestive SNP to candidat genes in Asian elephants.**

| Candidate Gene | Closest suggestive SNP   |      p-value |                  Distance from candidate gene | SNP chromosome           | Gene containing SNP | Gene coordinates        | Annotation                           |
| -------------- | ------------------------ | -----------: | --------------------------------------------: | ------------------------ | ------------------- | ----------------------- | ------------------------------------ |
| ENAM           | CM044024.1:87145701:G:T  | 1.846930e-08 |                           171,183 bp upstream | NC_064823.1 / CM044024.1 | LOC126077101        | 87,144,107–87,162,148   | G-rich RNA sequence binding factor 1 |
| AMBN           | CM044024.1:87145701:G:T  | 1.846930e-08 |                           208,589 bp upstream | NC_064823.1 / CM044024.1 | LOC126077101        | 87,144,107–87,162,148   | G-rich RNA sequence binding factor 1 |
| AMTN           | CM044024.1:87145701:G:T  | 1.846930e-08 |                           271,084 bp upstream | NC_064823.1 / CM044024.1 | LOC126077101        | 87,144,107–87,162,148   | G-rich RNA sequence binding factor 1 |
| ODAM           | CM044024.1:88318180:A:G  | 1.639486e-06 |                         476,845 bp downstream | NC_064823.1 / CM044024.1 | LOC126077138        | 88,320,890–88,369,850   | Sulfotransferase 1E1-like            |
| AMELX          | CM044047.1:168377782:C:T | 9.656607e-08 | 347,384 bp upstream (inside neighboring gene) | NC_064846.1 / CM044047.1 | LOC126069583        | 168,359,757–168,861,890 | Rho GTPase activating protein 6      |
| MEP1A          | CM044020.1:111267040:C:T | 1.321988e-05 | 412,908 bp upstream (inside neighboring gene) | NC_064819.1 / CM044020.1 | LOC126081574        | 110,947,996–111,303,545 | Regulator of calcineurin 2           |

Moreover, a SNP located on chromosome NC_064846.1 on gene lncRNA LOC126069593 (CM044047.1:170145312:G:T) is returned significant from GEMMA GWAS mixed-model (p = 1.56E-10). This region is detected 1.42Mb downstream of AMELX (amelogenin X-linked) candidate gene, on the same chromosome, previously involved in tusk morphology in African elephants (Campbell et al., 2021). No Bonferroni-significant SNP has been found directly in AMELX, but GEMMA GWAS signal is detected inside the same chromosomic bloc and might reflect a local association.

**Table of Closest Bonferroni SNP to AMELX candidat genes in Asian elephants.**

| Candidate Gene |   Closest Bonferroni SNP |      p-value | Distance from candidate gene | SNP chromosome           | Gene containing SNP | Gene coordinates        | Annotation             |
| -------------- | -----------------------: | -----------: | ---------------------------: | ------------------------ | ------------------- | ----------------------- | ---------------------- |
| AMELX          | CM044047.1:170145312:G:T | 1.563143e-10 |      1,415,894 bp downstream | NC_064846.1 / CM044047.1 | LOC126069593        | 170,107,747–170,332,119 | uncharacterized lncRNA |

**A functional annotation on this SNP and a LD analysis are needed to asses its potential link with AMELX region**.

The Campbell candidate gene analysis was expanded to include MEP1B and PLA2G7. Among the candidate genes, MEP1B emerged as one of the strongest local signals. In the mixed GWAS, MEP1B showed a stronger within-gene association than MEP1A (p = 5.34E-03 versus p = 1.52E-02) and comparable suggestive associations within larger windows.

## 4. GEMMA only-males GWAS mixed-model

The only-male GWAS is similar to the previous mixed-model GWAS. For males I analyzed 46 males with 20,405,248 SNPs, and a pve estimate of 0.99, meaning that the model explained almost all the observed variance. 31 SNPs returned significants for GWAS GEMMA only-males. This result is statistically consistent given that the sample size is reduced to 46 males while the previous GWAS was with 93 individuals. There is less statistic power, less observed recombination and a higher variance of the estimates. But it may also mean that these signals are robusts especially whether the same SNPs as the previous GWAS are returned. 

Population structure is taken into account and 3 PCs are kept as before. Kinship is also corrected. The equation of GEMMA GWAS model for male is :

>y=Wα+xβ+u+ϵ
>
>W = covariables (population + PC1 + PC2 + PC3)
>
>x = SNP tested
>
>u∼N(0,Kσ²g) = GRM
>
>K = kinship matrix

The script for the third Gemma model : [SEE 13.Prepare_GEMMA_male_PC3_inputs.R & 14.GEMMA_mixed_model_male_PC3.slurm]

### _4.1 Genomic Control Inflation Factor & Plots_

The lambda GC is equal to 0.93, meaning that the model is slightly conservative and well calibrated. 

<img width="630" height="730" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GEMMA_male_PC3_Manhattan_highlight.png" /> <img width="315" height="370" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/03.GWAS/GEMMA_male_PC3_QQ_pwald.png" />

The QQ plot shows a good calibration pattern, following the expected line at the beginning and then a strong deviation at the tail. The majority of the SNPs follow the expected but a small fraction presents a strong signal. For the Manhattan plot, a main signal is detected with a SNP (CM044020.1:85883031:A:T) with a pvalue at 1.25E-16 and a -log10(P) = 15.90. 

Script for the plots [SEE 15.GEMMA_male_Manhattan_PC3_QQ.R]

### _4.2 Bonferroni SNPs_

For the GEMMA mixed-model with only males and 3PCs, the ***Bonferroni threshold*** : 

> PBonf = 0.05 / Ntests  with Ntests = number of total variants tested, 20 405 248
> 
> PBonf = 2.45E-09
> 
> -log10(PBonf) = 8.61

So all SNPs with -log10(P) superior to 8.61 are genome-wide significants. **31 significants SNPs are superior to 8.61**. 
Each SNPs seem to be part of LD blocs by chromosomes. I am particularly interested by SNPs in chromosome 1. 

From 31 SNPs, 26 SNPs are inside coding genes and therefore have gene description and functions available from the gff annotation file. 

### _4.3 Annotation & Candidate genes_

For each Campbell candidate gene (i.e. AMBN, AMTN, ENAM, ODAM, AMELX, MEP1A, MEP1B and PLA2G7), the pipeline [SEE 17.Campbell_candidates_v2_all_and_male.sh] searches for the best GWAS SNP inside the gene and within ±500 kb, ±1 Mb and ±2 Mb windows. It also extracts all Bonferroni-significant SNPs falling in these regions.

Results reproduced the original mixed-model findings. In the male-only GWAS, the enamel-cluster signal persisted but was substantially weaker (best p = 8.1E-05). No Bonferroni-significant SNPs were detected within 2 Mb of any Campbell candidate gene opposite to AMELX associated SNP (CM044047.1:170145312:G:T; p = 1.56E-10) from GEMMA GWAS mixed-model. These results indicate that the Campbell candidate loci do not overlap the strongest male-only GWAS signals and suggest that the AMELX-associated signal observed in the mixed GWAS is not driven by males alone.

In the male-only GWAS, the strongest local association among the Campbell candidates was detected near MEP1B (p = 1.76E-05, approximately 53 kb from the gene). These results indicate that MEP1B deserves to be considered alongside MEP1A as a biologically relevant candidate locus for tusk variation in Asian elephants. **Overall, none of the Campbell candidate genes directly overlap the strongest male-only GWAS peaks, but several loci—particularly MEP1B and AMELX—remain promising targets for downstream LD, population genetic, and regulatory analyses**.



