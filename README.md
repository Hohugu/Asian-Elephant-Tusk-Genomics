# Selection and population differentiation scans

## 1. Introduction

This part investigates genomic differentiation between **tusked (TT)** and **tuskless (TX)** Asian elephants using population-genomic summary statistics. The aim is to identify genomic regions showing consistent differentiation between phenotype groups, independently of the GWAS association tests. 

The analyses focused on several complementary signals : 

- Genome-wide FST scans, to identify differences in allele frequencies between TT qnd TX groups.

- Observed heterozygosity (Ho) to detect regions with local reductions in within-group heterozygosity.

- Nucleotide diversity comparisons, to identify local reductions in genetic diversity.
  
- Tajima's D, to identify differences in the allele-frequency spectrum
  
- Local linkage disequilibrium (LD) and allele-frequency differences for regional refinement

- Candidate-region annotation for SNP and gene-level annotation

The objective is not to demonstrate selection directly, but to identify genomic regions showing convergent differentiation signals that can be prioritized for downstream biological interpretation. 

Three complementary analysis designs are considered : 
  **(i).** An initial genome-wide TT vs TX scan using the full available sample.
  **(ii).** A male-only analysis using all available males.
  **(iii).** A population-balanced male-only sensitivity analysis based on repeated 6 TT vs 6 TX resampling.

The distinction between these designs is important because the original phenotype groups were strongly imbalanced with respect to sex, whereas the male-only analyses remove this major source of confounding.

## 2. Initial all-individual genome-wide scan

The first population-genomic scan compared all available tusked and tuskless individuals. This analysis provided the initial genome-wide FST, Ho, pi and Tajima's D candidate regions.

However, the TT and TX groups were not balanced for sex: the tusked group was male, whereas most tuskless individuals were female. Consequently, some differentiation signals could reflect sex-linked structure rather than tusk phenotype alone.
For this reason, these results are retained as the initial exploratory genome-wide scan, but they are complemented by male-only analyses.

### 2.1 Genome-wide FST

Genome-wide FST was calculated in genomic windows to identify regions with elevated allele-frequency between tusked and tuskless elephants [SEE 01.FST_genomewide_TT_vs_TX.sh - 03.FST_make_windows.R]. FST scans identified genomic windows showing elevated differentiation between tusked and tuskless individuals [SEE 04.plot_FST_genomewide.R]. Then I summarized FST values in genomic windows, identified candidate FST-enriched regions, and annotated these regions with nearby genes [SEE 05.FST_merge_candidate_regions.R - 07.make_FST_candidate_regions_annotated.R].
These regions were used as the first layer of evidence for selection/differentiation candidate regions.

FST candidate regions were also compared with GWAS signals to test whether differentiated regions are located near GWAS-associated SNPs or candidate GWAS loci [SEE 08.compare_FST_vs_GWAS.R and 09.compare_FST_vs_GWAS_distance..R]. This step helped to distinct regions supported only by differentiation scans from regions also close to GWAS signals. These comparisons were later used during candidate-locus prioritization, but they are not equivalent to GWAS evidence.

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Genomewide_FST_1Mb_annotated.png" />
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig1_FST_50kb_GEMMA_Campbell_annotated.png" />

SNPs with the higher FST values are located in chromosome 1, 3, 16 and among sexual chromosomes. Most GEMMA Bonferroni SNPs did not fall within the strongest FST windows, but some GWAS signals were located near Campbell candidate genes and close to differentiated regions.
A high FST value usually means that allelic frequencies between tusked and tuskless individuals are differents. Thereby, if a GWAS GEMMA SNPs have a high FST value, then this GWAS SNPs is located inside a region where tusked and tuskless individuals are differents. Moreover, a Bonferroni SNPs near to Campbell candidate gene and inside a high FST window, would suggest spatial concordance between GWAS association, TT/TX differentiation and prior tooth/tusk candidate genes.

### 2.2 Observed Heterozygosity  

Observed heterozygosity was estimated to identify regions with reduced within-group genetic diversity [SEE 10.heterozygosity_TT_TX.sh - 15.plot_genomewide_delta_Ho.R].
The comparison is expressed as : **ΔHo = Ho_TX - Ho_TT** 

- **ΔHo < 0**, meaning that tuskless individuals display lower heterozygoty than tusked individuals.
- **ΔHo > 0**, tusked individuals have lower heterozygosity than tuskless individuals. A local decrease of the heterozygosity can be associated with recent selection, frequent haplotype in a group, or a local lost of diversity.
- For a potential selection in tuskless individuals, I should observe a high FST and a ΔHo < 0, concluding on a real difference in the associated region and that tuskless individuals carry less local diversity.

A local decrease in heterozygosity can be compatible with recent selection or the presence of a frequent local haplotype, but can also be produced by demographic structure, drift or technical effects It is therefore interpreted only as one component of a multi-metric signal

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/TT_TX_delta_Ho_50kb_genomewide.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig2_delta_Ho_50kb_GEMMA_Campbell_annotated.png" />

Similar to mean FST graphs, SNPs with significant heterozygosity values are located in chromosome 22, 32 and among sexual chromosomes. 

The final candidate selection/differentiation regions showed reduced heterozygosity in tuskless elephants.

### 2.3 Nucleotide diversity

Nucleotide diversity was calculated for each tusk types and summarized as : 
**Δπ = π_TX - π_TT** used to identify regions with one group showing a reduced local diversity 

- **Δπ <0**: tuskless individuals have a lower nucleotidic diversity.
- **Δπ >0**: tusked individuals have a lower nucleotidic diversity. 

[SEE 16.compare_pi_TT_TX.R - 18.annotate_reduced_pi_FST_overlap.R]. For a recent selection in TX, high FST, negatif Δπ and ΔHo are expected and would suggest that tuskless individuals carry a more homogeneous haplotype in the associated region. For a selection in tusked individuals, high FST and positif Δπ and ΔHo are wanted meaning that diversity is lesser in tusked than tuskless individuals. 
The interpretation of Δπ is directional. A negative Δπ indicates reduced local diversity in tuskless individuals, whereas a positive Δπ indicates reduced local diversity in tusked individuals.

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/TT_TX_delta_pi_50kb_genomewide.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig3_delta_pi_50kb_GEMMA_Campbell_annotated.png" />

Here I can observed that four regions are negatives for the three last metrics : chr 22, 26 and 32 and among sexual chromosomes. 

### 2.4 Tajima's D 

The Tajima's D was calculated in genomic windows for each phenotype group [SEE 19.TajimaD_TT_TX.sh]. Differences in Tajima's D is used to detect regions where the allele-frequency spectrum differed between tusked and tuskless elephants. I have intersected Tajima’s D differences with FST and pi signals to be able to compare with other metrics [SEE 20.compare_TajimaD_TT_TX.R and 21.intersect_FST_pi_TajimaD.R]. Tajima’s D differences were used as a third population-genomic evidence layer.

The comparison is expressed as : 
**ΔTajimaD = TajimaD_TX - TajimaD_TT**

- **ΔTajimaD < 0** : Tajima's D is lower in TX
- **ΔTajimaD > 0** : Tajima's D is lower in TT

[SEE 20.compare_TajimaD_TT_TX.R and 21.intersect_FST_pi_TajimaD.R]

About the Tajima's D, this metric compare two diversity forms : average diversity between sequences and the number of variants. A negatif Tajima's D can be associated with an excess of rare variants which can be suitable with recent or sweep selection, demographic expansion or purify selection. A positif Tajima's D means that there is an excess of variants with intermediate frequency, and can be associated with balancing selection, population structure, bottleneck event or population mixature.

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/TT_TX_delta_TajimaD_50kb_genomewide.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig4_delta_TajimaD_50kb_GEMMA_Campbell_annotated.png" />

**All the figures combining Bonferroni SNP GEMMA, Campbell candidate genes, and metrics were generated from [SEE 25bis.plot_genomewide_scans_with_GEMMA_Campbell.R]**.

### 2.5. Integrated Campbell candidate genes

For this part, I combined FST, nucleotide diversity and Tajima's D candidate regions with annotated convergent regions. Then I compared selection/differentiation signals with Campbell's candidate genes, to assess whether these regions are closed to Campbell's genes [SEE 22-29 scripts]. The script 25bis was computed during this step.

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Candidate_gene.png"/>

The SNPs and the region on chromosome 3 that I have circled here are closed to Campbell candidate's genes but also for one of them displayed a high FST value. This particular SNP can also be found in other metrics with under 0 values for Heterozygosity, nucleotide diversity and Tajima's D metrics and strenghten thereby its potential association with tusk phenotype, in addition to be relatively close with AMELX. Moreover, this SNP is over the Bonferroni threshold in GEMMA mixed-GWAS (female and male included). Because this SNP is significant in the all-sample GEMMA model but not recovered as the same lead SNP in the male-only analysis, this signal may be influenced by sex composition or sex-linked genetic structure. It should therefore be interpreted cautiously and not as direct evidence for a female-carried causal mutation.

Three final autosomal candidate regions were retained in the original all-individual analysis [SEE 42.make_final_selection_candidate_regions.R]:

| Region | Coordinates | Top SNP | Max FST | Δπ | ΔTajimaD | Ho direction | Evidence score | Genes |
|---|---|---|---:|---:|---:|---|---:|---|
| FST_region_18 | CM044020.1:77800001-77850000 | CM044020.1:77845197:C:T | 0.162061 | -0.00195076 | -2.422898 | TX lower | 4 | LOC126077071; LOC126077054 |
| FST_region_9 | CM044022.1:120700001-120750000 | CM044022.1:120714630:A:T | 0.139687 | -0.00145871 | -2.511417 | TX lower | 4 | no clear protein-coding gene |
| FST_region_56 | CM044021.1:137050001-137100000 | CM044021.1:137060324:A:G | 0.137149 | -0.00138704 | -1.731080 | TX lower | 4 | no clear protein-coding gene |

All three regions showed:
- elevated differentiation;
- reduced π in TX;
- lower Tajima's D in TX;
- lower observed heterozygosity in TX.
These regions are compatible with local differentiation or selection-related processes, but they are not proof of selection. In addition, the sex imbalance in the original dataset requires these signals to be interpreted cautiously.

## 3. Male-only genome-wide scan: 40 TT vs 6 TX

A second genome-wide scan was therefore performed using only males (40 tusked males (TT) & 6 tuskless males (TX)). This removes the major sex-composition confound present in the original analysis, but introduces a strong sample-size imbalance between phenotype groups. The male-only analyses include:
- genome-wide FST;
- observed heterozygosity;
- nucleotide diversity;
- Tajima's D;
- integrated multi-metric candidate windows;
- SNP- and gene-level annotation.
  
### 3.1 Genome-wide FST

A complete genome-wide FST analysis was repeated using only males. The final SNP-level FST table contained **21,295,223 SNPs**.
To remain consistent with the original pipeline, negative FST estimates were excluded before calculating candidate-window summaries. Fixed non-overlapping windows were then generated, and empirical upper-tail thresholds were used to identify differentiated regions.
Using 50-kb windows and the historical p99.9 threshold:
- 86 candidate windows were identified;
- these merged into 75 candidate regions.
All SNPs located in these candidate windows were extracted for annotation.
The resulting FST-only candidate annotation contained 1,741 SNPs. Most were successfully assigned to an intragenic, exonic or nearest-gene context. A small set of variants located on currently unresolved contigs remains to be checked during the final annotation-completeness audit.

### 3.2 Integrated male-only FST, Ho, π and Tajima's D scan

FST, observed heterozygosity, nucleotide diversity and Tajima's D were also integrated in the male-only dataset.
Under the exploratory strict 0.1% multi-metric framework, five FST-anchored regions were retained:

| Region | Coordinates | Strict 0.1% metrics | Broader 1% metrics | Main local annotation |
|---|---|---|---|---|
| PRIMARY_001 | CM044021.1:23200001-23250000 | FST + nucleotide diversity | FST + heterozygosity + nucleotide diversity | CDH18 |
| PRIMARY_002a | CM044021.1:200750001-200800000 | FST + Tajima's D | FST + heterozygosity + nucleotide diversity + Tajima's D | ZNF454 / MGAT2-like region |
| PRIMARY_002b | CM044021.1:200800001-200850000 | FST + nucleotide diversity + Tajima's D | FST + heterozygosity + nucleotide diversity + Tajima's D | ZNF454 / MGAT2-like region |
| PRIMARY_003 | CM044024.1:93500001-93550000 | FST + Tajima's D | FST + heterozygosity + nucleotide diversity + Tajima's D | TRNAK-UUU region |
| PRIMARY_004 | CM044024.1:164700001-164750000 | FST + Tajima's D | FST + nucleotide diversity + Tajima's D | SORCS2-like |
| PRIMARY_005 | CM044039.1:42550001-42600000 | FST + Tajima's D | FST + heterozygosity + Tajima's D | SEC61A1 / RPL9-like |

The five regions contained 2,828 unique SNPs, all of which were extracted and annotated.
The annotation includes SNP-level FST values, genomic context, gene assignment and distance to the assigned gene. Coding SNPs were present in CDH18 and SEC61A1.
These regions are useful as exploratory male-only candidates, but the 40 TT vs 6 TX sample-size imbalance can influence both allele-frequency and diversity estimates. This motivated the population-balanced 6v6 sensitivity analysis described below.

## 4. Population-balanced male-only sensitivity analysis

To evaluate the robustness of the male-only results to the 40 vs 6 imbalance, a repeated balanced design was implemented. All six tuskless males were retained in every replicate. For each replicate, six tusked males were sampled while preserving population representation:
- 3 TT males from India & 3 TT males from Myanmar randomly;
- 3 TX males from India;
- 3 TX males from Myanmar.
This 6 TT vs 6 TX design was repeated 100 times.
For each replicate, genome-wide FST, Ho, π and Tajima's D were recalculated. Candidate windows were then summarized across the 100 resamplings using both:
- the median signal across replicates;
- the recurrence of extreme windows across replicates.
This balanced sensitivity analysis is used to identify candidate regions that are less dependent on the original group-size imbalance.

### 4.1 Resampling design

The male-only dataset was reanalysed using 100 balanced 6 TT vs 6 TX replicates/iterations.
For every replicate:
- all 6 TX males were retained;
- 6 TT males were sampled (3 from India and 3 from Myanmar randomly);
- population composition was balanced between India and Myanmar.
Genome-wide FST, observed heterozygosity, nucleotide diversity and Tajima's D were recalculated independently in each replicate.

### 4.2 Consensus signal and recurrence

Two complementary quantities were used:
1. the median value across the 100 replicates, which describes the typical signal;
2. the proportion of replicates in which a window falls in an extreme tail, which measures recurrence.
These two quantities should not be conflated: a window can have a strong median signal without being extreme in most individual resamples, and vice versa.

The consensus empirical thresholds were:
| Metric | 1% threshold | 0.1% threshold |
|---|---:|---:|
| FST | 0.1732936880 | 0.2522454709 |
| ΔHo, lower tail | -0.1353420219 | -0.2152966895 |
| ΔHo, upper tail | 0.1791663332 | 0.2683106876 |
| Δπ, lower tail | -0.0009958934 | -0.0026413277 |
| Δπ, upper tail | 0.0009917412 | 0.0025367984 |
| ΔTajimaD, lower tail | -1.724626510 | -2.670842656 |
| ΔTajimaD, upper tail | 1.614918310 | 2.3503501607 |

### 4.3 Integrated candidate windows

Candidate windows were defined using multi-metric convergence at the 0.1% level.
The final strict set contained:
- 37 windows supported by at least two extreme metrics;
- 2 windows supported by at least three extreme metrics;
- 0 windows supported by all four metrics;
- 7 FST-anchored windows supported by FST plus at least one additional metric.

A separate recurrence criterion identified 10 windows supported by at least two recurrent extreme metrics in at least 50% of resamples. All 10 were already contained within the 37 strict windows.
One additional tail-contig window, JAMZQU010000059.1:100001-150000, was retained as an explicitly exploratory extreme-Ho signal.
The final annotation universe therefore contained 38 windows.

### 4.4 SNP-level extraction and annotation

All SNPs located within the 38 candidate windows were extracted from the full QC dataset. This produced exactly:
- 32,827 candidate-window SNPs;
- 32,827 annotated SNPs.

The annotation categories were:

| Annotation category | Number of SNPs |
|---|---:|
| Intergenic, nearest gene assigned | 27,751 |
| Non-coding intragenic, non-exonic | 4,207 |
| Intragenic exonic CDS | 697 |
| Intragenic exonic non-CDS | 172 |
| **Total** | **32,827** |

Among these variants, 869 SNPs were exonic.
The phrase candidate-window SNPs is important: these 32,827 SNPs are all variants located inside candidate windows. They are not all individually extreme for FST or another statistic.

### 4.5 Priority regional signals

Priority regional plots were generated from the strongest integrated signals.
Particularly notable regions include:

- CM044026.1:48.60-48.65 Mb: strict three-metric convergence, including FST, Ho and Tajima's D;
- CM044028.1:116.35-116.40 Mb: strict three-metric convergence involving Ho, π and Tajima's D;
- CM044020.1:114.70-114.75 Mb: FST-anchored and highly recurrent;
- several additional FST-anchored or recurrent multi-metric regions;
- the exploratory JAMZQU tail-contig Ho signal.
The current regional figures are considered working figures. Before final publication-level use, adjacent candidate windows should be merged consistently and SNP-level regional FST should be recalculated across the full local interval rather than only inside preselected candidate windows.

| Coordinates | Strict 0.1% metrics | Median FST | Median ΔHo | Median Δπ | Median ΔTajimaD |
|---|---|---:|---:|---:|---:|
| CM044026.1:48600001-48650000 | FST + Ho + TajimaD | 0.3883192 | -0.2468348 | -0.001022123 | -2.811294 |
| CM044020.1:114700001-114750000 | FST + TajimaD | 0.3841758 | -0.1619718 | -0.001229542 | -3.563275 |
| CM044039.1:42550001-42600000 | FST + TajimaD | 0.2784752 | 0.2029672 | 0.000748180 | 3.332765 |
| CM044027.1:84850001-84900000 | FST + TajimaD | 0.2782997 | -0.1110926 | -0.001201060 | -2.923185 |
| CM044026.1:55800001-55850000 | FST + nucleotide diversity | 0.2774050 | -0.2050883 | -0.009441868 | -2.127157 |
| CM044024.1:78450001-78500000 | FST + TajimaD | 0.2808902 | -0.1759868 | -0.000565607 | -3.216030 |
| CM044024.1:164700001-164750000 | FST + TajimaD | 0.2860438 | 0.1702142 | 0.001772594 | 2.377594 |

## 5. Local LD and allele-frequency refinement

This step refined candidate regions by identifying local allele-frequency patterns and LD structure within differentiated genomic intervals. For that I summarized local LD and allele-frequency differences within priority candidate regions [SEE 30-38 scripts]. The script and the methodology was insipired from : https://cloufield.github.io/GWASTutorial/19_ld/.

Main output:
local LD summary tables
50 kb LD window tables
priority-region variant tables
delta allele-frequency ranked variants
annotated top delta-AF variants

## 6. XY selection signal annotation

These scripts extracted and summarized XY-linked or XY-enriched selection signals and annotated top candidate regions.
The sex-linked signals were retained as an additional selection/differentiation evidence layer and were later used during rare-variant and final candidate-locus integration [SEE 39-41].

Regional plots on X chromosome were computed following the explanation of [GWASTutorial](https://cloufield.github.io/GWASTutorial/Visualization/#create-regional-plot).

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/AMELX_sex_linked_GWAS_peak_regional_GWAS_LD_peakLine.png"/>

AMELX gene is located on scaffold CM044047.1 between 168,725,166 and 168,729,418 bp. The best SNP directly located in the gene is CM044047.1:168726231:C:T with a pvalue = 0.10. No Bonferroni SNP was identified in AMELX or within a 100kb window. A suggestif signal exists near to the gene with a pvalue = 9.66E-08, approximately 347bp upstream from AMELX and located in LOC126069583/ARHGAP6-like. However, within 2Mb window around AMELX, only one Bonferroni SNP is returned with a pvalue = 1.56E-10, 1,4Mb (upstream) from AMELX, located in intragenic ncRNA LOC126069593. 

A whole sex-linked/Y-like scaffold GWAS/LD plot was generated for CM044048.1 sex-linked / Y-like scaffold to inspect the major differentiation signal on this scaffold.

Regional plots on Y chromosome : 
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/LOC126069858_GLRA3_like_XY_peak_regional_GWAS_LD_peakLine.png"/>

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/CM044048.1_Y_peak_whole_chr_GWAS_LD.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/CM044048.1_Y_peak_zoom_GWAS_LD.png"/>

lead FST = CM044048.1:14022478:C:T
position = 14.022 Mb

This signal is sex-linked / Y-like scaffold CM044048.1.

top GWAS = CM044048.1:14979125:G:T
position = 14.979 Mb

distance ≈ 0.96 Mb
LD between GWAS and FST peak = r² = 0.050

A chromosome-scale GWAS/LD plot was generated across CM044048.1 to inspect the major sex-linked differentiation signal. The strongest FST window on this scaffold was located at 14.00–14.05 Mb, whereas the strongest GWAS signal was located ~0.96 Mb away near LOC126069858. The weak LD between both lead SNPs suggests that the FST and GWAS peaks are distinct local signals rather than a single shared LD peak.

XY selection signal tables
XY top-region annotation tables


## 7. Final selection/differentiation candidate regions

The annotation framework distinguishes:
- coding CDS variants;
- exonic non-CDS variants;
- non-coding intragenic variants;
- intergenic variants with nearest-gene assignment.
For the current genome-wide analyses:
- the male-only integrated scan contains 2,828 annotated candidate-region SNPs;
- the male-only FST-only p99.9 scan contains 1,741 annotated candidate-region SNPs;
- the balanced 6v6 scan contains 32,827 annotated candidate-window SNPs.
The final annotation audit will explicitly verify, for every scan:
number of SNPs expected → number of unique SNPs extracted → number annotated → number missing → number unresolved
The goal is to ensure that no SNP returned by a genome-wide scan is lost before cross-analysis integration or biological interpretation.

All the previous metrics are then integrated into a final candidate-region table. If one region is supported by multiple statistics, this region is prioritized and should be consider as a strong region.

After integration of all the metrics together, I retained three candidate regions [SEE 42.make_final_selection_candidate_regions.R] : 
The three final selection/differentiation regions are not located on the main X/Y-associated peak.

All three final regions showed:

delta_pi < 0
delta_TajimaD < 0
TX lower heterozygosity
Evidence_score_with_Het = 4

**<ins>Resume table of the 3 identified regions from Genome-wide FST<ins>**

| Region        | Coordinates                    | Top SNP                  |  Max FST |    delta pi | delta TajimaD | Ho direction | Score | Genes                      |
| ------------- | ------------------------------ | ------------------------ | -------: | ----------: | ------------: | ------------ | ----: | -------------------------- |
| FST_region_18 | CM044020.1:77800001-77850000   | CM044020.1:77845197:C:T  | 0.162061 | -0.00195076 |     -2.422898 | TX lower     |     4 | LOC126077071; LOC126077054 |
| FST_region_9  | CM044022.1:120700001-120750000 | CM044022.1:120714630:A:T | 0.139687 | -0.00145871 |     -2.511417 | TX lower     |     4 | none clear                 |
| FST_region_56 | CM044021.1:137050001-137100000 | CM044021.1:137060324:A:G | 0.137149 | -0.00138704 |     -1.731080 | TX lower     |     4 | none clear                 |


These regions are compatible with selection or haplotypic differentiation between tusked and tuskless elephants. However, they do not prove selection. Genetic drift, residual population structure and technical artifacts cannot be fully excluded.

## 8. Conclusion

This part identified three final candidate regions of TT/TX genomic differentiation. These regions were supported by multiple population-genomic signals, including FST, nucleotide diversity, Tajima’s D and observed heterozygosity.

The results should be interpreted as evidence for candidate selection/differentiation regions, not as proof of selection and therefore should be interpretated with caution. These regions were later integrated with GWAS, rare-variant and functional annotation results in the final candidate-locus prioritization.

