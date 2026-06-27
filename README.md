# Selection and population differentiation scans

## 1. Introduction
This part investigates genomic differentiation between tusked and tuskless Asian elephants using population-genomic summary statistics. The aim is to identify genomic regions showing consistent differentiation between phenotype groups, independently of the GWAS association tests. 

The analyses focused on differents complementary signals : 

- Genome-wide FST scans.
  
- Nucleotide diversity comparisons.
  
- Tajima's D, to identify shifts in the allele-frequency spectrum
  
- Observed heterozygosity to detect regions with reduced within-group heterozygosity.

- Local LD

- Candidate-region annotation

These analyses are used to define candidate selection or reveal differentiation regions. 

## 2. FST scan and FST candidate regions

Genome-wide FST is calculated in genomic windows to identify regions with elevated allele-frequency between tusked and tuskless elephants [SEE 01.FST_genomewide_TT_vs_TX.sh - 03.FST_make_windows.R]. FST scans identified genomic windows showing elevated differentiation between tusked and tuskless individuals [SEE 04.plot_FST_genomewide.R]. Then I summarized FST values in genomic windows, identified candidate FST-enriched regions, and annotated these regions with nearby genes [SEE 05.FST_merge_candidate_regions.R - 07.make_FST_candidate_regions_annotated.R].
These regions were used as the first layer of evidence for selection/differentiation candidate regions.

FST candidate regions are then compared with GWAS signals, to test whether differentiated regions are located near GWAS-associated SNPs or candidate GWAS loci [SEE 08.compare_FST_vs_GWAS.R and 09.compare_FST_vs_GWAS_distance..R]. This step helped to distinct regions supported only by differentiation scans from regions also close to GWAS signals. These comparisons were later used during candidate-locus prioritization, but they are not equivalent to GWAS evidence.

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Genomewide_FST_1Mb_annotated.png" />
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig1_FST_50kb_GEMMA_Campbell_annotated.png" />

SNPs with the higher FST values are located in chromosome 1, 3, 11 and among sexual chromosomes. No GEMMA SNPs returned with high/significant FST value. 
A high FST value usually means that allelic frequencies between tusked and tuskless individuals are differents. Thereby, if a GWAS GEMMA SNPs have a high FST value, then this GWAS SNPs is located inside a region where tusked and tuskless individuals are differents. Morevoer, a Bonferroni SNPs near to Campbell candidate gene and inside a high FST window, would suggest spatial concordance between GWAS association, TT/TX differentiation and prior tooth/tusk candidate genes.

## 3. Heterozygosity scan 

Observed heterozygosity is estimated across candidate regions to test whether differentiated windows also showed reduced heterozygosity in one tusk type [SEE 10
heterozygosity_TT_TX.sh - 13.windowed_observed_heterozygosity.R]. 
Observed heterozygosity is used as an additional evidence layer. The final candidate selection/differentiation regions showed reduced heterozygosity in tuskless elephants [SEE 14.add_Het_to_integrated_candidate_table.R and 15.plot_genomewide_delta_Ho.R].

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/TT_TX_delta_Ho_50kb_genomewide.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig2_delta_Ho_50kb_GEMMA_Campbell_annotated.png" />

Similar to mean FST graphs, SNPs with significant heterozygosity values are located in chromosome 3, 11 and among sexual chromosomes. 

## 4. Nucleotide diversity scan

Nucleotide diversity is estimated for each tusk types and compared [SEE 16.compare_pi_TT_TX.R]. The difference Δπ is then used to identify regions with one group showing a reduced local diversity [SEE 17.compare_pi_FST_GWAS.R and 18.annotate_reduced_pi_FST_overlap.R]. 

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/TT_TX_delta_pi_50kb_genomewide.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig3_delta_pi_50kb_GEMMA_Campbell_annotated.png" />


## 5. Tajima's D scan

The Tajima's D is calculated in genomic windows for each phenotype group [SEE 19.TajimaD_TT_TX.sh]. Differences in Tajima's D is used to detect regions where the allele-frequency spectrum differed between tusked and tuskless elephants. I have intersected Tajima’s D differences with FST and pi signals to be able to compare with other metrics [SEE 20.compare_TajimaD_TT_TX.R and 21.intersect_FST_pi_TajimaD.R]. Tajima’s D differences were used as a third population-genomic evidence layer.

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/TT_TX_delta_TajimaD_50kb_genomewide.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/Fig4_delta_TajimaD_50kb_GEMMA_Campbell_annotated.png" />


## 6. Convergent candidate regions and Campbell candidate genes

For this part, I combined FST, nucleotide diversity and Tajima's D candidate regions with annotated convergent regions. Then I compared selection/differentiation signals with Campbell's candidate genes, to assess whether these regions are closed to Campbell's genes. [22-29]

## 7. Local LD and allele-frequency refinement

This step refined candidate regions by identifying local allele-frequency patterns and LD structure within differentiated genomic intervals. For that I summarized local LD and allele-frequency differences within priority candidate regions. [30-38] The script and the methodology was insipired from : https://cloufield.github.io/GWASTutorial/19_ld/.

Main output:
local LD summary tables
50 kb LD window tables
priority-region variant tables
delta allele-frequency ranked variants
annotated top delta-AF variants

## 8. XY selection signal annotation

These scripts extracted and summarized XY-linked or XY-enriched selection signals and annotated top candidate regions.
The sex-linked signals were retained as an additional selection/differentiation evidence layer and were later used during rare-variant and final candidate-locus integration [SEE 39-41].

The regional plots were computed following the explanation of [GWASTutorial](https://cloufield.github.io/GWASTutorial/Visualization/#create-regional-plot).

<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/AMELX_sex_linked_GWAS_peak_regional_GWAS_LD_peakLine.png"/>
<img width="900" height="700" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/04.Selection-Differentiation-scans/LOC126069858_GLRA3_like_XY_peak_regional_GWAS_LD_peakLine.png"/>

XY selection signal tables
XY top-region annotation tables


## 9. Final selection//differentiation candidate regions

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

## 10. Conclusion

This part identified three final candidate regions of TT/TX genomic differentiation. These regions were supported by multiple population-genomic signals, including FST, nucleotide diversity, Tajima’s D and observed heterozygosity.

The results should be interpreted as evidence for candidate selection/differentiation regions, not as proof of selection and therefore should be interpretated with caution. These regions were later integrated with GWAS, rare-variant and functional annotation results in the final candidate-locus prioritization.

