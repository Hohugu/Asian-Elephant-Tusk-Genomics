# Selection and population differentiation scans

## 1. Introduction
This part investigates genomic differentiation between tusked and tuskless Asian elephants using population-genomic summary statistics. The aim is to identify genomic regions showing consistent differentiation between phenotype groups, independently of the GWAS association tests. 

The analyses focused on 4 complementary signals : 

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

## 3. Heterozygosity scan 

Observed heterozygosity is estimated across candidate regions to test whether differentiated windows also showed reduced heterozygosity in one tusk type [SEE 10
heterozygosity_TT_TX.sh - 13.windowed_observed_heterozygosity.R]. 
Observed heterozygosity is used as an additional evidence layer. The final candidate selection/differentiation regions showed reduced heterozygosity in tuskless elephants [SEE 14.add_Het_to_integrated_candidate_table.R and 15.plot_genomewide_delta_Ho.R].

## 4. Nucleotide diversity scan

Nucleotide diversity is estimated for each tusk types and compared [SEE 16.compare_pi_TT_TX.R]. The difference Δπ is then used to identify regions with one group showing a reduced local diversity [SEE 17.compare_pi_FST_GWAS.R and 18.annotate_reduced_pi_FST_overlap.R]. 

## 5. Tajima's D scan

The Tajima's D is calculated in genomic windows for each phenotype group [SEE 19.TajimaD_TT_TX.sh]. Differences in Tajima's D is used to detect regions where the allele-frequency spectrum differed between tusked and tuskless elephants. I have intersected Tajima’s D differences with FST and pi signals to be able to compare with other metrics [SEE 20.compare_TajimaD_TT_TX.R and 21.intersect_FST_pi_TajimaD.R]. Tajima’s D differences were used as a third population-genomic evidence layer.

## 6. Convergent candidate regions and Campbell candidate genes


## 7. Local LD and allele-frequency refinement

## 8. XY selection signal annotation


## 9. Final selection//differentiation candidate regions

All the previous metrics are then integrated into a final candidate-region table. If one region is supported by multiple statistics, this region is prioritized and should be consider as a strong region.

After integration of all the metrics together, I retained three candidate regions : 




All three final regions showed:

delta_pi < 0
delta_TajimaD < 0
TX lower heterozygosity
Evidence_score_with_Het = 4

These regions are compatible with selection or haplotypic differentiation between tusked and tuskless elephants. However, they do not prove selection. Genetic drift, residual population structure and technical artifacts cannot be fully excluded.

## 10. Conclusion

This part identified three final candidate regions of TT/TX genomic differentiation. These regions were supported by multiple population-genomic signals, including FST, nucleotide diversity, Tajima’s D and observed heterozygosity.

The results should be interpreted as evidence for candidate selection/differentiation regions, not as proof of selection and therefore should be interpretated with caution. These regions were later integrated with GWAS, rare-variant and functional annotation results in the final candidate-locus prioritization.

