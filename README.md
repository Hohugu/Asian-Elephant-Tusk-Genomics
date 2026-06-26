# Selection and population differentiation scans

## 1. Introduction
This part investigates genomic differentiation between tusked and tuskless Asian elephants using population-genomic summary statistics. The aim is to identify genomic regions showing consistent differentiation between phenotype groups, independently of the GWAS association tests. 

The analyses focused on 4 complementary signals : 

- FST, to detect allele-frequency differentiation between tusked and tuskless males.
  
- Nucleotide diversity used to compare local genetic diversity between tusk types.
  
- Tajima's D, to identify shifts in the allele-frequency spectrum
  
- Observed heterozygosity to detect regions with reduced within-group heterozygosity. 

These analyses are used to define candidate selection or reveal differentiation regions. 

## 2. FST scan

Genome-wide FST is calculated in genomic windows to identify regions with elevated allele-frequency between tusked and tuskless elephants [SEE 01.FST_genomewide_TT_vs_TX.sh & 02.FST_postprocess_TT_vs_TX.sh].

## 3. Nucleotide diversity scan

Nucleotide diversity is estimated for each tusk types. The difference Δπ is then used to identify regions where one group showed reduced local diversity [SEE ....]. 

## 4. Tajima's D scan

The Tajima's D is calculated in genomic windows for each phenotype group. Differences in Tajima's D is used to detect regions where the allele-frequency spectrum differed between tusked and tuskless elephants [SEE ...]. 


## 5. Heterozygosity scan 

Observed heterozygosity was estimated across candidate regions to test whether differentiated windows also showed reduced heterozygosity in one tusk type [SEE ....]. 


## 6. Integrated candidate regions

All the previous metrics are then integrated into a final candidate-region table. If one region is supported by multiple statistics, this region is prioritized and should be consider as a strong region. 


After integration of all the metrics together, I retained three candidate regions : 


## 7. Conclusion

This three regions are compatible with selection or haplotypic differentiation between tusked and tuskless elephants. However, these results have to be interpreted cautionly. 

