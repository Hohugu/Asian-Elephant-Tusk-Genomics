# Rare-variant analysis

## 1. Introduction

This part investigates rare variants differing in allele frequency between tusked and tuskless Asian elephants. The goal was to identify rare alleles enriched in one phenotype group and to test whether these variants overlap with GWAS signals, selection/differentiation regions, sex-linked candidate regions or functional genomic annotations.

Rare variants were analyzed as an additional layer of evidence, complementary to GWAS and population-genomic scans. These analyses do not aim to prove causality, but to identify candidate rare variants or rare-variant regions that may contribute to tusk phenotype variation. https://cloufield.github.io/GWASTutorial/34_rare_variant/#variant-selection-and-annotation

**TX means tuskless and TT means tusked.**

## 2. Rare-variant allele-frequency scan

Rare variants were first extracted from the QC genotype dataset and compared between TT and TX elephants [SEE 43.rare_variant_AF_scan_TT_TX.sh and 44.summarize_rare_variants_TT_TX.R].

This scan identified rare variants showing strong allele-frequency differences between tusked and tuskless elephants. In the QC dataset, 7,635 rare differentiated variants were detected.

All rare differentiated variants from this QC scan were classified as TT-enriched rare variants. This means that these alleles were more frequent in tusked individuals than in tuskless individuals.

This result suggests that the rare-variant signal is not randomly distributed between phenotype groups, but it should be interpreted cautiously because rare-variant scans can be sensitive to filtering, missingness, sample size and population structure.

## 3. Raw VCF rare-variant scan

To avoid missing rare variants removed during GWAS QC filtering, I repeated the rare-variant scan using the raw genotyped VCF [SEE 45.rawVCF_rare_variant_AF_scan_TT_TX.sh and 46.summarize_rawVCF_rare_variants_TT_TX.R].

The raw VCF analysis included:

* 40 tusked elephants.
* 53 tuskless elephants.
* 26 423 639 merged biallelic variants.
* 8 008 rare differentiated variants.

Similar to the QC-based scan, all rare differentiated variants detected from the raw VCF were TT-enriched rare variants.

The raw VCF scan therefore confirmed the same global pattern observed in the QC dataset: rare differentiated alleles were predominantly enriched in tusked elephants.

## 4. Intersection with GWAS and selection signals

Rare differentiated variants were intersected with GWAS Bonferroni SNPs, convergent selection/differentiation regions and XY-linked candidate regions [SEE 47.intersect_rare_variants_with_selection_regions.R and 48.intersect_rawVCF_rare_variants_with_selection_regions.R].

Main results from the raw VCF rare-variant scan:

| Category                                       | Count |
| ---------------------------------------------- | ----: |
| Rare differentiated variants                   | 8 008 |
| In convergent selection/differentiation region |     1 |
| In XY score 3 region                           |    87 |
| Directly overlapping a GWAS Bonferroni SNP     |     0 |

No rare differentiated variant was identical to a GEMMA Bonferroni SNP. This indicates that the rare-variant signal and the common-variant GWAS signal are mostly distinct at the single-SNP level.

However, 87 rare differentiated variants were located within XY score 3 regions, suggesting that some rare variants occur in sex-linked or sex-enriched candidate regions.

## 5. Functional annotation of rare variants

Rare differentiated variants were annotated using gene and regulatory features derived from the reference GFF [SEE 57.annotate_rawVCF_rare_variants_regulatory_and_priority_genes.R].

This step classified rare variants according to their genomic context, including coding regions, intronic regions, promoter-proximal regions and intergenic regions.

This annotation was used to prioritize rare variants located near candidate genes, within regulatory-proximal regions or inside protein-coding genes.

However, GFF-based annotation alone cannot identify enhancers, transcription-factor binding sites or tissue-specific regulatory elements. Therefore, non-coding rare variants should be interpreted as candidate regulatory variants only, not as proven regulatory mutations.

## 6. Rare coding variants and predicted effects

Rare variants overlapping coding sequences were further analyzed to predict their potential effect on protein sequences [SEE 58.predict_rawVCF_rare_CDS_variant_effects.py and 59.make_rare_CDS_missense_unique_priority_table.R].

The coding-variant analysis found:

| Category                    | Count |
| --------------------------- | ----: |
| Rare CDS input variants     |    46 |
| VCF alleles found           |    44 |
| CDS transcript-overlap rows |    98 |
| Predicted effect rows       |    89 |
| Missense transcript rows    |    51 |
| Synonymous transcript rows  |    38 |
| Unique missense variants    |    24 |
| Unique synonymous variants  |    20 |

No stop-gained or stop-lost rare variant was detected.

The highest-priority rare missense variants were:

| Variant             | Gene         | Protein annotation                         | Effect | Priority score |
| ------------------- | ------------ | ------------------------------------------ | ------ | -------------: |
| CM044038.1:71818278 | LOC126062954 | carbonyl reductase [NADPH] 2-like          | S43P   |              5 |
| CM044047.1:97173410 | LOC126068753 | E3 ubiquitin-protein ligase DTX1-like      | G66S   |              5 |
| CM044035.1:51994427 | LOC126059344 | polycystic kidney disease 2-like 1 protein | M4V    |              4 |

These variants are interesting functional candidates, but none of the rare missense variants overlapped directly with a GWAS priority gene, convergent selection region or XY selection region.

## 7. Rare variants near GWAS priority genes

Rare differentiated variants were also searched near GWAS-prioritized genes [SEE 60.summarize_rare_variants_near_GWAS_priority_genes.R].

This step identified rare variants located near candidate GWAS loci, including the LOC126069858 / GLRA3-like region. Several rare variants near this region were located in promoter-proximal, intragenic non-exonic or nearby intergenic positions.

These results suggested that LOC126069858 may carry both common-variant GWAS evidence and rare-variant enrichment in the surrounding region.

However, this step was only a candidate-region screening. Formal burden and haplotype-like analyses of LOC126069858 were performed later in a separate section.

## 8. Interpretation

The rare-variant analysis identified thousands of rare variants with differentiated allele frequencies between tusked and tuskless elephants.

The main pattern was that rare differentiated variants were consistently enriched in tusked elephants. This was observed both in the QC genotype dataset and in the raw VCF scan.

Rare variants did not directly overlap GEMMA Bonferroni SNPs, suggesting that rare-variant signals and common-variant GWAS signals mostly represent different layers of genetic variation.

Some rare variants were located within XY-linked candidate regions or near GWAS-prioritized genes, especially around the LOC126069858 / GLRA3-like region. These results motivated further regional rare-variant burden and haplotype-like analyses.

## 9. Caution

Rare-variant enrichment should be interpreted carefully. Rare-variant scans can be influenced by sample size, genotype quality, missingness, read depth, local mapping quality and residual population structure.

The rare variants identified here are candidate variants or candidate regions. They should not be interpreted as causal mutations without additional validation.

In particular, the raw VCF analysis is useful for detecting rare variants that may have been removed by GWAS QC, but raw variants require cautious interpretation and should not replace QC-based association testing.

## 10. Conclusion

This part identified a set of rare differentiated variants between tusked and tuskless elephants. Most rare differentiated variants were enriched in tusked individuals, and no rare differentiated variant directly overlapped a GEMMA Bonferroni SNP.

The rare-variant results provided an additional candidate layer for downstream integration with GWAS, selection/differentiation scans and functional annotation. The strongest rare-variant regional follow-up was later focused on the LOC126069858 / GLRA3-like locus.
