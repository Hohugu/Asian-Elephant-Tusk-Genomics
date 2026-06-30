# # Rare-variant discovery and annotation

## 1. Introduction

This part investigates rare variants differing in allele frequency between tusked and tuskless Asian elephants. The goal was to identify rare alleles enriched in one phenotype group and to test whether these variants overlap with GWAS signals, selection/differentiation regions, sex-linked candidate regions or functional genomic annotations.

Rare variants were analyzed as an additional layer of evidence, complementary to GWAS and population-genomic scans. These analyses do not aim to prove causality, but to identify candidate rare variants or rare-variant regions that may contribute to tusk phenotype variation. This section does not perform genome-wide SKAT, SKAT-O or formal rare-variant burden testing. Although whole-genome sequencing is suitable for rare-variant discovery, the limited sample size makes formal genome-wide rare-variant association testing underpowered. Rare variants were therefore analyzed through allele-frequency enrichment, gene-context annotation and candidate-region prioritization. https://cloufield.github.io/GWASTutorial/34_rare_variant/#variant-selection-and-annotation

**TX and TT means tuskless and tusked.**

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

## 5. Gene-context annotation of rare variants

Rare differentiated variants were annotated according to their genomic context using GFF-derived genomic features [SEE 50.make_promoter_exon_CDS_beds_from_GFF.sh, 50b.convert_regulatory_beds_NC_to_CM.sh and 49.annotate_rawVCF_rare_variants_regulatory_and_priority_genes.R].

The annotation classified rare variants into several categories:

- CDS variants.
- Exon non-CDS variants.
- Intragenic non-exonic variants.
- Promoter-proximal variants, defined as variants located within 2 kb of a gene transcription start site.
- Intergenic variants.

This step was used to identify rare variants that may have functional or regulatory relevance, especially variants located near GWAS-prioritized genes or within sex-linked candidate regions.

A subset of rare variants near the LOC126069858 / GLRA3-like region was particularly interesting. Among the high-priority rare variants in this region, several were located in promoter-proximal or intragenic non-exonic contexts. For example, rare variants close to the LOC126069858 transcription start site included promoter-proximal variants located approximately 832 bp and 1,524 bp from the TSS, as well as intragenic non-exonic variants located less than 2 kb from the TSS.

These rare non-coding variants may indicate a regulatory or haplotypic signal around LOC126069858. However, this interpretation remains cautious. GFF-based annotation can identify genes, exons, CDS and promoter-proximal regions, but it cannot directly identify enhancers, transcription-factor binding sites, chromatin accessibility or tissue-specific regulatory elements.

Therefore, these variants should be interpreted as candidate regulatory or non-coding variants, not as proven regulatory mutations.

## 6. Rare coding variants and predicted effects

Rare variants overlapping coding sequences were further analyzed to predict their potential effect on protein sequences [SEE 51.predict_rawVCF_rare_CDS_variant_effects.py and 52.make_rare_CDS_missense_unique_priority_table.R].

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

Rare differentiated variants were also searched near GWAS-prioritized genes [SEE 53.summarize_rare_variants_near_GWAS_priority_genes.R].

Several rare variants near GWAS-prioritized genes were non-coding but located in regulatory-relevant contexts, including promoter-proximal or intragenic non-exonic regions. This was particularly clear around the LOC126069858 / GLRA3-like locus, where rare TT-enriched variants were found near the gene body and close to the predicted transcription start site.

This supports the hypothesis that rare-variant contribution at this locus may be non-coding, regulatory or haplotypic rather than protein-coding. However, this remains a candidate interpretation because no experimental regulatory data such as ATAC-seq, ChIP-seq or enhancer annotation were available.

| Variant | Context | Distance to TSS | AF_TT | AF_TX | Comment |
|---|---|---:|---:|---:|---|
| CM044048.1:15116231 | promoter-proximal | 832 bp | 0.1000 | 0.0094 | near LOC126069858 TSS |
| CM044048.1:15116923 | promoter-proximal | 1,524 bp | NA | NA | indel near LOC126069858 TSS |
| CM044048.1:15113612 | intragenic non-exonic | 1,787 bp | 0.1125 | 0 | TT-enriched, 9 TT carriers, 0 TX |

These variants were retained as candidate non-coding rare variants near LOC126069858. Their regional burden and haplotype-like structure are analyzed in a later section.

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
The rare-variant results provided an additional candidate layer for downstream integration with GWAS, selection/differentiation scans and functional annotation. In particular, the presence of rare promoter-proximal and intragenic non-coding variants near LOC126069858 motivated further regional burden and haplotype-like analyses.
