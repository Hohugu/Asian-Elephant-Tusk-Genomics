# Functional and regulatory-proximal annotation

## 1. Introduction

This part annotates GWAS Bonferroni SNPs according to their genomic and functional context. The goal was to determine whether associated SNPs were located in coding regions, exons, promoter-proximal regions, intragenic non-coding regions, or intergenic regions.

This analysis was performed after the GWAS and population-differentiation scans. It provides a functional and regulatory-proximal interpretation of GWAS signals, but it does not prove causality or regulatory mechanism.

The annotation was based on the reference genome FASTA and GFF files. Therefore, this analysis can identify genes, exons, CDS, transcript structure and promoter-proximal regions, but it cannot directly identify enhancers, transcription-factor binding sites, chromatin accessibility or tissue-specific regulatory elements.

## 2. GFF-derived annotation resources

Gene, exon, CDS and promoter-proximal annotation files were generated from the reference GFF [SEE 50.make_promoter_exon_CDS_beds_from_GFF.sh and 50b.convert_regulatory_beds_NC_to_CM.sh].

The resulting annotation resources included:

| Feature class                   |   Count |
| ------------------------------- | ------: |
| Genes                           |  27,210 |
| Promoter-proximal regions, 2 kb |  27,210 |
| Exons                           | 729,282 |
| CDS intervals                   | 637,957 |

Promoter-proximal regions were defined as regions located within 2 kb of a predicted gene transcription start site.

These files were used to classify GWAS Bonferroni SNPs into functional or genomic-context categories.

## 3. Full GFF annotation of GWAS Bonferroni SNPs

All GEMMA Bonferroni SNPs were annotated against the full GFF gene annotation [SEE 54.annotate_ALL_GWAS_Bonferroni_SNPs_fullGFF.R].

The annotation included both:

* GEMMA all-sample PC3 GWAS.
* GEMMA male-only PC3 GWAS.

The all-sample model contained 261 Bonferroni-significant SNPs, whereas the curated male-only model contained 23 Bonferroni-significant SNPs.

The full GFF annotation showed that many GWAS SNPs were located inside or near annotated genes, but a large fraction remained intergenic or non-coding.

| GWAS model     | Inside gene | Outside gene | Not annotated |
| -------------- | ----------: | -----------: | ------------: |
| GEMMA all PC3  |         117 |          137 |             7 |
| GEMMA male PC3 |           8 |           15 |             0 |

This result suggests that the GWAS architecture is not dominated by protein-coding substitutions. Instead, many significant SNPs are located in non-coding, intragenic or intergenic contexts.

## 4. Regulatory-proximal and gene-context classes

GWAS Bonferroni SNPs were classified into functional and regulatory-proximal categories [SEE 55.annotate_GWAS_Bonferroni_regulatory_classes.R].

The categories were:

* CDS.
* Exon non-CDS.
* Promoter-proximal, within 2 kb of a TSS.
* Intragenic non-exonic.
* Intergenic.

### GEMMA all-sample PC3

| Class                 | Count |
| --------------------- | ----: |
| Intergenic            |   142 |
| Intragenic non-exonic |   109 |
| Exon non-CDS          |     6 |
| Promoter-proximal     |     2 |
| CDS                   |     2 |

### GEMMA male-only PC3

| Class                 | Count |
| --------------------- | ----: |
| Intergenic            |    15 |
| Intragenic non-exonic |     7 |
| CDS                   |     1 |

Most Bonferroni SNPs were therefore non-CDS variants. This is compatible with a predominantly non-coding, regulatory-proximal or haplotypic architecture, but this should not be interpreted as proof of regulatory function.

## 5. Coding SNP effect prediction

GWAS Bonferroni SNPs overlapping CDS regions were analyzed to predict their possible coding effects.

Two GWAS CDS SNPs were retained:

| SNP                      | Gene                       | GWAS model | Effect     | Amino-acid change |
| ------------------------ | -------------------------- | ---------- | ---------- | ----------------- |
| CM044025.1:127341274:C:T | DOCK10-like / LOC126078159 | all + male | synonymous | P → P             |
| CM044026.1:117237421:A:G | OR5M3-like / LOC126079052  | all        | missense   | R260C             |

No stop-gained or stop-lost GWAS Bonferroni SNP was detected.

The OR5M3-like SNP is the strongest coding functional candidate because it is a predicted missense variant. The DOCK10-like SNP is present in both all-sample and male-only GWAS models, but its predicted coding effect is synonymous.

## 6. Promoter-proximal GWAS candidates

Two Bonferroni SNPs were located in promoter-proximal regions within 2 kb of a predicted transcription start site.

| SNP                     | Nearby gene                | Class             | Interpretation                               |
| ----------------------- | -------------------------- | ----------------- | -------------------------------------------- |
| CM044020.1:2362860:C:G  | SID1-like / LOC126072020   | promoter-proximal | candidate non-coding regulatory-proximal SNP |
| CM044020.1:14019031:A:G | SEC22A-like / LOC126072714 | promoter-proximal | candidate non-coding regulatory-proximal SNP |

These SNPs are interesting because promoter-proximal variants may influence gene expression. However, without experimental regulatory data, they should be interpreted only as candidate regulatory-proximal variants.

## 7. GWAS-selection overlap annotation

GWAS Bonferroni SNPs were also compared with population-differentiation and selection candidate regions [SEE 57.integrate_GWAS_Bonferroni_vs_selection.R and 59.annotate_GWAS_selection_overlap_SNPs_fullGFF.R].

Main result:

| GWAS model     | Bonferroni SNPs | SNPs overlapping selection/differentiation signals | Percentage |
| -------------- | --------------: | -------------------------------------------------: | ---------: |
| GEMMA all PC3  |             261 |                                                  8 |      3.07% |
| GEMMA male PC3 |              23 |                                                  0 |         0% |

The overlap was limited and occurred only in the all-sample GWAS model.

Notable annotated overlap SNPs included:

| SNP                        | Gene or region            | Annotation                               |
| -------------------------- | ------------------------- | ---------------------------------------- |
| CM044048.1:14979125:G:T    | LOC126069858 / GLRA3-like | intragenic non-exonic                    |
| CM044047.1:154448329:G:T   | LOC126068977 / PDK3       | intragenic                               |
| CM044033.1:6293379/6293387 | LOC126058128              | lncRNA region, approximately 311 kb away |
| JAMZQU010000039.1 variants | no clear gene             | no gene found on contig                  |

The LOC126069858 / GLRA3-like SNP is particularly important because it combines GWAS evidence, sex-linked differentiation evidence and later rare-variant support. However, this does not prove that this SNP is causal.

## 8. Functional priority candidates

The annotation results were integrated into a GWAS functional priority table [SEE 60.make_GWAS_Bonferroni_priority_table.R].

The highest-priority candidates included:

| SNP                      | Gene                       | Class                 | Main reason for priority                              |
| ------------------------ | -------------------------- | --------------------- | ----------------------------------------------------- |
| CM044026.1:117237421:A:G | OR5M3-like / LOC126079052  | CDS                   | missense variant                                      |
| CM044025.1:127341274:C:T | DOCK10-like / LOC126078159 | CDS                   | found in all-sample and male-only GWAS                |
| CM044020.1:2362860:C:G   | SID1-like / LOC126072020   | promoter-proximal     | near TSS                                              |
| CM044020.1:14019031:A:G  | SEC22A-like / LOC126072714 | promoter-proximal     | near TSS                                              |
| CM044047.1:154448329:G:T | PDK3 / LOC126068977        | intragenic            | GWAS-selection overlap                                |
| CM044048.1:14979125:G:T  | LOC126069858 / GLRA3-like  | intragenic non-exonic | GWAS-selection overlap and later rare-variant support |

These loci were retained for downstream candidate-locus integration.

## 9. Interpretation

The functional annotation showed that most GWAS Bonferroni SNPs were not located in CDS regions. Only a small number of GWAS SNPs were coding, and no stop-gained or stop-lost coding SNP was detected.

Several important GWAS SNPs were located in intragenic non-exonic, promoter-proximal or intergenic regions. This pattern is compatible with a non-coding, regulatory-proximal or haplotypic genetic architecture for tusk phenotype variation.

The strongest coding candidate was the OR5M3-like missense SNP. The strongest male-supported coding candidate was the DOCK10-like synonymous SNP. The LOC126069858 / GLRA3-like locus was not a coding-effect candidate, but it remained highly important because of its GWAS signal, selection overlap and later rare-variant evidence.

## 10. Caution

This analysis is based on GFF-derived annotation. Therefore, promoter-proximal and intragenic non-coding SNPs should not be interpreted as proven regulatory mutations.

The GFF annotation cannot identify:

* Enhancers.
* Transcription-factor binding sites.
* ATAC-seq peaks.
* ChIP-seq peaks.
* Tissue-specific regulatory activity.

Therefore, the correct interpretation is that many GWAS SNPs are candidate non-coding or regulatory-proximal variants. Functional validation would be required to demonstrate a true regulatory mechanism.

## 11. Conclusion

This part showed that the GWAS Bonferroni SNPs are mostly non-CDS variants, with only a small number of coding candidates. The results support a candidate architecture involving non-coding, regulatory-proximal and haplotypic signals rather than a simple coding loss-of-function mechanism.

The functional annotation prioritized several candidate loci, including OR5M3-like, DOCK10-like, SID1-like, SEC22A-like, PDK3 and LOC126069858 / GLRA3-like.
