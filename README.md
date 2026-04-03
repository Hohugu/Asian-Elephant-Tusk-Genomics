# WGS Variant Calling : Variant Calling with GATK | a Detailed NGS Analysis

## 1. Introduction

The aim of this protocle is to perform a series of steps to determine a set of genetic variants compiled in the ***VCF*** file. 

For that, I used GATK : **Genome Analysis ToolKit**

GATK is a standart tookit to analyse and identify variants in genome. GATK tools can be used individually or chained together into complete workflows. I followed the schematic picture in **README**  in the main branch. 

The reference for this workflow is called : _GCA.024166365.1_mEleMAX1_primary_haplotypes_genomic.fa_ 
This reference have been used for India individuals and Myanmar individuals. Whole genome sequencing have been done for 28 and 65 individuals respectively, for a total of 93 individuals sequenced.

For Myanmar individuals, the ***.fastq*** files are stored in online platform ALLAS from CSC and are not public. 
For India, the data comes from [_A. Khan et al., 2024_](https://doi.org/10.1016/j.cub.2024.08.062) _:"Divergence and serial colonization shape genetic variation and define conservation units in Asian elephants"_.
  > NCBI Bioproject link : _https://ncbi.nlm.nih.gov/bioproject/PRJNA1013751_

  > All the ***.fastq*** files can be reach at : [_https://www.ncbi.nlm.nih.gov/sra_](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR25983372)

There the aim of following Variant Calling with GATK protocol, is to generate ***.fastq*** files for India individuals and convert the ***.fastq*** files of Myanmar and Indian elephants to ***.bam*** files. Then the ***.bam*** files will be also convert after several steps into a final ***.vcf*** file to have Analysis-ready vcf file for Pre-GWAS step. 

For that, I followed those links: _www.youtube.com/watch?v=iHkiQvxyr5c_ & _https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-indels_

## 2. Detailed protocol

### 2.1 Data pre-processing & Variant discorery

For this first step, I started by downloading the ***.sra*** files of the Indian individuals stored in [NCBI  Bioproject link SRA](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR25983372), into Puhti server. Then I converted ***.sra*** files into ***.fastq*** files with **fastq-dump** fonction. All the individuals were paired-end sequenced, which means that individuals were sequenced in 5'3' and 3'5' directions. So the conversion returns 2 ***.fastq*** files for each individual [SEE 02.SRA_to_FASTq_slurm]. Then each ***.fastq*** files are controlled for quality with **fastqc** parameter [SEE 03.FASTQC MYA and SRR]. If any adapters were detected in the html returned files then ***.fastq*** files have to be trimmed. Here no adapters were found in the html generated files and very good quality has been returned. 

Once I get the aligned read ***.bam*** with **bwa mem**, I flagged the duplicate reads with ***MarkDuplicates Spark***. During sequencing process the same DNA fragments may be sequenced multiple times. Duplicate reads can arise during sample preparation step that is library construction during PCR. Duplicate reads are not informative and can be evidence for or against a variant, so they can be elimate. Once flagged, duplicate reads will be ignored for the rest of the downstream category tools. [SEE 04.Alignment_duplicate_BAM]

The last step of data pre-processing, is to recalibrate base quality scores. For that an algorithms is called to rely heavily on the quality scores that are assigned to individual bases in each sequencing read. The quality scores can tell us how much one particular base can be trust at one location. So if a base call with a low quality score, it means that this base is not sure. This quality score serve as evidence to decide on removing this base at this location.  
For this step, **Boostrapping** can be used: In the case that no known set of variants is available for the system studied (i.e Asian elephants), a first step can be to generate a raw ***.vcf*** file with standard filtration, but without using base quality score recalibration and filtering the variance to obtain high confidence set of variants and then using those variants as an input for the base quality score recalibration. This step is recommended but optional. As in this study, we lack of known variant and in a matter of time and HPC's space, I had to skip this step, but once the results done, it totally possible to re-run the process from ***.bam*** files and configurate **Boostrapping**, to the end compare with first results.

To optimize the space on HPC, I converted the ***.bam*** files into ***.cram*** files with a default compression level fixed at 5 [SEE 05.BAM_CRAM].

For the **Variant discovery** step, the parameter ***HaplotypeCaller*** is used from GATK, and returns a ***g.vcf*** file for each individuals [SEE 06.GVCF]. A step of verification of files is applied to ensure the good use of ***g.vcf*** files for further steps. ***HaplotypeCaller*** is used to call variance from the ***.cram*** reads. The choice of variant calling algorithm has to be dependent on certain critera. ***HaplotypeCaller*** can handle multiple samples but it is not recommended to use it when you are trying to analyse more than 100 samples at a time. In our study, I took the choice to apply ***HaplotypeCaller*** only for 1 individual at a time, to avoid bad scalability, batch effects and use large memory and CPU.
According to the [Broad Institute’s GATK documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller?utm_source=chatgpt.com), the recommended germline variant calling workflow is to run ***HaplotypeCaller*** in GVCF mode per individual to produce single-sample gVCFs, then consolidate those gVCFs using ***GenomicsDBImport*** (or CombineGVCFs) [SEE 07.GenomicsDBImport], and finally perform joint genotyping using ***GenotypeGVCFs*** to produce a multi-sample VCF [SEE 08.GenotypeGVCF].

This part of the protocole is resume in the following picture : 

  <img width="412" height="834" alt="image" src="https://github.com/Hohugu/Genomic-on-Asian-elephant-Tusk/blob/main/GenomicsDBIGenotypeVCFs.png" />

### 2.2 Filtering and Annotation

Once the ***.vcf*** file is generated, a new filtering step can be done. For that, they are two filtering way :
> - **Variant Quality Score Recalibration (VQSR)**
      VQSR is a sophisticated filtering technique applied on the variant callset.
      It uses machine learning to model the technical profile of variants in a training set.
      It uses that to filter out probable artificats from the callset.

> - **Hard filtering**
      Hard-filtering threshold for site-level manual filtering.
      A typical manual filtration is small cohort callsets (<30 samples)

The two techniques commonly uses and the aim is to look at the statistics associated with each variants and coming up with certain set ratios and then keep the variants that fall within the thresholds and filter out variants that fall either above or below the set threshold. One problem is to throwing good quality variance because one of their annotations look bad or keep bad variance.

VQSR needs good curated and quality training data from which the algorithm learns the annotation profile of goof quality vs. bad. But also the good curated training ressources are not always available, so it recommands large amount of variations sites to operate properly. If not possible, they are small scale experiments or such as targeted gene panels. 

The filter site-level and sample-levels can be found more in detail in [GATK training](gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels)

Following, a resume of the filter (source : https://www.youtube.com/watch?v=XZ8scaScfjw ) 

<img width="954" height="453" alt="image" src="https://github.com/user-attachments/assets/635ddfe0-7edc-431b-98b8-750097a8aeef" />

I didn't use this step of filtering because in the pre-GWAS step, initialization and filtering are performed to prepare the data for the GWAS [SEE Pre-GWAS branch].
