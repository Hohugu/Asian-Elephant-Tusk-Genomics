# WGS Variant Calling : Variant Calling with GATK | a Detailed NGS Analysis

## 1. Introduction

The aim of this protocle is to perform a series of steps to determine a set of genetic variants compiled in the ***VCF*** file. 

For that, I used GATK : **Genome Analysis ToolKit**

GATK is a standart tookit to analyse and identify variants in genome. GATK tools can be used individually or chained together into complete workflows. I followed the schematic picture in **README**  in the main branch. 

The reference for this workflow is called : _GCA.024166365.1_mEleMAX1_primary_haplotypes_genomic.fa_ 
This reference have been used for India individuals and Myanmar individuals. Whole genome sequencing have been done for 28 and 65 individuals respectively, for a total of 93 individuals sequenced.

For Myanmar individuals, all the ***.bam*** files and ***.vcf*** files were generated for all the individuals. 
For India, the data comes from _Anubhab et al., 2024 :"Divergence and serial colonization shape genetic variation and define conservation units in Asian elephants"_.
  > NCBI Bioproject link : _https://ncbi.nlm.nih.gov/bioproject/PRJNA1013751_

  > All the ***.fastq*** files can be reach at : _https://www.ncbi.nlm.nih.gov/sra_

There the aim of following Variant Calling with GATK protocol, is to generate ***.fastq*** files for India individuals and merge the final ***.bam*** and ***.vcf*** files to Myanmar files and have Analysis-ready vcf file for Pre-GWAS step. 

For that, I followed those links: _www.youtube.com/watch?v=iHkiQvxyr5c_ & _https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-indels_

## 2. Detailed protocol

### 2.1 Data pre-processing & Variant discorery

This step only concerns indian individuals.
For this first step, I started by downloading the ***.sra*** files stored in [NCBI  Bioproject link SRA](https://www.ncbi.nlm.nih.gov/sra), into Puhti server. Then I converted ***.sra*** files into ***.fastq*** files. All the individuals were paired-end sequenced, which means that individuals were sequenced in 5'3' and 3'5' directions. So the conversion returned two ***.fastq*** files for each individual [SEE 02.SCRIPT]. Then each ***.fastq*** files are controlled for quality with **fastqc** parameter. If any adapters were detected in the html returned files then ***.fastq*** files have to be trimmed. Here no adapters were found. 

Once I get the aligned read ***.bam***, I flagged the duplicate reads. During sequencing process the same DNA fragments may be sequenced multiple times. Duplicate reads can arise during sample preparation step that is library construction during PCR. Duplicate reads are not informative and can be evidence for or against a variant, so they can be elimate. Once flagged, duplicate reads will be ignored for the rest of the downstream category tools. To flag and eliminate duplicate reads, I used **MarkDuplicates Spark**.

The last step of data pre-processing, is to recalibrate base quality scores. For that an algorithms is called to rely heavily on the quality scores that are assigned to individual bases in each sequencing read. The quality scores can tell us how much one particular base can be trust at one location. So if a base call with a low quality score, it means that this base is not sure. This quality score serve as evidence to decide on removing this base at this location.  
For this step, I applied **Boostrapping**: In the case that no known set of variants is available for the system studied (i.e Asian elephants), a first step can be to generate a raw ***.vcf*** file filtered, but without using base quality score recalibration and filtering the variance to obtain high confidence set of variants and then using those variants as an input for the base quality score recalibration. 

For the **Variant discovery** step, the parameter _HaplotypeCaller_ is used from GATK, and returns a ***.vcf*** file. According the protocol, one file is generated with SNPs and another with indels, with their respective index file. _HaplotypeCaller_ is used to call variance from the reads. The choice of variant calling algorithm has to be dependent on certain critera. _HaplotypeCaller_ can handle multiple samples but it is not recommended to use it when you are trying to analyse more than 100 samples at a time. 

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
