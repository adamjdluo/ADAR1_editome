# ADAR1_editome
* This repository recorded original scripts for the data analysis of the article titled **Decoupling expression and editing preferences of ADAR1 p150 and p110 isoforms** (Sun et al, PNAS(2021), 118(12):e2021757118) [link](https://www.pnas.org/doi/abs/10.1073/pnas.2021757118). 

* Essential Resources
  + *STAR* aligner (ver 2.5.4b)
  + *Picard* (ver 2.18.1)
  + Genomic Analysis Toolkit (*GATK*, ver 4.0.8)
  + *annovar* (ver 2019-10-24 00:05:28 -0400) [link](https://annovar.openbioinformatics.org/en/latest/)
  + *bcftools* (ver 1.9)
  + *bam-readcount* (0.8.0-unstable-6-963acab-dirty) [link](https://github.com/genome/bam-readcount)
  + Reference genome sequences (fa): *BSGenome.Hsapiens.UCSC.hg19.fa*
  + Annotation table (gtf): *TxDb.Hsapiens.UCSC.hg19.knownGene.gtf*
  + annotation database constructed by annovar [link to the instruction](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)
  + Locations of Alu elements in human genome (hg19): *hg19_UCSC_Alu_env.RData (in this repository)*

* Workflow

This script is revised from RNAseq short variant discovery ([link](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-)) of GATK Best Practices Workflows, then followed by identifying and categorizing variations in a customized way. The steps were recorded in the following Rmd files

  + *align_and_variant_call_GATK_BestPrec_Rseq.Rmd*: Alignment, variant calling and annotation by following GATK Best Practice. Variant annotation was donw by annovar.
  + *cal_DP_AC_AF.Rmd*: calculate the total depth (DP), depth of minor allele (AC), and alleleic frequency (AF) of minor allele. 
  + *infoInt_catVar.Rmd*: integrate information for variants in each sample, identify candidates, and categorize variants
  + *report_rmd_v2.Rmd*: make report and visualization

Bioinformatics Resource Center, The Rockefeller University