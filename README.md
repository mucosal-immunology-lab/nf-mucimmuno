# Nextflow Workflows

As they are built and published, this repository will contain Nextflow workflows for processing of different data omic modalities.

Additionally, we may provide additional tools and code for further downstream processing, with the goal of standardising data analytic approaches within the Mucosal Immunology Lab.

- [Nextflow Workflows](#nextflow-workflows)
  - [Single-cell RNAseq FASTQ pre-processing](#single-cell-rnaseq-fastq-pre-processing)


## Single-cell RNAseq FASTQ pre-processing

[**nf-mucimmuno/scRNAseq**](./scRNAseq/) is a bioinformatics pipeline that can be used to run quality control steps and alignment to a host genome using STARsolo. It takes a samplesheet and FASTQ files as input, performs FastQC, trimming and alignment, and produces an output `.tar.gz` archive containing the collected outputs from STARsolo, ready for further processing downstream in R. MultiQC is run on the FastQC outputs both before and after TrimGalore! for visual inspection of sample quality &ndash; output `.html` files are collected in the results.

<div align="center">
<img src="./assets/images/nf-mucimmuno_scRNAseq-01.png" width=80%>
</div>