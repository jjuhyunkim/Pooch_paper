# Manuscript Analyses

This directory contains the scripts used to process data, perform statistical
analyses, train and evaluate models, and generate genome browser tracks
for the Pooch manuscript.

Plotting scripts and figure notebooks are organized separately in the
[figure index](../figures/README.md). Shared reference resources,
annotations, and chain files are documented in
[data/README.md](../data/README.md).

## Analysis Overview

| Directory | Description |
|---|---|
| [00_methylation_processing/](00_methylation_processing/) | Read alignment, BAM sorting/merging, haplotype-specific extraction, modkit pileup, and haplotype-to-CHM13 chain generation. The checked-in scripts do not include coordinate lifting or cohort-level CpG filtering. |
| [01_idmr_identification/](01_idmr_identification/) | Paired maternal/paternal differential methylation testing with DSS and DMR calling. |
| [02_model_training_validation/](02_model_training_validation/) | Grouped nested cross-validation, elastic-net model fitting, and comparison of parental-origin predictions with PatMat, including XIST-adjusted modeling and BED-based CpG exclusion/validation utilities. |
| [03_prediction_evaluation/](03_prediction_evaluation/) | XIST-adjusted chrX parental-origin prediction for WashU pedigree and Robertsonian translocation samples. |
| [04_chr3_chr17_chrX_analysis/](04_chr3_chr17_chrX_analysis/) | CpG-level XIST interaction analysis and BED methylation summarization into genomic bins. |
| [05_computational_benchmark/](05_computational_benchmark/) | Prediction CPU-time/peak-memory benchmarking and chr22 model training-size analysis using 10–80 training samples. |


## Analysis Documentation

Each analysis README summarizes:

- The purpose of the analysis.
- The methods and parameters visible in the checked-in scripts.
- Each script's inputs, outputs, and execution requirements where available.

The step-specific READMEs describe each workflow's scripts, inputs, outputs,
and execution requirements where available. Some workflows depend on data,
sample lists, reference files, or cluster modules that are not distributed
here; adapt paths and execution details to the local environment.

## Pooch Software

The Pooch software is maintained separately at
https://github.com/marbl/Pooch.