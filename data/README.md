# Data and Reference Resources

This directory contains reference resources, annotations, and supporting
files used in the manuscript analyses.

## Reference Genomes

| File | Assembly / version | Source |
|---|---|---|
| `chm13v2.0.fa.gz` | T2T-CHM13v2.0 | [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz) |


## Genomic Annotations

| File | Description | Genome assembly | Source |
|---|---|---|---|
| `Imprinted_DMR_List_V1.LiftedOverToT2T-CHM13v2Hs1.tsv` | PatMat iDMRs | T2T-CHM13v2.0 | [PatMat github](https://github.com/vahidAK/PatMat/raw/refs/heads/main/patmat/Imprinted_DMR_List_V1.LiftedOverToT2T-CHM13v2Hs1.tsv) |
| `putative_ICRs_35_65_final_N1488-chrX.gff3` | ICR candidates | GRCh38 | [humanicr](https://humanicr.org/) |
| `Sperm_with_sue_100-chrX.gff3` | Methylation fraction from Sperm | GRCh38 | [humanicr](https://humanicr.org/) |
| `Oocyte_100-chrX.gff3` | Methylation fraction from Oocyte | GRCh38 | [humanicr](https://humanicr.org/) |
| `chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz` | Gene annotations | T2T-CHM13v2.0 | [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz) |
| `chm13v2.0_censat_v2.1.bed` | A more comprehensive centromere/satellite repeat annotation | T2T-CHM13v2.0 | [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed) |

## Assemblies and sequencing data 

### Human Pangenome Reference Consortium (HPRC)
| Dataset | Description | Source |
|---|---|---|
| ONT sequencing data | WGS Oxford nanopore with methylation tag | [Human Pangenome Data Explorer](https://data.humanpangenome.org/raw-sequencing-data) |
| Assemblies | trio hic phased Verkko 2.2.1 | [Human Pangenome Data Explorer](https://data.humanpangenome.org/assemblies?filter=%5B%7B%22categoryKey%22%3A%22release%22%2C%22value%22%3A%5B%222%22%5D%7D%5D) |
| Sample meta data | including sex, population, biosample etc | [Human Pangenome Data Explorer](https://data.humanpangenome.org/samples) |

### Washu pedigree samples
| Dataset | Description | Source |
|---|---|---|
| Washu pedigree assemblies | Assemblies after panpatch with minor curation (v1.3.1) | [washu-pedigree github](https://github.com/biomonika/washu-pedigree) |
| ONT FASTQ from LCLs| CELL ULTRA-LONG | [HPP github](https://github.com/biomonika/HPP/tree/main/T2T-Pedigree-project%20) |
| HiFi FASTQ from Blood | HiFi Data (methylcalled — bam file with MM/ML tags) | [HPP github](https://github.com/biomonika/HPP/tree/main/T2T-Pedigree-project%20) |

### Robertsonian translocation (ROB) samples
| Dataset | Description | Source |
|---|---|---|
| ROB assemblies | NA | [dbgap phs003920.v1.p1](https://www.nature.com/articles/s41586-025-09540-8) |
| ONT FASTQ from fibroblasts | NA | [dbgap phs003920.v1.p1](https://www.nature.com/articles/s41586-025-09540-8) |


## Coordinate Conventions

- BED files use 0-based, half-open coordinates.
- The genome assembly associated with each coordinate file is listed above.
- Personalized assembly coordinates refer to the corresponding sample and haplotype.