# Pooch paper analyses

This repository contains the analysis scripts and figure-generation code accompanying the Pooch manuscript

The Pooch software is maintained separately at https://github.com/marbl/Pooch. This repository documents the configurations and analyses used to generate the results reported in the manuscript.

## UCSC Genome Browser Tracks

The mean CpG methylation differences between maternal and paternal haplotypes
(maternal minus paternal) and the iDMRs identified in this study can be viewed
in the [UCSC Genome Browser](UCSC_BROWSER_LINK).

## Repository organization

- `data/`: Input data documentation, reference annotations, and example datasets.
- `analysis/`: Data preparation, statistical analyses, prediction experiments, and benchmarking.
- `figures/`: Code for generating the manuscript figures.
- `results/`: Output organization and selected analysis results.

## Data availability

HPRC Sequencing data and assemblies are available from [Human Pangnome Data Explorer](https://data.humanpangenome.org/raw-sequencing-data).

Reference genomes and annotations are documented in `data/README.md`. Large input files and intermediate outputs are not included in this repository.

## Citation

If you use this code, please cite:
[TBD]

## License

See `LICENSE` for the terms governing the code in this repository. External datasets and annotations remain subject to their original licenses.