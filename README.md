# Pooch Manuscript Resources

This repository contains analysis workflows, figure-generation notebooks,
reference documentation, and selected results accompanying the Pooch
manuscript.

The Pooch software is maintained separately at
https://github.com/marbl/Pooch. This repository documents the analyses and
figure workflows used to generate the manuscript results.

## Documentation

- [Analysis workflows](analysis/README.md): data processing, iDMR
  identification, model training and validation, prediction evaluation,
  chromosome-specific analyses, and computational benchmarks.
- [Figure and notebook index](figures/README.md): figure-specific notebooks
  and descriptions.
- Figure documentation: [Figure 1](figures/figure1/README.md) ·
  [Figure 2](figures/figure2/README.md) ·
  [Figure 3](figures/figure3/README.md) ·
  [Figure 4](figures/figure4/README.md) ·
  [Figure 5](figures/figure5/README.md).
- [Data and reference resources](data/README.md).
- [Results and genome browser tracks](results/README.md).

## UCSC Genome Browser Tracks

The mean CpG methylation differences between maternal and paternal haplotypes
(maternal minus paternal) and the iDMRs identified in this study can be viewed
in the [UCSC Genome Browser](UCSC_BROWSER_LINK).

## Repository organization

- [`data/`](data/): Input data documentation and reference annotations.
- [`analysis/`](analysis/): Data preparation, statistical analyses,
  prediction experiments, and benchmarking.
- [`figures/`](figures/): Figure notebooks, figure-specific documentation,
  and figure-generation resources.
- [`results/`](results/): Track descriptions, supporting files, and links to
  selected analysis results.

## Data availability

HPRC sequencing data and assemblies are available from
[Human Pangenome Data Explorer](https://data.humanpangenome.org/raw-sequencing-data).

Reference genomes and annotations are documented in
[data/README.md](data/README.md). Large input files and intermediate outputs
are not included in this repository.

## Citation

If you use this code, please cite:
[TBD]

## License

See `LICENSE` for the terms governing the code in this repository. External datasets and annotations remain subject to their original licenses.