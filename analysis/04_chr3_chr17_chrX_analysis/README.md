# Chromosome 3, 17, and X Analyses

This directory contains two utilities for examining chromosome-level
methylation patterns and XIST-associated effects.

## Scripts

| Script | Purpose |
|---|---|
| [`fast_cpg_interaction.py`](scripts/fast_cpg_interaction.py) | Fits a paired linear model at each CpG: maternal-minus-paternal methylation difference versus sample XIST difference. Reports the haplotype intercept and XIST interaction slope with p-values and Benjamini-Hochberg FDR values. CpGs are processed in chunks. |
| [`bin_sum.py`](scripts/bin_sum.py) | Sums BED column 4 values into non-overlapping genomic bins (default 50 bp) and writes a four-column bedGraph-like TSV. |

## XIST interaction input/output

`fast_cpg_interaction.py` requires:

- A methylation TSV with `Sample`, `Label`, and one column per CpG position.
  Labels may be `Maternal`/`Paternal` or `1`/`2` (`1` maps to Paternal,
  `2` to Maternal).
- An XIST TSV with a sample column and numeric XIST difference column.

Example:

```bash
python scripts/fast_cpg_interaction.py \
  --methylation chrX_methylation.tsv --xist xist.tsv \
  --out chrX_interaction.tsv --chrom chrX
```

The default column names are `Sample`, `Label`, `sample`, and `XIST_diff`.
Options configure these names, the minimum number of paired samples
(default 5), processing chunk size (default 10,000 CpGs), and samples whose
maternal/paternal orientation must be reversed. For each CpG the output
includes the number of paired samples, haplotype coefficient and p-value,
interaction coefficient and p-value, and separate FDR-adjusted values.
XIST values must already use the corrected orientation for any samples
passed with `--swap-samples`.

## BED binning

`bin_sum.py` requires a BED-like input with at least four columns. It groups
records by chromosome and `floor(start / bin-size) * bin-size`, sums column
4, and writes `chrom`, bin start, bin end, and sum without a header:

```bash
python scripts/bin_sum.py --input methylation.bed \
  --output methylation.50bp.bedgraph --bin-size 50
```

The scripts require pandas, NumPy, SciPy, and statsmodels (for the
interaction analysis); `bin_sum.py` requires pandas.