# Pooch paper results

This directory contains track descriptions and the utility script used to
prepare the GRCh38 DMR track. Large track data files are hosted in the public
S3 location below.

- [Browse the Pooch 2026 tracks on S3](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/Pooch_2026/tracks/)

## Tracks

### Differential methylation regions (DMRs)

The DMR track reports the maternal-minus-paternal methylation difference.
Positive values indicate higher methylation on the maternal haplotype;
negative values indicate higher methylation on the paternal haplotype.
Regions are colored by parental direction and difference magnitude:

| Color | Meaning | Threshold |
|---|---|---|
| Red (`#FF0000`) | Strong maternal hypermethylation | `difference >= 0.15` |
| Pink (`#E78080`) | Moderate maternal hypermethylation | `0 < difference < 0.15` |
| Light blue (`#8484FC`) | Moderate paternal hypermethylation | `-0.15 < difference < 0` |
| Blue (`#0000FF`) | Strong paternal hypermethylation | `difference <= -0.15` |

The analysis used HPRC release v2 assemblies and Oxford Nanopore methylation
calls. Autosomal analyses included 123 samples; chromosome X analyses included
50 female samples balanced by inferred inactive-X parental origin. Differential
methylation was analyzed with DSS v2.58.0 in R v4.5.3 using a paired design
(`~ haplotype + sample`). CpG sites with more than 5% missing data were
excluded. DMRs were called at `P < 0.001`; neighboring regions within 1,000 bp
were merged when their direction of methylation difference agreed. The
CHM13v2.0 track contains 880 regions after this consolidation.

### Mean methylation difference

The BigWig tracks show the mean maternal-minus-paternal CpG methylation
difference. For each site, the difference was calculated within each sample
and then averaged across samples. Positive values indicate higher methylation
on the maternal haplotype; negative values indicate higher methylation on the
paternal haplotype.

Track descriptions:

- [`readme_iDMRs_T2T-CHM13v2.html`](readme_iDMRs_T2T-CHM13v2.html): CHM13v2.0 DMR track.
- [`readme_iDMRs_GRCh38.html`](readme_iDMRs_GRCh38.html): GRCh38 DMR track and post-processing.
- [`readme_mean_methylation_track.T2T-CHM13v2.html`](readme_mean_methylation_track.T2T-CHM13v2.html): CHM13v2.0 mean methylation track.
- [`readme_mean_methylation_track.GRCh38.html`](readme_mean_methylation_track.GRCh38.html): GRCh38 mean methylation track.

## Large track files on S3

The track files below are distributed from the S3 prefix at the top of this
README. Use the S3 listing to browse or download the objects. The listed DMR
counts are based on the prepared local files: 880 CHM13v2.0 regions and 870
GRCh38 regions.

| Track | S3 object |
|---|---|
| CHM13v2.0 DMR BED | `dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.bed.gz` |
| CHM13v2.0 DMR Tabix index | `dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.bed.gz.tbi` |
| GRCh38 DMR BED | `dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.mrg.coordinate_fix.bed.gz` |
| GRCh38 DMR Tabix index | `dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.mrg.coordinate_fix.bed.gz.tbi` |
| CHM13v2.0 mean methylation BigWig | `all_chr.missingrate_5.createDiffDatasetFromMethylDataset.flipped.mean.bw` |
| GRCh38 mean methylation BigWig | `all_chr.missingrate_5.createDiffDatasetFromMethylDataset.flipped.mean.grch38.bw` |

To download an object using the AWS CLI, substitute its name from the table:

```bash
aws s3 cp --no-sign-request \
  s3://human-pangenomics/publications/Pooch_2026/tracks/<object-name> .
```

The BED files contain nine tab-separated BED fields:

1. Chromosome.
2. Start (0-based, inclusive).
3. End (exclusive).
4. Hypermethylated haplotype (`Maternal` or `Paternal`).
5. Maternal-minus-paternal methylation difference.
6. Strand (`.`; not specified).
7. Thick start.
8. Thick end.
9. Item color (RGB hex value).

## Preparing the GRCh38 tracks

### DMR BED

```bash
CrossMap bed \
  chm13v2-grch38.chain \
  dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.bed.gz \
  dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.bed

python merge_bed_identical_annotation.py \
  dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.bed \
  dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.mrg.bed \
  --max-gap 10000

awk '{print $1,$2,$3,$4,$5,$6,$2,$3,$9}' OFS='\t' \
  dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.mrg.bed \
  > dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.mrg.coordinate_fix.bed

bgzip dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.mrg.coordinate_fix.bed
tabix -p bed dmrs.N05_filt.all_chr.1e3.mrg.forPlotting.grch38.mrg.coordinate_fix.bed.gz
```

CrossMap can split a source interval into multiple mapped fragments. Fragments originating from the same source interval were merged when they were located on the same chromosome and separated by no more than 10,000 bp. In the prepared outputs, 10 of the 880 source regions were
unmapped; the mapped fragments consolidated to 870 GRCh38 regions.

### Mean methylation BigWig

```bash
CrossMap bigwig \
  chm13v2-grch38.chain \
  all_chr.missingrate_5.createDiffDatasetFromMethylDataset.flipped.mean.bw \
  all_chr.missingrate_5.createDiffDatasetFromMethylDataset.flipped.mean.grch38.bw
```

## Local support files

| File | Description |
|---|---|
| `merge_bed_identical_annotation.py` | Merges nearby BED intervals with matching annotation columns. |
| `readme_iDMRs_GRCh38.html` | GRCh38 DMR track description. |
| `readme_iDMRs_T2T-CHM13v2.html` | CHM13v2.0 DMR track description. |
| `readme_mean_methylation_track.GRCh38.html` | GRCh38 mean methylation track description. |
| `readme_mean_methylation_track.T2T-CHM13v2.html` | CHM13v2.0 mean methylation track description. |
