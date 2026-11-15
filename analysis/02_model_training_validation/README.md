# Model Training and Validation

This step trains parental-origin classifiers from paired haplotype
methylation differences and evaluates them with sample-grouped cross-validation
and held-out samples. It also includes the manuscript comparison of
parental-origin predictions with PatMat; that comparison belongs to this
model validation stage, not to the cohort-specific prediction scripts in
[03_prediction_evaluation](../03_prediction_evaluation/).

## Scripts

| Script | Purpose |
|---|---|
| [`make_input.sh`](scripts/make_input.sh) | SLURM workflow to create a methylation dataset, filter CpGs by training-set missingness, and create the differential methylation input table used by the classifiers. |
| [`filter_missing_rate.py`](scripts/filter_missing_rate.py) | Keeps CpG columns whose missing-value count in the training samples is at or below a specified percentage. |
| [`elasticnet.nestedcv.py`](scripts/elasticnet.nestedcv.py) | Selects CpGs with an F-test/Benjamini-Hochberg FDR screen, tunes logistic elastic-net `C` and `l1_ratio` in grouped nested 5-fold cross-validation, fits the final model, and reports held-out performance. Takes a chromosome argument. |
| [`elasticnet_xist_adjusted.py`](scripts/elasticnet_xist_adjusted.py) | Fits a grouped nested-CV elastic-net classifier after estimating and removing CpG-specific XIST-associated components while preserving parental effects. Requires a chromosome and `--xist` TSV. |
| [`filter_positions.py`](scripts/filter_positions.py) | Removes wide-table position columns overlapping chromosome-specific BED intervals, optionally padded. |
| [`check_filtered_postiions.py`](scripts/check_filtered_postiions.py) | Checks whether the filtered position columns exactly match the original columns outside the BED intervals and padding. |

## Differential dataset preparation

`make_input.sh` describes the preprocessing sequence used before model
training:

1. Run the external `createMethylDataset.py` on a chromosome using a sample
   FOFN and sex-mapping table.
2. Run `filter_missing_rate.py <chromosome> <missing-rate-percent>` on the
   resulting `<chromosome>.createMethylDataset.out`.
3. Run the external `createDiffDatasetFromMethylDataset.py` on the filtered
   table to create
   `<chromosome>.missingrate_<percent>.createDiffDatasetFromMethylDataset.out`.

The filter reads the chromosome's methylation dataset and selects feature
columns based only on rows whose sample IDs appear in the training sample
list. A CpG is retained when its number of missing training values is no
greater than `training-row-count * missing-rate-percent / 100`. The first
three metadata columns are retained, and the selected columns are then
written for all rows in the input table. For chr19, HG02602 is removed before
filtering. The script uses fixed sample-list paths, which must be updated for
another environment.

The checked-in `make_input.sh` is a site-specific SLURM job: its array is
currently set to indices `0,8` (chr1 and chr9), it assumes
`filter_missing_rate.py` is in the working directory, and it calls external
`pofo_tools` scripts via absolute paths. Its inputs include the methylation
FOFN and sex-mapping file configured in the script. Adapt the array, paths,
resources, and external-tool locations before running elsewhere. The
classifier input described below is the final differential table from this
workflow.

## Model workflows

`elasticnet.nestedcv.py` reads
`<chromosome>.missingrate_5.createDiffDatasetFromMethylDataset.out` and
training/test sample lists. It reverses the signs of the three configured
swap samples, augments rows with the reversed parental orientation, and uses
`GroupKFold` so rows from the same sample stay in one fold. The pipeline
screens features with `SelectFdr(f_classif, alpha=0.05)` and tunes logistic
elastic-net `C` over 0.01–1000 and `l1_ratio` over 0–1. It selects by
negative Brier score, reports nested-CV AUC and test metrics, and saves the
model, grid search, and a one-row summary CSV.

`elasticnet_xist_adjusted.py` uses the same wide CpG-difference table shape:
the first three columns must be `Sample`, `LeftOpHapLabel`, and
`RightOpHapLabel`; remaining columns are numeric CpG positions. Its XIST
input must contain `Sample` and `XIST_diff`. The script makes sign-reversed
orientation rows, performs sample-grouped nested CV by default, tunes the
elastic-net grid by Brier score, and evaluates the held-out set. The
CpG-specific XIST adjustment and optional FDR feature screen are fitted
inside training folds. `--fast-final` skips outer nested CV. Additional
options configure sample lists, BED-based exclusion, caching, folds,
regularization grids, and sign-correction samples; see `--help`.

The XIST-adjusted workflow writes `<prefix>.metrics.csv`,
`.gridsearch.csv`, `.coefficients.tsv`, `.selected_terms.tsv`,
`.selected_cpgs.tsv`, `.adjustment_slopes.tsv`,
`.effective_raw_coefficients.tsv`, `.test.predictions.tsv`, and
`.model.pkl`; nested-CV predictions are also written when outer CV is run.

## Position filtering

Both position utilities accept BED intervals as 0-based, half-open
coordinates. For example:

```bash
python scripts/filter_positions.py \
  --table input.tsv --bed regions.bed --chrom chrX \
  --padding 100kb --out filtered.tsv

python scripts/check_filtered_postiions.py \
  --original input.tsv --filtered filtered.tsv --bed regions.bed \
  --chrom chrX --padding 100kb
```

The table's numeric CpG column names are interpreted as 0-based by default;
set `--position-base 1` for 1-based positions. The checker exits nonzero if
the filtered positions differ from the expected set.

## Data and environment

The model scripts refer to project-specific sample lists and input tables
that are not included here. The first script has fixed sample-list paths in
its source; update them for the local data layout. The XIST-adjusted script
allows these paths to be supplied with `--input`, `--train-samples`, and
`--test-samples`. Its current non-chrX default sample-list expression refers
to a commented-out `root` variable, so pass both sample-list options for
non-chrX runs. Python dependencies include pandas, NumPy, SciPy,
scikit-learn, and joblib; `cloudpickle` is used for model serialization.