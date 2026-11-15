# Prediction Evaluation

This step applies the XIST-adjusted chrX parental-origin model to the WashU
pedigree and Robertsonian translocation (ROB) methylation datasets. The two
checked-in shell scripts run the same three-stage workflow for separate
cohorts: build a model-filtered methylation dataset, create haplotype
differential features, then predict parental origin using the trained model
and cohort XIST values.

## Scripts

| Script | Cohort | Cohort inputs |
|---|---|---|
| [`prediction_washu.sh`](scripts/prediction_washu.sh) | WashU pedigree | `wash.fofn`, `sex.txt`, and `chrX_promoter_methylation_washu.tsv` |
| [`prediction_rob.sh`](scripts/prediction_rob.sh) | Robertsonian translocation | `rob.fofn`, `sex.txt`, and `chrX_promoter_methylation_rob.tsv` |

Each script calls the external Pooch `pofo_tools` utilities in this order:

1. `createMethylDataset_xist.py` builds `prediction_xist/chrX.createMethylDataset.out`
   using the cohort FOFN, sex mapping, and model CpG filter.
2. `createDiffDatasetFromMethylDataset.py` creates
   `prediction_xist/chrX.createDiffDatasetFromMethylDataset.out`, using
   `Hap1` and `Hap2` as the left and right haplotype labels.
3. `predict_parent_origin_xist.py` predicts from that differential table,
   the model, and cohort XIST values, writing `chrX.predictions.tsv`.

The intermediate datasets and logs are written under `prediction_xist/`;
the prediction TSV is written in the current working directory. Run the two
scripts in separate working directories to avoid overwriting the same
intermediate and prediction filenames.

## Execution requirements

These are site-specific scripts, not standalone commands: they use absolute
paths for `pofo_tools` and cohort data, and locate a model with
`realpath ../*pkl` relative to the current working directory. Update those
paths and ensure exactly the intended model file is found. The cohort files
and model are not included in this repository. For dataset context, see
[data documentation](../../data/README.md); the broader workflow is
summarized in the [analysis overview](../README.md).