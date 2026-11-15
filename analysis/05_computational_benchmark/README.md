# Computational Benchmark

This step contains two related but distinct evaluations:

1. Prediction resource usage (elapsed time and peak memory) across
   chromosomes.
2. The effect of training-set size on chr22 elastic-net model performance.

## Training-size analysis

The training-size experiment ran the same `elasticnet.nestedcv.py` workflow
from the separate `pooch_data/09.subsampling_analysis/` directory, with
`chr22` and each requested training size: 10, 15, 20, 25, 30, 40, 50, 60,
70, and 80 samples. The script uses the same fixed test set and repeats its
grouped nested 5-fold cross-validation, FDR feature selection, and
elastic-net hyperparameter search for each training size. It writes one
`elasticnet_selected_cpgs_chr22.<N>.csv` summary and a corresponding
grid-search result for each size, as well as fitted model files.

In the script, the size argument is used to keep the first N rows of the
training subset before the orientation-augmentation step. Thus, this is a
training-sample-count series, not a percentage series or a random sample
unless the input rows were randomized beforehand.

The included scripts are [`array.sh`](scripts/array.sh) and
[`elasticnet.nestedcv.py`](scripts/elasticnet.nestedcv.py). The array script
currently lists 5 as well as 10–80 in `props`, but sets `--array=0`, so one
job runs only the first value (5). To reproduce the reported 10–80 series,
set `props` to those ten sample counts and configure SLURM array indices
`0-9`.

## Prediction resource benchmark

The included HG002 workflow uses
[`benchmark.sh`](scripts/benchmark.sh)
to run `/usr/bin/time -v` and capture per-chromosome prediction resource
usage. It iterates through chr1–chr22 and chrX;
[`predict.autosome.sh`](scripts/predict.autosome.sh)
and
[`predict.chrX.sh`](scripts/predict.chrX.sh)
prepare methylation inputs and run the appropriate prediction command. This
is a prediction-time resource benchmark, distinct from the chr22
training-size experiment above.

The scripts depend on site-specific Pooch tools, data files, model files,
and absolute paths. Adapt these paths and the SLURM configuration before
reproduction; resource measurements may also depend on the cluster and
software environment.