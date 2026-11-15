# Figure 5: Pooch Snakemake Pipeline and Computational Resource Requirements

This figure outlines the Pooch prediction and custom model-training
pipelines, evaluates the effect of training sample size on prediction
performance, and summarizes computational resource usage by chromosome.

## Main Figure Panels

| Panel | Description | Script / notebook |
|---|---|---|
| Figure 5A | End-to-end Pooch prediction pipeline, showing computational steps and their inputs and outputs. | Pipeline diagram |
| Figure 5B | Pooch-build pipeline for training customized models. | Pipeline diagram |
| Figure 5C | Test accuracy and Brier score as a function of training sample size, using chromosome 22 as a test case. | [benchmark.ipynb](benchmark.ipynb) |
| Figure 5D | CPU time and peak memory usage for the Pooch prediction step by chromosome. | [benchmark.ipynb](benchmark.ipynb) |

## Notes

- In Figures 5A and 5B, nodes represent computational steps, and
  document symbols on the edges represent inputs and outputs.
- Figure 5C evaluates prediction performance using test accuracy and
  Brier score. The Brier score measures the mean squared difference
  between predicted probabilities and observed binary labels;
  lower values indicate better performance.
- In Figure 5D, chromosomes with the highest and lowest CPU time
  and memory usage are highlighted.