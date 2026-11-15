# Figure 2: Pooch Model Training and Parent-of-Origin Prediction

This figure illustrates Pooch model training and prediction, evaluates
classification performance, and characterizes informative CpGs selected
by the models. It also compares prediction accuracy across CpG subsets
and between Pooch and the previously published PatMat method.

## Main Figure Panels

| Panel | Description | Script / notebook |
|---|---|---|
| Figure 2A | Schematic of model training using CpG methylation differences between maternal and paternal haplotypes. | Schematic illustration |
| Figure 2B | Schematic of parent-of-origin prediction using methylation differences between haplotypes of unknown parental origin. | Schematic illustration |
| Figure 2C | Model accuracy from nested cross-validation on the training dataset and evaluation on the held-out test dataset. | [boxplot_prob.ipynb](boxplot_prob.ipynb) |
| Figure 2D | Predicted probabilities for the Mat − Pat configuration across held-out test samples, highlighting low-confidence and incorrect predictions. | [boxplot_prob.ipynb](boxplot_prob.ipynb) |
| Figure 2E | Chromosomal distribution of maternal-minus-paternal methylation differences, with point sizes indicating absolute model coefficients. | [CpG_Coefficient_distribution_peak_calling.ipynb](CpG_Coefficient_distribution_peak_calling.ipynb) |
| Figure 2F | Distribution of informative CpGs within and outside PatMat iDMRs and the iDMRs identified in this study. | [fisherTest.ipynb](fisherTest.ipynb) |
| Figure 2G | Test accuracy of models trained using all CpGs, CpGs within PatMat iDMRs, or CpGs at least 1 Mb outside PatMat iDMRs. | [compare_performance_between_models_with_flanking.ipynb](compare_performance_between_models_with_flanking.ipynb) |
| Figure 2H | Percentages of samples correctly and incorrectly classified by Pooch and PatMat. | [hprc_patmat.ipynb](hprc_patmat.ipynb) |

## Additional Notebooks

| Notebook | Description |
|---|---|
| [fisherTest.ipynb](fisherTest.ipynb) | Fisher's exact test analysis. |
| [missing_idmr.ipynb](missing_idmr.ipynb) | Investigation of missing iDMRs. |

## Notes

- The training dataset included 80 samples for autosomes and 46 samples
  for chromosome X.
- Prediction probabilities indicate whether the input haplotype
  configuration corresponds to Mat − Pat.
- In Figure 2C, red indicates nested cross-validation accuracy and blue
  indicates held-out test accuracy.
- In Figure 2D, samples with probabilities below 0.6 are labeled in gray,
  with those below 0.5 labeled in red.
- Informative CpGs are CpGs assigned nonzero coefficients by the models.
  Point sizes in Figure 2E represent absolute coefficient values.
- PatMat iDMRs and the PatMat method refer to the study associated with
  PMID: 36777186.