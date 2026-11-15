# Figure 3: Application of Pooch to Samples with Abnormal Karyotypes

This figure presents Pooch predictions for the T2T WashU pedigree,
the KOLF dataset across multiple cell types, and samples carrying
Robertsonian translocations. It also illustrates how parental origin
was evaluated separately for the two long arms of ROB chromosomes.

## Main Figure Panels

| Panel | Description | Script / notebook |
|---|---|---|
| Figure 3A | Pooch predictions for the T2T WashU pedigree samples. | [validation_samples_viz.ipynb](validation_samples_viz.ipynb) |
| Figure 3B | Pooch predictions for the KOLF dataset, including iPS cells, cortical neurons, NGN2-induced neurons, astrocytes, and microglia. | [validation_samples_viz.ipynb](validation_samples_viz.ipynb) |
| Figure 3C | Pooch predictions for samples carrying Robertsonian translocations. | [validation_samples_viz.ipynb](validation_samples_viz.ipynb) |
| Figure 3D | Schematic illustrating the application of Pooch to Robertsonian translocation chromosomes. | Schematic illustration |
| Figure 3E | Separate Pooch predictions for the two long arms forming each ROB chromosome. | [validation_samples_viz.ipynb](validation_samples_viz.ipynb) |

## Notes

- In Figures 3A–C, colors indicate the predicted probability that the
  input haplotype configuration corresponds to Mat − Pat.
- Values below 0.95 are highlighted and annotated.
- In Figure 3E, the left and right scores correspond to the two
  constituent long arms of each ROB chromosome, in the order shown.