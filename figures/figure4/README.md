# Figure 4: Detailed Characterization of Imprinting on Chromosomes 3, 17, and X

This figure examines methylation patterns and informative model coefficients
at selected loci on chromosomes 3, 17, and X. It further characterizes
parent-of-origin and XIST-associated methylation patterns at the DHRSX
and SRPK3 loci.

## Main Figure Panels

| Panel | Description | Script / notebook |
|---|---|---|
| Figure 4A | Genomic overviews of chromosomes 3 and 17 and promoter-associated DMRs, showing differential methylation, model coefficients, iDMRs, and parental methylation fractions. | IGV |
| Figure 4B | Genomic overview of chromosome X showing differential methylation, model coefficients, and iDMRs. | IGV |
| Figure 4C | Detailed views of the DHRSX intron and SRPK3 promoter, including parental and XIST effects and methylation patterns stratified by preferential X inactivation. | IGV |
| Figure 4D | Haplotype-specific methylation plotted against maternal-minus-paternal XIST promoter methylation differences. | [chrX_DHRSX_SRPK3.ipynb](chrX_DHRSX_SRPK3.ipynb) |
| Figure 4E | Paired comparisons of maternal and paternal methylation at the DHRSX iDMR and selected promoters, including SRPK3, XIST, and genes subject to or escaping XCI. | [chrX_DHRSX_SRPK3.ipynb](chrX_DHRSX_SRPK3.ipynb) |
| Figure 4F | Schematic illustrating the contributions of parental origin and random XCI to methylation patterns at DHRSX and SRPK3. | Schematic illustration |

## Additional Notebooks

| Notebook | Description |
|---|---|
| [dhrsx_xx_vs_xy.ipynb](dhrsx_xx_vs_xy.ipynb) | Comparison of DHRSX methylation between XX and XY samples. |
| [compare_chrX_models_using_washu_hprc.ipynb](compare_chrX_models_using_washu_hprc.ipynb) | Comparison of chromosome X models using WashU and HPRC samples. |

## Notes

- Differential methylation is calculated as maternal minus paternal
  methylation.
- Model coefficients were summed within non-overlapping 100-bp bins
  for visualization.
- Coefficient clusters contain at least three CpGs with coefficients
  of the same sign, with consecutive CpGs separated by no more than 500 bp.
- In Figures 4A and 4B, maternal and paternal methylation fractions
  are shown in red and blue, respectively.
- In Figure 4C, the parental-effect and XIST-effect tracks show the
  estimated contributions of the corresponding model terms to CpG
  methylation differences.
- Pat Xi and Mat Xi groups were defined by lower paternal and maternal
  XIST promoter methylation, respectively.
- In Figure 4D, negative XIST promoter methylation differences indicate
  preferential maternal X inactivation; positive differences indicate
  preferential paternal X inactivation.
- Comparisons in Figure 4E used paired t-tests.