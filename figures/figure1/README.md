# Figure 1: Categorization of DMRs Identified from 123 HPRC Trios, Including Chromosome X

This figure summarizes the genomic distribution, methylation differences,
and lengths of the iDMRs identified in this study. It also compares these
regions with previously reported PatMat iDMRs and presents representative
loci in the UCSC Genome Browser.

## Main Figure Panels

| Panel | Description | Script / notebook |
|---|---|---|
| Figure 1A | Genomic distribution of iDMRs identified in this study, colored by the magnitude and direction of maternal-minus-paternal methylation differences, alongside previously reported PatMat iDMRs. | [ideogram.ipynb](ideogram.ipynb) |
| Figure 1B | Comparison of absolute methylation differences between maternally and paternally hypermethylated iDMRs. | [dmr_investigation.ipynb](dmr_investigation.ipynb) |
| Figure 1C | Comparison of region lengths between maternally and paternally hypermethylated iDMRs. | [dmr_investigation.ipynb](dmr_investigation.ipynb) |
| Figure 1D | Overlap between iDMRs identified in this study and PatMat iDMRs, with chromosome X regions shown separately. | [dmr_investigation.ipynb](dmr_investigation.ipynb) |
| Figure 1E | Counts of iDMRs by methylation category and overlap with PatMat iDMRs, including matched, mixed-sign, and missing regions. | [dmr_investigation.ipynb](dmr_investigation.ipynb) |
| Figure 1F | Distribution of maternal-minus-paternal methylation differences, with classification thresholds at −0.15 and 0.15. | [dmr_investigation.ipynb](dmr_investigation.ipynb) |
| Figure 1G | Genome browser views of the KCNQ1OT1/KCNQ1 and DHRSX/LOC124905239 loci, showing mean CpG methylation differences, PatMat iDMRs, iDMRs identified in this study, and gene annotations. | [UCSC Genome Browser](UCSC_BROWSER_LINK) |

## Supporting Analysis

| Analysis | Description | Script / notebook |
|---|---|---|
| Chromosome X sample selection | Selection of 50 female samples balanced for preferential maternal or paternal X inactivation, inferred from haplotype-specific XIST promoter methylation. | [XIST_balance_sampling.ipynb](XIST_balance_sampling.ipynb) |

## Notes

- Methylation differences are calculated as maternal minus paternal
  methylation. Positive values indicate maternal hypermethylation;
  negative values indicate paternal hypermethylation.
- Strong methylation differences are defined as absolute differences
  ≥ 0.15; moderate differences are < 0.15.
- In Figure 1A, iDMRs were extended by 500 kb on each side for visualization.
  These extensions were not used for the region-length comparisons.
- P-values in Figures 1B and 1C were calculated using the Wilcoxon
  rank-sum test.
- PatMat iDMRs were obtained from the study associated with
  PMID: 36777186. Counts in parentheses in Figures 1D and 1E indicate
  unique PatMat iDMRs.
- Of the 48 PatMat iDMRs not recovered in the initial analysis, 10 showed
  significant parental-haplotype methylation differences in region-based
  paired t-tests, with directions concordant with previous reports.
- Gene annotations in Figure 1G use T2T-CHM13v2.0
  JHU RefSeq + Liftoff v5.3.