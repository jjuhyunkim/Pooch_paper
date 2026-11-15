# iDMR Identification

This step tests methylation differences between the two haplotypes across
samples and calls differentially methylated regions (DMRs) with DSS.

## Script

[`DSS.R`](scripts/DSS.R) takes one positional argument: a chromosome name
(for example, `chr1` or `chrX`). It loads `dss.fofn` from the working
directory. For `chrX`, it uses the sample IDs in `chrX.samples`; for other
chromosomes, it uses the unique sample IDs in the first column of `dss.fofn`.

For each sample, the script reads hap1 and hap2 methylation BED files in
CHM13 coordinates. It skips HG03942 on chr12 and HG02602 on chr19, and
swaps haplotype inputs for HG03579, HG02717, and HG02451 before assigning
the `mat` and `pat` labels. Samples with missing haplotype files are skipped.
Loci with zero coverage in more than 5% of the retained sample-haplotype
columns are removed. DSS then fits a paired model (`~ haplotype + pair`)
and tests the haplotype effect.

## Outputs

For each chromosome, the script writes:

- `dmls.N05_filt.<chromosome>.nCpG.csv`: retained and excluded CpG counts.
- `dmls.N05_filt.<chromosome>.all.csv`: per-CpG DSS test results.
- `dmrs.N05_filt.<chromosome>.5e2.csv`, `.1e3.csv`, and `.1e5.csv`: DMR
  calls at p-value thresholds 0.05, 0.001, and 0.00001.

All three DMR calls use `delta=0`, minimum length 50 bp, minimum 3 CpGs,
merge distance 100 bp, and at least 50% significant CpGs.

## Inputs and execution notes

The script expects files named
`winnowmap_ont-x-<sample>-verkkoV2.2.1-asmV4.F0x904.methyl.filt.hap1.chm13.bed`
and the corresponding `.hap2.chm13.bed` under
`/path/to/pooch/pooch_data/04.methyl_chm13/`. Replace that path and provide
the sample-list files and R packages `DSS`, `data.table`, and `bsseq` in the
execution environment. Outputs are written to the current working directory.

## Related analyses

- [Methylation processing](../00_methylation_processing/)
- [Model training and validation](../02_model_training_validation/)