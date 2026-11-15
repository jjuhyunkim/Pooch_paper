# Read Alignment and Methylation Processing

This directory contains the preprocessing steps used to prepare
haplotype-specific methylation data for the Pooch manuscript analyses.
The checked-in scripts align HiFi or ONT reads, sort and merge alignments,
separate haplotypes, call methylation with modkit, and generate
haplotype-to-CHM13 chain files.

## Analysis Summary

For each sample, read splits are aligned to the sample's Verkko v2.2.1
assembly using minimap2 or Winnowmap. `05-align.slurm` accepts HiFi or ONT
reads and removes unmapped, secondary, and supplementary alignments
(`samtools view -F 0x904`). The split BAMs are sorted/indexed and merged;
the merged alignment is then subset by the haplotype sequence IDs supplied
in BED files.

`09-modkitPileup.slurm` runs `modkit pileup` with the `traditional` preset
against the personalized assembly, then compresses and indexes the
bedMethyl output and creates haplotype-specific bedMethyl files. The chain
helper scripts set up the sample directory, run nf-LO/minimap2, and convert,
split, and invert the resulting alignments. For a new sample,
`sh.create_chm13_chain_for_sample` initially submits the haplotype 1 job;
it can resubmit either haplotype if an existing chain output is empty.

Coordinate conversion of the methylation calls to CHM13 and cohort-level
combination/missingness filtering are not implemented by the scripts in this
directory; those steps require other project workflows.

## Scripts

| Script | Description |
|---|---|
| [`05-align.slurm`](scripts/05-align.slurm) | Aligns specified read splits; arguments: sample ID, `hifi`/`ont`, `minimap2`/`winnowmap`, and optionally split numbers/ranges. |
| [`06-sortBam.slurm`](scripts/06-sortBam.slurm) | Sorts and indexes the split BAMs for the selected sample, read type, and aligner. |
| [`07-mergeSortedBams.slurm`](scripts/07-mergeSortedBams.slurm) | Merges sorted split BAMs into one indexed BAM (or renames the sole input BAM). |
| [`08-hapSplitMergedBam.slurm`](scripts/08-hapSplitMergedBam.slurm) | Subsets the merged BAM using per-haplotype sequence-ID BED files; supports `mat`/`pat` or `hap1`/`hap2`. |
| [`09-modkitPileup.slurm`](scripts/09-modkitPileup.slurm) | Calls bedMethyl data, compresses/indexes it, and subsets calls by the same haplotype BED files. |
| [`sh.create_chm13_chain_for_sample`](scripts/sh.create_chm13_chain_for_sample) | Sets up a sample directory and submits the initial haplotype 1 nf-LO job; it can resubmit either haplotype when an existing chain output is empty. Takes a sample ID. |
| [`sh.nflo_pipeline`](scripts/sh.nflo_pipeline) | Runs nf-LO and chain conversion for source FASTA, target FASTA, and output directory (three positional arguments). |

## Inputs and Outputs

| Stage | Main inputs | Main outputs |
|---|---|---|
| Alignment | `reads/{hifi,ont}/` FASTQ splits and `{sample}-verkkoV2.2.1-asmV4.fa.gz` | BAM per requested split in `alns/` |
| Sort/merge | Split BAMs | Sorted/indexed split BAMs and one merged BAM in `alns/` |
| Haplotype split | Merged BAM and `{sample}-verkkoV2.2.1-asmV4.{mat,pat,hap1,hap2}.bed` | Haplotype-specific BAMs |
| Methylation pileup | Merged BAM and uncompressed personalized reference FASTA | Compressed/indexed merged and haplotype-specific bedMethyl |
| Chain generation | CHM13 reference and haplotype FASTAs | Split/inverted chain files and intermediate PAF/chain outputs |

The SLURM scripts load site-specific software modules (including samtools
1.21, modkit 0.4.1, htslib 1.21, and bedtools 2.31.1). The chain scripts also
expect nf-LO, chain-tools, and cluster-specific resources. Update those
environment and data paths before running outside the original cluster.

## Related Analyses

- [iDMR identification](../01_idmr_identification/)
- [Model training and validation](../02_model_training_validation/)