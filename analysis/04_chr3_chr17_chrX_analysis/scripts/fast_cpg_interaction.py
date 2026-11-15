#!/usr/bin/env python3
"""Fast chromosome-wide paired maternal-minus-paternal CpG interaction analysis."""
import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--methylation', required=True, help='Wide TSV: Label, Sample, CpG position columns')
    p.add_argument('--xist', required=True, help='TSV with sample and XIST differential methylation')
    p.add_argument('--out', required=True, help='Output TSV')
    p.add_argument('--sample-col', default='Sample')
    p.add_argument('--label-col', default='Label')
    p.add_argument('--xist-sample-col', default='sample')
    p.add_argument('--xist-value-col', default='XIST_diff')
    p.add_argument('--min-pairs', type=int, default=5)
    p.add_argument('--chunk-size', type=int, default=10000)
    p.add_argument('--swap-samples', nargs='*', default=[], help='Samples whose maternal/paternal CpG orientation must be reversed; XIST values must already use corrected orientation')
    p.add_argument('--chrom', default='chrX')
    return p.parse_args()


def analyze_block(y, x, min_pairs):
    valid = np.isfinite(y) & np.isfinite(x[:, None])
    yy = np.where(valid, y, 0.0)
    w = valid.astype(float)
    xx = x[:, None]
    n = w.sum(0)
    sx = (w * xx).sum(0)
    sy = yy.sum(0)
    sxx = (w * xx * xx).sum(0)
    sxy = (yy * xx).sum(0)
    syy = (yy * yy).sum(0)
    centered_xx = sxx - np.divide(sx * sx, n, out=np.zeros_like(sxx), where=n>0)
    ok = (n >= max(min_pairs, 3)) & (centered_xx > 1e-12)
    intercept = np.full(y.shape[1], np.nan)
    slope = np.full(y.shape[1], np.nan)
    p_int = np.full(y.shape[1], np.nan)
    p_slope = np.full(y.shape[1], np.nan)
    slope[ok] = (sxy[ok] - sx[ok]*sy[ok]/n[ok]) / centered_xx[ok]
    intercept[ok] = (sy[ok] - slope[ok]*sx[ok]) / n[ok]
    rss = syy - intercept*sy - slope*sxy
    rss = np.maximum(rss, 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        var = rss / (n - 2)
        se_slope = np.sqrt(var / centered_xx)
        se_int = np.sqrt(var * (1/n + (sx/n)**2/centered_xx))
        t0 = intercept/se_int
        t1 = slope/se_slope
        p_int[ok] = 2*student_t.sf(np.abs(t0[ok]), n[ok]-2)
        p_slope[ok] = 2*student_t.sf(np.abs(t1[ok]), n[ok]-2)
    # Handle exactly zero residual variance (perfect fit).
    for beta, se, pp in [(intercept, se_int, p_int), (slope, se_slope, p_slope)]:
        perfect = ok & (se == 0)
        pp[perfect & (beta != 0)] = 0.0
        pp[perfect & (beta == 0)] = 1.0
    return n.astype(int), intercept, slope, p_int, p_slope


def main():
    a = parse_args()
    if a.chunk_size < 1 or a.min_pairs < 3:
        raise ValueError('--chunk-size must be >= 1 and --min-pairs must be >= 3')
    xdf = pd.read_csv(a.xist, sep='\t')
    if xdf[a.xist_sample_col].duplicated().any():
        raise ValueError('Duplicate sample IDs in XIST table')
    xist = pd.Series(pd.to_numeric(xdf[a.xist_value_col], errors='coerce').to_numpy(), index=xdf[a.xist_sample_col].astype(str))
    # Read only header first, then CpG columns in chunks to bound memory usage.
    columns = pd.read_csv(a.methylation, sep='\t', nrows=0).columns.tolist()
    if a.sample_col not in columns or a.label_col not in columns:
        raise ValueError('Methylation input must contain sample and label columns')
    cpg_cols = [c for c in columns if c not in (a.sample_col, a.label_col)]
    if not cpg_cols:
        raise ValueError('No CpG columns found')
    sample_labels = pd.read_csv(a.methylation, sep='\t', usecols=[a.sample_col, a.label_col], dtype=str)
    sample_labels[a.label_col] = sample_labels[a.label_col].replace({'1':'Paternal','2':'Maternal'})
    if sample_labels.duplicated([a.sample_col, a.label_col]).any():
        raise ValueError('Duplicate (Sample, Label) rows; aggregate them before analysis')
    m = sample_labels[sample_labels[a.label_col]=='Maternal'].set_index(a.sample_col)
    p = sample_labels[sample_labels[a.label_col]=='Paternal'].set_index(a.sample_col)
    samples = m.index.intersection(p.index).intersection(xist.index)
    if len(samples) < a.min_pairs:
        raise ValueError(f'Only {len(samples)} paired samples with XIST data')
    # Positions are unique because of the duplicate check above.
    maternal_row = sample_labels.index[(sample_labels[a.sample_col].isin(samples)) & (sample_labels[a.label_col]=='Maternal')]
    paternal_row = sample_labels.index[(sample_labels[a.sample_col].isin(samples)) & (sample_labels[a.label_col]=='Paternal')]
    m_rows = pd.Series(maternal_row.to_numpy(), index=sample_labels.loc[maternal_row,a.sample_col]).loc[samples].to_numpy()
    p_rows = pd.Series(paternal_row.to_numpy(), index=sample_labels.loc[paternal_row,a.sample_col]).loc[samples].to_numpy()
    x = xist.loc[samples].to_numpy(float)
    flip = np.isin(samples.to_numpy(), a.swap_samples)
    output = []
    for start in range(0, len(cpg_cols), a.chunk_size):
        cols = cpg_cols[start:start+a.chunk_size]
        block = pd.read_csv(a.methylation, sep='\t', usecols=[a.sample_col, a.label_col, *cols])
        # usecols does not guarantee requested order
        vals = block[cols].apply(pd.to_numeric, errors='coerce').to_numpy(float)
        y = vals[m_rows] - vals[p_rows]
        y[flip] *= -1
        n, b0, b1, p0, p1 = analyze_block(y, x, a.min_pairs)
        output.append(pd.DataFrame({'chr':a.chrom,'cpg':cols,'n_pairs':n,'haplotype_coef':b0,'interaction_coef':b1,'haplotype_pvalue':p0,'interaction_pvalue':p1}))
        print(f'Processed {min(start+len(cols),len(cpg_cols))}/{len(cpg_cols)} CpGs', file=sys.stderr)
    result = pd.concat(output, ignore_index=True)
    for term in ['haplotype','interaction']:
        pcol = f'{term}_pvalue'
        ok = np.isfinite(result[pcol].to_numpy(float))
        result[f'{term}_fdr'] = np.nan
        if ok.any():
            result.loc[ok, f'{term}_fdr'] = multipletests(result.loc[ok,pcol], method='fdr_bh')[1]
    result.to_csv(a.out, sep='\t', index=False)
    print(f'Saved {len(result)} CpGs to {a.out}; paired samples: {len(samples)}', file=sys.stderr)

if __name__ == '__main__':
    main()
