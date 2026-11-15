#!/usr/bin/env python3
"""Grouped nested CV: CpG-specific XIST adjustment, then logistic elastic net.

Within each training fold fit d_j = a_j*h + b_j*z by least squares,
where h=2*y-1 and z=XIST_diff. Then transform r_j=d_j-b_j*z.
Keep the parental component a_j*h; NEVER subtract the full fitted value.
Predict from adjusted CpGs only: eta=sum(beta_j*r_j), no intercept.
All classifier coefficients are penalized. No freely fitted XIST term.
Every adjusted feature reverses sign when d and z both reverse sign.

XIST TSV: Sample, XIST_diff (header required). Its difference direction MUST
match the input CpG differences BEFORE sample swaps. All feature columns after
the first three metadata columns must be numeric CpG differences.

Optional FDR screening uses adjusted CpGs, only within training folds.
Its p-values are a feature-screening heuristic, NOT inferential evidence:
adjustment uses labels and augmented rows are not independent replicates.
Default --fdr-alpha 1 skips screening.
No scaling or intercept is used, as in the original script.
"""

import argparse
import json
import time
import warnings
from pathlib import Path

from joblib import Parallel, delayed
from joblib.externals import cloudpickle
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy import special
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.feature_selection import SelectFdr, f_classif
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             brier_score_loss, classification_report,
                             confusion_matrix, roc_auc_score)
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.utils.validation import check_is_fitted


def log(message):
    print(message, flush=True)


def read_tsv(path, **kwargs):
    engine = kwargs.pop('engine')
    try:
        return pd.read_csv(path, sep='\t', engine=engine, **kwargs)
    except (ImportError, pd.errors.ParserError) as error:
        if engine != 'pyarrow':
            raise
        log(f'pyarrow failed; falling back to pandas C engine ({error})')
        return pd.read_csv(path, sep='\t', engine='c', **kwargs)


def file_fingerprint(path):
    if path is None:
        return None
    path = Path(path).resolve()
    stat = path.stat()
    return {'path': str(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def cache_key(args, chrom, input_path, train_path, test_path):
    return {
        'version': 1,
        'chromosome': chrom,
        'input': file_fingerprint(input_path),
        'xist': file_fingerprint(args.xist),
        'train_samples': file_fingerprint(train_path),
        'test_samples': file_fingerprint(test_path),
        'exclude_bed': file_fingerprint(args.exclude_bed),
        'cpg_coordinate_base': args.cpg_coordinate_base,
    }


def load_preprocessed_cache(path, expected_key, metadata):
    if not path or not Path(path).exists():
        return None
    t0 = time.perf_counter()
    with np.load(path, allow_pickle=False) as cache:
        actual_key = json.loads(str(cache['cache_key']))
        if actual_key != expected_key:
            log(f'Ignoring stale preprocessed cache: {path}')
            return None
        cpg_columns = cache['cpg_columns'].astype(str).tolist()
        row_metadata = pd.DataFrame(cache['row_metadata'].astype(str), columns=metadata)
        feature_values = cache['feature_values']
    log(
        f'Loaded preprocessed cache {path}: {len(row_metadata)} rows x '
        f'{feature_values.shape[1] - 1} CpGs in {time.perf_counter() - t0:.2f}s'
    )
    return row_metadata, cpg_columns, feature_values


def save_preprocessed_cache(path, key, row_metadata, cpg_columns, feature_values, metadata):
    if not path:
        return
    t0 = time.perf_counter()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        cache_key=json.dumps(key, sort_keys=True),
        row_metadata=np.asarray(row_metadata.loc[:, metadata].astype(str).to_numpy(), dtype=str),
        cpg_columns=np.asarray(cpg_columns, dtype=str),
        feature_values=feature_values,
    )
    log(f'Saved preprocessed cache {path} in {time.perf_counter() - t0:.2f}s')


class CpGXISTAdjuster(BaseEstimator, TransformerMixin):
    """Learn CpG-specific XIST slopes while preserving parental effects."""

    def __init__(self, alpha=1.0, cpg_names=None):
        self.alpha = alpha
        self.cpg_names = cpg_names

    def fit(self, X, y):
        if hasattr(X, 'columns'):
            self.feature_names_in_ = np.asarray(X.columns, dtype=object)
            self.cpg_names_ = np.array([c for c in X.columns if c != 'XIST_diff'])
            d = X.loc[:, self.cpg_names_].to_numpy(dtype=float)
            z = X['XIST_diff'].to_numpy(dtype=float)
        else:
            X = np.asarray(X, dtype=float)
            if self.cpg_names is None:
                raise ValueError('cpg_names are required when fitting CpGXISTAdjuster on a numpy array.')
            if X.ndim != 2 or X.shape[1] != len(self.cpg_names) + 1:
                raise ValueError('Array input must contain CpGs followed by one XIST_diff column.')
            self.cpg_names_ = np.asarray(self.cpg_names, dtype=object)
            self.feature_names_in_ = np.concatenate([self.cpg_names_, ['XIST_diff']])
            d = X[:, :len(self.cpg_names_)]
            z = X[:, -1]
        h = 2 * np.asarray(y, dtype=float) - 1
        if not np.array_equal(np.unique(h), [-1, 1]):
            raise ValueError('Adjustment requires binary labels with both classes.')
        if np.all(z == 0):
            self.parental_effects_ = h @ d / len(h)
            self.xist_slopes_ = np.zeros(d.shape[1])
        else:
            hh = float(h @ h)
            hz = float(h @ z)
            zz = float(z @ z)
            trace = hh + zz
            delta = np.hypot(hh - zz, 2 * hz)
            largest = 0.5 * (trace + delta)
            smallest = 0.5 * (trace - delta)
            rank_tol = np.finfo(float).eps * max(len(h), 2) * np.sqrt(largest)
            if smallest <= 0 or np.sqrt(smallest) <= rank_tol:
                raise ValueError('Cannot distinguish parental and XIST effects in this training fold: '
                                 'XIST_diff is proportional to the parental direction. '
                                 'More variation in parental-oriented XIST differences is needed.')
            if np.sqrt(largest / smallest) > 1e6:
                warnings.warn('Parental direction and XIST are nearly collinear; adjustment slopes may be unstable.')
            hd = h @ d
            zd = z @ d
            determinant = hh * zz - hz * hz
            if determinant <= np.finfo(float).eps * max(len(h), 2) * hh * zz:
                raise ValueError('Cannot distinguish parental and XIST effects in this training fold: '
                                 'XIST_diff is proportional to the parental direction. '
                                 'More variation in parental-oriented XIST differences is needed.')
            self.parental_effects_ = (zz * hd - hz * zd) / determinant
            self.xist_slopes_ = (hh * zd - hz * hd) / determinant
        adjusted = d - z[:, None] * self.xist_slopes_[None, :]
        self.selector_ = None
        if self.alpha < 1:
            self.selector_ = SelectFdr(f_classif, alpha=self.alpha)
            self.selector_.fit(adjusted, y)
            mask = self.selector_.get_support()
        else:
            mask = np.ones(len(self.cpg_names_), dtype=bool)
        self.selected_mask_ = mask
        self.selected_indices_ = np.flatnonzero(mask)
        self.selected_cpg_names_ = self.cpg_names_[mask]
        self.feature_names_out_ = self.selected_cpg_names_.copy()
        return self

    def transform(self, X):
        check_is_fitted(self, 'feature_names_out_')
        if hasattr(X, 'columns'):
            d = X.loc[:, self.selected_cpg_names_].to_numpy(dtype=float)
            z = X['XIST_diff'].to_numpy(dtype=float)[:, None]
        else:
            X = np.asarray(X, dtype=float)
            d = X[:, self.selected_indices_]
            z = X[:, -1:]
        return d - z * self.xist_slopes_[self.selected_mask_][None, :]

    def get_feature_names_out(self, input_features=None):
        check_is_fitted(self, 'feature_names_out_')
        return self.feature_names_out_.copy()


class XISTElasticNet(ClassifierMixin, BaseEstimator):
    """Convex logistic elastic net with every coefficient penalized.

    mean(logistic loss) + [l1_ratio*|beta|_1 +
        (1-l1_ratio)*||beta||^2/2] / (C*n)

    beta contains XIST-adjusted CpGs; all are penalized.
    No standalone XIST coefficient or intercept is fitted.
    Positive/negative splitting makes the L1 objective smooth with box bounds.
    C follows inverse-penalty scaling of sklearn's original logistic model.
    """

    def __init__(self, C=1.0, l1_ratio=0.5, max_iter=20000, tol=1e-7):
        self.C = C
        self.l1_ratio = l1_ratio
        self.max_iter = max_iter
        self.tol = tol

    def fit(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        if X.ndim != 2 or not np.isfinite(X).all():
            raise ValueError('Model input must be a finite 2D array.')
        if not np.array_equal(np.unique(y), [0, 1]):
            raise ValueError('Each training fold must contain both classes.')
        if self.C <= 0 or not 0 <= self.l1_ratio <= 1:
            raise ValueError('Require C > 0 and 0 <= l1_ratio <= 1.')
        n, p = X.shape
        k = p
        if p == 0:
            # Optional FDR can retain no CpGs: no-intercept baseline is p=0.5.
            self.coef_ = np.empty((1, 0))
            self.intercept_ = np.array([0.0])
            self.classes_ = np.array([0, 1])
            self.n_features_in_ = 0
            self.n_iter_ = np.array([0])
            return self
        lam = 1.0 / (self.C * n)
        l1 = lam * self.l1_ratio
        l2 = lam * (1 - self.l1_ratio)

        def objective(theta):
            if l1 > 0:
                beta = theta[:k] - theta[k:2*k]
            else:
                beta = theta
            z = X @ beta
            loss = np.mean(np.logaddexp(0, z) - y * z)
            loss += 0.5 * l2 * np.dot(beta, beta)
            residual = (special.expit(z) - y) / n
            grad_beta = X.T @ residual + l2 * beta
            if l1 > 0:
                loss += l1 * np.sum(theta[:2*k])
                gradient = np.r_[grad_beta + l1, -grad_beta + l1]
            else:
                gradient = grad_beta
            return loss, gradient

        size = 2*k if l1 > 0 else p
        bounds = [(0, None)] * (2*k) if l1 > 0 else None
        fit = minimize(
            objective, np.zeros(size), method='L-BFGS-B', jac=True,
            bounds=bounds,
            options={'maxiter': self.max_iter, 'gtol': self.tol,
                     'ftol': 1e-12, 'maxls': 50},
        )
        if not fit.success:
            raise RuntimeError(f'Optimizer failed: {fit.message}')
        beta = fit.x[:k] - fit.x[k:2*k] if l1 > 0 else fit.x[:k]
        self.coef_ = beta[None, :]
        self.intercept_ = np.array([0.0])
        self.classes_ = np.array([0, 1])
        self.n_features_in_ = p
        self.n_iter_ = np.array([fit.nit])
        return self

    def predict_proba(self, X):
        check_is_fitted(self, 'coef_')
        prob = special.expit(np.asarray(X, dtype=float) @ self.coef_[0])
        return np.column_stack([1 - prob, prob])

    def predict(self, X):
        return (self.predict_proba(X)[:, 1] >= 0.5).astype(int)


def read_samples(path):
    samples = pd.read_csv(path, header=None, sep='\t', dtype=str)[0]
    return set(samples.str.strip())


def augment_metadata(data):
    """Reverse haplotype labels for the sign-flipped feature rows."""
    flipped = data.copy()
    flipped[['LeftOpHapLabel', 'RightOpHapLabel']] = (
        data[['RightOpHapLabel', 'LeftOpHapLabel']].to_numpy()
    )
    return pd.concat([data, flipped], ignore_index=True)


def augment_features(X):
    """Reverse ALL signed differences."""
    return np.vstack([X, -X])


def fit_estimator(X, y, C, l1_ratio, fdr_alpha, max_iter, cpg_names, label='final'):
    log(f'{label}: fitting XIST adjustment/FDR on {len(y)} rows x {X.shape[1] - 1} CpGs')
    t0 = time.perf_counter()
    adjuster = CpGXISTAdjuster(alpha=fdr_alpha, cpg_names=tuple(cpg_names))
    X_adjusted = adjuster.fit_transform(X, y)
    log(
        f'{label}: adjustment kept {len(adjuster.selected_cpg_names_)}/{len(adjuster.cpg_names_)} '
        f'CpGs in {time.perf_counter() - t0:.2f}s'
    )
    log(
        f'{label}: fitting final elastic net C={C:g}, l1_ratio={l1_ratio:g} '
        f'on {X_adjusted.shape[1]} adjusted CpGs'
    )
    t0 = time.perf_counter()
    model = XISTElasticNet(C=C, l1_ratio=l1_ratio, max_iter=max_iter)
    model.fit(X_adjusted, y)
    log(f'{label}: final elastic net converged in {model.n_iter_[0]} iterations, {time.perf_counter() - t0:.2f}s')
    return Pipeline([('select', adjuster), ('model', model)])


def fit_and_score_model(X_tr, y_tr, X_va, y_va, params, max_iter):
    t0 = time.perf_counter()
    model = XISTElasticNet(
        C=params['C'],
        l1_ratio=params['l1_ratio'],
        max_iter=max_iter,
    )
    model.fit(X_tr, y_tr)
    fit_time = time.perf_counter() - t0

    t0 = time.perf_counter()
    prob = model.predict_proba(X_va)[:, 1]
    score_time = time.perf_counter() - t0
    neg_brier = -brier_score_loss(y_va, prob)
    accuracy = accuracy_score(y_va, prob >= 0.5)
    return (
        params,
        int(model.n_iter_[0]),
        fit_time,
        score_time,
        neg_brier,
        accuracy,
    )


def evaluate_parameter_grid(X, y, groups, C_grid, l1_grid, folds, fdr_alpha, max_iter, n_jobs, cpg_names, label):
    grid = [{'C': C, 'l1_ratio': l1_ratio} for C in C_grid for l1_ratio in l1_grid]
    neg_brier = np.zeros((len(grid), folds))
    accuracy = np.zeros((len(grid), folds))
    fit_time = np.zeros((len(grid), folds))
    score_time = np.zeros((len(grid), folds))
    adjustment_time = np.zeros(folds)
    inner = GroupKFold(n_splits=folds)
    log(
        f'{label}: evaluating {len(grid)} C/l1 candidates with {folds}-fold inner CV '
        f'({n_jobs} job{"s" if n_jobs != 1 else ""})'
    )

    for split, (tr, va) in enumerate(inner.split(X, y, groups)):
        log(
            f'{label}: inner fold {split + 1}/{folds}: fitting XIST adjustment/FDR '
            f'on {len(tr)} train rows, validating {len(va)} rows'
        )
        t0 = time.perf_counter()
        adjuster = CpGXISTAdjuster(alpha=fdr_alpha, cpg_names=tuple(cpg_names))
        X_tr = adjuster.fit_transform(X[tr], y[tr])
        X_va = adjuster.transform(X[va])
        adjustment_time[split] = time.perf_counter() - t0
        log(
            f'{label}: inner fold {split + 1}/{folds}: kept '
            f'{len(adjuster.selected_cpg_names_)}/{len(adjuster.cpg_names_)} CpGs '
            f'after adjustment/FDR in {adjustment_time[split]:.2f}s'
        )

        scores = Parallel(n_jobs=n_jobs, prefer='threads', require='sharedmem')(
            delayed(fit_and_score_model)(
                X_tr, y[tr], X_va, y[va], params, max_iter,
            )
            for params in grid
        )
        for i, (params, n_iter, fit, score, brier, acc) in enumerate(scores):
            fit_time[i, split] = fit
            score_time[i, split] = score
            neg_brier[i, split] = brier
            accuracy[i, split] = acc
            log(
                f'{label}: inner fold {split + 1}/{folds}: '
                f'C={params["C"]:g}, l1={params["l1_ratio"]:g}, '
                f'iter={n_iter}, neg_brier={brier:.6g}, accuracy={acc:.4f}, '
                f'fit={fit:.2f}s, score={score:.2f}s'
            )

    mean_neg_brier = neg_brier.mean(axis=1)
    mean_accuracy = accuracy.mean(axis=1)
    order = np.argsort(-mean_neg_brier, kind='mergesort')
    ranks = np.empty(len(grid), dtype=int)
    ranks[order] = np.arange(1, len(grid) + 1)
    accuracy_order = np.argsort(-mean_accuracy, kind='mergesort')
    accuracy_ranks = np.empty(len(grid), dtype=int)
    accuracy_ranks[accuracy_order] = np.arange(1, len(grid) + 1)

    rows = []
    for i, params in enumerate(grid):
        row = {
            'mean_fit_time': fit_time[i].mean(),
            'std_fit_time': fit_time[i].std(),
            'mean_score_time': score_time[i].mean(),
            'std_score_time': score_time[i].std(),
            'mean_adjustment_time': adjustment_time.mean(),
            'std_adjustment_time': adjustment_time.std(),
            'param_model__C': params['C'],
            'param_model__l1_ratio': params['l1_ratio'],
            'params': {'model__C': params['C'], 'model__l1_ratio': params['l1_ratio']},
            'mean_test_neg_brier_score': mean_neg_brier[i],
            'std_test_neg_brier_score': neg_brier[i].std(),
            'rank_test_neg_brier_score': ranks[i],
            'mean_test_accuracy': mean_accuracy[i],
            'std_test_accuracy': accuracy[i].std(),
            'rank_test_accuracy': accuracy_ranks[i],
        }
        for split in range(folds):
            row[f'split{split}_test_neg_brier_score'] = neg_brier[i, split]
            row[f'split{split}_test_accuracy'] = accuracy[i, split]
        rows.append(row)

    results = pd.DataFrame(rows)
    best = grid[int(order[0])]
    log(
        f'{label}: best C={best["C"]:g}, l1_ratio={best["l1_ratio"]:g}, '
        f'mean_neg_brier={mean_neg_brier[int(order[0])]:.6g}, '
        f'mean_accuracy={mean_accuracy[int(order[0])]:.4f}'
    )
    return best, results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('chromosome')
    parser.add_argument('--xist', required=True, help='TSV: Sample, XIST_diff')
    parser.add_argument('--input', help='CpG difference TSV; original filename by default')
    parser.add_argument('--train-samples', help='One sample ID per line')
    parser.add_argument('--test-samples', help='One sample ID per line')
    parser.add_argument('--output-prefix')
    parser.add_argument('--n-jobs', type=int, default=1)
    parser.add_argument('--folds', type=int, default=5)
    parser.add_argument('--fdr-alpha', type=float, default=1.0,
                        help='Optional heuristic FDR screen on adjusted CpGs; default 1 disables it')
    parser.add_argument('--exclude-bed', help='Optional BED of CpGs to exclude, e.g. XIST promoter')
    parser.add_argument('--cpg-coordinate-base', type=int, choices=[0, 1], default=0,
                        help='Coordinate base of numeric CpG column names for BED exclusion [0]')
    parser.add_argument('--cache-npz', help='Optional preprocessed numeric matrix cache for repeated runs')
    parser.add_argument('--refresh-cache', action='store_true',
                        help='Rebuild --cache-npz even if a matching cache already exists')
    parser.add_argument('--read-engine', choices=['c', 'pyarrow'], default='c',
                        help='pandas TSV parser engine for the wide CpG and XIST TSVs [c]')
    parser.add_argument('--fast-final', action='store_true',
                        help='Skip outer nested CV and fit only the final tuned model')
    parser.add_argument('--max-iter', type=int, default=20000)
    parser.add_argument('--C-grid', nargs='+', type=float,
                        default=[0.01, 0.1, 1, 10, 100, 1000])
    parser.add_argument('--l1-grid', nargs='+', type=float,
                        default=[0, 0.25, 0.5, 0.75, 1])
    parser.add_argument('--swap-samples', nargs='*',
                        default=['HG02451', 'HG02717', 'HG03579'],
                        help='Sign-correct mislabeled input rows, as in original script')
    args = parser.parse_args()
    if args.folds < 2 or not 0 < args.fdr_alpha <= 1:
        parser.error('Require folds >= 2 and 0 < fdr-alpha <= 1.')

    total_start = time.perf_counter()
    chrom = args.chromosome
    #root = Path('/path/to/pooch/model_validation/model_juhyun/test_data')
    default_train = ('balance_n46.training.sample.list'
                     if chrom == 'chrX' else root / 'trainSet.samplelist')
    default_test = ('balance_n46.test.sample.list'
                    if chrom == 'chrX' else root / 'testSet.samplelist')
    train_path = args.train_samples or default_train
    test_path = args.test_samples or default_test
    train_ids = read_samples(train_path)
    test_ids = read_samples(test_path)
    if train_ids & test_ids:
        raise ValueError(f'Train/test samples overlap: {sorted(train_ids & test_ids)}')
    input_path = args.input or f'{chrom}.missingrate_5.createDiffDatasetFromMethylDataset.out'
    metadata = ['Sample', 'LeftOpHapLabel', 'RightOpHapLabel']
    key = cache_key(args, chrom, input_path, train_path, test_path)
    cached = None if args.refresh_cache else load_preprocessed_cache(args.cache_npz, key, metadata)

    if cached:
        data, cpg_columns, feature_values = cached
    else:
        log(f'Reading CpG differences from {input_path}')
        t0 = time.perf_counter()
        data = read_tsv(
            input_path, engine=args.read_engine,
            dtype={'Sample': str, 'LeftOpHapLabel': str, 'RightOpHapLabel': str},
        )
        log(f'Read {len(data)} rows x {data.shape[1]} columns in {time.perf_counter() - t0:.2f}s')
        if list(data.columns[:3]) != metadata:
            raise ValueError(f'First three columns must be {metadata}')
        if 'XIST_diff' in data.columns:
            raise ValueError('Input already has XIST_diff; remove it before merging --xist.')
        cpg_columns = list(data.columns[3:])
        if args.exclude_bed:
            log(f'Applying BED exclusion from {args.exclude_bed} to {len(cpg_columns)} CpGs')
            t0 = time.perf_counter()
            # BED intervals are always zero-based, end-exclusive. Only this chromosome
            # is used; numeric input column names are interpreted using the base flag.
            regions = pd.read_csv(args.exclude_bed, sep=r'\s+', header=None,
                                  usecols=[0, 1, 2], comment='#')
            regions = regions.loc[regions[0] == chrom]
            positions = pd.to_numeric(pd.Index(cpg_columns), errors='raise').to_numpy()
            positions = positions - args.cpg_coordinate_base
            excluded = np.zeros(len(cpg_columns), dtype=bool)
            for _, start, end in regions.itertuples(index=False, name=None):
                if int(start) < 0 or int(end) <= int(start):
                    raise ValueError('Invalid BED interval.')
                excluded |= (positions >= int(start)) & (positions < int(end))
            drop_columns = [c for c, drop in zip(cpg_columns, excluded) if drop]
            data = data.drop(columns=drop_columns)
            cpg_columns = [c for c, drop in zip(cpg_columns, excluded) if not drop]
            log(
                f'Excluded {len(drop_columns)} CpGs using {args.exclude_bed}; '
                f'{len(cpg_columns)} remain ({time.perf_counter() - t0:.2f}s)'
            )
        if not cpg_columns:
            raise ValueError('No CpG columns in input.')
        # Restrict validation to samples actually used, with the original split.
        data = data.loc[data.Sample.isin(train_ids | test_ids)].copy()
        log(f'Retained {len(data)} rows after train/test sample filtering')
        if data.Sample.duplicated().any():
            raise ValueError('Training input must contain one row per sample, before augmentation.')
        cov = read_tsv(args.xist, engine=args.read_engine, usecols=['Sample', 'XIST_diff'])
        data = data.merge(cov, on='Sample', how='left', validate='many_to_one', sort=False)
        features = cpg_columns + ['XIST_diff']
        log(f'Materializing {len(features)} numeric model input columns')
        t0 = time.perf_counter()
        feature_values = data.loc[:, features].to_numpy(dtype=float, copy=True)
        log(f'Model input matrix materialized in {time.perf_counter() - t0:.2f}s')
        data = data.loc[:, metadata].copy()
        save_preprocessed_cache(args.cache_npz, key, data, cpg_columns, feature_values, metadata)

    features = cpg_columns + ['XIST_diff']
    finite = np.isfinite(feature_values).all(axis=1)
    if not finite.all():
        raise ValueError(f'Missing/nonfinite CpG or XIST values: {data.loc[~finite, "Sample"].unique().tolist()}')
    for left, right in zip(data.LeftOpHapLabel, data.RightOpHapLabel):
        if {left, right} != {'Maternal', 'Paternal'}:
            raise ValueError('Each row must have one Maternal and one Paternal label.')
    swap = data.Sample.isin(args.swap_samples)
    feature_values[swap.to_numpy()] *= -1
    # Original sample correction changes signs, not labels.
    train_rows = data.Sample.isin(train_ids).to_numpy()
    test_rows = data.Sample.isin(test_ids).to_numpy()
    train = augment_metadata(data.loc[train_rows, metadata].reset_index(drop=True))
    test = augment_metadata(data.loc[test_rows, metadata].reset_index(drop=True))
    X = augment_features(feature_values[train_rows])
    Xt = augment_features(feature_values[test_rows])
    if train.empty or test.empty:
        raise ValueError('Training or test set is empty; check sample lists.')
    y = (train.RightOpHapLabel == 'Paternal').astype(int).to_numpy()
    groups = train.Sample.to_numpy()
    yt = (test.RightOpHapLabel == 'Paternal').astype(int).to_numpy()
    log(
        f'Prepared {len(train)} augmented training rows from {len(np.unique(groups))} samples '
        f'and {len(test)} augmented test rows with {len(cpg_columns)} CpGs'
    )
    if len(np.unique(groups)) < args.folds:
        raise ValueError('Not enough distinct training samples for outer CV.')
    outer = GroupKFold(n_splits=args.folds)
    splits = list(outer.split(X, y, groups))
    if any(len(np.unique(groups[idx])) < args.folds for idx, _ in splits):
        raise ValueError('Not enough samples for inner CV; lower --folds.')

    cv_prob = None
    if args.fast_final:
        log('Fast-final mode: skipping outer nested CV')
    else:
        cv_prob = np.zeros(len(y))
        for fold, (tr, va) in enumerate(splits, 1):
            log(
                f'Outer fold {fold}/{args.folds}: '
                f'{len(tr)} train rows / {len(np.unique(groups[tr]))} samples, '
                f'{len(va)} validation rows / {len(np.unique(groups[va]))} samples'
            )
            best_params, _ = evaluate_parameter_grid(
                X[tr], y[tr], groups[tr], args.C_grid, args.l1_grid,
                args.folds, args.fdr_alpha, args.max_iter, args.n_jobs, cpg_columns,
                f'outer {fold}/{args.folds}',
            )
            best_fold = fit_estimator(
                X[tr], y[tr], best_params['C'], best_params['l1_ratio'],
                args.fdr_alpha, args.max_iter, cpg_columns, f'outer {fold}/{args.folds} refit',
            )
            cv_prob[va] = best_fold.predict_proba(X[va])[:, 1]
            log(f'Outer fold {fold}/{args.folds}: best {{"model__C": {best_params["C"]}, "model__l1_ratio": {best_params["l1_ratio"]}}}')

    log('Fitting final model on all training rows')
    best_params, gridsearch = evaluate_parameter_grid(
        X, y, groups, args.C_grid, args.l1_grid, args.folds,
        args.fdr_alpha, args.max_iter, args.n_jobs, cpg_columns, 'final grid',
    )
    best = fit_estimator(
        X, y, best_params['C'], best_params['l1_ratio'],
        args.fdr_alpha, args.max_iter, cpg_columns, 'final refit',
    )
    log('Scoring final model on test rows')
    names = best.named_steps['select'].get_feature_names_out()
    coefficients = best.named_steps['model'].coef_[0]
    coef_table = pd.DataFrame({'feature': names, 'coefficient': coefficients})
    adjuster = best.named_steps['select']
    selected_cpg_names = adjuster.selected_cpg_names_
    k = len(selected_cpg_names)
    coef_table['feature_type'] = 'XIST_adjusted_CpG'
    coef_table['CpG'] = selected_cpg_names
    coef_table['penalized'] = True
    nonzero = coef_table.loc[coef_table.coefficient.abs() > 1e-7].copy()
    nonzero = nonzero.iloc[np.argsort(-nonzero.coefficient.abs().to_numpy())]
    all_betas = np.zeros(len(adjuster.cpg_names_))
    all_betas[adjuster.selected_mask_] = coefficients
    cpg_summary = pd.DataFrame({
        'CpG': adjuster.cpg_names_,
        'parental_effect_a': adjuster.parental_effects_,
        'XIST_slope_b': adjuster.xist_slopes_,
        'passed_screen': adjuster.selected_mask_,
        'elasticnet_beta': all_betas,
    })
    active = cpg_summary.elasticnet_beta.abs() > 1e-7
    selected_cpgs = cpg_summary.loc[active]
    induced_xist_coef = -float(np.dot(all_betas, adjuster.xist_slopes_))
    prob = best.predict_proba(Xt)[:, 1]
    pred = (prob >= 0.5).astype(int)
    metrics = {
        'chromosome': chrom,
        'original_CpG': len(cpg_columns),
        'total_CpG': k,
        'selected_CpG': len(selected_cpgs),
        'total_penalized_terms': k,
        'nonzero_penalized_terms': len(nonzero),
        'adjusted_CpG_FDR_screen_enabled': args.fdr_alpha < 1,
        'XIST_free_standalone_term': False,
        'induced_raw_XIST_coefficient': induced_xist_coef,
        'adjustment': 'd_j - b_j*XIST_diff; slopes fitted with parental direction in each training fold',
        'all_model_terms_penalized': True,
        'FDR_is_inferential': False,
        'outer_nested_cv_enabled': not args.fast_final,
        'cv_auc': None if cv_prob is None else roc_auc_score(y, cv_prob),
        'cv_brier': None if cv_prob is None else brier_score_loss(y, cv_prob),
        'best_c': best.named_steps['model'].C,
        'best_l1_ratio': best.named_steps['model'].l1_ratio,
        'Test_AUC': roc_auc_score(yt, prob),
        'Test_Accuracy': accuracy_score(yt, pred),
        'Test_Balanced_Accuracy': balanced_accuracy_score(yt, pred),
        'Test_Brier_Score': brier_score_loss(yt, prob),
    }
    prefix = args.output_prefix or f'elasticnet_xist_adjusted_{chrom}'
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    log(f'Writing outputs to {prefix}.*')
    coef_table.to_csv(f'{prefix}.coefficients.tsv', sep='\t', index=False)
    nonzero.to_csv(f'{prefix}.selected_terms.tsv', sep='\t', index=False)
    selected_cpgs.to_csv(f'{prefix}.selected_cpgs.tsv', sep='\t', index=False)
    cpg_summary.to_csv(f'{prefix}.adjustment_slopes.tsv', sep='\t', index=False)
    # An audit of the equivalent linear score on RAW inputs; not a separate fit.
    pd.DataFrame({'feature': adjuster.cpg_names_.tolist() + ['XIST_diff'],
                  'effective_coefficient': np.r_[all_betas, induced_xist_coef]}).to_csv(
        f'{prefix}.effective_raw_coefficients.tsv', sep='\t', index=False)
    pd.DataFrame([metrics]).to_csv(f'{prefix}.metrics.csv', index=False)
    gridsearch.to_csv(f'{prefix}.gridsearch.csv', index=False)
    prediction_outputs = [('test', test, yt, prob)]
    if cv_prob is not None:
        prediction_outputs.insert(0, ('nested_cv', train, y, cv_prob))
    for label, frame, truth, probabilities in prediction_outputs:
        predictions = frame[metadata].copy()
        predictions['y_true'] = truth
        predictions['p_right_paternal'] = probabilities
        predictions.to_csv(f'{prefix}.{label}.predictions.tsv', sep='\t', index=False)
    # cloudpickle preserves these custom classes when this file runs as __main__.
    with open(f'{prefix}.model.pkl', 'wb') as handle:
        cloudpickle.dump({'pipeline': best, 'input_features': features,
                          'metrics': metrics, 'arguments': vars(args)}, handle)
    log(json.dumps(metrics, indent=2))
    log('Test confusion matrix:\n' + str(confusion_matrix(yt, pred)))
    log(classification_report(yt, pred, zero_division=0))
    log(f'Saved outputs: {prefix}.*')
    log(f'Total runtime: {time.perf_counter() - total_start:.2f}s')


if __name__ == '__main__':
    main()
