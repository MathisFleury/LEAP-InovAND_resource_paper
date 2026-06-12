#!/usr/bin/env python3
"""
Sensitivity analysis #1 for the cluster-stratified fMRI analyses:
restrict to participants with >= 6 min of post-censoring quality data
AND mean_fd < 1 mm.

Mirrors the primary 01_generate_cluster_fmri_inputs.py:
  - Cluster C{k} (Autism) vs pooled NT — Welch's t + Cohen's d, FDR
  - Pairwise (C1 vs C2, C1 vs C3, C2 vs C3) — Welch's t + Cohen's d, FDR

Outputs live under outputs/figures/sensitivity/ for this section, with
the *_6min suffix.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# Pull the shared QC helpers from section 6 to avoid duplication.
_THIS_DIR = Path(__file__).parent
_SHARED_UTILS_DIR = _THIS_DIR.parents[1] / '6_functional_analysis' / 'scripts'
sys.path.insert(0, str(_SHARED_UTILS_DIR))
from _sensitivity_utils import (  # noqa: E402
    apply_sensitivity_filters,
    cohens_d,
    load_connectivity_and_clusters,
)

_SECTION_DIR = _THIS_DIR.parent
OUT_DIR = _SECTION_DIR / 'outputs' / 'figures' / 'sensitivity'
R_INPUT_DIR = OUT_DIR / 'r_input_files'
OUT_DIR.mkdir(parents=True, exist_ok=True)
R_INPUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_MINUTES = 6.0
MAX_MEAN_FD = 1.0
SUFFIX = '6min'

CLUSTERS = ['C1', 'C2', 'C3']
PAIRWISE_PAIRS = [('C1', 'C2'), ('C1', 'C3'), ('C2', 'C3')]


def _cohens_d_ci(d: float, na: int, nb: int) -> tuple[float, float]:
    if not np.isfinite(d) or na < 2 or nb < 2:
        return (np.nan, np.nan)
    se = np.sqrt((na + nb) / (na * nb) + d ** 2 / (2.0 * (na + nb)))
    z = 1.959963984540054
    return (float(d - z * se), float(d + z * se))


def _fdr(p: pd.Series) -> np.ndarray:
    mask = ~p.isna()
    out = np.full(len(p), np.nan)
    if mask.sum() > 0:
        out[mask] = multipletests(p[mask], method='fdr_bh')[1]
    return out


def run_cluster_vs_nt(df: pd.DataFrame, conn_cols: list[str]) -> dict:
    nt = df[df['PopulationS1'] == 'NT']
    print(f'\nPooled NT reference: {len(nt)}')
    sig_counts = {}

    for cluster in CLUSTERS:
        autism = df[(df['Cluster'] == cluster) & (df['PopulationS1'] == 'Autism')]
        print(f'  Autism {cluster}: {len(autism)}')
        if len(autism) < 5:
            print(f'  Skipping {cluster} (< 5 participants)')
            continue

        rows = []
        for feat in conn_cols:
            a = autism[feat].dropna().values
            b = nt[feat].dropna().values
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                t_stat, p_val = ttest_ind(a, b, equal_var=False)
            except Exception:
                continue
            d = cohens_d(a, b)
            d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
            rows.append({
                'region': feat, 't_stat': t_stat, 'p_val': p_val,
                'cohens_d': d, 'cohens_d_lo': d_lo, 'cohens_d_hi': d_hi,
                'autism_mean': float(np.mean(a)), 'td_mean': float(np.mean(b)),
                'autism_n': int(len(a)), 'td_n': int(len(b)),
            })
        if not rows:
            continue
        df_r = pd.DataFrame(rows)
        df_r['p_fdr'] = _fdr(df_r['p_val'])
        df_r = df_r[['region', 't_stat', 'p_val', 'cohens_d',
                     'cohens_d_lo', 'cohens_d_hi', 'p_fdr',
                     'autism_mean', 'td_mean', 'autism_n', 'td_n']]
        n_sig = int((df_r['p_fdr'] < 0.05).sum())
        sig_counts[cluster] = n_sig
        print(f'  Significant (FDR < 0.05): {n_sig} / {len(df_r)}')

        df_r.to_csv(OUT_DIR / f'cluster_{cluster}_connectivity_autism_vs_td_{SUFFIX}.csv', index=False)

        sig_set = df_r[df_r['p_fdr'] < 0.05][
            ['region', 't_stat', 'cohens_d', 'cohens_d_lo', 'cohens_d_hi', 'p_fdr']
        ].rename(columns={'region': 'feature'})
        sig_set.to_csv(OUT_DIR / f'cluster_{cluster}_vs_NT_significant_results_{SUFFIX}.csv', index=False)

        df_r_ready = df_r[['region', 't_stat', 'cohens_d', 'cohens_d_lo', 'cohens_d_hi', 'p_fdr']
                          ].rename(columns={'t_stat': 'mean_value'})
        df_r_ready.to_csv(R_INPUT_DIR / f'cluster_{cluster}_connectivity_for_r_{SUFFIX}.csv', index=False)

    return sig_counts


def run_pairwise(df: pd.DataFrame, conn_cols: list[str]) -> dict:
    autism = df[df['PopulationS1'] == 'Autism']
    print(f'\nAutism-only sample for pairwise: {len(autism)}')
    sig_counts = {}
    for c1, c2 in PAIRWISE_PAIRS:
        g1 = autism[autism['Cluster'] == c1]
        g2 = autism[autism['Cluster'] == c2]
        print(f'  {c1} ({len(g1)}) vs {c2} ({len(g2)})')
        if len(g1) < 5 or len(g2) < 5:
            print(f'  Skipping {c1}/{c2}')
            continue
        rows = []
        for feat in conn_cols:
            a = g1[feat].dropna().values
            b = g2[feat].dropna().values
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                t_stat, p_val = ttest_ind(a, b, equal_var=False)
            except Exception:
                continue
            d = cohens_d(a, b)
            d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
            rows.append({
                'region': feat, 't_stat': t_stat, 'p_val': p_val,
                'cohens_d': d, 'cohens_d_lo': d_lo, 'cohens_d_hi': d_hi,
                f'{c1}_mean': float(np.mean(a)), f'{c2}_mean': float(np.mean(b)),
                f'{c1}_n': int(len(a)), f'{c2}_n': int(len(b)),
            })
        if not rows:
            continue
        df_pw = pd.DataFrame(rows)
        df_pw['p_fdr'] = _fdr(df_pw['p_val'])
        n_sig = int((df_pw['p_fdr'] < 0.05).sum())
        sig_counts[f'{c1}_vs_{c2}'] = n_sig
        print(f'  Significant (FDR < 0.05): {n_sig} / {len(df_pw)}')
        df_pw.to_csv(OUT_DIR / f'pairwise_{c1}_vs_{c2}_connectivity_{SUFFIX}.csv', index=False)

    return sig_counts


def main() -> int:
    print('=' * 60)
    print('CLUSTER fMRI — SENSITIVITY (6-min cutoff + FD < 1 mm)')
    print('=' * 60)

    df, conn_cols = load_connectivity_and_clusters()
    print(f'Initial merged sample: {len(df)}')
    print(f'Connectivity features: {len(conn_cols)}')

    df_filt, counts = apply_sensitivity_filters(df, min_minutes=MIN_MINUTES, max_mean_fd=MAX_MEAN_FD)
    print('\nFilter trace:', counts)
    print('PopulationS1 counts:', df_filt['PopulationS1'].value_counts().to_dict())
    autism_only = df_filt[df_filt['PopulationS1'] == 'Autism']
    print('Autism cluster counts:', autism_only['Cluster'].value_counts().to_dict())

    cluster_sigs = run_cluster_vs_nt(df_filt, conn_cols)
    pairwise_sigs = run_pairwise(df_filt, conn_cols)

    summary_rows = []
    for cluster, n_sig in cluster_sigs.items():
        n_autism = int(((df_filt['Cluster'] == cluster) & (df_filt['PopulationS1'] == 'Autism')).sum())
        n_nt = int((df_filt['PopulationS1'] == 'NT').sum())
        summary_rows.append({
            'analysis': f'cluster_{cluster}_vs_NT_{SUFFIX}',
            'n_autism': n_autism, 'n_nt': n_nt,
            'n_significant_fdr05': n_sig,
            'min_minutes_quality_data': MIN_MINUTES, 'max_mean_fd': MAX_MEAN_FD,
        })
    for pair, n_sig in pairwise_sigs.items():
        c1, c2 = pair.split('_vs_')
        n1 = int(((df_filt['Cluster'] == c1) & (df_filt['PopulationS1'] == 'Autism')).sum())
        n2 = int(((df_filt['Cluster'] == c2) & (df_filt['PopulationS1'] == 'Autism')).sum())
        summary_rows.append({
            'analysis': f'pairwise_{pair}_{SUFFIX}',
            'n_autism': n1, 'n_nt': n2,
            'n_significant_fdr05': n_sig,
            'min_minutes_quality_data': MIN_MINUTES, 'max_mean_fd': MAX_MEAN_FD,
        })

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            OUT_DIR / f'cluster_fmri_sensitivity_{SUFFIX}_summary.csv', index=False
        )
    print(f'\nOutputs in {OUT_DIR}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
