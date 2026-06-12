#!/usr/bin/env python3
"""
Sensitivity analysis #2 for the cluster-stratified fMRI analyses:
keep the original sample (with mean_fd < 1 mm) and fit OLS

    connectivity_ij ~ group + mean_fd + minutes_quality_data

per feature for:
  - Cluster Ck (Autism) vs pooled NT  (group_dummy = 1 for Ck)
  - Pairwise Ck vs Cl within Autism   (group_dummy = 1 for Ck)

Reports t-stat, beta, FDR-corrected p-value of the group dummy.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_THIS_DIR = Path(__file__).parent
_SHARED_UTILS_DIR = _THIS_DIR.parents[1] / '6_functional_analysis' / 'scripts'
sys.path.insert(0, str(_SHARED_UTILS_DIR))
from _sensitivity_utils import (  # noqa: E402
    apply_sensitivity_filters,
    load_connectivity_and_clusters,
    ols_per_feature,
)

_SECTION_DIR = _THIS_DIR.parent
OUT_DIR = _SECTION_DIR / 'outputs' / 'figures' / 'sensitivity'
R_INPUT_DIR = OUT_DIR / 'r_input_files'
OUT_DIR.mkdir(parents=True, exist_ok=True)
R_INPUT_DIR.mkdir(parents=True, exist_ok=True)

MAX_MEAN_FD = 1.0
COVARIATES = ['mean_fd', 'minutes_quality_data']
SUFFIX = 'covadj'

CLUSTERS = ['C1', 'C2', 'C3']
PAIRWISE_PAIRS = [('C1', 'C2'), ('C1', 'C3'), ('C2', 'C3')]


def _tag_group(df: pd.DataFrame, autism_mask: pd.Series, nt_mask: pd.Series, autism_label: str, nt_label: str) -> pd.DataFrame:
    df = df.copy()
    df['_group'] = np.where(autism_mask, autism_label,
                            np.where(nt_mask, nt_label, np.nan))
    return df[df['_group'].notna()].copy()


def run_cluster_vs_nt(df: pd.DataFrame, conn_cols: list[str]) -> dict:
    nt_mask = df['PopulationS1'] == 'NT'
    print(f'\nPooled NT reference: {nt_mask.sum()}')
    sig_counts = {}
    for cluster in CLUSTERS:
        autism_mask = (df['Cluster'] == cluster) & (df['PopulationS1'] == 'Autism')
        n_aut = int(autism_mask.sum())
        print(f'  Autism {cluster}: {n_aut}')
        if n_aut < 5:
            continue
        sub = _tag_group(df, autism_mask, nt_mask, autism_label=cluster, nt_label='NT')
        results = ols_per_feature(
            sub, conn_cols,
            group_col='_group', group_a=cluster, group_b='NT',
            covariates=COVARIATES,
        )
        if len(results) == 0:
            continue
        results = results.rename(columns={
            f'n_{cluster}': 'autism_n',
            'n_NT': 'td_n',
            'p_value': 'p_val',
            'p_value_fdr': 'p_fdr',
        })
        n_sig = int((results['p_fdr'] < 0.05).sum())
        sig_counts[cluster] = n_sig
        print(f'  Significant (FDR<0.05, group effect): {n_sig} / {len(results)}')

        results.to_csv(OUT_DIR / f'cluster_{cluster}_connectivity_autism_vs_td_{SUFFIX}.csv', index=False)

        sig_set = results[results['p_fdr'] < 0.05][['feature', 't_stat', 'beta', 'p_fdr']]
        sig_set.to_csv(OUT_DIR / f'cluster_{cluster}_vs_NT_significant_results_{SUFFIX}.csv', index=False)

        r_ready = results[['feature', 't_stat', 'beta', 'p_fdr']].rename(
            columns={'feature': 'region', 't_stat': 'mean_value'}
        )
        r_ready.to_csv(R_INPUT_DIR / f'cluster_{cluster}_connectivity_for_r_{SUFFIX}.csv', index=False)
    return sig_counts


def run_pairwise(df: pd.DataFrame, conn_cols: list[str]) -> dict:
    autism = df[df['PopulationS1'] == 'Autism']
    print(f'\nAutism-only sample for pairwise: {len(autism)}')
    sig_counts = {}
    for c1, c2 in PAIRWISE_PAIRS:
        sub = autism[autism['Cluster'].isin([c1, c2])].copy()
        n1 = int((sub['Cluster'] == c1).sum())
        n2 = int((sub['Cluster'] == c2).sum())
        print(f'  {c1} ({n1}) vs {c2} ({n2})')
        if n1 < 5 or n2 < 5:
            continue
        results = ols_per_feature(
            sub, conn_cols,
            group_col='Cluster', group_a=c1, group_b=c2,
            covariates=COVARIATES,
        )
        if len(results) == 0:
            continue
        results = results.rename(columns={'p_value': 'p_val', 'p_value_fdr': 'p_fdr'})
        n_sig = int((results['p_fdr'] < 0.05).sum())
        sig_counts[f'{c1}_vs_{c2}'] = n_sig
        print(f'  Significant (FDR<0.05, group effect): {n_sig} / {len(results)}')
        results.to_csv(OUT_DIR / f'pairwise_{c1}_vs_{c2}_connectivity_{SUFFIX}.csv', index=False)
    return sig_counts


def main() -> int:
    print('=' * 60)
    print('CLUSTER fMRI — COVARIATE-ADJUSTED (FD + minutes of quality data)')
    print('=' * 60)

    df, conn_cols = load_connectivity_and_clusters()
    print(f'Initial merged sample: {len(df)}')
    print(f'Connectivity features: {len(conn_cols)}')

    df_filt, counts = apply_sensitivity_filters(df, min_minutes=None, max_mean_fd=MAX_MEAN_FD)
    print('\nFilter trace:', counts)
    before = len(df_filt)
    df_filt = df_filt.dropna(subset=COVARIATES)
    print(f'Dropped for missing covariates: {before - len(df_filt)}')
    print('PopulationS1 counts:', df_filt['PopulationS1'].value_counts().to_dict())
    autism_only = df_filt[df_filt['PopulationS1'] == 'Autism']
    print('Autism cluster counts:', autism_only['Cluster'].value_counts().to_dict())

    cluster_sigs = run_cluster_vs_nt(df_filt, conn_cols)
    pairwise_sigs = run_pairwise(df_filt, conn_cols)

    summary_rows = []
    for cluster, n_sig in cluster_sigs.items():
        n_aut = int(((df_filt['Cluster'] == cluster) & (df_filt['PopulationS1'] == 'Autism')).sum())
        n_nt = int((df_filt['PopulationS1'] == 'NT').sum())
        summary_rows.append({
            'analysis': f'cluster_{cluster}_vs_NT_{SUFFIX}',
            'n_autism': n_aut, 'n_nt': n_nt,
            'n_significant_fdr05': n_sig,
            'max_mean_fd': MAX_MEAN_FD,
            'covariates': '+'.join(COVARIATES),
        })
    for pair, n_sig in pairwise_sigs.items():
        c1, c2 = pair.split('_vs_')
        n1 = int(((df_filt['Cluster'] == c1) & (df_filt['PopulationS1'] == 'Autism')).sum())
        n2 = int(((df_filt['Cluster'] == c2) & (df_filt['PopulationS1'] == 'Autism')).sum())
        summary_rows.append({
            'analysis': f'pairwise_{pair}_{SUFFIX}',
            'n_autism': n1, 'n_nt': n2,
            'n_significant_fdr05': n_sig,
            'max_mean_fd': MAX_MEAN_FD,
            'covariates': '+'.join(COVARIATES),
        })
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            OUT_DIR / f'cluster_fmri_sensitivity_{SUFFIX}_summary.csv', index=False
        )
    print(f'\nOutputs in {OUT_DIR}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
