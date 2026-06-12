#!/usr/bin/env python3
"""
Sensitivity analysis #2 for Autism vs NT connectivity:
keep the original sample (with mean_fd < 1 mm only) and fit OLS

    connectivity_ij ~ group + mean_fd + minutes_quality_data

per feature. We report the t-statistic and FDR-corrected p-value on the
group dummy (group=1 for Autism, 0 for NT).

This addresses the reviewer comment that mean FD and total amount of
quality fMRI data should be included as covariates.

All outputs land under outputs/figures/sensitivity/ and outputs/tables/sensitivity/.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(_SCRIPT_DIR))
from _sensitivity_utils import (  # noqa: E402
    ATLAS_FILE,
    add_source_target_cols,
    apply_sensitivity_filters,
    load_connectivity_and_clusters,
    ols_per_feature,
    split_hyper_hypo,
)

_SECTION_DIR = _SCRIPT_DIR.parent
OUT_FIG = _SECTION_DIR / 'outputs' / 'figures' / 'sensitivity'
OUT_TAB = _SECTION_DIR / 'outputs' / 'tables' / 'sensitivity'
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

MAX_MEAN_FD = 1.0
COVARIATES = ['mean_fd', 'minutes_quality_data']
SUFFIX = 'covadj'


def main() -> int:
    print('=' * 60)
    print('AUTISM vs NT — COVARIATE-ADJUSTED (FD + minutes of quality data)')
    print('=' * 60)

    df, conn_cols = load_connectivity_and_clusters()
    print(f'Initial merged sample: {len(df)}')
    print(f'Connectivity features: {len(conn_cols)}')

    df = df[df['PopulationS1'].isin(['Autism', 'NT'])].copy()
    print(f'After restricting to Autism/NT: {len(df)} '
          f'({df["PopulationS1"].value_counts().to_dict()})')

    df_filt, counts = apply_sensitivity_filters(df, min_minutes=None, max_mean_fd=MAX_MEAN_FD)
    print('\nSensitivity filter trace:')
    for k, v in counts.items():
        print(f'  {k}: {v}')

    # Subjects also need non-NaN values on the covariates, drop here:
    before = len(df_filt)
    df_filt = df_filt.dropna(subset=COVARIATES)
    print(f'Dropped for missing covariates: {before - len(df_filt)}')
    print('Final group sizes:')
    print(' ', df_filt['PopulationS1'].value_counts().to_dict())

    if df_filt['PopulationS1'].value_counts().min() < 5:
        print('ERROR: at least one group has < 5 participants after filtering.')
        return 1

    df_t = ols_per_feature(
        df_filt, conn_cols,
        group_col='PopulationS1', group_a='Autism', group_b='NT',
        covariates=COVARIATES,
    )
    if len(df_t) == 0:
        print('ERROR: no valid OLS results.')
        return 1

    results = split_hyper_hypo(df_t)
    df_sig = results['significant_results']
    if results['fallback_used']:
        print('Note: no FDR-survivors; reporting uncorrected p<0.05.')
    print(f'Significant connections (group effect): {len(df_sig)} '
          f'(hyper={len(results["hyper_results"])}, hypo={len(results["hypo_results"])})')

    # Save tables
    df_t.to_csv(OUT_TAB / f'autism_vs_td_{SUFFIX}_all_results.csv', index=False)
    df_sig.to_csv(OUT_TAB / f'autism_vs_td_{SUFFIX}_significant_results.csv', index=False)

    # Hyper / hypo full data + R-ready counts
    for df_sub, label in [
        (results['hyper_results'], f'autism_vs_td_hyperconnectivity_{SUFFIX}'),
        (results['hypo_results'], f'autism_vs_td_hypoconnectivity_{SUFFIX}'),
    ]:
        if len(df_sub) == 0:
            continue
        df_sub_st = add_source_target_cols(df_sub)
        df_sub_st[['source', 'target', 't_stat', 'beta', 'p_value', 'p_value_fdr']].to_csv(
            OUT_FIG / f'{label}_full_data.csv', index=False
        )
        counts_src = df_sub_st['source'].value_counts().reset_index()
        if 'count' in counts_src.columns:
            counts_src = counts_src.rename(columns={'count': 't_stat'})
        elif 'source' in counts_src.columns:
            counts_src = counts_src.rename(columns={'source': 't_stat', 'index': 'region'})
        else:
            counts_src.columns = ['region', 't_stat']
        if 'region' not in counts_src.columns:
            counts_src = counts_src.rename(columns={counts_src.columns[0]: 'region'})
        counts_src['mean_value'] = counts_src['t_stat']
        counts_src.to_csv(OUT_FIG / f'{label}_for_r.csv', index=False)

    # Network breakdown
    if ATLAS_FILE.exists():
        df_atlas = pd.read_csv(ATLAS_FILE, sep='\t')
        for df_sub, suffix in [
            (results['significant_results'], 'network_analysis'),
            (results['hyper_results'], 'hyperconnectivity_network_analysis'),
            (results['hypo_results'], 'hypoconnectivity_network_analysis'),
        ]:
            if len(df_sub) == 0:
                continue
            df_sub_st = add_source_target_cols(df_sub)
            merged = df_sub_st.merge(
                df_atlas[['label', 'network_label']],
                left_on='source', right_on='label', how='left',
            )
            net_counts = (
                merged.groupby('network_label').size().reset_index(name='count')
                .sort_values('count', ascending=False)
            )
            net_counts.to_csv(OUT_FIG / f'autism_vs_td_{suffix}_{SUFFIX}.csv', index=False)
    else:
        print(f'Atlas file not found, skipping network breakdown: {ATLAS_FILE}')

    summary = pd.DataFrame([{
        'analysis': f'autism_vs_td_{SUFFIX}',
        'max_mean_fd': MAX_MEAN_FD,
        'covariates': '+'.join(COVARIATES),
        'n_autism': int((df_filt['PopulationS1'] == 'Autism').sum()),
        'n_nt': int((df_filt['PopulationS1'] == 'NT').sum()),
        'total_significant': len(df_sig),
        'hyperconnectivity': len(results['hyper_results']),
        'hypoconnectivity': len(results['hypo_results']),
        'fdr_threshold_used': not results['fallback_used'],
    }])
    summary.to_csv(OUT_TAB / f'autism_vs_td_{SUFFIX}_summary.csv', index=False)
    print(f'\nWrote outputs to:\n  {OUT_FIG}\n  {OUT_TAB}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
