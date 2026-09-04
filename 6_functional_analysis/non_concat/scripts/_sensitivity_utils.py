"""
Shared helpers for fMRI connectivity sensitivity analyses.

Centralises the QC merge, the connectivity-column filter, the per-feature
Welch t-test, and the per-feature OLS model with FD + minutes-of-quality-data
covariates. Used by both Autism-vs-NT (section 6) and cluster-level (section 7)
sensitivity scripts.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

EEG_MRI_RESULTS = Path('/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results')
DF_CLUSTERS_FILE = EEG_MRI_RESULTS / 'dataset_paper' / 'dataframes' / 'df_clusters_complete_kmeans.csv'
# Default to this section's harmonised connectivity; the concatenated-runs
# pipeline (6_functional_analysis/concat) overrides via the FMRI_CONN_FILE env.
FMRI_CONN_FILE = Path(os.environ.get(
    'FMRI_CONN_FILE',
    Path(__file__).resolve().parents[1] / 'preprocessing' / 'outputs' / 'df_conn_cohort_norm.csv'))
ATLAS_FILE = Path(
    '/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/'
    'SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv'
)

QC_TABLE_PATH = Path(os.environ.get(
    'FMRI_QC_TABLE',
    '/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/'
    '6_functional_analysis/non_concat/outputs/tables/fmri_qc_per_subject.csv'))

QC_COLS_KEEP = [
    'mean_fd',
    'mean_fd_post_censoring',
    'num_retained_volumes',
    'num_censored_volumes',
    'TR_sec',
    'minutes_quality_data',
]


def filter_connectivity_cols(cols: list[str]) -> list[str]:
    """Cerebellar filter applied in the primary analyses."""
    filt = [c for c in cols if c.startswith('con_')]
    filt = [c for c in filt if 'Cerebellar' not in c or 'con_Cerebellar_Region1' in c]
    filt = [c for c in filt if 'Cerebellar_Region10' not in c]
    return filt


def load_qc_table() -> pd.DataFrame:
    if not QC_TABLE_PATH.exists():
        raise FileNotFoundError(
            f'QC table not found at {QC_TABLE_PATH}. '
            'Run 00_build_fmri_qc_table.py first.'
        )
    qc = pd.read_csv(QC_TABLE_PATH)
    qc['ID'] = qc['ID'].astype(str)
    return qc


def load_connectivity_and_clusters() -> tuple[pd.DataFrame, list[str]]:
    """
    Load the harmonised connectivity dataframe, merge cluster labels +
    QC metrics, rename TD->NT, and return (df_merged, connectivity_cols).
    Equivalent to the loaders used in 01_autism_td_connectivity_analysis.py
    and 01_generate_cluster_fmri_inputs.py, but augmented with QC.
    """
    df_conn = pd.read_csv(FMRI_CONN_FILE, low_memory=False)
    df_conn['ID'] = df_conn['ID'].astype(str)
    df_conn = df_conn.drop(columns=['PopulationS1', 'Population1'], errors='ignore')
    # Drop any pre-existing QC columns so the QC merge below provides the
    # canonical, recomputed values without _x/_y suffix collisions.
    df_conn = df_conn.drop(columns=[c for c in QC_COLS_KEEP if c in df_conn.columns],
                           errors='ignore')

    df_clusters = pd.read_csv(DF_CLUSTERS_FILE, low_memory=False)
    df_clusters['ID'] = df_clusters['ID'].astype(str)

    df = df_conn.merge(
        df_clusters[['ID', 'Cluster', 'PopulationS1', 'Population1']],
        on='ID',
        how='inner',
    )

    qc = load_qc_table()
    df = df.merge(qc[['ID'] + QC_COLS_KEEP], on='ID', how='left')

    df['PopulationS1'] = df['PopulationS1'].replace('TD', 'NT')

    conn_cols = filter_connectivity_cols(df.columns.tolist())
    return df, conn_cols


def apply_sensitivity_filters(
    df: pd.DataFrame,
    min_minutes: float | None = None,
    max_mean_fd: float | None = None,
) -> tuple[pd.DataFrame, dict]:
    """
    Return df filtered by minutes_quality_data >= min_minutes and
    mean_fd < max_mean_fd. Subjects with NaN QC values are dropped only
    when the corresponding filter is active.
    Also returns a dict summarising how many subjects were dropped at each step.
    """
    n0 = len(df)
    df = df.copy()
    counts = {'n_initial': n0}

    if min_minutes is not None:
        before = len(df)
        df = df[df['minutes_quality_data'].notna() & (df['minutes_quality_data'] >= min_minutes)]
        counts[f'dropped_minutes_lt_{min_minutes:g}'] = before - len(df)

    if max_mean_fd is not None:
        before = len(df)
        df = df[df['mean_fd'].notna() & (df['mean_fd'] < max_mean_fd)]
        counts[f'dropped_mean_fd_ge_{max_mean_fd:g}'] = before - len(df)

    counts['n_remaining'] = len(df)
    return df, counts


def welch_per_feature(
    df: pd.DataFrame,
    feature_cols: list[str],
    group_col: str,
    group_a: str,
    group_b: str,
) -> pd.DataFrame:
    """Two-sided Welch's t (group_a vs group_b) per feature, BH-FDR over all features."""
    a_df = df[df[group_col] == group_a]
    b_df = df[df[group_col] == group_b]
    rows = []
    for feat in feature_cols:
        a = a_df[feat].dropna().values
        b = b_df[feat].dropna().values
        if len(a) < 3 or len(b) < 3:
            continue
        try:
            t_stat, p_val = ttest_ind(a, b, equal_var=False, nan_policy='omit')
        except Exception:
            continue
        rows.append({
            'feature': feat,
            't_stat': float(t_stat),
            'p_value': float(p_val),
            f'n_{group_a}': int(len(a)),
            f'n_{group_b}': int(len(b)),
        })
    if not rows:
        return pd.DataFrame(columns=['feature', 't_stat', 'p_value', 'p_value_fdr',
                                     f'n_{group_a}', f'n_{group_b}'])
    out = pd.DataFrame(rows)
    out['p_value_fdr'] = multipletests(out['p_value'].values, method='fdr_bh')[1]
    return out


def ols_per_feature(
    df: pd.DataFrame,
    feature_cols: list[str],
    group_col: str,
    group_a: str,
    group_b: str,
    covariates: list[str],
) -> pd.DataFrame:
    """
    For each feature fit OLS y ~ group_dummy + covariates, return the
    t-statistic and p-value of the group dummy. group_dummy = 1 for group_a,
    0 for group_b. BH-FDR over all features.
    Subjects with missing covariate values are dropped per feature.
    """
    sub = df[df[group_col].isin([group_a, group_b])].copy()
    sub['group_dummy'] = (sub[group_col] == group_a).astype(int)

    rows = []
    for feat in feature_cols:
        cols_needed = ['group_dummy', feat] + covariates
        data = sub[cols_needed].dropna()
        if data['group_dummy'].nunique() < 2:
            continue
        n_a = int(data['group_dummy'].sum())
        n_b = int((1 - data['group_dummy']).sum())
        if n_a < 3 or n_b < 3:
            continue
        X = sm.add_constant(data[['group_dummy'] + covariates])
        y = data[feat]
        try:
            model = sm.OLS(y, X, missing='drop').fit()
            t_stat = float(model.tvalues['group_dummy'])
            p_val = float(model.pvalues['group_dummy'])
            beta = float(model.params['group_dummy'])
        except Exception:
            continue
        rows.append({
            'feature': feat,
            't_stat': t_stat,
            'p_value': p_val,
            'beta': beta,
            f'n_{group_a}': n_a,
            f'n_{group_b}': n_b,
        })
    if not rows:
        return pd.DataFrame(columns=['feature', 't_stat', 'p_value', 'beta', 'p_value_fdr',
                                     f'n_{group_a}', f'n_{group_b}'])
    out = pd.DataFrame(rows)
    out['p_value_fdr'] = multipletests(out['p_value'].values, method='fdr_bh')[1]
    return out


def split_hyper_hypo(df_t_results: pd.DataFrame, alpha: float = 0.05) -> dict:
    """Match the primary script's significant/hyper/hypo split with FDR-then-uncorrected fallback."""
    df_sig = df_t_results[df_t_results['p_value_fdr'] < alpha].copy()
    fallback_used = False
    if len(df_sig) == 0:
        df_sig = df_t_results[df_t_results['p_value'] < alpha].copy()
        fallback_used = True
    return {
        'all_results': df_t_results,
        'significant_results': df_sig,
        'hyper_results': df_sig[df_sig['t_stat'] > 0].copy(),
        'hypo_results': df_sig[df_sig['t_stat'] < 0].copy(),
        'fallback_used': fallback_used,
    }


def add_source_target_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Split 'con_<source>/<target>' into 'source' + 'target' columns."""
    df = df.copy()
    df['source'] = df['feature'].str.split('/').str[0].str.replace('con_', '', regex=False)
    df['target'] = df['feature'].str.split('/').str[1]
    return df


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float('nan')
    pooled = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    if pooled == 0:
        return float('nan')
    return float((np.mean(a) - np.mean(b)) / pooled)
