#!/usr/bin/env python3
"""
Cluster-stratified fMRI Connectivity Analysis

For each autism cluster (C1, C2, C3):
  - Compares functional connectivity against a pooled TD reference
  - Welch's t-tests per connectivity pair with FDR correction
  - Cohen's d effect sizes
  - Pairwise cluster comparisons (autism-only)
  - Saves R-ready input files for Schaefer atlas visualisation

Adapted from eeg_mri-pipeline — self-contained, no sys.path modifications.
"""

import os
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent

OUTPUT_DIR = _SECTION_DIR / 'outputs' / 'tables'  # pure-table step: no plots
R_INPUT_DIR = OUTPUT_DIR / 'r_input_files'

EEG_MRI_RESULTS = Path('/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results')
DF_CLUSTERS_FILE = EEG_MRI_RESULTS / 'dataset_paper' / 'dataframes' / 'df_clusters_complete_kmeans.csv'
# Default to the manuscript's primary regime (concat/nogsr_concat); override via
# FMRI_CONN_FILE (same env var 06_cluster_nbs.py / _sensitivity_utils honour).
FMRI_CONN_FILE = Path(os.environ.get(
    'FMRI_CONN_FILE',
    str(_SECTION_DIR.parent / '6_functional_analysis' / 'concat' / 'preprocessing'
        / 'outputs' / 'nogsr_concat' / 'df_conn_cohort_norm_nogsr_concat.csv')))
LEAP_QC_FILE     = EEG_MRI_RESULTS / 'anat_qc' / 'dataframes' / 'subjects_to_remove_LEAP.txt'
ATLAS_FILE       = (
    '/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/'
    'SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv'
)

CLUSTERS = ['C1', 'C2', 'C3']
PAIRWISE_PAIRS = [('C1', 'C2'), ('C1', 'C3'), ('C2', 'C3')]


# =============================================================================
# HELPERS
# =============================================================================

def _cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Pooled-SD Cohen's d (Welch-style: no assumption of equal n)."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    pooled_sd = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    if pooled_sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled_sd


def _cohens_d_ci(d: float, na: int, nb: int) -> tuple[float, float]:
    """Analytic 95 % CI for Cohen's d (Hedges & Olkin 1985)."""
    if not np.isfinite(d) or na < 2 or nb < 2:
        return (np.nan, np.nan)
    se = np.sqrt((na + nb) / (na * nb) + d ** 2 / (2.0 * (na + nb)))
    z = 1.959963984540054  # qnorm(0.975)
    return (float(d - z * se), float(d + z * se))


def _apply_cerebellar_filter(cols: list[str]) -> list[str]:
    """
    Keep only con_* columns; within those remove all Cerebellar regions
    EXCEPT con_Cerebellar_Region1, and also remove con_Cerebellar_Region10.
    """
    filtered = [c for c in cols if c.startswith('con_')]
    filtered = [
        c for c in filtered
        if 'Cerebellar' not in c or 'con_Cerebellar_Region1' in c
    ]
    filtered = [c for c in filtered if 'Cerebellar_Region10' not in c]
    return filtered


def _fdr_correct(p_values: pd.Series) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction; handles NaN gracefully."""
    mask = ~p_values.isna()
    p_fdr = np.full(len(p_values), np.nan)
    if mask.sum() > 0:
        p_fdr[mask] = multipletests(p_values[mask], method='fdr_bh')[1]
    return p_fdr


# =============================================================================
# DATA LOADING
# =============================================================================

def load_and_prepare_data():
    """
    Returns
    -------
    df_merged : pd.DataFrame
        Inner join of fMRI connectivity + cluster labels, LEAP-QC filtered.
    connectivity_cols : list[str]
        Cerebellar-filtered list of connectivity feature column names.
    """
    print("=" * 60)
    print("LOADING DATA")
    print("=" * 60)

    # -- Clusters --
    if not DF_CLUSTERS_FILE.exists():
        raise FileNotFoundError(f"Cluster file not found: {DF_CLUSTERS_FILE}")
    df_clusters = pd.read_csv(DF_CLUSTERS_FILE, low_memory=False)
    df_clusters['ID'] = df_clusters['ID'].astype(str)
    print(f"Loaded cluster data: {len(df_clusters)} participants")

    # -- fMRI connectivity --
    if not FMRI_CONN_FILE.exists():
        raise FileNotFoundError(f"fMRI connectivity file not found: {FMRI_CONN_FILE}")
    df_fmri = pd.read_csv(FMRI_CONN_FILE)
    df_fmri['ID'] = df_fmri['ID'].astype(str)
    print(f"Loaded fMRI connectivity data: {len(df_fmri)} participants")

    # -- LEAP QC filtering --
    if LEAP_QC_FILE.exists():
        leap_qc = pd.read_csv(LEAP_QC_FILE, header=None)
        leap_qc[0] = leap_qc[0].astype(str)
        ids_to_remove = set(leap_qc[0])
        before = len(df_fmri)
        df_fmri = df_fmri[~df_fmri['ID'].isin(ids_to_remove)].copy()
        print(f"LEAP QC: removed {before - len(df_fmri)} participants ({len(df_fmri)} remaining)")
    else:
        print(f"Warning: LEAP QC file not found, skipping QC filter: {LEAP_QC_FILE}")

    # -- Merge fMRI with cluster labels --
    # Diagnosis comes from the curated clinical merge (population_group) carried
    # in the preprocessing output; only the Cluster labels are taken from the
    # cluster file (paper k-means).
    df_merged = df_fmri.merge(
        df_clusters[['ID', 'Cluster']],
        on='ID',
        how='inner',
    )
    print(f"After inner merge: {len(df_merged)} participants")

    # Curated diagnosis (TD -> NT for display).
    df_merged['PopulationS1'] = df_merged['population_group'].replace('TD', 'NT')

    print(f"Population breakdown: {df_merged['PopulationS1'].value_counts().to_dict()}")
    if 'Cluster' in df_merged.columns:
        autism_only = df_merged[df_merged['PopulationS1'] == 'Autism']
        print(f"Cluster breakdown (Autism): {autism_only['Cluster'].value_counts().to_dict()}")

    # -- Connectivity column selection + cerebellar filter --
    raw_conn_cols = [c for c in df_merged.columns if c.startswith('con_')]
    connectivity_cols = _apply_cerebellar_filter(raw_conn_cols)
    print(f"Connectivity features: {len(raw_conn_cols)} raw → {len(connectivity_cols)} after cerebellar filter")

    return df_merged, connectivity_cols


# =============================================================================
# CLUSTER vs TD ANALYSIS
# =============================================================================

def run_cluster_vs_td_analysis(df_merged: pd.DataFrame, connectivity_cols: list[str], output_dir: Path):
    """
    For each cluster (C1, C2, C3):
      - autism_data  = Autism subjects assigned to that cluster
      - td_data      = ALL TD subjects (pooled across clusters)
      - Welch's t-test + Cohen's d per connectivity feature
      - FDR correction (Benjamini-Hochberg)
      - Saves full results CSV and R-ready input CSV
    """
    print("\n" + "=" * 60)
    print("CLUSTER vs POOLED-TD CONNECTIVITY ANALYSIS")
    print("=" * 60)

    td_data = df_merged[df_merged['PopulationS1'] == 'NT'].copy()
    print(f"Pooled NT reference: {len(td_data)} participants")

    r_input_dir = output_dir / 'r_input_files'
    r_input_dir.mkdir(parents=True, exist_ok=True)

    summary = {}

    for cluster in CLUSTERS:
        print(f"\n--- Cluster {cluster} ---")
        autism_data = df_merged[
            (df_merged['Cluster'] == cluster) &
            (df_merged['PopulationS1'] == 'Autism')
        ].copy()
        print(f"  Autism {cluster}: {len(autism_data)} participants")

        if len(autism_data) < 5:
            print(f"  Skipping {cluster}: fewer than 5 autism participants")
            continue

        rows = []
        for feature in connectivity_cols:
            a = autism_data[feature].dropna().values
            b = td_data[feature].dropna().values
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                t_stat, p_val = ttest_ind(a, b, equal_var=False)
                d = _cohens_d(a, b)
                d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
                rows.append({
                    'region':      feature,
                    't_stat':      t_stat,
                    'p_val':       p_val,
                    'cohens_d':    d,
                    'cohens_d_lo': d_lo,
                    'cohens_d_hi': d_hi,
                    'autism_mean': float(np.mean(a)),
                    'td_mean':     float(np.mean(b)),
                    'autism_n':    int(len(a)),
                    'td_n':        int(len(b)),
                })
            except Exception:
                continue

        if not rows:
            print(f"  No valid test results for {cluster}")
            continue

        df_results = pd.DataFrame(rows)
        df_results['p_fdr'] = _fdr_correct(df_results['p_val'])

        # Reorder columns
        df_results = df_results[
            ['region', 't_stat', 'p_val', 'cohens_d',
             'cohens_d_lo', 'cohens_d_hi', 'p_fdr',
             'autism_mean', 'td_mean', 'autism_n', 'td_n']
        ]

        n_sig = int((df_results['p_fdr'] < 0.05).sum())
        print(f"  Significant connections (FDR < 0.05): {n_sig} / {len(df_results)}")
        summary[cluster] = n_sig

        # Full results
        full_out = output_dir / f'cluster_{cluster}_connectivity_autism_vs_td.csv'
        df_results.to_csv(full_out, index=False)
        print(f"  Saved: {full_out.name}")

        # Significant results for network matrix (R script): include both
        # the t-statistic (for backward compatibility) and Cohen's d.
        df_sig = df_results[df_results['p_fdr'] < 0.05][
            ['region', 't_stat', 'cohens_d', 'cohens_d_lo', 'cohens_d_hi', 'p_fdr']
        ].copy()
        df_sig = df_sig.rename(columns={'region': 'feature'})
        sig_out = output_dir / f'cluster_{cluster}_vs_NT_significant_results.csv'
        df_sig.to_csv(sig_out, index=False)
        print(f"  Saved significant results: {sig_out.name}")

        # R-ready input: region, mean_value (kept = t_stat for compatibility),
        # cohens_d (with 95 % CI) — used by the Schaefer brain map.
        df_r = df_results[
            ['region', 't_stat', 'cohens_d', 'cohens_d_lo', 'cohens_d_hi', 'p_fdr']
        ].copy()
        df_r = df_r.rename(columns={'t_stat': 'mean_value'})
        r_out = r_input_dir / f'cluster_{cluster}_connectivity_for_r.csv'
        df_r.to_csv(r_out, index=False)
        print(f"  Saved R input: {r_out.name}")

    return summary


# =============================================================================
# PAIRWISE CLUSTER ANALYSIS
# =============================================================================

def run_pairwise_cluster_analysis(df_merged: pd.DataFrame, connectivity_cols: list[str], output_dir: Path):
    """
    Pairwise Welch's t-tests between autism subjects of C1 vs C2, C1 vs C3, C2 vs C3.
    FDR correction applied across all features for each pair.
    Saves pairwise_{c1}_vs_{c2}_connectivity.csv.
    """
    print("\n" + "=" * 60)
    print("PAIRWISE CLUSTER CONNECTIVITY ANALYSIS (AUTISM ONLY)")
    print("=" * 60)

    autism_df = df_merged[df_merged['PopulationS1'] == 'Autism'].copy()

    for c1, c2 in PAIRWISE_PAIRS:
        print(f"\n--- {c1} vs {c2} ---")
        grp1 = autism_df[autism_df['Cluster'] == c1]
        grp2 = autism_df[autism_df['Cluster'] == c2]
        print(f"  {c1}: {len(grp1)} | {c2}: {len(grp2)}")

        if len(grp1) < 5 or len(grp2) < 5:
            print(f"  Skipping: one group has fewer than 5 participants")
            continue

        rows = []
        for feature in connectivity_cols:
            a = grp1[feature].dropna().values
            b = grp2[feature].dropna().values
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                t_stat, p_val = ttest_ind(a, b, equal_var=False)
                d = _cohens_d(a, b)
                d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
                rows.append({
                    'region':    feature,
                    't_stat':    t_stat,
                    'p_val':     p_val,
                    'cohens_d':  d,
                    'cohens_d_lo': d_lo,
                    'cohens_d_hi': d_hi,
                    f'{c1}_mean': float(np.mean(a)),
                    f'{c2}_mean': float(np.mean(b)),
                    f'{c1}_n':   int(len(a)),
                    f'{c2}_n':   int(len(b)),
                })
            except Exception:
                continue

        if not rows:
            print(f"  No valid test results for {c1} vs {c2}")
            continue

        df_pw = pd.DataFrame(rows)
        df_pw['p_fdr'] = _fdr_correct(df_pw['p_val'])

        n_sig = int((df_pw['p_fdr'] < 0.05).sum())
        print(f"  Significant (FDR < 0.05): {n_sig} / {len(df_pw)}")

        out_file = output_dir / f'pairwise_{c1}_vs_{c2}_connectivity.csv'
        df_pw.to_csv(out_file, index=False)
        print(f"  Saved: {out_file.name}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    R_INPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("CLUSTER fMRI INPUT GENERATION")
    print("=" * 60)
    print(f"Output directory : {OUTPUT_DIR}")
    print(f"R input directory: {R_INPUT_DIR}")

    df_merged, connectivity_cols = load_and_prepare_data()

    # Summary counts per cluster
    print("\n--- Participant counts ---")
    for cluster in CLUSTERS:
        n = int(((df_merged['Cluster'] == cluster) & (df_merged['PopulationS1'] == 'Autism')).sum())
        print(f"  Autism {cluster}: {n}")
    n_td = int((df_merged['PopulationS1'] == 'NT').sum())
    print(f"  NT (pooled): {n_td}")

    sig_summary = run_cluster_vs_td_analysis(df_merged, connectivity_cols, OUTPUT_DIR)
    run_pairwise_cluster_analysis(df_merged, connectivity_cols, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("SUMMARY — significant connections (FDR < 0.05) per cluster")
    print("=" * 60)
    for cluster, n_sig in sig_summary.items():
        print(f"  {cluster}: {n_sig}")

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
