#!/usr/bin/env python3
"""
EEG Alpha Peak Cluster Analysis

Analyzes individual alpha frequency comparing Autism clusters vs TD:
- Loads combined corrected alpha peak data merged with cluster assignments
- Age + sex regression and z-score normalization using TD as reference
- Cluster vs TD statistical comparisons (t-test or Mann-Whitney, FDR-corrected)
- Pairwise cluster comparisons
- Violin plots with significance brackets

Adapted from eeg_mri-pipeline/analysis/figures_papers/eeg_cluster_analysis/run_eeg_alpha_peak.py
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent  # 9_cluster_eeg_analysis/

FIG_DIR = _SECTION_DIR / 'outputs' / 'figures'

EEG_MRI_RESULTS = Path('/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results')
DATA_DIR = EEG_MRI_RESULTS / 'dataset_paper' / 'dataframes'
DF_CLUSTERS_FILE = DATA_DIR / 'df_clusters_complete_kmeans.csv'
ALPHA_PEAK_FILE = DATA_DIR / 'Alpha_peak_combined_corrected.csv'
POWER_SPECTRUM_FILE = DATA_DIR / 'Power_spectrum_combined_corrected.csv'

# =============================================================================
# COLORS / ORDER
# =============================================================================
CLUSTER_COLORS = {'C1': '#7A8B47', 'C2': '#ff9fa0', 'C3': '#e7ba52', 'NT': '#C1C2BC'}
CLUSTER_ORDER = ['C1', 'C2', 'C3', 'NT']


# =============================================================================
# DATA LOADING
# =============================================================================

def load_alpha_peak_with_clusters() -> pd.DataFrame:
    """Load and merge alpha peak data with cluster assignments.

    Returns a DataFrame with columns:
        ID, Cluster, PopulationS1, age_yrs, Sex, alpha_peak_corrected
    filtered to Autism (with cluster) and TD subjects only.
    """
    clusters = pd.read_csv(DF_CLUSTERS_FILE, low_memory=False)
    alpha_peak_combined = pd.read_csv(ALPHA_PEAK_FILE)

    print(f"Loaded {len(alpha_peak_combined)} subjects from combined alpha peak file")

    # Align ID types for merging
    alpha_peak_combined['ID'] = alpha_peak_combined['ID'].astype(str)
    clusters_subset = clusters[['ID', 'Cluster', 'PopulationS1', 'age_yrs', 'Sex', 'control_status']].copy()
    clusters_subset['ID'] = clusters_subset['ID'].astype(str)

    # Inner join: keep subjects present in both files
    alpha_peak = alpha_peak_combined.merge(clusters_subset, on='ID', how='inner',
                                           suffixes=('', '_clusters'))

    # Resolve duplicated columns from merge: prefer the clusters file values where available
    for col in ['PopulationS1', 'age_yrs', 'Sex', 'control_status']:
        col_clusters = f'{col}_clusters'
        if col_clusters in alpha_peak.columns:
            alpha_peak[col] = alpha_peak[col_clusters].combine_first(alpha_peak[col])
            alpha_peak.drop(columns=[col_clusters], inplace=True)

    # Rename population label TD -> NT for display
    alpha_peak['PopulationS1'] = alpha_peak['PopulationS1'].replace('TD', 'NT')

    # Keep only needed columns
    keep_cols = ['ID', 'Cluster', 'PopulationS1', 'age_yrs', 'Sex',
                 'alpha_peak_corrected', 'control_status']
    keep_cols = [c for c in keep_cols if c in alpha_peak.columns]
    alpha_peak = alpha_peak[keep_cols].copy()

    # Quality filters
    alpha_peak = alpha_peak.dropna(subset=['age_yrs'])
    alpha_peak = alpha_peak[alpha_peak['alpha_peak_corrected'].notna()]
    alpha_peak = alpha_peak[alpha_peak['alpha_peak_corrected'] != 0]

    # Population filter: Autism with cluster OR NT
    autism_mask = (alpha_peak['PopulationS1'] == 'Autism') & alpha_peak['Cluster'].notna()
    td_mask = alpha_peak['PopulationS1'] == 'NT'
    alpha_peak = alpha_peak[autism_mask | td_mask].reset_index(drop=True)

    try:
        counts = alpha_peak['Cluster'].fillna('NT').value_counts(dropna=False).to_dict()
        print(f"Subjects after filtering — {counts}")
        pop_counts = alpha_peak['PopulationS1'].value_counts(dropna=False).to_dict()
        print(f"By PopulationS1: {pop_counts}")
    except Exception as exc:
        print(f"Warning: could not compute counts: {exc}")

    return alpha_peak


# =============================================================================
# REGRESSION & Z-SCORE
# =============================================================================

def regress_and_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Regress alpha_peak_corrected on age, age^2, sex; z-score residuals using TD mean/std.

    Parameters
    ----------
    df : DataFrame with columns age_yrs, Sex, alpha_peak_corrected, control_status

    Returns
    -------
    df with alpha_peak_corrected replaced by z-scored residuals.
    """
    df = df.copy()
    df['age_yrs_sq'] = df['age_yrs'] ** 2
    X = df[['age_yrs', 'age_yrs_sq', 'Sex']].copy()
    X['Sex'] = X['Sex'].fillna(0)
    y = df['alpha_peak_corrected']

    model = LinearRegression().fit(X, y)
    residuals = y - model.predict(X)

    # Z-score using TD (control) group statistics
    td_mask = (
        df['control_status'].str.contains('control', case=False, na=False)
        if 'control_status' in df.columns
        else df['PopulationS1'] == 'NT'
    )
    td_residuals = residuals[td_mask]
    if len(td_residuals) < 5:
        print("Warning: fewer than 5 TD subjects for z-scoring; using all subjects.")
        td_residuals = residuals

    mu = td_residuals.mean()
    sd = td_residuals.std()
    if sd == 0 or np.isnan(sd):
        print("Warning: TD std is 0 or NaN; skipping z-scoring.")
        df['alpha_peak_corrected'] = residuals
    else:
        df['alpha_peak_corrected'] = (residuals - mu) / sd

    df.drop(columns=['age_yrs_sq'], inplace=True)
    return df


# =============================================================================
# STATISTICS
# =============================================================================

def autism_clusters_vs_td_stats(df: pd.DataFrame, value_col: str,
                                 order: tuple = ('C1', 'C2', 'C3')) -> pd.DataFrame:
    """Compute cluster vs TD and pairwise cluster comparisons with FDR correction.

    For each cluster vs TD: Shapiro test → t-test or Mann-Whitney.
    FDR correction (BH) across all cluster vs TD tests.
    Pairwise cluster comparisons with separate FDR correction.

    Parameters
    ----------
    df        : DataFrame containing Cluster, PopulationS1, and value_col.
    value_col : Name of the column to compare.
    order     : Tuple of cluster names to include.

    Returns
    -------
    DataFrame with columns:
        feature, test_type, comparison, cluster_n, td_n, cluster_mean, td_mean,
        t_or_u, p_uncorrected, p_fdr, significant
    """
    autism_data = df[df['PopulationS1'] == 'Autism'].dropna(subset=['Cluster'])
    td_data = df[df['PopulationS1'] == 'NT']
    td_values = td_data[value_col].dropna()

    rows = []
    cluster_vs_td_buffer = []

    # ---- Cluster vs TD -------------------------------------------------------
    for cluster in order:
        cluster_values = autism_data[autism_data['Cluster'] == cluster][value_col].dropna()
        if len(cluster_values) < 3 or len(td_values) < 3:
            continue

        # Normality test
        cluster_normal = stats.shapiro(cluster_values)[1] > 0.05 if len(cluster_values) >= 3 else False
        td_normal = stats.shapiro(td_values)[1] > 0.05 if len(td_values) >= 3 else False
        test_type = 't-test' if (cluster_normal and td_normal) else 'Mann-Whitney'

        if test_type == 't-test':
            t_val, p = stats.ttest_ind(cluster_values, td_values, equal_var=False)
        else:
            t_val, p = stats.mannwhitneyu(cluster_values, td_values, alternative='two-sided')

        cluster_vs_td_buffer.append({
            'feature': value_col,
            'test_type': test_type,
            'comparison': f'{cluster} vs NT',
            'cluster_n': len(cluster_values),
            'td_n': len(td_values),
            'cluster_mean': float(cluster_values.mean()),
            'td_mean': float(td_values.mean()),
            't_or_u': float(t_val),
            'p_uncorrected': float(p),
        })

    # FDR correction for cluster vs TD
    if cluster_vs_td_buffer:
        pvals = [r['p_uncorrected'] for r in cluster_vs_td_buffer]
        _, pvals_fdr, _, _ = smm.multipletests(pvals, method='fdr_bh', alpha=0.05)
        for rec, p_fdr in zip(cluster_vs_td_buffer, pvals_fdr):
            rec['p_fdr'] = float(p_fdr)
            rec['significant'] = bool(p_fdr < 0.05)
            rows.append(rec)

    # ---- Pairwise cluster comparisons ----------------------------------------
    cluster_arrays = [
        autism_data[autism_data['Cluster'] == c][value_col].dropna()
        for c in order
        if len(autism_data[autism_data['Cluster'] == c]) >= 3
    ]
    valid_clusters = [
        c for c in order
        if len(autism_data[autism_data['Cluster'] == c]) >= 3
    ]

    if len(cluster_arrays) >= 2:
        # Global normality check to decide ANOVA vs Kruskal
        all_normal = all(
            stats.shapiro(arr)[1] > 0.05 if len(arr) >= 3 else False
            for arr in cluster_arrays
        )
        pw_test_type = 'ANOVA' if all_normal else 'Kruskal'

        pair_buffer = []
        for c1, c2 in combinations(valid_clusters, 2):
            v1 = autism_data[autism_data['Cluster'] == c1][value_col].dropna()
            v2 = autism_data[autism_data['Cluster'] == c2][value_col].dropna()
            if len(v1) < 3 or len(v2) < 3:
                continue
            if pw_test_type == 'ANOVA':
                t_val, p = stats.ttest_ind(v1, v2, equal_var=False)
            else:
                t_val, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
            pair_buffer.append({
                'feature': value_col,
                'test_type': pw_test_type,
                'comparison': f'{c1} vs {c2}',
                'cluster_n': len(v1),
                'td_n': len(v2),
                'cluster_mean': float(v1.mean()),
                'td_mean': float(v2.mean()),
                't_or_u': float(t_val),
                'p_uncorrected': float(p),
            })

        if pair_buffer:
            pvals_pw = [r['p_uncorrected'] for r in pair_buffer]
            rej_pw, pvals_pw_fdr, _, _ = smm.multipletests(pvals_pw, method='fdr_bh', alpha=0.05)
            for rec, p_fdr, rejected in zip(pair_buffer, pvals_pw_fdr, rej_pw):
                rec['p_fdr'] = float(p_fdr)
                rec['significant'] = bool(rejected)
                rows.append(rec)

    return pd.DataFrame(rows)


# =============================================================================
# PLOTTING
# =============================================================================

def plot_alpha_peak_clusters_vs_td(df: pd.DataFrame, stats_df: pd.DataFrame,
                                    order: tuple = ('C1', 'C2', 'C3')) -> None:
    """Violin + strip + pointplot of alpha peak by cluster and TD.

    Adds significance brackets for FDR-significant cluster vs TD comparisons.
    Saves figure to FIG_DIR/eeg_alpha_peak_clusters_vs_td_violin.pdf.
    """
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    autism_data = df[df['PopulationS1'] == 'Autism'].dropna(subset=['Cluster'])
    td_data = df[df['PopulationS1'] == 'NT'].copy()
    td_data['Cluster'] = 'NT'

    plot_data = pd.concat([autism_data, td_data], ignore_index=True)
    plot_order = list(order) + ['NT']

    fig, ax = plt.subplots(1, 1, figsize=(8, 10))

    sns.violinplot(
        y='alpha_peak_corrected', x='Cluster', data=plot_data,
        inner=None, linewidth=0, order=plot_order, palette=CLUSTER_COLORS, ax=ax,
    )
    sns.stripplot(
        y='alpha_peak_corrected', x='Cluster', data=plot_data,
        color='black', size=2, jitter=True, alpha=0.5, order=plot_order, ax=ax,
    )
    sns.pointplot(
        y='alpha_peak_corrected', x='Cluster', data=plot_data,
        color='#8A0201', join=False, errorbar=('ci', 95), scale=1.2,
        order=plot_order, ax=ax,
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set_facecolor('white')

    counts = plot_data['Cluster'].value_counts()
    x_labels = [f"{c} (n={counts.get(c, 0)})" for c in plot_order]
    ax.set_xticklabels(x_labels, rotation=15, ha='right')
    ax.set_ylabel('alpha peak (z-score)', fontsize=18)
    ax.set_xlabel('')
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Significance brackets for cluster vs TD (FDR-corrected)
    if stats_df is not None and not stats_df.empty:
        cluster_vs_td = stats_df[stats_df['comparison'].str.contains(' vs NT', na=False)]
        y_vals = plot_data['alpha_peak_corrected'].dropna()
        y_max = y_vals.max()
        y_range = y_vals.max() - y_vals.min()
        y_step = y_range * 0.08
        pos_map = {name: i for i, name in enumerate(plot_order)}

        bracket_y = y_max + y_step
        for _, row in cluster_vs_td.iterrows():
            if not row.get('significant', False):
                continue
            cluster_label = row['comparison'].split(' vs NT')[0]
            if cluster_label not in pos_map or 'NT' not in pos_map:
                continue
            x1 = pos_map[cluster_label]
            x2 = pos_map['NT']
            p_fdr = row['p_fdr']

            # Bracket
            ax.plot(
                [x1, x1, x2, x2],
                [bracket_y, bracket_y + y_step * 0.3, bracket_y + y_step * 0.3, bracket_y],
                color='black', linewidth=1.2,
            )
            # P-value text
            if p_fdr < 0.001:
                p_text = f'{p_fdr:.0e}'
            else:
                p_text = f'{p_fdr:.3f}'
            ax.text(
                (x1 + x2) / 2, bracket_y + y_step * 0.4, p_text,
                ha='center', va='bottom', fontsize=11, fontweight='bold',
            )
            bracket_y += y_step * 1.2

    out_path = FIG_DIR / 'eeg_alpha_peak_clusters_vs_td_violin.pdf'
    fig.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load data
    alpha_peak = load_alpha_peak_with_clusters()
    print(f"Total subjects for analysis: {len(alpha_peak)}")

    # 2. Regression + z-score (optional — data is already combat-corrected;
    #    this step additionally removes residual age/sex effects and
    #    normalises to the TD group scale)
    alpha_peak = regress_and_zscore(alpha_peak)

    # 3. Statistics
    order = ('C1', 'C2', 'C3')
    stats_df = autism_clusters_vs_td_stats(alpha_peak, 'alpha_peak_corrected', order)
    print(stats_df[['comparison', 'test_type', 'p_uncorrected', 'p_fdr', 'significant']].to_string())

    # Save stats
    stats_out = FIG_DIR / 'alpha_peak_clusters_vs_nt_stats.csv'
    stats_df.to_csv(stats_out, index=False)
    print(f"Stats saved: {stats_out}")

    # 4. Plot
    plot_alpha_peak_clusters_vs_td(alpha_peak, stats_df, order)

    print(f"\nAll outputs saved to: {FIG_DIR}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
