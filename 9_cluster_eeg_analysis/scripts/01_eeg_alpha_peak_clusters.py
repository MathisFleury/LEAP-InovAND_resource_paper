#!/usr/bin/env python3
"""
EEG Alpha Peak — Autism k-means clusters vs NT.

Consumes the curated, already-corrected alpha peak built by
  8_eeg_analysis/preprocessing/build_alpha_peak_corrected.py
(raw $IMG5 rebuild + curated demographics + curated k-means Cluster +
full-sample regress/z-score). This script does NOT correct again — it filters
to Autism-with-cluster + NT, then runs cluster-vs-NT and pairwise stats + plot.
"""

import sys
from pathlib import Path
import pandas as pd
from scipy import stats
from itertools import combinations
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# PATHS  (k-means + curated only)
# =============================================================================
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent  # 9_cluster_eeg_analysis/
sys.path.insert(0, str(_SECTION_DIR.parent / "8_eeg_analysis" / "preprocessing"))
from config import ALPHA_PEAK_CORRECTED  # noqa: E402

FIG_DIR = _SECTION_DIR / 'outputs' / 'figures'
TABLES_DIR = _SECTION_DIR / 'outputs' / 'tables'

# =============================================================================
# COLORS / ORDER
# =============================================================================
CLUSTER_COLORS = {'C1': '#7A8B47', 'C2': '#ff9fa0', 'C3': '#e7ba52',
                  'NT': '#C1C2BC', 'IDD': '#D8A4CB'}
CLUSTER_ORDER = ['C1', 'C2', 'C3', 'NT', 'IDD']


# =============================================================================
# DATA LOADING
# =============================================================================

def load_alpha_peak_with_clusters() -> pd.DataFrame:
    """Load the corrected alpha peak; keep Autism-with-cluster + NT subjects.

    Cluster / PopulationS1 / demographics are already curated k-means in the file.
    """
    alpha_peak = pd.read_csv(ALPHA_PEAK_CORRECTED)
    print(f"Loaded {len(alpha_peak)} subjects from {ALPHA_PEAK_CORRECTED.name}")
    alpha_peak['ID'] = alpha_peak['ID'].astype(str)
    alpha_peak = alpha_peak[alpha_peak['alpha_peak_corrected'].notna()]

    # Population filter: Autism split by cluster (C1/C2/C3), plus NT (reference)
    # and IDD, each kept as its own category.
    autism_mask = (alpha_peak['PopulationS1'] == 'Autism') & alpha_peak['Cluster'].notna()
    td_mask = alpha_peak['PopulationS1'] == 'NT'
    idd_mask = alpha_peak['PopulationS1'] == 'IDD'
    alpha_peak = alpha_peak[autism_mask | td_mask | idd_mask].reset_index(drop=True)

    autism_by_cluster = (alpha_peak.loc[alpha_peak['PopulationS1'] == 'Autism', 'Cluster']
                         .value_counts().to_dict())
    n_nt = int((alpha_peak['PopulationS1'] == 'NT').sum())
    n_idd = int((alpha_peak['PopulationS1'] == 'IDD').sum())
    print(f"Autism by cluster: {autism_by_cluster} | NT (reference): {n_nt} | IDD: {n_idd}")
    return alpha_peak


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

    # ---- Group vs NT (clusters C1/C2/C3 and IDD) -----------------------------
    groups_vs_nt = [(c, autism_data[autism_data['Cluster'] == c][value_col].dropna())
                    for c in order]
    groups_vs_nt.append(('IDD', df[df['PopulationS1'] == 'IDD'][value_col].dropna()))
    for group_label, group_values in groups_vs_nt:
        if len(group_values) < 3 or len(td_values) < 3:
            continue

        # Normality test
        group_normal = stats.shapiro(group_values)[1] > 0.05 if len(group_values) >= 3 else False
        td_normal = stats.shapiro(td_values)[1] > 0.05 if len(td_values) >= 3 else False
        test_type = 't-test' if (group_normal and td_normal) else 'Mann-Whitney'

        if test_type == 't-test':
            t_val, p = stats.ttest_ind(group_values, td_values, equal_var=False)
        else:
            t_val, p = stats.mannwhitneyu(group_values, td_values, alternative='two-sided')

        cluster_vs_td_buffer.append({
            'feature': value_col,
            'test_type': test_type,
            'comparison': f'{group_label} vs NT',
            'cluster_n': len(group_values),
            'td_n': len(td_values),
            'cluster_mean': float(group_values.mean()),
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
    idd_data = df[df['PopulationS1'] == 'IDD'].copy()
    idd_data['Cluster'] = 'IDD'

    plot_data = pd.concat([autism_data, td_data, idd_data], ignore_index=True)
    plot_order = list(order) + ['NT', 'IDD']

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

        def _draw_bracket(x1, x2, y, p_fdr):
            ax.plot(
                [x1, x1, x2, x2],
                [y, y + y_step * 0.3, y + y_step * 0.3, y],
                color='black', linewidth=1.2,
            )
            p_text = f'{p_fdr:.0e}' if p_fdr < 0.001 else f'{p_fdr:.3f}'
            ax.text(
                (x1 + x2) / 2, y + y_step * 0.4, p_text,
                ha='center', va='bottom', fontsize=11, fontweight='bold',
            )

        bracket_y = y_max + y_step
        # Group vs NT (C1/C2/C3 and IDD)
        for _, row in cluster_vs_td.iterrows():
            if not row.get('significant', False):
                continue
            group_label = row['comparison'].split(' vs NT')[0]
            if group_label not in pos_map or 'NT' not in pos_map:
                continue
            _draw_bracket(pos_map[group_label], pos_map['NT'], bracket_y, row['p_fdr'])
            bracket_y += y_step * 1.2

        # Pairwise inter-cluster comparisons (C1 vs C2, etc.), stacked above
        pairwise = stats_df[~stats_df['comparison'].str.contains(' vs NT', na=False)]
        for _, row in pairwise.iterrows():
            if not row.get('significant', False):
                continue
            c1, c2 = [s.strip() for s in row['comparison'].split(' vs ')]
            if c1 not in pos_map or c2 not in pos_map:
                continue
            _draw_bracket(pos_map[c1], pos_map[c2], bracket_y, row['p_fdr'])
            bracket_y += y_step * 1.2

    out_path = FIG_DIR / 'eeg_alpha_peak_clusters_vs_td_violin.pdf'
    fig.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:
    print(f"k-means clusters vs NT  →  {FIG_DIR}")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load already-corrected data (curated k-means)
    alpha_peak = load_alpha_peak_with_clusters()
    print(f"Total subjects for analysis: {len(alpha_peak)}")

    # 2. Statistics
    order = ('C1', 'C2', 'C3')
    stats_df = autism_clusters_vs_td_stats(alpha_peak, 'alpha_peak_corrected', order)
    print(stats_df[['comparison', 'test_type', 'p_uncorrected', 'p_fdr', 'significant']].to_string())

    # Save stats
    stats_out = TABLES_DIR / 'alpha_peak_clusters_vs_nt_stats.csv'
    stats_df.to_csv(stats_out, index=False)
    print(f"Stats saved: {stats_out}")

    # 4. Plot
    plot_alpha_peak_clusters_vs_td(alpha_peak, stats_df, order)

    print(f"\nAll outputs saved to: {FIG_DIR}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
