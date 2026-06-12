#!/usr/bin/env python3
"""
EEG Power Bands Cluster Analysis

Analyzes EEG absolute power spectra comparing Autism clusters vs TD:
- Loads combined corrected power spectrum data merged with cluster assignments
- Reshapes to wide format (one row per subject, one column per region/side/band)
- Per-feature age regression and z-score normalization using TD group
- Cluster vs TD statistical comparisons (t-test or Mann-Whitney, FDR-corrected)
- Grid violin plots per frequency band
- Synthesis heatmap (% significant per band per cluster)

Adapted from eeg_mri-pipeline/analysis/figures_papers/eeg_cluster_analysis/run_eeg_power_bands.py
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
STATS_DIR = FIG_DIR / 'stats'
FDR_DIR = FIG_DIR / 'fdr_corrected'
PLOTS_DIR = FIG_DIR / 'plots'

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

def load_power_with_clusters() -> pd.DataFrame:
    """Load and merge power spectrum data with cluster assignments.

    The combined power spectrum file (Power_spectrum_combined_corrected.csv)
    contains columns: ID, Subject_ID, cohort, Channel, Brain_region, Brain_side,
    Eye_status, Freq_band, Quant_status, PSD, PSD_corrected, control_status,
    age_yrs, Sex, PopulationS1, Population1.

    Cluster labels are taken from df_clusters_complete_kmeans.csv (inner join on ID).
    Only Autism subjects with a Cluster assignment and TD subjects are retained.

    Returns
    -------
    Long-format DataFrame merged with cluster metadata.
    """
    clusters = pd.read_csv(DF_CLUSTERS_FILE, low_memory=False)
    power = pd.read_csv(POWER_SPECTRUM_FILE, low_memory=False)

    print(f"Loaded {len(power)} rows from power spectrum file")

    # Align ID types for merging
    power['ID'] = power['ID'].astype(str)
    clusters_subset = clusters[['ID', 'Cluster', 'PopulationS1', 'age_yrs',
                                 'Sex', 'control_status']].copy()
    clusters_subset['ID'] = clusters_subset['ID'].astype(str)

    # Inner join
    power = power.merge(clusters_subset, on='ID', how='inner',
                        suffixes=('', '_clusters'))

    # Resolve duplicated columns: prefer clusters file values
    for col in ['PopulationS1', 'age_yrs', 'Sex', 'control_status']:
        col_clusters = f'{col}_clusters'
        if col_clusters in power.columns:
            power[col] = power[col_clusters].combine_first(power[col])
            power.drop(columns=[col_clusters], inplace=True)

    # Rename population label TD -> NT for display
    power['PopulationS1'] = power['PopulationS1'].replace('TD', 'NT')

    # Quality filter
    power = power.dropna(subset=['age_yrs'])

    # Population filter: Autism with cluster OR NT
    autism_mask = (power['PopulationS1'] == 'Autism') & power['Cluster'].notna()
    td_mask = power['PopulationS1'] == 'NT'
    power = power[autism_mask | td_mask].reset_index(drop=True)

    try:
        unique_subjects = power[['ID', 'PopulationS1', 'Cluster']].drop_duplicates(subset=['ID'])
        pop_counts = unique_subjects['PopulationS1'].value_counts(dropna=False).to_dict()
        cluster_counts = unique_subjects['Cluster'].fillna('NT').value_counts(dropna=False).to_dict()
        print(f"Unique subjects by PopulationS1: {pop_counts}")
        print(f"Unique subjects by Cluster: {cluster_counts}")
    except Exception as exc:
        print(f"Warning: could not compute counts: {exc}")

    return power


# =============================================================================
# RESHAPE
# =============================================================================

def reshape_power(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape long-format power to wide format with one row per subject.

    Filters:
    - Brain_region in ['Central', 'Frontal', 'Posterior'] (if column present)
    - Brain_side in ['Left', 'Right', 'Mid'] (if column present)
    - Quant_status == 'Absolute' (if column present)
    - Eye_status == 'eyeo' (if column present)

    Aggregates PSD by median per (ID, Brain_region, Brain_side, Freq_band).
    Pivots to wide format: columns named 'con_{Brain_region}_{Brain_side}_{Freq_band}'.

    Metadata columns (Cluster, PopulationS1, age_yrs, Sex, control_status)
    are preserved via a first-value aggregation and merged back.

    Returns
    -------
    Wide-format DataFrame indexed by ID with power features and metadata.
    """
    filt = df.copy()

    # Apply available filters
    if 'Brain_region' in filt.columns:
        filt = filt[filt['Brain_region'].isin(['Central', 'Frontal', 'Posterior'])]
    if 'Brain_side' in filt.columns:
        filt = filt[filt['Brain_side'].isin(['Left', 'Right', 'Mid'])]
    if 'Quant_status' in filt.columns:
        filt = filt[filt['Quant_status'] == 'Absolute']
    if 'Eye_status' in filt.columns:
        filt = filt[filt['Eye_status'] == 'eyeo']

    if filt.empty:
        print("Warning: no rows remain after reshape_power filters.")
        return pd.DataFrame()

    # Determine the band column name ('Freq_band' or 'band')
    band_col = 'Freq_band' if 'Freq_band' in filt.columns else 'band'

    # Aggregate PSD by median per subject/region/side/band
    group_cols = ['ID', 'Brain_region', 'Brain_side', band_col]
    psd_col = 'PSD'  # use raw absolute PSD (already combat-corrected in this file)
    agg = filt.groupby(group_cols, as_index=False).agg(PSD=(psd_col, 'median'))

    # Pivot to wide format
    piv = agg.pivot(index='ID', columns=['Brain_region', 'Brain_side', band_col], values='PSD')
    piv.columns = ['con_' + '_'.join(col) for col in piv.columns]
    piv = piv.reset_index()

    # Merge metadata back
    metadata_cols = ['ID', 'Cluster', 'PopulationS1', 'age_yrs', 'Sex', 'control_status']
    metadata_cols = [c for c in metadata_cols if c in filt.columns]
    metadata = filt.groupby('ID')[metadata_cols[1:]].first().reset_index()
    piv = piv.merge(metadata, on='ID', how='left')

    print(f"Reshaped power: {len(piv)} subjects × {len(piv.columns)} columns")
    return piv


# =============================================================================
# REGRESSION & Z-SCORE
# =============================================================================

def regress_and_zscore_features(df: pd.DataFrame, features: list) -> pd.DataFrame:
    """For each feature, regress on age and age^2, then z-score residuals using TD group.

    TD group is identified by control_status containing 'control' (case-insensitive).
    Adds columns named '{feature}_corrected'.

    Parameters
    ----------
    df       : Wide-format DataFrame with age_yrs, control_status, and feature columns.
    features : List of feature column names (e.g., con_Central_Left_alpha).

    Returns
    -------
    df with additional '{feature}_corrected' columns.
    """
    df = df.copy()
    df['_age_sq'] = df['age_yrs'] ** 2

    # TD mask for z-scoring reference
    if 'control_status' in df.columns:
        td_mask = df['control_status'].str.contains('control', case=False, na=False)
    else:
        td_mask = df['PopulationS1'] == 'NT'

    for feat in features:
        valid_mask = ~df[feat].isna()
        if valid_mask.sum() < 10:
            continue

        X = df.loc[valid_mask, ['age_yrs', '_age_sq']]
        y = df.loc[valid_mask, feat]

        if y.isna().all() or X.isna().any().any():
            continue

        model = LinearRegression().fit(X, y)
        residuals = y - model.predict(X)

        # Store residuals temporarily to compute TD statistics
        df.loc[valid_mask, f'_res_{feat}'] = residuals

        control_residuals = df.loc[valid_mask & td_mask, f'_res_{feat}'].dropna()
        if len(control_residuals) < 5:
            df.drop(columns=[f'_res_{feat}'], inplace=True)
            continue

        mu = control_residuals.mean()
        sd = control_residuals.std()
        if sd == 0 or np.isnan(sd):
            df.drop(columns=[f'_res_{feat}'], inplace=True)
            continue

        df.loc[valid_mask, f'{feat}_corrected'] = (residuals - mu) / sd
        df.drop(columns=[f'_res_{feat}'], inplace=True)

    df.drop(columns=['_age_sq'], inplace=True)
    return df


# =============================================================================
# STATISTICS
# =============================================================================

def clusters_vs_td_stats(df: pd.DataFrame, corrected_cols: list,
                          order: tuple = ('C1', 'C2', 'C3')) -> pd.DataFrame:
    """Compute cluster vs TD statistics with global FDR correction.

    For each (feature, cluster) pair:
        - Shapiro normality test → t-test or Mann-Whitney
    Global FDR correction (BH) across all (feature × cluster) tests.

    Parameters
    ----------
    df             : Wide-format DataFrame with Cluster, PopulationS1, and feature columns.
    corrected_cols : List of z-scored feature column names.
    order          : Tuple of cluster names.

    Returns
    -------
    DataFrame with columns:
        feature, test_type, comparison, cluster_n, td_n, cluster_mean, td_mean,
        t_or_u, p_uncorrected, p_fdr, significant
    """
    autism_data = df[df['PopulationS1'] == 'Autism'].dropna(subset=['Cluster'])
    td_data = df[df['PopulationS1'] == 'NT']

    all_rows = []

    for feat in corrected_cols:
        td_values = td_data[feat].dropna()
        for cluster in order:
            cluster_values = autism_data[autism_data['Cluster'] == cluster][feat].dropna()
            if len(cluster_values) < 3 or len(td_values) < 3:
                continue

            cluster_normal = stats.shapiro(cluster_values)[1] > 0.05 if len(cluster_values) >= 3 else False
            td_normal = stats.shapiro(td_values)[1] > 0.05 if len(td_values) >= 3 else False
            test_type = 't-test' if (cluster_normal and td_normal) else 'Mann-Whitney'

            if test_type == 't-test':
                t_val, p = stats.ttest_ind(cluster_values, td_values, equal_var=False)
            else:
                t_val, p = stats.mannwhitneyu(cluster_values, td_values, alternative='two-sided')

            all_rows.append({
                'feature': feat,
                'test_type': test_type,
                'comparison': f'{cluster} vs NT',
                'cluster_n': int(len(cluster_values)),
                'td_n': int(len(td_values)),
                'cluster_mean': float(cluster_values.mean()),
                'td_mean': float(td_values.mean()),
                't_or_u': float(t_val),
                'p_uncorrected': float(p),
            })

    if not all_rows:
        return pd.DataFrame()

    stats_df = pd.DataFrame(all_rows)

    # Global FDR correction
    pvals = stats_df['p_uncorrected'].values
    _, pvals_fdr, _, _ = smm.multipletests(pvals, method='fdr_bh', alpha=0.05)
    stats_df['p_fdr'] = pvals_fdr
    stats_df['significant'] = pvals_fdr < 0.05

    return stats_df


# =============================================================================
# PLOTTING
# =============================================================================

def plot_power_bands_clusters(df: pd.DataFrame, band: str, corrected_cols: list,
                               stats_df: pd.DataFrame,
                               order: tuple = ('C1', 'C2', 'C3')) -> None:
    """Create a grid of violin plots for a given frequency band, one subplot per feature.

    Each subplot shows violin + strip + pointplot for all clusters and TD.
    Significant (FDR) cluster vs TD comparisons are marked with a '*' annotation.
    Saves to PLOTS_DIR/power_{band}_clusters_vs_td.pdf.

    Parameters
    ----------
    df            : Wide-format DataFrame.
    band          : Frequency band string (e.g., 'alpha').
    corrected_cols: List of z-scored feature column names for this band.
    stats_df      : Statistics DataFrame returned by clusters_vs_td_stats.
    order         : Tuple of cluster names.
    """
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    if not corrected_cols:
        print(f"No corrected columns for band {band}, skipping plot.")
        return

    autism_data = df[df['PopulationS1'] == 'Autism'].dropna(subset=['Cluster'])
    td_data = df[df['PopulationS1'] == 'NT'].copy()
    td_data['Cluster'] = 'NT'
    plot_data = pd.concat([autism_data, td_data], ignore_index=True)
    plot_order = list(order) + ['NT']

    n_features = len(corrected_cols)
    n_cols = min(3, n_features)
    n_rows = int(np.ceil(n_features / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(5 * n_cols, 6 * n_rows),
                              squeeze=False)
    fig.suptitle(f'EEG {band.capitalize()} Power — Autism Clusters vs NT', fontsize=14, y=1.01)

    for ax_idx, feat in enumerate(corrected_cols):
        row_idx = ax_idx // n_cols
        col_idx = ax_idx % n_cols
        ax = axes[row_idx][col_idx]

        feat_data = plot_data[['Cluster', feat]].copy().dropna(subset=[feat])

        sns.violinplot(
            y=feat, x='Cluster', data=feat_data,
            inner=None, linewidth=0, order=plot_order,
            palette=CLUSTER_COLORS, ax=ax,
        )
        sns.stripplot(
            y=feat, x='Cluster', data=feat_data,
            color='black', size=1.5, jitter=True, alpha=0.4,
            order=plot_order, ax=ax,
        )
        sns.pointplot(
            y=feat, x='Cluster', data=feat_data,
            color='#8A0201', join=False, errorbar=('ci', 95),
            scale=0.8, order=plot_order, ax=ax,
        )

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.axhline(0, color='grey', linewidth=0.8, linestyle='--', alpha=0.5)

        # Derive a readable title from feature name: con_Frontal_Left_alpha → Frontal Left
        parts = feat.replace('_corrected', '').split('_')
        roi_label = ' '.join(parts[1:-1]) if len(parts) > 2 else feat
        ax.set_title(roi_label, fontsize=9)
        ax.set_xlabel('')
        ax.set_ylabel('z-score', fontsize=8)
        ax.tick_params(axis='x', labelsize=7, rotation=30)
        ax.tick_params(axis='y', labelsize=8)

        # Mark significant clusters with '*'
        if stats_df is not None and not stats_df.empty:
            feat_stats = stats_df[
                (stats_df['feature'] == feat) &
                (stats_df['significant'] == True) &
                (stats_df['comparison'].str.contains(' vs NT', na=False))
            ]
            pos_map = {name: i for i, name in enumerate(plot_order)}
            y_vals = feat_data[feat].dropna()
            if not y_vals.empty:
                y_annot = y_vals.max() + (y_vals.max() - y_vals.min()) * 0.05
                for _, srow in feat_stats.iterrows():
                    cluster_label = srow['comparison'].split(' vs NT')[0]
                    if cluster_label in pos_map:
                        ax.text(
                            pos_map[cluster_label], y_annot, '*',
                            ha='center', va='bottom', fontsize=12, color='black',
                        )

    # Hide unused subplots
    for ax_idx in range(n_features, n_rows * n_cols):
        row_idx = ax_idx // n_cols
        col_idx = ax_idx % n_cols
        axes[row_idx][col_idx].set_visible(False)

    fig.tight_layout()
    out_path = PLOTS_DIR / f'power_{band}_clusters_vs_nt.pdf'
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


# =============================================================================
# SYNTHESIS HEATMAP
# =============================================================================

def create_synthesis_heatmap(bands: list, order: tuple = ('C1', 'C2', 'C3')) -> None:
    """Read per-band stats CSVs and create a heatmap of % significant features.

    Reads STATS_DIR/power_{band}_clusters_vs_td_stats.csv for each band.
    Saves heatmap to FIG_DIR/power_synthesis_heatmap.pdf and summary CSV to
    FDR_DIR/power_synthesis_summary.csv.
    """
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    FDR_DIR.mkdir(parents=True, exist_ok=True)

    all_band_stats = []
    for band in bands:
        stats_file = STATS_DIR / f'power_{band}_clusters_vs_nt_stats.csv'
        if stats_file.exists():
            band_stats = pd.read_csv(stats_file)
            band_stats['band'] = band
            all_band_stats.append(band_stats)

    if not all_band_stats:
        print("No per-band stats found for synthesis heatmap.")
        return

    combined = pd.concat(all_band_stats, ignore_index=True)

    summary_rows = []
    for band in bands:
        bd = combined[combined['band'] == band]
        if bd.empty:
            continue
        for cluster in order:
            cd = bd[bd['comparison'] == f'{cluster} vs NT']
            n_sig = int(cd['significant'].sum()) if 'significant' in cd.columns else 0
            n_total = len(cd)
            summary_rows.append({
                'band': band,
                'cluster': cluster,
                'n_significant_features': n_sig,
                'n_total_features': n_total,
                'percentage_significant': (n_sig / n_total * 100) if n_total > 0 else 0.0,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = FDR_DIR / 'power_synthesis_summary.csv'
    summary_df.to_csv(summary_csv, index=False)
    print(f"Synthesis summary saved: {summary_csv}")

    if summary_df.empty:
        return

    heatmap_data = summary_df.pivot(index='band', columns='cluster',
                                     values='percentage_significant')
    # Reorder columns to match CLUSTER_ORDER
    heatmap_cols = [c for c in order if c in heatmap_data.columns]
    heatmap_data = heatmap_data[heatmap_cols]

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(
        heatmap_data, annot=True, fmt='.1f', cmap='Reds',
        cbar_kws={'label': '% Significant Features (FDR)'},
        ax=ax, linewidths=0.5,
    )
    ax.set_title('EEG Power: % Significant Features\n(Autism Clusters vs NT, FDR-corrected)', fontsize=12)
    ax.set_xlabel('Autism Clusters', fontsize=11)
    ax.set_ylabel('Frequency Band', fontsize=11)
    fig.tight_layout()
    heatmap_out = FIG_DIR / 'power_synthesis_heatmap.pdf'
    fig.savefig(heatmap_out, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Synthesis heatmap saved: {heatmap_out}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:
    # Create output directories
    for d in [FIG_DIR, STATS_DIR, FDR_DIR, PLOTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # 1. Load and merge power data
    power = load_power_with_clusters()
    print(f"Total rows loaded: {len(power)}")

    # 2. Reshape to wide format
    power_wide = reshape_power(power)
    if power_wide.empty:
        print("ERROR: reshape_power returned empty DataFrame. Aborting.")
        return 1

    # 3. Identify all band-feature columns
    all_con_cols = [c for c in power_wide.columns if c.startswith('con_')]
    if not all_con_cols:
        print("ERROR: no power features found after reshape. Aborting.")
        return 1

    bands = sorted({c.split('_')[-1] for c in all_con_cols})
    print(f"Frequency bands detected: {bands}")

    all_band_stats = []

    for band in bands:
        print(f"\n{'='*60}")
        print(f"Processing band: {band.upper()}")
        print(f"{'='*60}")

        band_feats = [c for c in all_con_cols if c.endswith(f'_{band}')]
        if not band_feats:
            print(f"No features for band {band}, skipping.")
            continue

        # 4. Regress and z-score features
        df_band = regress_and_zscore_features(power_wide, band_feats)
        corrected_cols = [f'{c}_corrected' for c in band_feats
                          if f'{c}_corrected' in df_band.columns]
        if not corrected_cols:
            print(f"No corrected features for band {band} after regression, skipping.")
            continue

        print(f"  {len(corrected_cols)} corrected features for band {band}")

        # 5. Statistics with global FDR correction
        stats_df = clusters_vs_td_stats(df_band, corrected_cols)
        if stats_df.empty:
            print(f"  No stats computed for band {band}.")
            continue

        n_sig = int(stats_df['significant'].sum())
        n_total = len(stats_df)
        print(f"  {n_sig}/{n_total} comparisons significant (FDR-corrected p < 0.05)")

        # Save per-band stats
        stats_out = STATS_DIR / f'power_{band}_clusters_vs_nt_stats.csv'
        stats_df.to_csv(stats_out, index=False)
        print(f"  Stats saved: {stats_out}")

        # Save FDR summary for this band
        fdr_out = FDR_DIR / f'power_{band}_fdr_corrected.csv'
        stats_df.to_csv(fdr_out, index=False)
        print(f"  FDR results saved: {fdr_out}")

        all_band_stats.append(stats_df.assign(band=band))

        # 6. Plots
        plot_power_bands_clusters(df_band, band, corrected_cols, stats_df)

    # 7. Global FDR re-correction across all bands (optional synthesis)
    if all_band_stats:
        combined_stats = pd.concat(all_band_stats, ignore_index=True)
        all_pvals = combined_stats['p_uncorrected'].values
        _, pvals_global_fdr, _, _ = smm.multipletests(all_pvals, method='fdr_bh', alpha=0.05)
        combined_stats['p_fdr_global'] = pvals_global_fdr
        combined_stats['significant_global'] = pvals_global_fdr < 0.05
        global_out = FDR_DIR / 'power_all_bands_global_fdr.csv'
        combined_stats.to_csv(global_out, index=False)
        n_global_sig = int(combined_stats['significant_global'].sum())
        print(f"\nGlobal FDR (all bands): {n_global_sig}/{len(combined_stats)} significant")
        print(f"Global FDR saved: {global_out}")

    # 8. Synthesis heatmap
    create_synthesis_heatmap(bands)

    print(f"\nAll outputs saved to: {FIG_DIR}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
