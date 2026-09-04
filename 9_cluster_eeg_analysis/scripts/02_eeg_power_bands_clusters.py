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
_PREP = _SECTION_DIR.parent / "8_eeg_analysis" / "preprocessing"  # shared preprocessing
sys.path.insert(0, str(_PREP))
from step_05_regress import regress  # noqa: E402
from step_06_zscore import zscore    # noqa: E402
from config import POWER_COMBINED    # noqa: E402
from curated_map import attach_curated_demographics  # noqa: E402

FIG_DIR = _SECTION_DIR / 'outputs' / 'figures'
TABLES_DIR = _SECTION_DIR / 'outputs' / 'tables'
STATS_DIR = TABLES_DIR / 'stats'
FDR_DIR = TABLES_DIR / 'fdr_corrected'
PLOTS_DIR = FIG_DIR / 'plots'

# =============================================================================
# COLORS / ORDER
# =============================================================================
CLUSTER_COLORS = {'C1': '#7A8B47', 'C2': '#ff9fa0', 'C3': '#e7ba52',
                  'NT': '#C1C2BC', 'IDD': '#D8A4CB'}
CLUSTER_ORDER = ['C1', 'C2', 'C3', 'NT', 'IDD']


# =============================================================================
# DATA LOADING
# =============================================================================

def load_power_with_clusters() -> pd.DataFrame:
    """Load and merge power spectrum data with cluster assignments.

    The combined power spectrum file (Power_spectrum_combined_corrected.csv)
    contains columns: ID, Subject_ID, cohort, Channel, Brain_region, Brain_side,
    Eye_status, Freq_band, Quant_status, PSD, PSD_corrected, control_status,
    age_yrs, Sex, PopulationS1, Population1.

    Cluster labels + demographics are the curated k-means values (attach_curated_demographics).
    Only Autism subjects with a Cluster assignment and TD subjects are retained.

    Returns
    -------
    Long-format DataFrame merged with cluster metadata.
    """
    power = pd.read_csv(POWER_COMBINED, low_memory=False)
    print(f"Loaded {len(power)} rows from power spectrum file")

    # Curated demographics + k-means Cluster (curated is the only clinical source)
    power = attach_curated_demographics(power)

    # Quality filter
    power = power.dropna(subset=['age_yrs'])

    # Population filter: Autism split by cluster (C1/C2/C3), plus NT (reference)
    # and IDD, each kept as its own category.
    autism_mask = (power['PopulationS1'] == 'Autism') & power['Cluster'].notna()
    td_mask = power['PopulationS1'] == 'NT'
    idd_mask = power['PopulationS1'] == 'IDD'
    power = power[autism_mask | td_mask | idd_mask].reset_index(drop=True)

    uniq = power[['ID', 'PopulationS1', 'Cluster']].drop_duplicates(subset=['ID'])
    autism_by_cluster = (uniq.loc[uniq['PopulationS1'] == 'Autism', 'Cluster']
                         .value_counts().to_dict())
    n_nt = int((uniq['PopulationS1'] == 'NT').sum())
    n_idd = int((uniq['PopulationS1'] == 'IDD').sum())
    print(f"Unique subjects — Autism by cluster: {autism_by_cluster} | "
          f"NT (reference): {n_nt} | IDD: {n_idd}")

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
    """Regress age/age²/sex out of each feature and z-score on the FULL sample.

    No NT / control reference — see preprocessing/ (mirrors the anat z-scoring
    pipeline, DX_COL=None). Adds columns named '{feature}_corrected'.

    Parameters
    ----------
    df       : Wide-format DataFrame with age_yrs, Sex, and feature columns.
    features : List of feature column names (e.g., con_Central_Left_alpha).

    Returns
    -------
    df with additional '{feature}_corrected' columns.
    """
    df = df.copy()
    corrected = [f'{f}_corrected' for f in features]
    for f in features:
        df[f'{f}_corrected'] = df[f]
    df = regress(df, corrected)
    df = zscore(df, corrected)
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
    idd_data = df[df['PopulationS1'] == 'IDD']

    def _test(a, b):
        a_normal = stats.shapiro(a)[1] > 0.05 if len(a) >= 3 else False
        b_normal = stats.shapiro(b)[1] > 0.05 if len(b) >= 3 else False
        if a_normal and b_normal:
            return 't-test', *stats.ttest_ind(a, b, equal_var=False)
        return 'Mann-Whitney', *stats.mannwhitneyu(a, b, alternative='two-sided')

    vs_nt_rows, pw_rows = [], []

    for feat in corrected_cols:
        td_values = td_data[feat].dropna()
        # ---- Group vs NT (clusters C1/C2/C3 and IDD) ----
        groups = [(c, autism_data[autism_data['Cluster'] == c][feat].dropna()) for c in order]
        groups.append(('IDD', idd_data[feat].dropna()))
        for label, vals in groups:
            if len(vals) < 3 or len(td_values) < 3:
                continue
            test_type, t_val, p = _test(vals, td_values)
            vs_nt_rows.append({
                'feature': feat, 'test_type': test_type, 'comparison': f'{label} vs NT',
                'cluster_n': int(len(vals)), 'td_n': int(len(td_values)),
                'cluster_mean': float(vals.mean()), 'td_mean': float(td_values.mean()),
                't_or_u': float(t_val), 'p_uncorrected': float(p),
            })
        # ---- Pairwise inter-cluster (C1 vs C2, C1 vs C3, C2 vs C3) ----
        for c1, c2 in combinations(order, 2):
            v1 = autism_data[autism_data['Cluster'] == c1][feat].dropna()
            v2 = autism_data[autism_data['Cluster'] == c2][feat].dropna()
            if len(v1) < 3 or len(v2) < 3:
                continue
            test_type, t_val, p = _test(v1, v2)
            pw_rows.append({
                'feature': feat, 'test_type': test_type, 'comparison': f'{c1} vs {c2}',
                'cluster_n': int(len(v1)), 'td_n': int(len(v2)),
                'cluster_mean': float(v1.mean()), 'td_mean': float(v2.mean()),
                't_or_u': float(t_val), 'p_uncorrected': float(p),
            })

    # FDR-BH within each family separately (vs-NT global, pairwise global).
    frames = []
    for rows in (vs_nt_rows, pw_rows):
        if rows:
            sdf = pd.DataFrame(rows)
            _, pf, _, _ = smm.multipletests(sdf['p_uncorrected'].values, method='fdr_bh', alpha=0.05)
            sdf['p_fdr'] = pf
            sdf['significant'] = pf < 0.05
            frames.append(sdf)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


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
    idd_data = df[df['PopulationS1'] == 'IDD'].copy()
    idd_data['Cluster'] = 'IDD'
    plot_data = pd.concat([autism_data, td_data, idd_data], ignore_index=True)
    plot_order = list(order) + ['NT', 'IDD']

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

        # Annotate significant results: '*' over each group for group-vs-NT,
        # brackets with p-value for pairwise inter-cluster comparisons.
        if stats_df is not None and not stats_df.empty:
            pos_map = {name: i for i, name in enumerate(plot_order)}
            y_vals = feat_data[feat].dropna()
            if not y_vals.empty:
                yr = (y_vals.max() - y_vals.min()) or 1.0
                # group vs NT stars
                feat_stats = stats_df[
                    (stats_df['feature'] == feat) & (stats_df['significant']) &
                    (stats_df['comparison'].str.contains(' vs NT', na=False))
                ]
                y_annot = y_vals.max() + yr * 0.05
                for _, srow in feat_stats.iterrows():
                    group_label = srow['comparison'].split(' vs NT')[0]
                    if group_label in pos_map:
                        ax.text(pos_map[group_label], y_annot, '*',
                                ha='center', va='bottom', fontsize=12, color='black')
                # pairwise inter-cluster brackets (significant only), stacked above
                feat_pw = stats_df[
                    (stats_df['feature'] == feat) & (stats_df['significant']) &
                    (~stats_df['comparison'].str.contains(' vs NT', na=False))
                ]
                ystep = yr * 0.09
                by = y_vals.max() + yr * 0.13
                for _, prow in feat_pw.iterrows():
                    c1, c2 = [s.strip() for s in prow['comparison'].split(' vs ')]
                    if c1 not in pos_map or c2 not in pos_map:
                        continue
                    x1, x2 = pos_map[c1], pos_map[c2]
                    ax.plot([x1, x1, x2, x2],
                            [by, by + ystep * 0.3, by + ystep * 0.3, by],
                            color='black', linewidth=1.0)
                    pf = prow['p_fdr']
                    ptext = f'{pf:.0e}' if pf < 0.001 else f'{pf:.3f}'
                    ax.text((x1 + x2) / 2, by + ystep * 0.35, ptext,
                            ha='center', va='bottom', fontsize=7, fontweight='bold')
                    by += ystep * 1.3

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
    print(f"k-means clusters vs NT  →  {FIG_DIR}")
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
