#!/usr/bin/env python3
"""
Cluster MRI Input Generation

Generates t-statistic CSV files comparing each EEG cluster (C1, C2, C3) against
pooled TD controls across cortical thickness, surface area, and subcortical volume
metrics. Also produces pairwise cluster comparisons.

Outputs written to: ../outputs/figures/r_input_files/
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--nt-c1", action="store_true",
    help="Restrict the NT control pool to NT subjects whose own Cluster == 'C1'. "
         "Mirrors the _nt_c1 variant used by 3_cluster_genetic_analysis scripts.",
)
_args, _ = _parser.parse_known_args()
NT_C1_ONLY = _args.nt_c1

# =============================================================================
# PATH CONSTANTS
# =============================================================================
_SCRIPT_DIR  = Path(__file__).parent
_LIB_DIR     = _SCRIPT_DIR.parent.parent.parent.parent  # parent of the repo (legacy/5_.../scripts -> up 4)
DATA_IMAG2GENET_DIR = str(_LIB_DIR / 'imaging2genet')
DF_CLUSTERS_FILE    = str(_LIB_DIR / 'eeg_mri-pipeline' / 'results' / 'dataset_paper' / 'dataframes' / 'df_clusters_complete_kmeans.csv')
MRI_FILE            = DATA_IMAG2GENET_DIR + '/0_input/dataframes/MRI_ANAT_INOVAND_LEAP_COMBAT.tsv'
_SECTION_DIR = _SCRIPT_DIR.parent
OUTPUT_DIR   = _SECTION_DIR / 'outputs' / 'figures' / 'r_input_files'

# =============================================================================
# ATLAS DICTIONARIES
# =============================================================================
data_dk = {
    'lh_bankssts_thickness': 'bankssts',
    'lh_caudalanteriorcingulate_thickness': 'caudalanteriorcingulate',
    'lh_caudalmiddlefrontal_thickness': 'caudalmiddlefrontal',
    'lh_cuneus_thickness': 'cuneus',
    'lh_entorhinal_thickness': 'entorhinal',
    'lh_fusiform_thickness': 'fusiform',
    'lh_inferiorparietal_thickness': 'inferiorparietal',
    'lh_inferiortemporal_thickness': 'inferiortemporal',
    'lh_isthmuscingulate_thickness': 'isthmuscingulate',
    'lh_lateraloccipital_thickness': 'lateraloccipital',
    'lh_lateralorbitofrontal_thickness': 'lateralorbitofrontal',
    'lh_lingual_thickness': 'lingual',
    'lh_medialorbitofrontal_thickness': 'medialorbitofrontal',
    'lh_middletemporal_thickness': 'middletemporal',
    'lh_parahippocampal_thickness': 'parahippocampal',
    'lh_paracentral_thickness': 'paracentral',
    'lh_parsopercularis_thickness': 'parsopercularis',
    'lh_parsorbitalis_thickness': 'parsorbitalis',
    'lh_parstriangularis_thickness': 'parstriangularis',
    'lh_pericalcarine_thickness': 'pericalcarine',
    'lh_postcentral_thickness': 'postcentral',
    'lh_posteriorcingulate_thickness': 'posteriorcingulate',
    'lh_precentral_thickness': 'precentral',
    'lh_precuneus_thickness': 'precuneus',
    'lh_rostralanteriorcingulate_thickness': 'rostralanteriorcingulate',
    'lh_rostralmiddlefrontal_thickness': 'rostralmiddlefrontal',
    'lh_superiorfrontal_thickness': 'superiorfrontal',
    'lh_superiorparietal_thickness': 'superiorparietal',
    'lh_superiortemporal_thickness': 'superiortemporal',
    'lh_supramarginal_thickness': 'supramarginal',
    'lh_frontalpole_thickness': 'frontalpole',
    'lh_temporalpole_thickness': 'temporalpole',
    'lh_transversetemporal_thickness': 'transversetemporal',
    'lh_insula_thickness': 'insula',
    'rh_bankssts_thickness': 'bankssts',
    'rh_caudalanteriorcingulate_thickness': 'caudalanteriorcingulate',
    'rh_caudalmiddlefrontal_thickness': 'caudalmiddlefrontal',
    'rh_cuneus_thickness': 'cuneus',
    'rh_entorhinal_thickness': 'entorhinal',
    'rh_fusiform_thickness': 'fusiform',
    'rh_inferiorparietal_thickness': 'inferiorparietal',
    'rh_inferiortemporal_thickness': 'inferiortemporal',
    'rh_isthmuscingulate_thickness': 'isthmuscingulate',
    'rh_lateraloccipital_thickness': 'lateraloccipital',
    'rh_lateralorbitofrontal_thickness': 'lateralorbitofrontal',
    'rh_lingual_thickness': 'lingual',
    'rh_medialorbitofrontal_thickness': 'medialorbitofrontal',
    'rh_middletemporal_thickness': 'middletemporal',
    'rh_parahippocampal_thickness': 'parahippocampal',
    'rh_paracentral_thickness': 'paracentral',
    'rh_parsopercularis_thickness': 'parsopercularis',
    'rh_parsorbitalis_thickness': 'parsorbitalis',
    'rh_parstriangularis_thickness': 'parstriangularis',
    'rh_pericalcarine_thickness': 'pericalcarine',
    'rh_postcentral_thickness': 'postcentral',
    'rh_posteriorcingulate_thickness': 'posteriorcingulate',
    'rh_precentral_thickness': 'precentral',
    'rh_precuneus_thickness': 'precuneus',
    'rh_rostralanteriorcingulate_thickness': 'rostralanteriorcingulate',
    'rh_rostralmiddlefrontal_thickness': 'rostralmiddlefrontal',
    'rh_superiorfrontal_thickness': 'superiorfrontal',
    'rh_superiorparietal_thickness': 'superiorparietal',
    'rh_superiortemporal_thickness': 'superiortemporal',
    'rh_supramarginal_thickness': 'supramarginal',
    'rh_frontalpole_thickness': 'frontalpole',
    'rh_temporalpole_thickness': 'temporalpole',
    'rh_transversetemporal_thickness': 'transversetemporal',
    'rh_insula_thickness': 'insula',
}

data_aseg = {
    'Left-Lateral-Ventricle': 'Left-Lateral-Ventricle',
    'Left-Thalamus': 'Left-Thalamus',
    'Left-Caudate': 'Left-Caudate',
    'Left-Putamen': 'Left-Putamen',
    'Left-Pallidum': 'Left-Pallidum',
    '3rd-Ventricle': '3rd-Ventricle',
    '4th-Ventricle': '4th-Ventricle',
    'Brain-Stem': 'Brain-Stem',
    'Left-Hippocampus': 'Left-Hippocampus',
    'Left-Amygdala': 'Left-Amygdala',
    'Left-VentralDC': 'Left-VentralDC',
    'Right-Lateral-Ventricle': 'Right-Lateral-Ventricle',
    'Right-Thalamus': 'Right-Thalamus',
    'Right-Caudate': 'Right-Caudate',
    'Right-Putamen': 'Right-Putamen',
    'Right-Pallidum': 'Right-Pallidum',
    'Right-Hippocampus': 'Right-Hippocampus',
    'Right-Amygdala': 'Right-Amygdala',
    'Right-VentralDC': 'Right-VentralDC',
    'CC_Posterior': 'CC_Posterior',
    'CC_Mid_Posterior': 'CC_Mid_Posterior',
    'CC_Central': 'CC_Central',
    'CC_Mid_Anterior': 'CC_Mid_Anterior',
    'CC_Anterior': 'CC_Anterior',
    'Right-Cerebellum-Cortex': 'Right-Cerebellum-Cortex',
    'Right-Cerebellum-White-Matter': 'Right-Cerebellum-White-Matter',
    'Left-Cerebellum-Cortex': 'Left-Cerebellum-Cortex',
    'Left-Cerebellum-White-Matter': 'Left-Cerebellum-White-Matter',
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_and_prepare_data():
    """Load cluster assignments and MRI data, merge, apply QC, return df_merged."""
    print("Loading cluster assignments...")
    df_clusters = pd.read_csv(DF_CLUSTERS_FILE)
    df_clusters['ID'] = df_clusters['ID'].astype(str)

    print("Loading MRI data...")
    df_mri = pd.read_csv(MRI_FILE, sep='\t')

    # Fill ID from MRI_ID if ID column is missing or all NaN
    if 'ID' not in df_mri.columns or df_mri['ID'].isna().all():
        df_mri['ID'] = df_mri['MRI_ID']
    df_mri['ID'] = df_mri['ID'].astype(str)

    # Drop conflicting columns from MRI to avoid suffix conflicts
    conflict_cols = ['control_status', 'PopulationS1', 'Population1', 'Cluster']
    df_mri = df_mri.drop(columns=[c for c in conflict_cols if c in df_mri.columns])

    print("Merging datasets on ID (inner join)...")
    cluster_cols = ['ID', 'Cluster', 'PopulationS1', 'Population1']
    available_cluster_cols = [c for c in cluster_cols if c in df_clusters.columns]
    df_merged = pd.merge(
        df_clusters[available_cluster_cols],
        df_mri,
        on='ID',
        how='inner'
    )
    print(f"  Merged shape: {df_merged.shape}")

    # QC is handled upstream (during FreeSurfer ComBat/regression z-scoring) —
    # the per-script LEAP-only blacklist filter has been removed.

    # Keep clustered subjects AND the NT pool (autism-only clustering mode
    # leaves NT subjects without a cluster label but they're still needed
    # for cluster-vs-NT statistics).
    df_merged = df_merged[df_merged['Cluster'].notna() |
                          (df_merged['PopulationS1'].isin(['NT', 'TD']))]
    print(f"  Final shape after dropping NaN clusters (NT kept): {df_merged.shape}")

    # Rename population label TD -> NT for display
    df_merged['PopulationS1'] = df_merged['PopulationS1'].replace('TD', 'NT')

    return df_merged


# =============================================================================
# COLUMN NAME TRANSFORMS
# =============================================================================

def transform_name(name, metric):
    """
    Convert a column name to a region label suitable for ggseg:
      - lh_/rh_ prefix is retained; _{metric} suffix is stripped.
      - Legacy _left/_right suffix style is converted to lh_/rh_ prefix.
      - All other names are returned unchanged.
    """
    if name.endswith('_left'):
        base = name[: -len('_left')]
        base = base.replace(f'_{metric}', '')
        return f'lh_{base}'
    if name.endswith('_right'):
        base = name[: -len('_right')]
        base = base.replace(f'_{metric}', '')
        return f'rh_{base}'
    if name.startswith('lh_') or name.startswith('rh_'):
        suffix = f'_{metric}'
        if name.endswith(suffix):
            return name[: -len(suffix)]
        return name
    return name


def transform_column_names(df, metric):
    """Rename all columns in df using transform_name."""
    return df.rename(columns={col: transform_name(col, metric) for col in df.columns})


# =============================================================================
# STATISTICAL HELPERS
# =============================================================================

def _cohens_d(a, b):
    """Pooled-std Cohen's d (a relative to b)."""
    n_a, n_b = len(a), len(b)
    if n_a < 2 or n_b < 2:
        return np.nan
    pooled_std = np.sqrt(((n_a - 1) * np.var(a, ddof=1) + (n_b - 1) * np.var(b, ddof=1)) / (n_a + n_b - 2))
    if pooled_std == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled_std


def _cohens_d_ci(d, n_a, n_b, alpha=0.05):
    """Analytic 95 % CI for Cohen's d (Hedges & Olkin 1985 approximation).

    SE(d) = sqrt((n1+n2) / (n1*n2)  +  d² / (2*(n1+n2)))
    Returns (lo, hi). Returns (nan, nan) if d is nan or n insufficient.
    """
    if not np.isfinite(d) or n_a < 2 or n_b < 2:
        return (np.nan, np.nan)
    se = np.sqrt((n_a + n_b) / (n_a * n_b) + d ** 2 / (2.0 * (n_a + n_b)))
    z = 1.959963984540054  # qnorm(0.975)
    return (float(d - z * se), float(d + z * se))


def _run_tests(autism_data, td_data, columns, label_map):
    """
    Run Welch's t-tests for each column, compute Cohen's d (+ 95 % CI), apply FDR.

    Returns a DataFrame with columns: label, t_stat, p_val,
    cohens_d, cohens_d_lo, cohens_d_hi, n_a, n_b, p_fdr.
    """
    rows = []
    for col in columns:
        if col not in autism_data.columns or col not in td_data.columns:
            continue
        a = autism_data[col].dropna().values
        b = td_data[col].dropna().values
        if len(a) < 2 or len(b) < 2:
            continue
        t_stat, p_val = ttest_ind(a, b, equal_var=False)
        d = _cohens_d(a, b)
        d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
        label = label_map.get(col, col)
        rows.append({
            'label': label,
            't_stat': t_stat,
            'p_val': p_val,
            'cohens_d': d,
            'cohens_d_lo': d_lo,
            'cohens_d_hi': d_hi,
            'n_a': int(len(a)),
            'n_b': int(len(b)),
        })

    if not rows:
        return pd.DataFrame(columns=['label', 't_stat', 'p_val',
                                     'cohens_d', 'cohens_d_lo', 'cohens_d_hi',
                                     'n_a', 'n_b', 'p_fdr'])

    result_df = pd.DataFrame(rows)
    _, p_fdr, _, _ = multipletests(result_df['p_val'].values, method='fdr_bh')
    result_df['p_fdr'] = p_fdr
    return result_df


# =============================================================================
# CLUSTER vs TD FILES
# =============================================================================

def create_cluster_vs_td_files(df_merged, output_dir):
    """
    For each cluster (C1, C2, C3) and each metric (thickness, area, volume):
      - Compare autism subjects in that cluster against all TD subjects (pooled).
      - Save t-stat CSV.
    Also saves cluster_counts_mri.csv.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    clusters = ['C1', 'C2', 'C3']

    # Identify metric column sets
    thickness_cols = [c for c in df_merged.columns
                      if c.endswith('_thickness') and (c.startswith('lh_') or c.startswith('rh_'))]
    area_cols      = [c for c in df_merged.columns
                      if c.endswith('_area') and (c.startswith('lh_') or c.startswith('rh_'))]
    volume_cols    = [c for c in df_merged.columns if c in data_aseg]

    metrics = [
        ('thickness', 'dk',   thickness_cols, data_dk),
        ('area',      'dk',   area_cols,      data_dk),
        ('volume',    'aseg', volume_cols,    data_aseg),
    ]

    td_data = df_merged[df_merged['PopulationS1'] == 'NT']
    if NT_C1_ONLY:
        before = len(td_data)
        td_data = td_data[td_data['Cluster'] == 'C1']
        print(f"  [--nt-c1] NT pool restricted to NT with Cluster == 'C1' "
              f"({before} → {len(td_data)})")
    print(f"  NT pool size: {len(td_data)}")

    # Cluster counts
    counts = []
    for cluster in clusters:
        cluster_df = df_merged[df_merged['Cluster'] == cluster]
        autism_df  = cluster_df[cluster_df['PopulationS1'] == 'Autism']
        counts.append({'Cluster': cluster, 'N_autism': len(autism_df), 'N_total': len(cluster_df)})
    counts_df = pd.DataFrame(counts)
    counts_df.to_csv(output_dir / 'cluster_counts_mri.csv', index=False)
    print(f"  Saved: cluster_counts_mri.csv")

    for cluster in clusters:
        cluster_df  = df_merged[df_merged['Cluster'] == cluster]
        autism_data = cluster_df[cluster_df['PopulationS1'] == 'Autism']
        print(f"\n  Cluster {cluster}: {len(autism_data)} autism subjects")

        for metric_short, atlas_type, cols, label_map in metrics:
            if not cols:
                print(f"    No {metric_short} columns found — skipping")
                continue

            # Build label map restricted to columns that exist in data
            # For area, reuse data_dk but strip _thickness from values (labels are region names, not col names)
            if metric_short == 'area':
                area_label_map = {
                    col: data_dk.get(col.replace('_area', '_thickness'), col.replace('_area', '').lstrip('lh_').lstrip('rh_'))
                    for col in cols
                }
                # simpler: strip metric suffix from col name
                area_label_map = {col: transform_name(col, metric_short) for col in cols}
                result_df = _run_tests(autism_data, td_data, cols, area_label_map)
            elif metric_short == 'thickness':
                thickness_label_map = {col: transform_name(col, metric_short) for col in cols}
                result_df = _run_tests(autism_data, td_data, cols, thickness_label_map)
            else:
                result_df = _run_tests(autism_data, td_data, cols, label_map)

            if result_df.empty:
                print(f"    No results for {metric_short} — skipping")
                continue

            fname = f't_stat_cluster_{cluster}_{metric_short}_{atlas_type}_mri_cluster_vs_td.csv'
            result_df.to_csv(output_dir / fname, index=False)
            sig = (result_df['p_fdr'] < 0.05).sum()
            print(f"    Saved: {fname}  (n_regions={len(result_df)}, n_fdr_sig={sig})")


# =============================================================================
# PAIRWISE CLUSTER FILES
# =============================================================================

def create_pairwise_cluster_files(df_merged, output_dir):
    """
    Pairwise t-tests between autism subjects in each pair of clusters
    (C1vsC2, C1vsC3, C2vsC3) for all metrics, with FDR correction.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs = [('C1', 'C2'), ('C1', 'C3'), ('C2', 'C3')]

    thickness_cols = [c for c in df_merged.columns
                      if c.endswith('_thickness') and (c.startswith('lh_') or c.startswith('rh_'))]
    area_cols      = [c for c in df_merged.columns
                      if c.endswith('_area') and (c.startswith('lh_') or c.startswith('rh_'))]
    volume_cols    = [c for c in df_merged.columns if c in data_aseg]

    metrics = [
        ('thickness', 'dk',   thickness_cols, {col: transform_name(col, 'thickness') for col in thickness_cols}),
        ('area',      'dk',   area_cols,      {col: transform_name(col, 'area')      for col in area_cols}),
        ('volume',    'aseg', volume_cols,    data_aseg),
    ]

    for c1, c2 in pairs:
        autism_c1 = df_merged[(df_merged['Cluster'] == c1) & (df_merged['PopulationS1'] == 'Autism')]
        autism_c2 = df_merged[(df_merged['Cluster'] == c2) & (df_merged['PopulationS1'] == 'Autism')]
        print(f"\n  Pairwise {c1} vs {c2}: n={len(autism_c1)} vs n={len(autism_c2)}")

        for metric_short, atlas, cols, label_map in metrics:
            if not cols:
                continue
            result_df = _run_tests(autism_c1, autism_c2, cols, label_map)
            if result_df.empty:
                continue
            fname = f'pairwise_{c1}_vs_{c2}_{metric_short}_{atlas}.csv'
            result_df.to_csv(output_dir / fname, index=False)
            sig = (result_df['p_fdr'] < 0.05).sum()
            print(f"    Saved: {fname}  (n_regions={len(result_df)}, n_fdr_sig={sig})")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("Cluster MRI Input Generation")
    print("=" * 60)

    print("\n[1/3] Loading and preparing data...")
    df_merged = load_and_prepare_data()

    print(f"\n[2/3] Creating cluster vs NT comparison files...")
    print(f"      Output directory: {OUTPUT_DIR}")
    create_cluster_vs_td_files(df_merged, OUTPUT_DIR)

    print(f"\n[3/3] Creating pairwise cluster comparison files...")
    create_pairwise_cluster_files(df_merged, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Cluster MRI input generation complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
