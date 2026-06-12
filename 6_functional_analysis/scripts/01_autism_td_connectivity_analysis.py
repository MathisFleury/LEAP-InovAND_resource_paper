#!/usr/bin/env python3
"""
Autism vs TD fMRI Connectivity Analysis

Compares functional connectivity patterns between Autism and TD groups:
- Welch's t-tests per connectivity pair
- FDR correction
- Hyper/hypo connectivity separation
- Network-based analysis using Schaefer atlas

Adapted from eeg_mri-pipeline/analysis/figures_papers/fmri_autism_td_analysis/autism_td_connectivity_analysis.py
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_SECTION_DIR = os.path.dirname(_SCRIPT_DIR)  # 4_functional_analysis/

EEG_MRI_RESULTS = '/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results'
DF_CLUSTERS_FILE = os.path.join(EEG_MRI_RESULTS, 'dataset_paper', 'dataframes', 'df_clusters_complete_kmeans.csv')
ATLAS_FILE = '/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv'

# Default: read from the refactored preprocessing pipeline, write under
# outputs/figures/revised/. Override either path via env var to compare
# pipelines without overwriting prior runs:
#   AUTISM_TD_FMRI_CSV   — input connectivity CSV
#   AUTISM_TD_OUTPUT_DIR — output folder for CSVs + intermediate tables
FMRI_FILE = os.environ.get(
    'AUTISM_TD_FMRI_CSV',
    os.path.join(_SECTION_DIR, 'preprocessing', 'outputs', 'df_conn_cohort_norm.csv'),
)
OUTPUT_DIR = os.environ.get(
    'AUTISM_TD_OUTPUT_DIR',
    os.path.join(_SECTION_DIR, 'outputs', 'figures', 'revised'),
)
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_and_prepare_data():
    print("Loading fMRI connectivity data...")

    if not os.path.exists(DF_CLUSTERS_FILE):
        raise FileNotFoundError(f"Cluster file not found: {DF_CLUSTERS_FILE}")
    df_clusters = pd.read_csv(DF_CLUSTERS_FILE)
    print(f"Loaded cluster data: {len(df_clusters)} participants")

    if not os.path.exists(FMRI_FILE):
        raise FileNotFoundError(f"fMRI connectivity file not found: {FMRI_FILE}")
    df_fmri = pd.read_csv(FMRI_FILE)
    print(f"Loaded fMRI connectivity data: {len(df_fmri)} participants")

    df_fmri = df_fmri.drop(columns=['PopulationS1', 'Population1'], errors='ignore')
    df_fmri['ID'] = df_fmri['ID'].astype(str)
    df_clusters['ID'] = df_clusters['ID'].astype(str)
    df_merged = df_fmri.merge(df_clusters[['ID', 'Cluster', 'PopulationS1', 'Population1']], on='ID', how='inner')

    not_merged = df_fmri[~df_fmri['ID'].isin(set(df_merged['ID']))]
    if len(not_merged) > 0:
        print(f"Warning: {len(not_merged)} participants not merged due to missing cluster/population info")

    # Rename population label TD -> NT for display
    df_merged['PopulationS1'] = df_merged['PopulationS1'].replace('TD', 'NT')

    df_autism_td = df_merged[df_merged['PopulationS1'].isin(['Autism', 'NT'])].copy()
    print(f"Final Autism vs NT dataset: {len(df_autism_td)} participants")
    print(f"Populations: {df_autism_td['PopulationS1'].value_counts().to_dict()}")
    return df_autism_td


def get_connectivity_columns(df):
    conn_cols = [col for col in df.columns if '/' in col and col.startswith('con_')]
    print(f"Found {len(conn_cols)} connectivity features before filtering")
    conn_cols_filtered = [col for col in conn_cols if 'Cerebellar' not in col or 'con_Cerebellar_Region1' in col]
    conn_cols_filtered = [col for col in conn_cols_filtered if 'Cerebellar_Region10' not in col]
    print(f"After cerebellar filtering: {len(conn_cols_filtered)} connectivity features")
    return conn_cols_filtered


def perform_autism_vs_td_analysis(df_autism_td, conn_cols):
    print("\nPerforming Autism vs NT analysis...")

    df_autism = df_autism_td[df_autism_td['PopulationS1'] == 'Autism'].copy()
    df_td = df_autism_td[df_autism_td['PopulationS1'] == 'NT'].copy()
    print(f"Autism group: {len(df_autism)} participants")
    print(f"NT group: {len(df_td)} participants")

    if len(df_autism) < 5 or len(df_td) < 5:
        print(f"Warning: One or both groups have < 5 participants")
        return None

    t_results = {}
    for feature in conn_cols:
        try:
            autism_data = df_autism[feature].dropna()
            td_data = df_td[feature].dropna()
            if len(autism_data) < 3 or len(td_data) < 3:
                continue
            t_stat, p_val = ttest_ind(autism_data, td_data, equal_var=False, nan_policy='omit')
            t_results[feature] = {'t_stat': t_stat, 'p_value': p_val, 'feature': feature}
        except Exception as e:
            continue

    df_t_results = pd.DataFrame(t_results).T
    if len(df_t_results) == 0:
        print("No valid results found")
        return None

    df_t_results['p_value_fdr'] = multipletests(df_t_results['p_value'], method='fdr_bh')[1]
    df_sig = df_t_results[df_t_results['p_value_fdr'] < 0.05].copy()

    if len(df_sig) == 0:
        print("No significant results with FDR < 0.05, trying uncorrected p < 0.05...")
        df_sig = df_t_results[df_t_results['p_value'] < 0.05].copy()
        print(f"Significant connections (uncorrected p < 0.05): {len(df_sig)}")

    df_hyper = df_sig[df_sig['t_stat'] > 0].copy()
    df_hypo = df_sig[df_sig['t_stat'] < 0].copy()

    print(f"Final significant connections: {len(df_sig)}")
    print(f"Hyperconnectivity (Autism > NT): {len(df_hyper)}")
    print(f"Hypoconnectivity (Autism < NT): {len(df_hypo)}")

    return {'all_results': df_t_results, 'significant_results': df_sig, 'hyper_results': df_hyper, 'hypo_results': df_hypo}


def create_connectivity_matrices(results):
    if results is None:
        return
    df_sig = results['significant_results']
    if len(df_sig) == 0:
        print("No significant results, skipping matrix creation")
        return

    df_sig = df_sig.copy()
    df_sig['source'] = df_sig['feature'].str.split('/').str[0].str.replace('con_', '')
    df_sig['target'] = df_sig['feature'].str.split('/').str[1]

    df_hyper = df_sig[df_sig['t_stat'] > 0].copy()
    df_hypo = df_sig[df_sig['t_stat'] < 0].copy()

    for df_sub, label in [(df_hyper, 'autism_vs_td_hyperconnectivity'), (df_hypo, 'autism_vs_td_hypoconnectivity')]:
        # Always write a fresh file (empty or not) so downstream plot scripts
        # never read a stale CSV from a previous run with different QC.
        df_r_data = df_sub[['source', 'target', 't_stat', 'p_value', 'p_value_fdr']].copy()
        df_r_data.to_csv(os.path.join(OUTPUT_DIR, f"{label}_full_data.csv"), index=False)
        print(f"Saved: {label}_full_data.csv ({len(df_sub)} edges)")

        source_counts = (df_sub['source'].value_counts().reset_index()
                         if len(df_sub) > 0
                         else pd.DataFrame({'region': [], 't_stat': []}))
        source_counts.columns = ['region', 't_stat']
        source_counts['mean_value'] = source_counts['t_stat']
        source_counts.to_csv(os.path.join(OUTPUT_DIR, f"{label}_for_r.csv"), index=False)


def create_network_based_analysis(df_autism_td, conn_cols, results):
    print("\nCreating network-based analysis...")
    if results is None:
        return

    if not os.path.exists(ATLAS_FILE):
        print(f"Warning: Atlas file not found: {ATLAS_FILE}")
        return

    df_atlas = pd.read_csv(ATLAS_FILE, sep='\t')
    df_sig = results['significant_results']
    df_hyper = results['hyper_results']
    df_hypo = results['hypo_results']

    for df, suffix in [(df_sig, 'network_analysis'), (df_hyper, 'hyperconnectivity_network_analysis'), (df_hypo, 'hypoconnectivity_network_analysis')]:
        output_file = os.path.join(OUTPUT_DIR, f"autism_vs_td_{suffix}.csv")
        if len(df) > 0:
            df = df.copy()
            df['source'] = df['feature'].str.split('/').str[0].str.replace('con_', '')
            df['target'] = df['feature'].str.split('/').str[1]
            df_merged = df.merge(df_atlas[['label', 'network_label']], left_on='source', right_on='label', how='left')
            network_counts = df_merged.groupby('network_label').size().reset_index(name='count')
            network_counts = network_counts.sort_values('count', ascending=False)
        else:
            network_counts = pd.DataFrame({'network_label': [], 'count': []})
        network_counts.to_csv(output_file, index=False)
        print(f"Saved network analysis: {output_file} ({len(network_counts)} rows)")


def create_summary_statistics(results):
    if results is None:
        return None

    df_sig = results['significant_results']
    df_hyper = results['hyper_results']
    df_hypo = results['hypo_results']

    summary_data = {
        'comparison': 'Autism vs NT',
        'total_significant': len(df_sig),
        'hyperconnectivity': len(df_hyper),
        'hypoconnectivity': len(df_hypo),
        'hyperconnectivity_pct': len(df_hyper) / len(df_sig) * 100 if len(df_sig) > 0 else 0,
        'hypoconnectivity_pct': len(df_hypo) / len(df_sig) * 100 if len(df_sig) > 0 else 0
    }

    df_summary = pd.DataFrame([summary_data])
    output_file = os.path.join(OUTPUT_DIR, "autism_td_connectivity_analysis_summary.csv")
    df_summary.to_csv(output_file, index=False)
    print(f"Saved summary: {output_file}")

    print("\n" + "=" * 60)
    print("AUTISM VS NT CONNECTIVITY ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"Total significant connections: {summary_data['total_significant']}")
    print(f"Hyperconnectivity (Autism > NT): {summary_data['hyperconnectivity']} ({summary_data['hyperconnectivity_pct']:.1f}%)")
    print(f"Hypoconnectivity (Autism < NT): {summary_data['hypoconnectivity']} ({summary_data['hypoconnectivity_pct']:.1f}%)")

    return df_summary


def main():
    print("=" * 60)
    print("AUTISM VS NT CONNECTIVITY ANALYSIS")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}")

    df_autism_td = load_and_prepare_data()
    conn_cols = get_connectivity_columns(df_autism_td)
    results = perform_autism_vs_td_analysis(df_autism_td, conn_cols)
    create_connectivity_matrices(results)
    create_network_based_analysis(df_autism_td, conn_cols, results)
    create_summary_statistics(results)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print(f"Results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
