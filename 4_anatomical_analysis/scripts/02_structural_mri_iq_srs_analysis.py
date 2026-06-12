#!/usr/bin/env python3
"""
Structural MRI-IQ-SRS Correlation Analysis for Autistic Individuals

Correlates structural MRI features with IQ, SRS, and Vineland scores for autism individuals.

Adapted from eeg_mri-pipeline/analysis/figures_papers/mri_autism_td_analysis/structural_mri_iq_srs_analysis.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from pathlib import Path
from scipy import stats
from scipy.stats import pearsonr
import warnings

warnings.filterwarnings('ignore')

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_SECTION_DIR = os.path.dirname(_SCRIPT_DIR)  # 3_anatomical_analysis/

OUTPUT_DIR = os.path.join(_SECTION_DIR, 'outputs', 'figures')
R_INPUT_DIR = os.path.join(OUTPUT_DIR, 'r_input_files')
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
Path(R_INPUT_DIR).mkdir(parents=True, exist_ok=True)

EEG_MRI_RESULTS = '/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results'
DATA_IMAG2GENET_DIR = '/Users/mfleury/POSTDOC/LIBRAIRY/imaging2genet'
DF_CLUSTERS_FILE = os.path.join(EEG_MRI_RESULTS, 'dataset_paper', 'dataframes', 'df_clusters_complete_kmeans.csv')

POPULATION_COLORS = {
    'Autism': '#8991FA',
    'Relatives': '#9AD5D3',
    'IDD': '#D8A4CB',
    'TD': '#C1C2BC',
    'Autism with IDD': '#324095',
    'Autism without IDD': '#5CAEE1'
}


def load_and_prepare_data():
    print("Loading structural MRI dataset...")

    df_individuals = pd.read_csv(DF_CLUSTERS_FILE)
    df_individuals['ID'] = df_individuals['ID'].astype(str)

    df_autism = df_individuals[df_individuals['PopulationS1'] == 'Autism'].copy()

    df_autism['total_IQ'] = pd.to_numeric(df_autism['total_IQ'], errors='coerce').fillna(df_autism['performance_IQ'])
    df_autism['SRS_t_score'] = pd.to_numeric(df_autism['SRS_tscore'], errors='coerce')
    df_autism['age_yrs'] = pd.to_numeric(df_autism['age_yrs'], errors='coerce')
    df_autism['vabsabcabc_standard'] = pd.to_numeric(df_autism['vabsabcabc_standard'], errors='coerce')

    df_autism = df_autism[
        (df_autism['total_IQ'] != 999) &
        (df_autism['SRS_t_score'] != 999) &
        (df_autism['total_IQ'] > 0) &
        (df_autism['SRS_t_score'] > 0) &
        (df_autism['vabsabcabc_standard'] > 0) &
        (df_autism['vabsabcabc_standard'] != 999)
    ].copy()

    print(f"Loaded {len(df_autism)} autism individuals with valid IQ, SRS, and Vineland data")
    print(f"IQ range: {df_autism['total_IQ'].min():.0f} to {df_autism['total_IQ'].max():.0f}")
    print(f"SRS range: {df_autism['SRS_t_score'].min():.0f} to {df_autism['SRS_t_score'].max():.0f}")

    print("Loading structural MRI data...")
    df_mri = pd.read_csv(os.path.join(DATA_IMAG2GENET_DIR, '0_input', 'dataframes', 'MRI_ANAT_INOVAND_LEAP_COMBAT.tsv'), sep='\t')
    df_mri['ID'] = df_mri['ID'].fillna(df_mri['MRI_ID'])
    df_mri['ID'] = df_mri['ID'].astype(str)

    thickness_cols = [col for col in df_mri.columns if 'thickness' in col.lower()]
    area_cols = [col for col in df_mri.columns if 'area' in col.lower() and 'thickness' not in col.lower()]
    volume_cols = [col for col in df_mri.columns if 'volume' in col.lower()]
    subcortical_cols = [col for col in df_autism.columns if any(region in col.lower() for region in
        ['hippocampus', 'amygdala', 'accumbens', 'caudate', 'putamen', 'pallidum', 'thalamus',
         'ventraldc', 'choroid', 'cc_', 'cerebellum', 'ventricle', 'brain-stem'])
        and not any(exclude in col.lower() for exclude in ['pgs_', 'pegs_', 'brain_biton'])]

    print(f"Found {len(thickness_cols)} thickness, {len(area_cols)} area, {len(volume_cols)} volume, {len(subcortical_cols)} subcortical features")

    mri_columns_to_exclude = thickness_cols + area_cols + volume_cols + subcortical_cols
    df_autism = df_autism.drop(columns=mri_columns_to_exclude, errors='ignore')

    df_merged = df_autism.merge(df_mri, on='ID', how='inner')

    df_merged['total_IQ'] = df_merged['total_IQ_x']
    df_merged['age_yrs'] = df_merged['age_yrs_x']

    vineland_cols = ['vabsabcabc_standard']
    for col in vineland_cols:
        if col + '_x' in df_merged.columns:
            df_merged[col] = df_merged[col + '_x']

    df_merged.mri_features = {
        'thickness': thickness_cols,
        'area': area_cols,
        'volume': volume_cols,
        'subcortical': subcortical_cols
    }

    print(f"Merged dataset: {len(df_merged)} autism individuals with MRI data")
    return df_merged


def prepare_mri_features(df):
    print("Preparing MRI features...")
    if hasattr(df, 'mri_features'):
        return df.mri_features

    thickness_cols = [col for col in df.columns if 'thickness' in col.lower()]
    area_cols = [col for col in df.columns if 'area' in col.lower() and 'thickness' not in col.lower()]
    volume_cols = [col for col in df.columns if 'volume' in col.lower()]
    return {'thickness': thickness_cols, 'area': area_cols, 'volume': volume_cols}


def perform_correlations(df, mri_features):
    print("Performing MRI-IQ/SRS/Vineland correlations...")

    correlation_results = {}

    for feature_type, feature_cols in mri_features.items():
        print(f"\nAnalyzing {feature_type} features...")

        iq_correlations = []
        srs_correlations = []
        vineland_correlations = []

        for col in feature_cols:
            merged_col = col
            if merged_col not in df.columns:
                continue

            atlas_region_name = col

            df_IQ = df[[merged_col, 'total_IQ']].dropna()
            if len(df_IQ) > 1:
                iq_corr, iq_p = pearsonr(df_IQ[merged_col], df_IQ['total_IQ'])
                iq_correlations.append({'region': atlas_region_name, 'correlation': iq_corr, 'p_value': iq_p, 'n': len(df_IQ)})

            df_SRS = df[[merged_col, 'SRS_t_score']].dropna()
            if len(df_SRS) > 1:
                srs_corr, srs_p = pearsonr(df_SRS[merged_col], df_SRS['SRS_t_score'])
                srs_correlations.append({'region': atlas_region_name, 'correlation': srs_corr, 'p_value': srs_p, 'n': len(df_SRS)})

            df_Vineland = df[[merged_col, 'vabsabcabc_standard']].dropna()
            if len(df_Vineland) > 1:
                vineland_corr, vineland_p = pearsonr(df_Vineland[merged_col], df_Vineland['vabsabcabc_standard'])
                vineland_correlations.append({'region': atlas_region_name, 'correlation': vineland_corr, 'p_value': vineland_p, 'n': len(df_Vineland)})

        iq_df = pd.DataFrame(iq_correlations)
        srs_df = pd.DataFrame(srs_correlations)
        vineland_df = pd.DataFrame(vineland_correlations)

        from statsmodels.stats.multitest import multipletests

        if len(iq_df) > 0:
            _, iq_df['p_value_fdr'], _, _ = multipletests(iq_df['p_value'], method='fdr_bh')
            iq_df['significant'] = iq_df['p_value_fdr'] < 0.05

        if len(srs_df) > 0:
            _, srs_df['p_value_fdr'], _, _ = multipletests(srs_df['p_value'], method='fdr_bh')
            srs_df['significant'] = srs_df['p_value_fdr'] < 0.05

        if len(vineland_df) > 0:
            _, vineland_df['p_value_fdr'], _, _ = multipletests(vineland_df['p_value'], method='fdr_bh')
            vineland_df['significant'] = vineland_df['p_value_fdr'] < 0.05

        correlation_results[feature_type] = {'iq': iq_df, 'srs': srs_df, 'vineland': vineland_df}

        print(f"  IQ: {len(iq_df)} regions, {iq_df['significant'].sum() if len(iq_df) > 0 else 0} significant")
        print(f"  SRS: {len(srs_df)} regions, {srs_df['significant'].sum() if len(srs_df) > 0 else 0} significant")
        print(f"  Vineland: {len(vineland_df)} regions, {vineland_df['significant'].sum() if len(vineland_df) > 0 else 0} significant")

    return correlation_results


def create_r_input_files(correlation_results, df):
    print("Creating R input files...")
    Path(R_INPUT_DIR).mkdir(parents=True, exist_ok=True)

    for feature_type, results in correlation_results.items():
        results['iq'].to_csv(os.path.join(R_INPUT_DIR, f'autism_{feature_type}_iq_correlations.csv'), index=False)
        results['srs'].to_csv(os.path.join(R_INPUT_DIR, f'autism_{feature_type}_srs_correlations.csv'), index=False)
        results['vineland'].to_csv(os.path.join(R_INPUT_DIR, f'autism_{feature_type}_vineland_correlations.csv'), index=False)

    summary_stats = {
        'Metric': ['N', 'Mean Age', 'SD Age', 'Mean IQ', 'SD IQ', 'Mean SRS', 'SD SRS', 'Mean Vineland', 'SD Vineland'],
        'Value': [
            len(df),
            round(df['age_yrs'].mean(), 2), round(df['age_yrs'].std(), 2),
            round(df['total_IQ'].mean(), 2), round(df['total_IQ'].std(), 2),
            round(df['SRS_t_score'].mean(), 2), round(df['SRS_t_score'].std(), 2),
            round(df['vabsabcabc_standard'].mean(), 2), round(df['vabsabcabc_standard'].std(), 2)
        ]
    }
    pd.DataFrame(summary_stats).to_csv(os.path.join(R_INPUT_DIR, 'autism_structural_mri_summary_statistics.csv'), index=False)
    print(f"R input files saved to: {R_INPUT_DIR}")
    return R_INPUT_DIR


def create_scatterplots(df, mri_features):
    print("Creating specific scatterplots...")

    key_regions = {
        'thickness': ['lh_pericalcarine_thickness', 'lh_superiortemporal_thickness'],
    }

    scatterplots = []

    for feature_type, regions in key_regions.items():
        for region in regions:
            merged_col = region + '_y' if region + '_y' in df.columns else (region + '_x' if region + '_x' in df.columns else region)

            if merged_col not in df.columns:
                continue

            for score_col, score_label, color in [
                ('total_IQ', 'Total IQ', '#2E6B8A'),
                ('SRS_t_score', 'SRS T Score', '#228B22'),
                ('vabsabcabc_standard', 'Vineland Total Score', '#FF8C00')
            ]:
                fig, ax = plt.subplots(figsize=(10, 8))
                df_score = df[[merged_col, score_col]].dropna()
                ax.scatter(df_score[merged_col], df_score[score_col], alpha=0.6, color=color, s=50)
                z = np.polyfit(df_score[merged_col], df_score[score_col], 1)
                p_fit = np.poly1d(z)
                ax.plot(df_score[merged_col], p_fit(df_score[merged_col]), color='#B30000', linewidth=2)
                corr, p_val = pearsonr(df_score[merged_col], df_score[score_col])
                ax.set_xlabel(f'{region.replace("_", " ").title()}')
                ax.set_ylabel(score_label)
                ax.set_title(f'{region.replace("_", " ").title()} vs {score_label}\nr = {corr:.3f}, p = {p_val:.3f}, n = {len(df_score)}')
                ax.grid(True, alpha=0.3)
                output_file = os.path.join(OUTPUT_DIR, f'autism_{region}_vs_{score_col}_scatterplot.pdf')
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()
                scatterplots.append(output_file)

    print(f"Created {len(scatterplots)} scatterplots")
    return scatterplots


def main():
    print("=" * 60)
    print("Structural MRI-IQ-SRS Correlation Analysis for Autism")
    print("=" * 60)

    df = load_and_prepare_data()
    mri_features = prepare_mri_features(df)
    correlation_results = perform_correlations(df, mri_features)
    create_r_input_files(correlation_results, df)
    create_scatterplots(df, mri_features)

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print(f"Final sample size: {len(df)} autism individuals")
    print(f"R input files: {R_INPUT_DIR}")


if __name__ == "__main__":
    main()
