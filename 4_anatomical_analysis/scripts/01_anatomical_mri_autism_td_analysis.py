#!/usr/bin/env python3
"""
Anatomical MRI Autism vs TD Analysis

Compares autism vs TD groups for:
- Cortical thickness (Desikan-Killiany atlas)
- Surface area (Desikan-Killiany atlas)
- Subcortical volumes (ASEG atlas)

Adapted from eeg_mri-pipeline/analysis/figures_papers/mri_autism_td_analysis/anatomical_mri_autism_td_analysis.py
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
from scipy.stats import ttest_ind
from sklearn.linear_model import LinearRegression
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from statsmodels.stats.multitest import multipletests

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

# Input data paths (referencing eeg_mri-pipeline)
EEG_MRI_RESULTS = '/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results'
DATA_IMAG2GENET_DIR = '/Users/mfleury/POSTDOC/LIBRAIRY/imaging2genet'

# =============================================================================
# ATLAS DICTIONARIES (from eeg_mri-pipeline/scripts/config/global_config.py)
# =============================================================================
data_dk = {'bankssts_left': -2, 'caudalanteriorcingulate_left': -2,
           'caudalmiddlefrontal_left': -2, 'cuneus_left': -2,
           'entorhinal_left': -2, 'fusiform_left': -2,
           'inferiorparietal_left': -2, 'inferiortemporal_left': -2,
           'isthmuscingulate_left': -2, 'lateraloccipital_left': -2,
           'lateralorbitofrontal_left': -2, 'lingual_left': -2,
           'medialorbitofrontal_left': -2, 'middletemporal_left': -2,
           'parahippocampal_left': -2, 'paracentral_left': -2,
           'parsopercularis_left': -2, 'parsorbitalis_left': -2,
           'parstriangularis_left': -2, 'pericalcarine_left': -2,
           'postcentral_left': -2, 'posteriorcingulate_left': -2,
           'precentral_left': -2, 'precuneus_left': -2,
           'rostralanteriorcingulate_left': -2, 'rostralmiddlefrontal_left': 3,
           'superiorfrontal_left': -2, 'superiorparietal_left': -2,
           'superiortemporal_left': -2, 'supramarginal_left': -2,
           'frontalpole_left': -2, 'temporalpole_left': -2,
           'transversetemporal_left': -2, 'insula_left': -2,
           'bankssts_right': -2, 'caudalanteriorcingulate_right': -2,
           'caudalmiddlefrontal_right': -2, 'cuneus_right': -2,
           'entorhinal_right': -2, 'fusiform_right': -2,
           'inferiorparietal_right': -2, 'inferiortemporal_right': -2,
           'isthmuscingulate_right': -2, 'lateraloccipital_right': -2,
           'lateralorbitofrontal_right': -2, 'lingual_right': -2,
           'medialorbitofrontal_right': -2, 'middletemporal_right': -2,
           'parahippocampal_right': -2, 'paracentral_right': -2,
           'parsopercularis_right': -2, 'parsorbitalis_right': -2,
           'parstriangularis_right': -2, 'pericalcarine_right': -2,
           'postcentral_right': -2, 'posteriorcingulate_right': -2,
           'precentral_right': -2, 'precuneus_right': -2,
           'rostralanteriorcingulate_right': -2, 'rostralmiddlefrontal_right': -2,
           'superiorfrontal_right': -2, 'superiorparietal_right': -2,
           'superiortemporal_right': -2, 'supramarginal_right': -2,
           'frontalpole_right': -2, 'temporalpole_right': -2,
           'transversetemporal_right': -2, 'insula_right': -2}

data_aseg = {'Left-Lateral-Ventricle': 12289.6, 'Left-Thalamus': 8158.3,
             'Left-Caudate': 3463.3, 'Left-Putamen': 4265.3,
             'Left-Pallidum': 1620.9, '3rd-Ventricle': 1635.6,
             '4th-Ventricle': 1115.6, 'Brain-Stem': 20393.1,
             'Left-Hippocampus': 4307.1, 'Left-Amygdala': 1512.4,
             'Left-VentralDC': 3871.9, 'Right-Lateral-Ventricle': 9530.7,
             'Right-Thalamus': 7613.9, 'Right-Caudate': 3605.0,
             'Right-Putamen': 4316.9, 'Right-Pallidum': 1696.9,
             'Right-Hippocampus': 4516.9, 'Right-Amygdala': 1781.7,
             'Right-VentralDC': 4056.8, 'CC_Posterior': 939.8,
             'CC_Mid_Posterior': 347.4, 'CC_Central': 399.0,
             'CC_Mid_Anterior': 443.7, 'CC_Anterior': 804.8,
             'Right-Cerebellum-Cortex': 21554.0, 'Right-Cerebellum-White-Matter': 2155.4,
             'Left-Cerebellum-Cortex': 21554.0, 'Left-Cerebellum-White-Matter': 2155.4}


def transform_name(name, metric):
    if name.startswith('lh_'):
        return name.replace(f'_{metric}', '')
    elif name.startswith('rh_'):
        return name.replace(f'_{metric}', '')
    else:
        return name


def transform_column_names(df, metric):
    df.columns = [transform_name(col, metric) for col in df.columns]
    return df


def cohens_d(group1, group2):
    n1 = len(group1)
    n2 = len(group2)
    if n1 < 2 or n2 < 2:
        return float('nan')
    mean1 = np.mean(group1)
    mean2 = np.mean(group2)
    std1 = np.std(group1, ddof=1)
    std2 = np.std(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return float('nan')
    return (mean1 - mean2) / pooled_std


def load_and_prepare_data():
    print("Loading structural MRI dataset...")

    df_mri = pd.read_csv(os.path.join(DATA_IMAG2GENET_DIR, '0_input', 'dataframes', 'MRI_ANAT_INOVAND_LEAP_COMBAT.tsv'), sep='\t')
    df_mri['ID'] = df_mri['ID'].fillna(df_mri['MRI_ID'])
    df_mri['ID'] = df_mri['ID'].astype(str)

    # QC is handled upstream (during FreeSurfer ComBat/regression z-scoring) —
    # the per-script INOVAND and LEAP blacklist filters have been removed.

    df_mri = df_mri[df_mri['PopulationS1'].isin(['TD', 'Autism'])].copy()

    # Rename population label TD -> NT for display
    df_mri['PopulationS1'] = df_mri['PopulationS1'].replace('TD', 'NT')

    print(f"Final dataset size: {len(df_mri)} participants")
    print(f"NT: {len(df_mri[df_mri['PopulationS1'] == 'NT'])}")
    print(f"Autism: {len(df_mri[df_mri['PopulationS1'] == 'Autism'])}")

    return df_mri


def compute_t_statistics(df, metric_name, atlas_regions, atlas_type='dk'):
    print(f"\nComputing t-statistics for {metric_name} ({atlas_type} atlas)...")

    df_control = df[df['PopulationS1'] == 'NT'].copy()
    df_autism = df[df['PopulationS1'] == 'Autism'].copy()

    print(f"NT subjects: {len(df_control)}")
    print(f"Autism subjects: {len(df_autism)}")

    t_stat_results = []

    if atlas_type == 'dk':
        excluded = ['lh_MeanThickness', 'rh_MeanThickness', 'lh_WhiteSurfArea', 'rh_WhiteSurfArea']
        dk_cols = [col for col in df_control.columns
                   if col.startswith(('lh_', 'rh_')) and col not in excluded]

        print(f"Found {len(dk_cols)} DK atlas columns to process")

        for col in dk_cols:
            if col in df_control.columns and col in df_autism.columns:
                control_values = df_control[col].dropna()
                autism_values = df_autism[col].dropna()
                if len(control_values) > 1 and len(autism_values) > 1:
                    t_stat, p_val = ttest_ind(autism_values, control_values, equal_var=False)
                    cohens_d_val = cohens_d(autism_values, control_values)

                    if col.startswith('lh_'):
                        region = col.replace('lh_', '')
                        label = f'{region}_left'
                    elif col.startswith('rh_'):
                        region = col.replace('rh_', '')
                        label = f'{region}_right'
                    else:
                        label = col

                    t_stat_results.append({'label': label, 't_stat': t_stat, 'p_val': p_val, 'cohens_d': cohens_d_val})
    else:
        for label in atlas_regions.keys():
            if label in df_control.columns and label in df_autism.columns:
                control_values = df_control[label].dropna()
                autism_values = df_autism[label].dropna()
                if len(control_values) > 1 and len(autism_values) > 1:
                    t_stat, p_val = ttest_ind(autism_values, control_values, equal_var=False)
                    cohens_d_val = cohens_d(autism_values, control_values)
                else:
                    t_stat, p_val = float('nan'), float('nan')
                    cohens_d_val = float('nan')
                t_stat_results.append({'label': label, 't_stat': t_stat, 'p_val': p_val, 'cohens_d': cohens_d_val})

    t_stat_df = pd.DataFrame(t_stat_results)

    valid_p_mask = t_stat_df['p_val'].notna()
    if valid_p_mask.sum() > 0:
        p_fdr_values = np.full(len(t_stat_df), np.nan)
        p_fdr_values[valid_p_mask] = multipletests(t_stat_df.loc[valid_p_mask, 'p_val'], method='bonferroni')[1]
        t_stat_df['p_fdr'] = p_fdr_values
    else:
        t_stat_df['p_fdr'] = np.nan

    output_file = os.path.join(R_INPUT_DIR, f't_stat_anat_{atlas_type}_{metric_name}_mri_autism_vs_control.csv')

    if len(t_stat_df) > 0:
        t_stat_df[['label', 't_stat', 'p_val', 'cohens_d', 'p_fdr']].to_csv(output_file, index=False, header=False)
        print(f"Saved t-statistics to: {output_file}")
        print(f"Total regions: {len(t_stat_df)}")
        print(f"Valid t-statistics: {t_stat_df['t_stat'].notna().sum()}")
    else:
        print(f"Warning: No regions found for {metric_name} ({atlas_type} atlas)")
        with open(output_file, 'w') as f:
            f.write("label,t_stat,p_val,cohens_d,p_fdr\n")

    return t_stat_df


def compute_cohort_specific_statistics(df, metric_name, atlas_regions, atlas_type='dk'):
    print(f"\nComputing cohort-specific t-statistics for {metric_name}...")

    cohort_output_dir = os.path.join(R_INPUT_DIR, 'cohort_specific')
    Path(cohort_output_dir).mkdir(parents=True, exist_ok=True)

    for cohort_name in ['INOVAND', 'LEAP']:
        print(f"\nProcessing {cohort_name} cohort...")
        df_cohort = df[df['cohort'] == cohort_name].copy()
        if len(df_cohort) == 0:
            print(f"Warning: No data found for {cohort_name} cohort")
            continue

        df_control = df_cohort[df_cohort['PopulationS1'] == 'NT'].copy()
        df_autism = df_cohort[df_cohort['PopulationS1'] == 'Autism'].copy()

        print(f"  NT subjects: {len(df_control)}")
        print(f"  Autism subjects: {len(df_autism)}")

        if len(df_control) < 2 or len(df_autism) < 2:
            print(f"  Warning: Insufficient sample size for {cohort_name}")
            continue

        t_stat_results = []

        if atlas_type == 'dk':
            excluded = ['lh_MeanThickness', 'rh_MeanThickness', 'lh_WhiteSurfArea', 'rh_WhiteSurfArea']
            dk_cols = [col for col in df_control.columns
                       if col.startswith(('lh_', 'rh_')) and col not in excluded]
            for col in dk_cols:
                if col in df_control.columns and col in df_autism.columns:
                    control_values = df_control[col].dropna()
                    autism_values = df_autism[col].dropna()
                    if len(control_values) > 1 and len(autism_values) > 1:
                        t_stat, p_val = ttest_ind(autism_values, control_values, equal_var=False)
                        cohens_d_val = cohens_d(autism_values, control_values)
                        if col.startswith('lh_'):
                            label = f'{col.replace("lh_", "")}_left'
                        elif col.startswith('rh_'):
                            label = f'{col.replace("rh_", "")}_right'
                        else:
                            label = col
                        t_stat_results.append({'label': label, 't_stat': t_stat, 'p_val': p_val, 'cohens_d': cohens_d_val})
        else:
            for label in atlas_regions.keys():
                if label in df_control.columns and label in df_autism.columns:
                    control_values = df_control[label].dropna()
                    autism_values = df_autism[label].dropna()
                    if len(control_values) > 1 and len(autism_values) > 1:
                        t_stat, p_val = ttest_ind(autism_values, control_values, equal_var=False)
                        cohens_d_val = cohens_d(autism_values, control_values)
                    else:
                        t_stat, p_val = float('nan'), float('nan')
                        cohens_d_val = float('nan')
                    t_stat_results.append({'label': label, 't_stat': t_stat, 'p_val': p_val, 'cohens_d': cohens_d_val})

        t_stat_df = pd.DataFrame(t_stat_results)

        valid_p_mask = t_stat_df['p_val'].notna()
        if valid_p_mask.sum() > 0:
            p_fdr_values = np.full(len(t_stat_df), np.nan)
            p_fdr_values[valid_p_mask] = multipletests(t_stat_df.loc[valid_p_mask, 'p_val'], method='fdr_bh')[1]
            t_stat_df['p_fdr'] = p_fdr_values
        else:
            t_stat_df['p_fdr'] = np.nan

        output_file = os.path.join(cohort_output_dir, f't_stat_anat_{atlas_type}_{metric_name}_mri_autism_vs_control_{cohort_name}.csv')
        if len(t_stat_df) > 0:
            t_stat_df[['label', 't_stat', 'p_val', 'cohens_d', 'p_fdr']].to_csv(output_file, index=False, header=False)
            print(f"  Saved {cohort_name} t-statistics to: {output_file}")


def main():
    print("=" * 60)
    print("Anatomical MRI Autism vs NT Analysis")
    print("=" * 60)

    df_mri = load_and_prepare_data()

    metrics = {
        'thickness': {'atlas': data_dk, 'atlas_type': 'dk'},
        'area': {'atlas': data_dk, 'atlas_type': 'dk'},
        'volume': {'atlas': data_aseg, 'atlas_type': 'aseg'}
    }

    for metric_name, config in metrics.items():
        print(f"\n{'=' * 60}")
        print(f"Processing {metric_name.upper()}")
        print(f"{'=' * 60}")

        df_metric = df_mri.copy()

        if config['atlas_type'] == 'dk':
            metric_cols = [col for col in df_metric.columns
                           if col.endswith(f'_{metric_name}') and col.startswith(('lh_', 'rh_'))]
            non_metric_cols = [col for col in df_metric.columns
                               if not col.endswith(('_thickness', '_area', '_volume')) or
                               col in ['lh_MeanThickness', 'rh_MeanThickness', 'lh_WhiteSurfArea', 'rh_WhiteSurfArea']]
            cols_to_keep = list(set(metric_cols + non_metric_cols))
            df_metric = df_metric[cols_to_keep].copy()

        df_metric = transform_column_names(df_metric, metric_name)

        compute_t_statistics(df_metric, metric_name, config['atlas'], atlas_type=config['atlas_type'])

        if 'cohort' in df_metric.columns:
            compute_cohort_specific_statistics(df_metric, metric_name, config['atlas'], config['atlas_type'])

    print(f"\n{'=' * 60}")
    print("Analysis Complete!")
    print(f"{'=' * 60}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"R input files directory: {R_INPUT_DIR}")
    print("\nNext step: Run the R visualization script:")
    print("  Rscript scripts/03_plot_anatomical_mri_brain_visualizations.R")


if __name__ == "__main__":
    main()
