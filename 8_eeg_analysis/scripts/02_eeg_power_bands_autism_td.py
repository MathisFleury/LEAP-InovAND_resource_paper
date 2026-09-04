#!/usr/bin/env python3
"""
EEG Multi-Band Power Autism vs TD Analysis

Analyzes EEG power spectra (absolute and relative) comparing Autism vs TD:
- Per-band statistical comparisons (t-test or Mann-Whitney)
- FDR correction across all features
- Covariate effect analysis (age, sex, cohort)
- Per-feature violin plots

Adapted from eeg_mri-pipeline/analysis/figures_papers/eeg_autism_td_analysis/run_eeg_power_bands_autism_td.py
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent  # 8_eeg_analysis/
sys.path.insert(0, str(_SECTION_DIR / "preprocessing"))
from step_05_regress import regress  # noqa: E402
from step_06_zscore import zscore    # noqa: E402
from config import POWER_COMBINED    # noqa: E402
from curated_map import attach_curated_demographics  # noqa: E402

FIG_DIR = _SECTION_DIR / 'outputs' / 'figures'
TABLES_DIR = _SECTION_DIR / 'outputs' / 'tables'
STATS_DIR = TABLES_DIR / "stats"
FDR_DIR = TABLES_DIR / "fdr_corrected"
COVARIATE_DIR = TABLES_DIR / "covariate_effects"
PLOTS_DIR = FIG_DIR / "plots"


def reshape_absolute_power(df_eeg: pd.DataFrame) -> pd.DataFrame:
    filt = df_eeg.loc[
        (df_eeg['Brain_region'].isin(['Central', 'Frontal', 'Posterior'])) &
        (df_eeg['Brain_side'].isin(['Left', 'Right', 'Mid'])) &
        (df_eeg['Quant_status'] == 'Absolute') &
        (df_eeg['Eye_status'] == 'eyeo')
    ]
    metadata_cols = [col for col in ['PopulationS1', 'age_yrs', 'Sex', 'cohort'] if col in filt.columns]
    if metadata_cols:
        metadata = filt.groupby('ID')[metadata_cols].first().reset_index()
    psd_col = 'PSD_corrected' if 'PSD_corrected' in filt.columns else 'PSD'
    agg = filt.groupby(['ID', 'Brain_region', 'Brain_side', 'Freq_band'], as_index=False).agg(PSD=(psd_col, 'median'))
    piv = agg.pivot(index=['ID'], columns=['Brain_region', 'Brain_side', 'Freq_band'], values='PSD')
    piv.columns = ['abs_' + '_'.join(col) for col in piv.columns]
    piv = piv.reset_index()
    if metadata_cols:
        piv = piv.merge(metadata, on='ID', how='left')
    return piv


def reshape_relative_power(df_eeg: pd.DataFrame) -> pd.DataFrame:
    filt = df_eeg.loc[
        (df_eeg['Brain_region'].isin(['Central', 'Frontal', 'Posterior'])) &
        (df_eeg['Brain_side'].isin(['Left', 'Right', 'Mid'])) &
        (df_eeg['Quant_status'] == 'Relative') &
        (df_eeg['Eye_status'] == 'eyeo')
    ]
    metadata_cols = [col for col in ['PopulationS1', 'age_yrs', 'Sex', 'cohort'] if col in filt.columns]
    if metadata_cols:
        metadata = filt.groupby('ID')[metadata_cols].first().reset_index()
    psd_col = 'PSD_corrected' if 'PSD_corrected' in filt.columns else 'PSD'
    agg = filt.groupby(['ID', 'Brain_region', 'Brain_side', 'Freq_band'], as_index=False).agg(PSD=(psd_col, 'median'))
    piv = agg.pivot(index=['ID'], columns=['Brain_region', 'Brain_side', 'Freq_band'], values='PSD')
    piv.columns = ['rel_' + '_'.join(col) for col in piv.columns]
    piv = piv.reset_index()
    if metadata_cols:
        piv = piv.merge(metadata, on='ID', how='left')
    return piv


def load_power_autism_td() -> pd.DataFrame:
    # Frozen long-format power (LEAP has no raw power on $IMG5), with curated
    # demographics + k-means Cluster re-attached (curated is the only clinical source).
    print(f"Loading combined corrected power spectrum from: {POWER_COMBINED}")
    power = pd.read_csv(POWER_COMBINED, low_memory=False)
    power = attach_curated_demographics(power)

    # Filter to valid corrected values
    if 'PSD_corrected' in power.columns:
        power = power[power['PSD_corrected'].notna()]
        power = power[power['PSD_corrected'] != 0]

    power = power.dropna(subset=['age_yrs'])

    # Print sample counts
    try:
        power_unique = power[['ID', 'PopulationS1']].drop_duplicates(subset=['ID'])
        print("Power spectrum counts by PopulationS1 (unique IDs):", power_unique["PopulationS1"].value_counts(dropna=False).to_dict())
    except Exception as e:
        print(f"Warning: could not compute counts: {e}")

    return power


def analyze_covariate_effects(df: pd.DataFrame, features: list) -> pd.DataFrame:
    df = df.copy()
    df['age_yrs_sq'] = df['age_yrs'] ** 2
    all_results = []

    for feat in features:
        valid_mask = ~df[feat].isna()
        if valid_mask.sum() < 20:
            continue
        df_subset = df.loc[valid_mask].copy()
        base_covars = ['age_yrs', 'age_yrs_sq']
        if 'Sex' in df_subset.columns:
            base_covars.append('Sex')
        valid_X_mask = ~df_subset[base_covars].isna().any(axis=1)
        df_subset = df_subset.loc[valid_X_mask]
        y = df_subset[feat].values

        X_list = [df_subset['age_yrs'].values, df_subset['age_yrs_sq'].values]
        predictor_names = ['age_yrs', 'age_yrs_sq']
        if 'Sex' in df_subset.columns:
            X_list.append(df_subset['Sex'].values)
            predictor_names.append('Sex')
        if 'cohort' in df_subset.columns:
            X_list.append((df_subset['cohort'].values == 'INOVAND').astype(float))
            predictor_names.append('cohort_INOVAND')

        X = np.column_stack(X_list)
        X = pd.DataFrame(X, columns=predictor_names, index=df_subset.index)

        if len(y) < 20 or np.isnan(y).all() or X.isna().any().any():
            continue

        X_const = sm.add_constant(X)
        try:
            model = sm.OLS(y, X_const).fit()
            residuals = y - model.predict(X_const)
            residual_sd = np.std(residuals, ddof=model.df_resid)
            for predictor in X.columns:
                if predictor in model.params.index:
                    coef = model.params[predictor]
                    pval = model.pvalues[predictor]
                    ci_lower, ci_upper = model.conf_int().loc[predictor]
                    t_stat = model.tvalues[predictor]
                    cohens_d_val = coef / residual_sd if residual_sd > 0 else np.nan
                    all_results.append({
                        'feature': feat, 'predictor': predictor, 'coefficient': coef,
                        'cohens_d': cohens_d_val, 'p_value': pval,
                        'ci_lower': ci_lower, 'ci_upper': ci_upper,
                        't_statistic': t_stat, 'significant': pval < 0.05,
                        'r_squared': model.rsquared, 'n_samples': len(y)
                    })
        except Exception as e:
            continue

    return pd.DataFrame(all_results)


def remove_outliers(df: pd.DataFrame, value_cols: list, threshold: float = 10.0) -> pd.DataFrame:
    outlier_mask = pd.Series([False] * len(df), index=df.index)
    for value_col in value_cols:
        if value_col in df.columns:
            outlier_mask = outlier_mask | (df[value_col] > threshold)
    return df.loc[~outlier_mask].copy()


def autism_vs_td_stats(df: pd.DataFrame, value_cols: list, apply_fdr: bool = True) -> pd.DataFrame:
    df = remove_outliers(df, value_cols, threshold=10.0)
    autism_data_df = df[df["PopulationS1"] == "Autism"]
    td_data_df = df[df["PopulationS1"] == "NT"]

    all_rows = []
    p_values = []

    for value_col in value_cols:
        autism_values = autism_data_df[value_col].dropna()
        td_values = td_data_df[value_col].dropna()
        if len(autism_values) < 3 or len(td_values) < 3:
            continue
        autism_normal = stats.shapiro(autism_values)[1] > 0.05 if len(autism_values) >= 3 else False
        td_normal = stats.shapiro(td_values)[1] > 0.05 if len(td_values) >= 3 else False
        test_type = "t-test" if (autism_normal and td_normal) else "Mann-Whitney"
        if test_type == "t-test":
            t_val, p = stats.ttest_ind(autism_values, td_values, equal_var=False)
        else:
            t_val, p = stats.mannwhitneyu(autism_values, td_values, alternative="two-sided")
        all_rows.append({
            "feature": value_col, "test_type": test_type, "comparison": "Autism vs NT",
            "autism_n": len(autism_values), "td_n": len(td_values),
            "autism_mean": autism_values.mean(), "td_mean": td_values.mean(),
            "autism_std": autism_values.std(), "td_std": td_values.std(),
            "t_or_u": t_val, "p_uncorrected": p, "significant": p < 0.05,
        })
        p_values.append(p)

    stats_df = pd.DataFrame(all_rows)
    if apply_fdr and len(p_values) > 0:
        rejected, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
        stats_df['p_fdr_corrected'] = p_corrected
        stats_df['significant_fdr'] = rejected
    else:
        stats_df['p_fdr_corrected'] = stats_df.get('p_uncorrected', np.nan)
        stats_df['significant_fdr'] = stats_df.get('significant', False)
    return stats_df


def extract_band_and_roi(feature_name: str):
    feature_clean = feature_name.replace('abs_', '').replace('rel_', '')
    parts = feature_clean.split('_')
    band = parts[-1] if len(parts) > 0 else None
    roi = '_'.join(parts[:-1]) if len(parts) > 1 else None
    return roi, band


def save_fdr_results_by_band_roi(stats_df: pd.DataFrame, power_type: str, cohort_name: str = None):
    if stats_df.empty:
        return
    stats_df = stats_df.copy()
    stats_df['ROI'] = stats_df['feature'].apply(lambda x: extract_band_and_roi(x)[0])
    stats_df['Band'] = stats_df['feature'].apply(lambda x: extract_band_and_roi(x)[1])
    output_cols = ['Band', 'ROI', 'feature', 'test_type', 'autism_n', 'td_n',
                   'autism_mean', 'td_mean', 'autism_std', 'td_std',
                   't_or_u', 'p_uncorrected', 'p_fdr_corrected', 'significant', 'significant_fdr']
    output_cols = [col for col in output_cols if col in stats_df.columns]
    fdr_results = stats_df[output_cols].sort_values(['Band', 'ROI', 'p_fdr_corrected'])
    suffix = f'_{cohort_name}' if cohort_name else ''
    filename = f'{power_type}_power_fdr_corrected_by_band_roi{suffix}.csv'
    FDR_DIR.mkdir(parents=True, exist_ok=True)
    output_path = FDR_DIR / filename
    fdr_results.to_csv(output_path, index=False)
    print(f"  FDR results saved: {output_path}")
    print(f"    Significant (FDR): {fdr_results['significant_fdr'].sum() if 'significant_fdr' in fdr_results.columns else 'N/A'}")
    return fdr_results


def plot_per_feature_autism_vs_td(df: pd.DataFrame, corrected_cols: list, band: str):
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    df = remove_outliers(df, corrected_cols, threshold=10.0)

    autism_data = df[df["PopulationS1"] == "Autism"].copy()
    td_data = df[df["PopulationS1"] == "NT"].copy()
    autism_data["Group"] = "Autism"
    td_data["Group"] = "NT"
    plot_data = pd.concat([autism_data, td_data], ignore_index=True)

    colors = {"Autism": "#8991FA", "NT": "#C1C2BC"}
    plot_order = ["Autism", "NT"]

    for col in corrected_cols:
        fig, ax = plt.subplots(figsize=(8, 10))
        sns.violinplot(data=plot_data, x='Group', y=col, order=plot_order, palette=colors, inner=None, ax=ax)
        np.random.seed(42)
        for i, group in enumerate(plot_order):
            group_data = plot_data[plot_data['Group'] == group][col].dropna()
            if len(group_data) > 0:
                x_jitter = np.clip(np.random.normal(i, 0.15, size=len(group_data)), i - 0.4, i + 0.4)
                ax.scatter(x_jitter, group_data.values, color='black', alpha=0.6, s=25,
                           edgecolors='white', linewidths=0.5, zorder=10)
        ax.set_xlabel('')
        ax.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.5)
        ax.set_ylabel(col.replace('_', ' ').replace('abs ', '').replace('rel ', ''), fontsize=12)
        ax.set_title(f'{col.replace("_", " ").replace("abs ", "").replace("rel ", "")} - {band}', fontsize=14)
        safe_feat = col.replace('/', '_').replace(' ', '_')
        fig.savefig(PLOTS_DIR / f'violinplot_{safe_feat}_power_{band}_autism_vs_td.png', bbox_inches='tight', dpi=300)
        plt.close(fig)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    power = load_power_autism_td()
    print(f"Loaded data: {len(power)} rows")

    power_types = [
        ('absolute', reshape_absolute_power, 'abs_'),
        ('relative', reshape_relative_power, 'rel_')
    ]

    for power_type, reshape_func, prefix in power_types:
        print(f"\n=== Processing {power_type.upper()} POWER (COMBINED) ===")

        power_reshaped = reshape_func(power)
        autism_mask = power_reshaped["PopulationS1"] == "Autism"
        td_mask = power_reshaped["PopulationS1"] == "NT"
        power_reshaped = power_reshaped[autism_mask | td_mask].dropna(subset=['age_yrs'])

        all_cols = [c for c in power_reshaped.columns if c.startswith(prefix)]
        if not all_cols:
            print(f'No {power_type} power features found.')
            continue

        # Regress age/age²/sex and z-score on the FULL sample (no NT reference) —
        # see preprocessing/ (mirrors the anat z-scoring pipeline, DX_COL=None).
        power_reshaped = regress(power_reshaped, all_cols)
        power_reshaped = zscore(power_reshaped, all_cols)

        bands = sorted({c.split('_')[-1] for c in all_cols})
        print(f"Found {len(all_cols)} {power_type} power features across {len(bands)} bands")

        # Covariate effects
        print(f"\nAnalyzing covariate effects for all {power_type} power features...")
        covariate_results = analyze_covariate_effects(power_reshaped, all_cols)
        if not covariate_results.empty:
            COVARIATE_DIR.mkdir(parents=True, exist_ok=True)
            covariate_file = COVARIATE_DIR / f'{power_type}_covariate_effects_all_features.csv'
            covariate_results.to_csv(covariate_file, index=False)
            print(f"Covariate effects saved: {covariate_file}")

        # Collect stats across all bands
        all_stats_list = []
        for band in bands:
            band_feats = [c for c in all_cols if c.endswith(f'_{band}')]
            if not band_feats:
                continue
            print(f"  Using already-corrected features for {band} band")
            df_band = power_reshaped.copy()
            stats_df = autism_vs_td_stats(df_band, band_feats, apply_fdr=False)
            stats_df['band'] = band
            all_stats_list.append(stats_df)
            STATS_DIR.mkdir(parents=True, exist_ok=True)
            stats_df.to_csv(STATS_DIR / f'{power_type}_power_{band}_autism_vs_td_stats.csv', index=False)
            plot_per_feature_autism_vs_td(df_band, band_feats, f'{power_type}_{band}')

        # FDR correction across all bands
        if all_stats_list:
            all_stats_combined = pd.concat(all_stats_list, ignore_index=True)
            p_values = all_stats_combined['p_uncorrected'].values
            rejected, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
            all_stats_combined['p_fdr_corrected'] = p_corrected
            all_stats_combined['significant_fdr'] = rejected
            print(f"\n  FDR correction applied to {len(all_stats_combined)} features")
            print(f"  Significant (FDR-corrected): {all_stats_combined['significant_fdr'].sum()}")
            save_fdr_results_by_band_roi(all_stats_combined, power_type, cohort_name=None)

        # Cohort-specific analysis
        if 'cohort' in power.columns:
            for cohort_name in ['LEAP', 'INOVAND']:
                print(f"\n{'='*60}")
                print(f"=== Processing {cohort_name} {power_type.upper()} POWER (COHORT-SPECIFIC) ===")
                power_cohort = power[power['cohort'] == cohort_name].copy()
                if len(power_cohort) == 0:
                    print(f"Warning: No data found for {cohort_name}")
                    continue
                power_reshaped_cohort = reshape_func(power_cohort)
                autism_m = power_reshaped_cohort["PopulationS1"] == "Autism"
                td_m = power_reshaped_cohort["PopulationS1"] == "NT"
                power_reshaped_cohort = power_reshaped_cohort[autism_m | td_m].dropna(subset=['age_yrs'])
                cohort_cols = [c for c in power_reshaped_cohort.columns if c.startswith(prefix)]
                if not cohort_cols:
                    continue
                # Regress + full-sample z-score within the cohort (no NT reference).
                power_reshaped_cohort = regress(power_reshaped_cohort, cohort_cols)
                power_reshaped_cohort = zscore(power_reshaped_cohort, cohort_cols)
                cohort_bands = sorted({c.split('_')[-1] for c in cohort_cols})
                cohort_stats_list = []
                for band in cohort_bands:
                    band_feats = [c for c in cohort_cols if c.endswith(f'_{band}')]
                    if not band_feats:
                        continue
                    stats_df = autism_vs_td_stats(power_reshaped_cohort, band_feats, apply_fdr=False)
                    stats_df['band'] = band
                    cohort_stats_list.append(stats_df)
                    stats_df.to_csv(STATS_DIR / f'{power_type}_power_{band}_autism_vs_td_stats_{cohort_name}.csv', index=False)
                    plot_per_feature_autism_vs_td(power_reshaped_cohort, band_feats, f'{power_type}_{band}_{cohort_name}')
                if cohort_stats_list:
                    cohort_combined = pd.concat(cohort_stats_list, ignore_index=True)
                    rej, p_corr, _, _ = multipletests(cohort_combined['p_uncorrected'].values, alpha=0.05, method='fdr_bh')
                    cohort_combined['p_fdr_corrected'] = p_corr
                    cohort_combined['significant_fdr'] = rej
                    print(f"  {cohort_name} significant (FDR): {cohort_combined['significant_fdr'].sum()}")
                    save_fdr_results_by_band_roi(cohort_combined, power_type, cohort_name=cohort_name)

    print(f"\n=== ANALYSIS COMPLETED ===")
    print(f"Results saved to: {FIG_DIR}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
