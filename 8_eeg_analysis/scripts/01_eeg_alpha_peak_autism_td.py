#!/usr/bin/env python3
"""
EEG Alpha Peak Autism vs TD Analysis

Analyzes individual alpha frequency comparing Autism vs TD:
- Age + sex regression
- Z-score normalization using TD as reference
- Violin plots with individual data points

Adapted from eeg_mri-pipeline/analysis/figures_papers/eeg_autism_td_analysis/run_eeg_alpha_peak_autism_td.py
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent  # 5_eeg_analysis/

FIG_DIR = _SECTION_DIR / 'outputs' / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

EEG_MRI_RESULTS = Path('/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results')
DATA_DIR = EEG_MRI_RESULTS / 'dataset_paper' / 'dataframes'
DF_CLUSTERS_FILE = DATA_DIR / 'df_clusters_complete_kmeans.csv'


def load_alpha_peak_autism_td():
    clusters = pd.read_csv(DF_CLUSTERS_FILE)

    alpha_peak_combined = pd.read_csv(DATA_DIR / "Alpha_peak_combined_corrected.csv")
    print(f"Loaded {len(alpha_peak_combined)} subjects from combined alpha peak file")

    alpha_peak_combined["ID"] = alpha_peak_combined["ID"].astype(str)

    clusters_subset = clusters[["ID", "Sex", "genetics_score", "Population1"]].copy()
    clusters_subset["ID"] = clusters_subset["ID"].astype(str)

    alpha_peak = alpha_peak_combined.merge(clusters_subset, on="ID", how="inner")

    # Rename population label TD -> NT for display
    alpha_peak['PopulationS1'] = alpha_peak['PopulationS1'].replace('TD', 'NT')

    alpha_peak = alpha_peak.dropna(subset=["age_yrs"])
    alpha_peak = alpha_peak[alpha_peak["alpha_peak_corrected"].notna()]
    alpha_peak = alpha_peak[alpha_peak["alpha_peak_corrected"] != 0]

    try:
        s1_counts = alpha_peak["PopulationS1"].value_counts(dropna=False).to_dict()
        print("Alpha peak counts by PopulationS1:", s1_counts)
    except Exception as e:
        print(f"Warning: could not compute counts: {e}")

    return alpha_peak


def regress_and_zscore(alpha_peak: pd.DataFrame) -> pd.DataFrame:
    alpha_peak = alpha_peak.copy()
    alpha_peak["age_yrs_sq"] = alpha_peak["age_yrs"] ** 2
    X = alpha_peak[["age_yrs", "Sex", "age_yrs_sq"]].copy()
    X['Sex'] = X['Sex'].fillna(0)
    y = alpha_peak["alpha_peak_corrected"]
    model = LinearRegression().fit(X, y)
    residuals = y - model.predict(X)

    controls = alpha_peak[alpha_peak["PopulationS1"] == "NT"]
    mu = controls["alpha_peak_residuals"].mean() if "alpha_peak_residuals" in controls.columns else residuals[alpha_peak["PopulationS1"] == "NT"].mean()
    sd_val = residuals[alpha_peak["PopulationS1"] == "NT"].std()

    alpha_peak["alpha_peak_corrected"] = (residuals - mu) / sd_val
    return alpha_peak


def autism_vs_td_stats(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    autism_data = df[df["PopulationS1"] == "Autism"][value_col].dropna()
    td_data = df[df["PopulationS1"] == "NT"][value_col].dropna()

    if len(autism_data) < 3 or len(td_data) < 3:
        return pd.DataFrame()

    autism_normal = stats.shapiro(autism_data)[1] > 0.05 if len(autism_data) >= 3 else False
    td_normal = stats.shapiro(td_data)[1] > 0.05 if len(td_data) >= 3 else False
    test_type = "t-test" if (autism_normal and td_normal) else "Mann-Whitney"

    if test_type == "t-test":
        t_val, p = stats.ttest_ind(autism_data, td_data, equal_var=False)
    else:
        t_val, p = stats.mannwhitneyu(autism_data, td_data, alternative="two-sided")

    return pd.DataFrame([{
        "feature": value_col, "test_type": test_type, "comparison": "Autism vs NT",
        "autism_n": len(autism_data), "td_n": len(td_data),
        "autism_mean": autism_data.mean(), "td_mean": td_data.mean(),
        "autism_std": autism_data.std(), "td_std": td_data.std(),
        "t_or_u": t_val, "p_uncorrected": p, "significant": p < 0.05,
    }])


def plot_alpha_peak_autism_vs_td(df: pd.DataFrame):
    autism_data = df[df["PopulationS1"] == "Autism"].copy()
    td_data = df[df["PopulationS1"] == "NT"].copy()
    autism_data["Group"] = "Autism"
    td_data["Group"] = "NT"
    plot_data = pd.concat([autism_data, td_data], ignore_index=True)

    fig, ax = plt.subplots(1, 1, figsize=(6, 10))
    colors = {"Autism": "#8991FA", "NT": "#C1C2BC"}
    plot_order = ["Autism", "NT"]

    sns.violinplot(y="alpha_peak_corrected", x="Group", data=plot_data, inner=None,
                   linewidth=0, order=plot_order, palette=colors)
    sns.stripplot(y="alpha_peak_corrected", x="Group", data=plot_data, color="black", size=2,
                  jitter=True, alpha=0.5, order=plot_order)
    sns.pointplot(y="alpha_peak_corrected", x="Group", data=plot_data, color="#8A0201",
                  join=False, errorbar=("ci", 95), scale=1.2, order=plot_order)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_facecolor('white')

    counts = plot_data['Group'].value_counts()
    x_labels = [f"{g} (n={counts.get(g, 0)})" for g in plot_order]
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("alpha peak (z-score)", fontsize=18)
    ax.set_xlabel('')
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=14)

    fig.savefig(FIG_DIR / 'eeg_alpha_peak_autism_vs_td_violin.pdf', dpi=600, bbox_inches='tight')
    plt.close(fig)
    print("Saved: eeg_alpha_peak_autism_vs_td_violin.pdf")


def plot_alpha_peak_population1(df: pd.DataFrame):
    relevant_populations = ["Autism", "ID", "NT", "Relatives"]
    plot_data = df[df["PopulationS1"].isin(relevant_populations)].copy()
    if plot_data.empty:
        print("No data available for PopulationS1 plot")
        return

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    colors = {"Autism": "#5CAEE1", "Autism with ID": "#324095", "IDD": "#D8A4CB", "NT": "#C1C2BC", "Relatives": "#9AD5D3"}
    plot_data["PopulationS1"] = plot_data["PopulationS1"].replace({"ID": "IDD"})
    plot_order = ["NT", "Autism", "IDD", "Relatives"]

    sns.violinplot(y="alpha_peak_corrected", x="PopulationS1", data=plot_data, inner=None,
                   linewidth=0, order=plot_order, palette=colors)
    sns.stripplot(y="alpha_peak_corrected", x="PopulationS1", data=plot_data, color="black", size=2,
                  jitter=True, alpha=0.5, order=plot_order)
    sns.pointplot(y="alpha_peak_corrected", x="PopulationS1", data=plot_data, color="#8A0201",
                  join=False, errorbar=("ci", 95), scale=1.2, order=plot_order)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_facecolor('white')

    counts = plot_data['PopulationS1'].value_counts()
    x_labels = [f"{pop} (n={counts.get(pop, 0)})" for pop in plot_order if pop in counts.index]
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    ax.set_ylabel("alpha peak (z-score)", fontsize=18)
    ax.set_xlabel('')
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    fig.savefig(FIG_DIR / 'eeg_alpha_peak_population1_violin.pdf', dpi=600, bbox_inches='tight')
    plt.close(fig)
    print("Saved: eeg_alpha_peak_population1_violin.pdf")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    alpha_peak = load_alpha_peak_autism_td()

    # Save merged dataframe
    alpha_peak.to_csv(DATA_DIR / 'df_alpha_peak_autism_td.csv', index=False)

    # Stats
    stats_df = autism_vs_td_stats(alpha_peak, "alpha_peak_corrected")
    stats_df.to_csv(FIG_DIR / 'alpha_peak_autism_vs_td_stats.csv', index=False)
    print(f"Stats saved. Significant: {stats_df['significant'].sum() if len(stats_df) > 0 else 0}")

    # Plots
    plot_alpha_peak_autism_vs_td(alpha_peak)
    plot_alpha_peak_population1(alpha_peak)

    print(f"\nAll outputs saved to: {FIG_DIR}")


if __name__ == "__main__":
    sys.exit(main())
