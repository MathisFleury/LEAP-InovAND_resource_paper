#!/usr/bin/env python3
"""
EEG Alpha Peak — Autism vs NT.

Consumes the curated, already-corrected alpha peak built by
  preprocessing/build_alpha_peak_corrected.py
(raw $IMG5 rebuild + curated demographics + full-sample regress/z-score).
This script does NOT correct again — it only runs stats and plots.
"""

import sys
from pathlib import Path
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent  # 8_eeg_analysis/
sys.path.insert(0, str(_SECTION_DIR / "preprocessing"))
from config import ALPHA_PEAK_CORRECTED  # noqa: E402

FIG_DIR = _SECTION_DIR / 'outputs' / 'figures'
TABLES_DIR = _SECTION_DIR / 'outputs' / 'tables'
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)


def load_alpha_peak_autism_td():
    alpha_peak = pd.read_csv(ALPHA_PEAK_CORRECTED)
    print(f"Loaded {len(alpha_peak)} subjects from {ALPHA_PEAK_CORRECTED.name}")
    alpha_peak["ID"] = alpha_peak["ID"].astype(str)
    alpha_peak = alpha_peak[alpha_peak["alpha_peak_corrected"].notna()]
    print("Alpha peak counts by PopulationS1:",
          alpha_peak["PopulationS1"].value_counts(dropna=False).to_dict())
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
    relevant_populations = ["Autism", "IDD", "NT", "Relatives"]
    plot_data = df[df["PopulationS1"].isin(relevant_populations)].copy()
    if plot_data.empty:
        print("No data available for PopulationS1 plot")
        return

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    colors = {"Autism": "#5CAEE1", "Autism with ID": "#324095", "IDD": "#D8A4CB", "NT": "#C1C2BC", "Relatives": "#9AD5D3"}
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

    alpha_peak = load_alpha_peak_autism_td()  # already corrected in preprocessing

    # Stats
    stats_df = autism_vs_td_stats(alpha_peak, "alpha_peak_corrected")
    stats_df.to_csv(TABLES_DIR / 'alpha_peak_autism_vs_td_stats.csv', index=False)
    print(f"Stats saved. Significant: {stats_df['significant'].sum() if len(stats_df) > 0 else 0}")

    # Plots
    plot_alpha_peak_autism_vs_td(alpha_peak)
    plot_alpha_peak_population1(alpha_peak)

    print(f"\nAll outputs saved to: {FIG_DIR}")


if __name__ == "__main__":
    sys.exit(main())
