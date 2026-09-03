#!/usr/bin/env python3
# =============================================================================
# 09 - Clinical scores, LEAP vs INOVAND, Autism vs NT
# =============================================================================
# Updates supplementary_figure_clinical_scores_leap_inovand_by_population.pdf
# (originally built in eeg_mri-pipeline/analysis/figures_papers/... from the
# frozen df_clusters_complete.csv) with the current CURATED clinical data.
#
# 13 clinical measures x {INOVAND, LEAP} x {Autism, NT}, one violin panel per
# measure, same ragged 4/3/3/3 row layout and score wording as the original
# figure so it drops in as a like-for-like replacement.
#
# Input:  LEAP_CLINICAL_CURATED (t1) + INOVAND_CLINICAL_CURATED, participants
#         only, Autism/NT only (matches the original figure's scope).
# Output: figures/clinical_scores_leap_inovand_by_population.pdf
# =============================================================================
import os
import sys

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import LEAP_CLINICAL_CURATED, INOVAND_CLINICAL_CURATED, FIGURES_DIR, COLORS

# column -> display label, in the exact panel order of the original figure
SCORE_LABELS = {
    "age_yrs": "Age (years)",
    "total_IQ": "Full Scale IQ",
    "verbal_IQ": "Verbal IQ",
    "performance_IQ": "Non Verbal IQ",
    "SRS_tscore": "SRS-2 t-scores",
    "ssp_total": "Short Sensory Profile",
    "RBS-R_total": "RBS-R Total",
    "ADI_communication": "ADI Communication",
    "ADI_rrb": "ADI RRB",
    "ADOS_total_css": "ADOS Total CSS",
    "vabsdscoresc_dss": "VABS-II Communication",
    "vabsdscoresd_dss": "VABS-II Daily Living Skills",
    "vabsdscoress_dss": "VABS-II Socialization",
}
ROW_SPLITS = [4, 3, 3, 3]  # panels per row -- matches the original figure's layout
GRID_COLS = 12             # LCM(4, 3): equal panel width whether a row holds 4 or 3

COHORT_ORDER = ["INOVAND", "LEAP"]
POP_ORDER = ["Autism", "NT"]
SENTINELS = [999, 998, 777]


def load_data():
    leap = pd.read_csv(LEAP_CLINICAL_CURATED, sep="\t", low_memory=False)
    leap["cohort"] = "LEAP"
    inovand = pd.read_csv(INOVAND_CLINICAL_CURATED, sep="\t", low_memory=False)
    inovand["cohort"] = "INOVAND"
    df = pd.concat([leap, inovand], ignore_index=True)
    df = df[df["relation_to_proposant"] == "participant"]
    df = df[df["population_group"].isin(POP_ORDER)].copy()

    for col in SCORE_LABELS:
        df[col] = pd.to_numeric(df[col], errors="coerce").replace(SENTINELS, np.nan)
    # same outlier clips as the original figure (data-entry artifacts, not real values)
    df.loc[df["age_yrs"] >= 70, "age_yrs"] = np.nan
    df.loc[df["ADI_communication"] >= 60, "ADI_communication"] = np.nan
    df.loc[df["ADOS_total_css"] >= 10, "ADOS_total_css"] = np.nan
    return df


def plot_panel(ax, df, score, label):
    sub = df.dropna(subset=[score])
    sns.violinplot(
        x="cohort", y=score, hue="population_group", data=sub, ax=ax,
        order=COHORT_ORDER, hue_order=POP_ORDER,
        palette={"Autism": COLORS["Autism"], "NT": COLORS["TD"]},
        inner="box", linewidth=1, cut=0,
    )
    if ax.legend_ is not None:
        ax.legend_.remove()

    counts = sub.groupby(["cohort", "population_group"]).size()
    xlabels = [
        f"{cohort}\n(Autism: n={counts.get((cohort, 'Autism'), 0)}, "
        f"NT: n={counts.get((cohort, 'NT'), 0)})"
        for cohort in COHORT_ORDER
    ]
    ax.set_xticks(range(len(COHORT_ORDER)))
    ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(label, fontsize=10, fontweight="bold", pad=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot(df):
    fig = plt.figure(figsize=(16, 13))
    gs = gridspec.GridSpec(len(ROW_SPLITS), GRID_COLS, figure=fig,
                            left=0.04, right=0.98, top=0.93, bottom=0.04,
                            hspace=0.55, wspace=0.35)

    scores = list(SCORE_LABELS.items())
    i = 0
    for r, n_panels in enumerate(ROW_SPLITS):
        span = GRID_COLS // n_panels
        for c in range(n_panels):
            score, label = scores[i]
            i += 1
            ax = fig.add_subplot(gs[r, c * span:(c + 1) * span])
            plot_panel(ax, df, score, label)

    handles = [Patch(facecolor=COLORS["Autism"], label="Autism"),
               Patch(facecolor=COLORS["TD"], label="NT")]
    fig.legend(handles=handles, title="Population", loc="upper right",
               bbox_to_anchor=(0.98, 0.985), fontsize=9, title_fontsize=10, frameon=True)
    return fig


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    df = load_data()
    fig = plot(df)
    out = os.path.join(FIGURES_DIR, "clinical_scores_leap_inovand_by_population.pdf")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
