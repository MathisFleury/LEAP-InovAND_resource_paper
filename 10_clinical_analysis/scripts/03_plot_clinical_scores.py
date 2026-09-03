#!/usr/bin/env python3
# =============================================================================
# 03 - Clinical scores by cluster (paper panel C)
# =============================================================================
# Cluster-wise distributions of four clinical scores:
#   full-scale IQ, SRS-2 total t-score, SSP total score, VABS-II ABC score.
# Violins per group in NT, C1, C2, C3, IDD order (NT far left, IDD far right),
# styled like psychomotor_milestones_combined.pdf. Two horizontal dashed lines
# mark the whole-cohort medians of autistic people (with or without IDD; blue)
# and neurotypical people (grey). Pairwise Mann-Whitney-Wilcoxon between groups
# is annotated with stacked brackets; the full table is written to tables/.
#
# Input:  individuals_metrics.tsv (scores)  +  df_clusters_complete_kmeans.csv
# Output: figures/clinical_scores_by_cluster.pdf
#         tables/clinical_scores_mww.csv
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    CLUSTERS_CURATED_FILE, FIGURES_DIR, TABLES_DIR,
    PALETTE_CLUSTERS, ORDER_CLUSTERS, CLINICAL_SCORES,
    COLOR_AUTISM_MEDIAN, COLOR_NT_MEDIAN,
)


def load_data():
    """Clinical scores + k-means cluster labels from the CURATED cluster table
    (curated clinical IQ/SRS/SSP/VABS already merged with the Cluster label).
    Builds a Cluster column spanning NT, C1, C2, C3, IDD."""
    df = pd.read_csv(CLUSTERS_CURATED_FILE, low_memory=False)
    df = df.drop_duplicates("ID").replace({999: np.nan, 998: np.nan})
    # curated already labels NT/IDD (not TD/ID), but normalise defensively.
    df["PopulationS1"] = df["PopulationS1"].replace({"TD": "NT", "ID": "IDD"})

    # NT/IDD participants get their own group; autism keeps its C1/C2/C3 label
    df.loc[df["PopulationS1"] == "NT", "Cluster"] = "NT"
    df.loc[df["PopulationS1"] == "IDD", "Cluster"] = "IDD"
    return df


def compute_stats(df, col, groups):
    """Pairwise Mann-Whitney-Wilcoxon between groups for one score."""
    data = {g: df.loc[df["Cluster"] == g, col].dropna().values for g in groups}
    data = {g: v for g, v in data.items() if len(v) >= 3}
    rows, grps = [], list(data)
    for i in range(len(grps)):
        for j in range(i + 1, len(grps)):
            g1, g2 = grps[i], grps[j]
            u, p = mannwhitneyu(data[g1], data[g2], alternative="two-sided")
            rows.append({
                "score": col, "group1": g1, "n1": len(data[g1]), "median1": np.median(data[g1]),
                "group2": g2, "n2": len(data[g2]), "median2": np.median(data[g2]),
                "U_statistic": u, "pvalue": p,
            })
    return rows


def plot():
    df = load_data()
    groups = [g for g in ORDER_CLUSTERS if (df["Cluster"] == g).sum() >= 2]
    scores = list(CLINICAL_SCORES)

    is_autism = df["PopulationS1"] == "Autism"
    is_nt = df["PopulationS1"] == "NT"

    # --- publication sizing (≤ 9 cm wide) ---
    TITLE_FS, SMALL_FS = 13, 7
    CM = 1 / 2.54
    NR = 200 - 20  # nominal score range, for proportional offsets/steps

    # Tight x geometry: violins adjacent within a measure, small gap between measures
    intra, gap = 1.0, 1.2
    block = len(groups) * intra + gap
    xpos = {(col, g): j * intra + i * block
            for i, col in enumerate(scores) for j, g in enumerate(groups)}

    # 18 cm wide (double-column); height leaves room for the significant-pairwise
    # connector stacks (up to 10 per measure)
    fig, ax = plt.subplots(1, 1, figsize=(18 * CM, 13 * CM))
    all_stats = []

    # --- violins ---
    for i, col in enumerate(scores):
        for g in groups:
            vals = df.loc[df["Cluster"] == g, col].dropna().values
            if len(vals) < 2:
                continue
            x = xpos[(col, g)]
            vp = ax.violinplot([vals], positions=[x], widths=0.9 * intra,
                               showmeans=False, showmedians=False)
            vp["bodies"][0].set_facecolor(PALETTE_CLUSTERS[g])
            vp["bodies"][0].set_alpha(0.7)
            vp["bodies"][0].set_edgecolor("black")
            vp["bodies"][0].set_linewidth(1.1)
            for key in ("cmaxes", "cmins", "cbars"):
                if key in vp:
                    vp[key].set_visible(False)
            q25, q50, q75 = np.percentile(vals, [25, 50, 75])
            ax.plot([x, x], [q25, q75], "k-", lw=0.8)
            ax.plot(x, q50, "ko", ms=2.5)

    # --- per-measure autism/NT median segments + n= labels + MWW brackets ---
    for i, col in enumerate(scores):
        xs = [xpos[(col, g)] for g in groups]
        x0, x1 = min(xs) - 0.5 * intra, max(xs) + 0.5 * intra


        present = {g: v for g in groups
                   if len((v := df.loc[df["Cluster"] == g, col].dropna().values)) >= 2}
        for g, v in present.items():
            n_color = "black" if g == "NT" else PALETTE_CLUSTERS[g]
            ax.text(xpos[(col, g)], np.percentile(v, 98) + NR * 0.015,
                    f"n={len(v)}", ha="center", va="bottom", fontsize=SMALL_FS,
                    fontweight="bold", color=n_color)

        # pairwise MWW connectors (significant only), stacked above the group
        stats = compute_stats(df, col, groups)
        all_stats.extend(stats)
        sig = [s for s in stats if s["pvalue"] < 0.05 and s["group1"] in present and s["group2"] in present]
        sig.sort(key=lambda s: abs(xpos[(col, s["group1"])] - xpos[(col, s["group2"])]))
        local_max = max(np.percentile(v, 98) for v in present.values()) if present else 200
        for k, s in enumerate(sig):
            xa, xb = xpos[(col, s["group1"])], xpos[(col, s["group2"])]
            y = local_max + NR * 0.05 * (k + 1)
            ax.plot([xa, xb], [y, y], "k-", lw=0.7)  # straight connector
            ax.text((xa + xb) / 2, y, f"{s['pvalue']:.1e}", ha="center", va="bottom", fontsize=SMALL_FS)

        autism_med = df.loc[is_autism, col].dropna().median()
        nt_med = df.loc[is_nt, col].dropna().median()
        if pd.notna(autism_med):
            ax.plot([x0, x1], [autism_med, autism_med], ls="--", color=COLOR_AUTISM_MEDIAN, lw=1.8, zorder=3)
        if pd.notna(nt_med):
            ax.plot([x0, x1], [nt_med, nt_med], ls="--", color=COLOR_NT_MEDIAN, lw=1.8, zorder=3)

    # measure titles centred below each group on the x-axis
    short = {"IQ": "Measured IQ \n values", "SRS_tscore": "SRS-2 \n t-score", "ssp_total": "SSP total", "vabsabcabc_standard": "VABS-II \n ABC score"}
    ax.set_xticks([(len(groups) - 1) / 2 * intra + i * block for i in range(len(scores))])
    ax.set_xticklabels([short.get(c, c) for c in scores], fontweight="bold", fontsize=SMALL_FS)
    ax.set_ylabel("Scale score", fontsize=SMALL_FS + 1)
    ax.tick_params(axis="y", labelsize=SMALL_FS)
    # group labels: horizontal, just above the bottom axis line, over each violin
    for (col, g), x in xpos.items():
        ax.text(x, 15 + NR * 0.012, g, ha="center", va="bottom", fontsize=SMALL_FS,
                # color=("black" if g == "NT" else PALETTE_CLUSTERS[g]), 
                rotation=0,
                bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.75))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # thicker x/y axis lines + ticks
    for sp in ("bottom", "left"):
        ax.spines[sp].set_linewidth(1.6)
    ax.tick_params(axis="both", width=1.6)
    # tick labels stop at 200; axis extends above only to fit the stat connectors
    ax.set_ylim(5, max(200, ax.get_ylim()[1]))
    ax.set_yticks(np.arange(25, 201, 25))

    # compact legend below the axis: group patches + reference-line entries
    handles = [Patch(facecolor=PALETTE_CLUSTERS[g], alpha=0.7, edgecolor="black", label=g) for g in groups]
    handles += [
        Line2D([0], [0], ls="--", lw=1.8, color=COLOR_AUTISM_MEDIAN, label="Autistic median"),
        Line2D([0], [0], ls="--", lw=1.8, color=COLOR_NT_MEDIAN, label="NT median"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06),
              ncol=4, fontsize=SMALL_FS, frameon=False, handlelength=1.2,
              columnspacing=1.0, handletextpad=0.4)

    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "clinical_scores_by_cluster.pdf")
    plt.savefig(out, dpi=500, format="pdf", transparent=True, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

    df_stats = pd.DataFrame(all_stats)
    tpath = os.path.join(TABLES_DIR, "clinical_scores_mww.csv")
    df_stats.to_csv(tpath, index=False)
    print(f"  Saved: {tpath}")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(TABLES_DIR, exist_ok=True)
    plot()


if __name__ == "__main__":
    main()
