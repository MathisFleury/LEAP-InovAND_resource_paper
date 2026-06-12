#!/usr/bin/env python3
# =============================================================================
# 02 - Verbal vs Non-verbal Status by Cluster
# =============================================================================
# Stacked bar chart showing the proportion of verbal vs non-verbal participants
# across clinical clusters, with significance brackets (Fisher exact, FDR)
# and Unknown/NA counts.
#
# Input:  df_verbalornot.tsv         (verbal/non-verbal classification)
#         df_clusters_complete_kmeans.csv   (cluster assignments)
# Output: figures/  verbal status bar plot with significance brackets
#         tables/   chi-squared statistics + counts + pairwise + cohort breakdown
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    VERBAL_TSV, CLUSTERS_FILE,
    FIGURES_DIR, TABLES_DIR,
    PALETTE_CLUSTERS, ORDER_CLUSTERS,
    VERBAL_LABELS, PALETTE_VERBAL,
)

CATEGORIES = ["Verbal", "Speech impairment/non verbal"]


def load_and_merge(verbal_path, cluster_path):
    """Load verbal status and cluster assignments, merge on ID."""
    df_verbal = pd.read_csv(verbal_path, sep="\t")
    df_clusters = pd.read_csv(cluster_path, low_memory=False)

    df_verbal["ID"] = df_verbal["ID"].astype(str)
    df_clusters["ID"] = df_clusters["ID"].astype(str)

    df = df_verbal.merge(
        df_clusters[["ID", "Cluster", "PopulationS1", "cohort"]],
        on="ID",
        how="inner",
    )
    df.loc[df["PopulationS1"] == "TD", "Cluster"] = "NT"
    df.loc[df["PopulationS1"] == "ID", "Cluster"] = "IDD"

    print(f"Loaded {len(df_verbal)} verbal rows, {len(df_clusters)} cluster rows")
    print(f"Merged: {len(df)} participants")
    return df


def compute_stats(df_known, groups_with_data):
    """Chi-squared omnibus + pairwise Fisher exact with FDR correction."""
    contingency = []
    for g in groups_with_data:
        grp = df_known[df_known["Cluster"] == g]
        row = [len(grp[grp["verbal_ornot1"] == cat]) for cat in CATEGORIES]
        contingency.append(row)
    contingency = np.array(contingency)

    chi2, p_chi2, dof, _ = chi2_contingency(contingency)

    pairwise = []
    for i in range(len(groups_with_data)):
        for j in range(i + 1, len(groups_with_data)):
            g1, g2 = groups_with_data[i], groups_with_data[j]
            table_2x2 = contingency[[i, j], :]
            odds_ratio, p_fisher = fisher_exact(table_2x2)
            pairwise.append({
                "group1": g1, "group2": g2,
                "odds_ratio": odds_ratio, "pvalue": p_fisher,
            })

    df_pairwise = pd.DataFrame(pairwise)
    if len(df_pairwise) > 0:
        _, pvals_corr, _, _ = multipletests(df_pairwise["pvalue"], method="fdr_bh")
        df_pairwise["pvalue_fdr"] = pvals_corr

    rows = []
    for g in groups_with_data:
        grp = df_known[df_known["Cluster"] == g]
        total = len(grp)
        for cat in CATEGORIES:
            n = len(grp[grp["verbal_ornot1"] == cat])
            rows.append({
                "cluster": g, "status": cat,
                "n": n, "pct": n / total * 100 if total > 0 else 0,
                "total": total,
            })
    df_counts = pd.DataFrame(rows)
    df_counts["chi2"] = chi2
    df_counts["chi2_p"] = p_chi2
    df_counts["chi2_dof"] = dof

    return df_counts, df_pairwise


def compute_cohort_breakdown(df, groups_with_data):
    """Proportion of verbal/non-verbal per cluster, split by cohort (LEAP vs InovAND).
    Returns a DataFrame with % computed on known (non-Unknown) participants only."""
    known = df[df["verbal_ornot1"].isin(CATEGORIES)].copy()
    rows = []
    for g in groups_with_data:
        for cohort in ["LEAP", "INOVAND"]:
            sub = known[(known["Cluster"] == g) & (known["cohort"] == cohort)]
            total = len(sub)
            for cat in CATEGORIES:
                n = len(sub[sub["verbal_ornot1"] == cat])
                rows.append({
                    "cluster": g, "cohort": cohort, "status": cat,
                    "n": n, "pct": n / total * 100 if total > 0 else 0,
                    "total_known": total,
                })
    return pd.DataFrame(rows)


def plot_verbal_proportions(df, groups, df_pairwise, output_path):
    """Stacked bar chart with significance brackets and NA counts."""
    df_known = df[df["verbal_ornot1"].isin(CATEGORIES)].copy()

    if len(df_known) == 0:
        print("  No known verbal status data to plot.")
        return

    # Compute counts per group (known + unknown)
    counts = {}
    for g in groups:
        grp_all = df[df["Cluster"] == g]
        grp_known = df_known[df_known["Cluster"] == g]
        vc = grp_known["verbal_ornot1"].value_counts()
        n_known = len(grp_known)
        n_unknown = len(grp_all) - n_known
        counts[g] = {cat: vc.get(cat, 0) for cat in CATEGORIES}
        counts[g]["total"] = n_known
        counts[g]["unknown"] = n_unknown

    groups_with_data = [g for g in groups if counts.get(g, {}).get("total", 0) > 0]
    if len(groups_with_data) == 0:
        print("  No groups with known verbal status data.")
        return

    x_pos = {g: i for i, g in enumerate(groups_with_data)}
    x = np.arange(len(groups_with_data))
    width = 0.6

    fig, ax = plt.subplots(figsize=(5, 5.5))

    # Stacked bars
    bottom = np.zeros(len(groups_with_data))
    for cat in CATEGORIES:
        proportions = np.array([
            counts[g][cat] / counts[g]["total"] * 100 if counts[g]["total"] > 0 else 0
            for g in groups_with_data
        ])
        ax.bar(x, proportions, width, bottom=bottom,
               label=VERBAL_LABELS.get(cat, cat),
               color=PALETTE_VERBAL[cat], edgecolor="black", linewidth=0.6)
        bottom += proportions

    # Count annotations inside bars
    for i, g in enumerate(groups_with_data):
        total = counts[g]["total"]
        cum = 0
        for cat in CATEGORIES:
            prop = counts[g][cat] / total * 100 if total > 0 else 0
            if prop > 8:
                ax.text(i, cum + prop / 2,
                        f"{counts[g][cat]}\n({prop:.0f}%)",
                        ha="center", va="center", fontsize=8, fontweight="bold",
                        color="white" if cat == "Speech impairment/non verbal" else "black")
            cum += prop

    # x labels with n known + n unknown
    ax.set_xticks(x)
    xlabels = [f"{g}\n(n={counts[g]['total']}, NA={counts[g]['unknown']})"
               for g in groups_with_data]
    ax.set_xticklabels(xlabels, fontweight="bold", fontsize=9)

    ax.set_ylabel("Proportion (%)", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Significance brackets
    if df_pairwise is not None and len(df_pairwise) > 0:
        sig = df_pairwise[df_pairwise["pvalue_fdr"] < 0.05].copy()
        if len(sig) > 0:
            y_base = 105
            step = 7
            for k, (_, row) in enumerate(sig.iterrows()):
                g1, g2 = row["group1"], row["group2"]
                if g1 not in x_pos or g2 not in x_pos:
                    continue
                x1, x2 = x_pos[g1], x_pos[g2]
                y_line = y_base + k * step

                ax.plot([x1, x2], [y_line, y_line], "k-", linewidth=1.5)
                ax.text((x1 + x2) / 2, y_line + step * 0.15,
                        f"{row['pvalue_fdr']:.2e}",
                        ha="center", va="bottom", fontsize=8, fontweight="bold")

            ax.set_ylim(0, y_base + len(sig) * step + 8)
        else:
            ax.set_ylim(0, 105)
    else:
        ax.set_ylim(0, 105)

    ax.legend(loc="upper right", frameon=True, fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=500, format="pdf",
                transparent=True, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("=" * 60)
    print("VERBAL STATUS BY CLUSTER")
    print("=" * 60)

    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(TABLES_DIR, exist_ok=True)

    # Load and merge
    df = load_and_merge(VERBAL_TSV, CLUSTERS_FILE)

    # Report distribution
    print("\nVerbal status distribution (all merged):")
    print(df["verbal_ornot1"].value_counts().to_string())

    # Groups present
    groups = [g for g in ORDER_CLUSTERS if g in df["Cluster"].unique()]
    print(f"\nGroups present: {groups}")

    # Filter to known status
    known = df[df["verbal_ornot1"].isin(CATEGORIES)]
    groups_with_data = [g for g in groups if len(known[known["Cluster"] == g]) > 0]
    print(f"Participants with known verbal status: {len(known)}")
    for g in groups_with_data:
        grp_all = df[df["Cluster"] == g]
        n_known = len(known[known["Cluster"] == g])
        n_unknown = len(grp_all) - n_known
        print(f"  {g}: known={n_known}, unknown/NA={n_unknown}")

    # Statistics
    print("\nComputing statistics...")
    df_counts, df_pairwise = compute_stats(known, groups_with_data)

    if len(df_counts) > 0:
        counts_path = os.path.join(TABLES_DIR, "verbal_status_counts.csv")
        df_counts.to_csv(counts_path, index=False)
        print(f"  Saved: {counts_path}")

        chi2_val = df_counts["chi2"].iloc[0]
        chi2_p = df_counts["chi2_p"].iloc[0]
        print(f"\n--- Statistical Summary ---")
        print(f"  Chi-squared = {chi2_val:.2f}, p = {chi2_p:.4g}")

    if len(df_pairwise) > 0:
        pairwise_path = os.path.join(TABLES_DIR, "verbal_status_pairwise.csv")
        df_pairwise.to_csv(pairwise_path, index=False)
        print(f"  Saved: {pairwise_path}")

        sig = df_pairwise[df_pairwise["pvalue_fdr"] < 0.05]
        for _, row in sig.iterrows():
            print(f"    {row['group1']} vs {row['group2']}: "
                  f"OR = {row['odds_ratio']:.2f}, p_fdr = {row['pvalue_fdr']:.4g}")

    # Cohort breakdown
    print("\n--- Cohort Breakdown (% of known, non-NA) ---")
    df_cohort = compute_cohort_breakdown(df, groups_with_data)
    if len(df_cohort) > 0:
        cohort_path = os.path.join(TABLES_DIR, "verbal_status_cohort_breakdown.csv")
        df_cohort.to_csv(cohort_path, index=False)
        print(f"  Saved: {cohort_path}")

        for g in groups_with_data:
            print(f"\n  {g}:")
            for cohort in ["LEAP", "INOVAND"]:
                sub = df_cohort[(df_cohort["cluster"] == g) & (df_cohort["cohort"] == cohort)]
                if sub.empty or sub["total_known"].iloc[0] == 0:
                    print(f"    {cohort}: no known data")
                    continue
                total = sub["total_known"].iloc[0]
                for _, row in sub.iterrows():
                    print(f"    {cohort}: {row['status']} = {row['n']}/{total} ({row['pct']:.1f}%)")

    # Plot (with significance brackets + NA counts)
    print("\nCreating verbal status plot...")
    plot_verbal_proportions(
        df, groups, df_pairwise,
        os.path.join(FIGURES_DIR, "verbal_status_by_cluster.pdf"),
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
