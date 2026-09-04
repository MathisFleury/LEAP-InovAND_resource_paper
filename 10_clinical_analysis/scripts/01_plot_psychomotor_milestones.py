#!/usr/bin/env python3
# =============================================================================
# 01 - Psychomotor Milestones by Cluster
# =============================================================================
# Violin plots of age at psychomotor milestones (sitting, walking, first words,
# first phrases) across clinical clusters (NT, C1, C2, C3, IDD).
#
# Clustering methodology (curated priority regime): K-means (n_init=25,
# random_state=42, k=3) on standardised IQ x SRS-2 (IQ = total_IQ fillna
# performance_IQ), with LEAP wave backfill T1->T2->T3. Written by
# 1_clustering/scripts/1_run_clustering.py -> CLUSTERS_CURATED_FILE.
# NT and IDD are assigned from the broad curated clinical rosters
# (population_group), NOT from the k-means fit — those rosters are not gated on
# IQ+SRS completeness, so controls with milestone data but no IQ/SRS are kept.
# Ward/GMM are retained only as sensitivity analyses (see CLAUDE.md /
# 6_method_comparison.py).
#
# Input:  concat_fonda_FIRST_ACQUISITIONS_251112.csv  (InovAND first acquisitions)
#         individuals_metrics_with_clusters_curated.csv       (curated C1/C2/C3 labels)
#         INOVAND/LEAP_clinical_curated*.tsv                  (NT/IDD membership)
# Output: figures/  milestone violin plots (individual + combined)
#         tables/   Kruskal-Wallis + pairwise statistics
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    FONDA_CSV, CLUSTERS_CURATED_FILE, LEAP_ADIR_CSV, LEAP_TO_FONDA,
    INOVAND_CLINICAL_CURATED, LEAP_CLINICAL_CURATED,
    FIGURES_DIR, TABLES_DIR,
    PALETTE_CLUSTERS, ORDER_CLUSTERS,
    MILESTONE_COLS, MILESTONE_LABELS,
)


def build_roster(clusters_path, inovand_clinical, leap_clinical):
    """ID -> (PopulationS1, Cluster) roster for the milestones merge.

    NT/IDD membership comes from the broad curated clinical rosters
    (population_group), which are NOT gated on IQ+SRS completeness. The curated
    k-means clustering table only supplies the autistic C1/C2/C3 labels — using
    it to define NT/IDD would silently drop controls that lack IQ/SRS but do
    have milestone data (that gating cost ~80 NT).
    """
    frames = []
    for f in (inovand_clinical, leap_clinical):
        if os.path.exists(f):
            d = pd.read_csv(f, sep="\t", low_memory=False)[["ID", "population_group"]]
            frames.append(d)
        else:
            print(f"WARNING: broad clinical roster not found: {f}")
    roster = pd.concat(frames, ignore_index=True).dropna(subset=["ID"])
    roster["ID"] = roster["ID"].astype(str)
    roster = roster.rename(columns={"population_group": "PopulationS1"})
    roster = roster.drop_duplicates("ID")

    clusters = pd.read_csv(clusters_path, low_memory=False)
    clusters["ID"] = clusters["ID"].astype(str)
    roster = roster.merge(
        clusters[["ID", "Cluster"]].drop_duplicates("ID"), on="ID", how="left",
    )
    print(f"Roster: {len(roster)} participants "
          f"(NT={(roster['PopulationS1']=='NT').sum()}, "
          f"IDD={(roster['PopulationS1']=='IDD').sum()}, "
          f"clustered={roster['Cluster'].notna().sum()})")
    return roster


def load_and_merge(fonda_path, df_clusters):
    """Load FONDA (InovAND) milestones + LEAP ADIR milestones, merge with the
    roster (ID, PopulationS1, Cluster) built by build_roster().

    Uses CODE_MERGE (not CODE_PATIENT) because 161 AURD participants
    have CODE_PATIENT=NaN but a valid CODE_MERGE.
    LEAP participants come from a separate ADIR CSV with different column
    names; their milestone values are mapped to the FONDA column names
    (AGEPAS, AGEMOTS, AGEPHRAS).
    """
    df_fonda = pd.read_csv(fonda_path, sep=";")
    df_fonda["CODE_MERGE"] = df_fonda["CODE_MERGE"].astype(str)
    df_clusters["ID"] = df_clusters["ID"].astype(str)

    # --- InovAND merge (FONDA) ---
    df = df_fonda.merge(
        df_clusters[["ID", "Cluster", "PopulationS1"]],
        left_on="CODE_MERGE",
        right_on="ID",
        how="inner",
    )
    print(f"Loaded {len(df_fonda)} FONDA rows, {len(df_clusters)} cluster rows")
    print(f"InovAND merge (FONDA): {len(df)} participants")

    # --- LEAP merge (ADIR early development) ---
    if os.path.exists(LEAP_ADIR_CSV):
        df_leap = pd.read_csv(LEAP_ADIR_CSV)
        df_leap["Study_ID"] = df_leap["Study_ID"].astype(str)

        # Keep only the milestone columns + Study_ID
        leap_cols = [c for c in LEAP_TO_FONDA if c in df_leap.columns] + ["Study_ID"]
        df_leap_sub = df_leap[leap_cols].copy()

        # Rename to FONDA column names
        rename = {k: v for k, v in LEAP_TO_FONDA.items() if k in df_leap_sub.columns}
        df_leap_sub = df_leap_sub.rename(columns=rename)

        # Merge with cluster assignments
        df_leap_merged = df_leap_sub.merge(
            df_clusters[["ID", "Cluster", "PopulationS1"]],
            left_on="Study_ID", right_on="ID", how="inner",
        )
        print(f"LEAP ADIR: {len(df_leap_merged)} participants with cluster data")

        # Add LEAP participants that are NOT already in df
        existing_ids = set(df["ID"].unique())
        df_leap_new = df_leap_merged[~df_leap_merged["ID"].isin(existing_ids)].copy()

        if len(df_leap_new) > 0:
            # Align columns: add missing FONDA columns as NaN
            for col in df.columns:
                if col not in df_leap_new.columns:
                    df_leap_new[col] = np.nan
            df = pd.concat([df, df_leap_new[df.columns]], ignore_index=True)
            print(f"Added {len(df_leap_new)} new LEAP participants")

        # Update existing LEAP participants (fill NaN milestone values)
        leap_in_existing = df_leap_merged[df_leap_merged["ID"].isin(existing_ids)]
        if len(leap_in_existing) > 0:
            fonda_cols = list(LEAP_TO_FONDA.values())
            leap_dict = leap_in_existing.set_index("ID")[fonda_cols].to_dict("index")
            for idx in df[df["ID"].isin(leap_dict)].index:
                pid = df.loc[idx, "ID"]
                for col in fonda_cols:
                    if pd.isna(df.loc[idx, col]) and col in leap_dict[pid]:
                        val = leap_dict[pid][col]
                        if pd.notna(val):
                            df.loc[idx, col] = val
            print(f"Updated {len(leap_in_existing)} existing participants with LEAP data")
    else:
        print(f"LEAP ADIR file not found: {LEAP_ADIR_CSV}")

    # Assign neurotypical -> NT, intellectual-disability -> IDD (from PopulationS1,
    # not the k-means fit). Curated PopulationS1 uses NT/IDD; frozen uses TD/ID.
    df.loc[df["PopulationS1"].isin(["TD", "NT"]), "Cluster"] = "NT"
    df.loc[df["PopulationS1"].isin(["ID", "IDD"]), "Cluster"] = "IDD"

    print(f"Final merged: {len(df)} participants (InovAND + LEAP)")
    return df


def clean_milestone(series):
    """Return cleaned numeric milestone values (months), dropping invalids.

    Note: value 9 is a valid age in months for these columns (e.g. sitting
    at 9 months).  The code 9='unknown' only applies to binary yes/no columns
    (RETPSY, RETLANG …), not to continuous age columns.
    """
    s = pd.to_numeric(series, errors="coerce")
    s = s.replace({999: np.nan, 998: np.nan, 777: np.nan, 0: np.nan})
    s[s > 100] = np.nan
    return s


def compute_stats(df, col, groups):
    """Kruskal-Wallis omnibus + pairwise Mann-Whitney with FDR correction."""
    group_data = {g: df.loc[df["Cluster"] == g, col].dropna().values for g in groups}
    group_data = {g: v for g, v in group_data.items() if len(v) >= 3}

    if len(group_data) < 2:
        return pd.DataFrame()

    # Omnibus
    stat_kw, p_kw = kruskal(*group_data.values())

    # Pairwise
    pairs = []
    grp_list = list(group_data.keys())
    for i in range(len(grp_list)):
        for j in range(i + 1, len(grp_list)):
            g1, g2 = grp_list[i], grp_list[j]
            stat_u, p_u = mannwhitneyu(group_data[g1], group_data[g2], alternative="two-sided")
            pairs.append({
                "milestone": col,
                "group1": g1, "n1": len(group_data[g1]),
                "median1": np.median(group_data[g1]),
                "group2": g2, "n2": len(group_data[g2]),
                "median2": np.median(group_data[g2]),
                "U_statistic": stat_u,
                "pvalue": p_u,
            })

    df_stats = pd.DataFrame(pairs)
    if len(df_stats) > 0:
        _, pvals_corr, _, _ = multipletests(df_stats["pvalue"], method="fdr_bh")
        df_stats["pvalue_fdr"] = pvals_corr
    df_stats["kruskal_wallis_H"] = stat_kw
    df_stats["kruskal_wallis_p"] = p_kw
    return df_stats


def plot_single_milestone(df, col, label, groups, palette, output_path):
    """Violin plot for a single milestone across clusters."""
    fig, ax = plt.subplots(figsize=(4.5, 5))

    plot_data = []
    for g in groups:
        vals = df.loc[df["Cluster"] == g, col].dropna().values
        if len(vals) > 0:
            plot_data.append((g, vals))

    if len(plot_data) == 0:
        plt.close()
        return

    positions = list(range(len(plot_data)))
    for i, (g, vals) in enumerate(plot_data):
        vp = ax.violinplot([vals], positions=[i], widths=0.7,
                           showmeans=False, showmedians=False)
        vp["bodies"][0].set_facecolor(palette[g])
        vp["bodies"][0].set_alpha(0.75)
        vp["bodies"][0].set_edgecolor("black")
        vp["bodies"][0].set_linewidth(0.6)
        for key in ("cmaxes", "cmins", "cbars"):
            if key in vp:
                vp[key].set_visible(False)

        # Quartile lines + median dot
        q25, q50, q75 = np.percentile(vals, [25, 50, 75])
        ax.plot([i, i], [q25, q75], "k-", linewidth=1.2)
        ax.plot(i, q50, "ko", markersize=5, zorder=5)

        # Sample size
        ax.text(i, ax.get_ylim()[1] * 0.02 + np.max(vals),
                f"n={len(vals)}", ha="center", va="bottom",
                fontsize=8, fontweight="bold")

    ax.set_xticks(positions)
    ax.set_xticklabels([g for g, _ in plot_data], fontweight="bold")
    ax.set_ylabel(label, fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=500, format="pdf",
                transparent=True, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_combined_milestones(df, milestone_cols, labels, groups, palette, output_path):
    """Consolidated violin plot with all milestones on one axis, significance
    brackets and p-values.  Matches the style of
    fonda_first_acquisitions_autism_td.pdf from the eeg_mri-pipeline."""

    n_measures = len(milestone_cols)
    n_groups = len(groups)
    spacing = n_groups + 2  # gap between milestone groups

    fig, ax = plt.subplots(1, 1, figsize=(max(14, n_measures * 3.5), 6))
    colors = [palette[g] for g in groups]

    # --- draw violins --------------------------------------------------------
    for i, col in enumerate(milestone_cols):
        for j, group in enumerate(groups):
            vals = df.loc[df["Cluster"] == group, col].dropna().values
            if len(vals) < 2:
                continue
            x_pos = j + i * spacing

            vp = ax.violinplot(
                [vals], positions=[x_pos], widths=0.8,
                showmeans=False, showmedians=False,
            )
            vp["bodies"][0].set_facecolor(colors[j])
            vp["bodies"][0].set_alpha(0.7)
            for key in ("cmaxes", "cmins", "cbars"):
                if key in vp:
                    vp[key].set_visible(False)

            # IQR line + median dot
            q25, q50, q75 = np.percentile(vals, [25, 50, 75])
            ax.plot([x_pos, x_pos], [q25, q75], "k-", linewidth=1)
            ax.plot(x_pos, q50, "ko", markersize=5)

            # n= label
            y_top = np.percentile(vals, 98)
            ax.text(x_pos, y_top + 1, f"n={len(vals)}", ha="center",
                    va="bottom", fontweight="bold", fontsize=8)

    # --- significance brackets -----------------------------------------------
    # first pass: collect all y-tops per milestone for bracket positioning
    y_min_all, y_max_all = ax.get_ylim()

    for i, col in enumerate(milestone_cols):
        # gather data per group for this milestone
        group_data = {}
        for j, group in enumerate(groups):
            vals = df.loc[df["Cluster"] == group, col].dropna().values
            if len(vals) >= 2:
                group_data[group] = (j + i * spacing, vals)

        if len(group_data) < 2:
            continue

        # pairwise Mann-Whitney
        pairs = []
        grp_list = list(group_data.keys())
        for a in range(len(grp_list)):
            for b in range(a + 1, len(grp_list)):
                g1, g2 = grp_list[a], grp_list[b]
                _, p = mannwhitneyu(group_data[g1][1], group_data[g2][1],
                                    alternative="two-sided")
                if p < 0.05:
                    pairs.append((g1, g2, p))

        # draw brackets for significant pairs
        # start bracket y above the tallest violin in this milestone
        local_max = max(np.max(v) for _, v in group_data.values())
        y_base = local_max + (y_max_all - y_min_all) * 0.08
        step = (y_max_all - y_min_all) * 0.06

        for k, (g1, g2, p) in enumerate(pairs):
            x1 = group_data[g1][0]
            x2 = group_data[g2][0]
            y_line = y_base + k * step

            ax.plot([x1, x2], [y_line, y_line], "k-", linewidth=1.5)
            ax.text((x1 + x2) / 2, y_line + step * 0.15, f"{p:.2e}",
                    ha="center", va="bottom", fontsize=8, fontweight="bold")

    # --- axes ----------------------------------------------------------------
    # x-ticks at centre of each milestone group
    tick_positions = [(n_groups - 1) / 2 + i * spacing for i in range(n_measures)]
    tick_labels = [labels.get(c, c).replace("\n", " ") for c in milestone_cols]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontweight="bold", fontsize=10)
    ax.set_ylabel("Age (months)", fontsize=12)
    ax.set_title("Psychomotor Milestones by Cluster", fontsize=16, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # expand y to fit brackets
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.25)

    # legend
    legend_handles = [Patch(facecolor=palette[g], alpha=0.7, label=g) for g in groups]
    ax.legend(handles=legend_handles, loc="upper right")

    plt.tight_layout()
    plt.savefig(output_path, dpi=500, format="pdf",
                transparent=True, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("=" * 60)
    print("PSYCHOMOTOR MILESTONES BY CLUSTER")
    print("=" * 60)

    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(TABLES_DIR, exist_ok=True)

    # Load and merge (NT/IDD from broad clinical rosters, C1/C2/C3 from curated k-means)
    roster = build_roster(CLUSTERS_CURATED_FILE, INOVAND_CLINICAL_CURATED, LEAP_CLINICAL_CURATED)
    df = load_and_merge(FONDA_CSV, roster)

    # Clean milestone columns
    for col in MILESTONE_COLS:
        if col in df.columns:
            df[col] = clean_milestone(df[col])

    available_cols = [c for c in MILESTONE_COLS if c in df.columns]
    if not available_cols:
        print("ERROR: No milestone columns found in FONDA data.")
        return

    # Determine which groups are present
    groups = [g for g in ORDER_CLUSTERS if g in df["Cluster"].unique()]
    print(f"Groups present: {groups}")
    for g in groups:
        print(f"  {g}: {(df['Cluster'] == g).sum()} participants")

    # Individual violin plots
    print("\nCreating individual milestone plots...")
    for col in available_cols:
        label = MILESTONE_LABELS.get(col, col)
        out_path = os.path.join(FIGURES_DIR, f"{col.lower()}_by_cluster.pdf")
        plot_single_milestone(df, col, label, groups, PALETTE_CLUSTERS, out_path)

    # Combined 2x2 plot
    print("\nCreating combined milestone plot...")
    plot_combined_milestones(
        df, available_cols, MILESTONE_LABELS, groups, PALETTE_CLUSTERS,
        os.path.join(FIGURES_DIR, "psychomotor_milestones_combined.pdf"),
    )

    # Statistics
    print("\nComputing statistics...")
    all_stats = []
    for col in available_cols:
        stats_df = compute_stats(df, col, groups)
        if len(stats_df) > 0:
            all_stats.append(stats_df)

    if all_stats:
        df_stats = pd.concat(all_stats, ignore_index=True)
        stats_path = os.path.join(TABLES_DIR, "psychomotor_milestones_stats.csv")
        df_stats.to_csv(stats_path, index=False)
        print(f"  Saved: {stats_path}")

        # Print summary
        print("\n--- Statistical Summary ---")
        for col in available_cols:
            sub = df_stats[df_stats["milestone"] == col]
            if len(sub) > 0:
                kw_p = sub["kruskal_wallis_p"].iloc[0]
                print(f"  {MILESTONE_LABELS.get(col, col).replace(chr(10), ' ')}: "
                      f"Kruskal-Wallis p = {kw_p:.4g}")
                sig = sub[sub["pvalue_fdr"] < 0.05]
                for _, row in sig.iterrows():
                    print(f"    {row['group1']} vs {row['group2']}: "
                          f"p_fdr = {row['pvalue_fdr']:.4g}")

    print("\nDone.")


if __name__ == "__main__":
    main()
