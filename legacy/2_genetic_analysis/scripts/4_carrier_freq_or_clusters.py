#!/usr/bin/env python3
# =============================================================================
# Carrier Frequencies, Odds Ratios & PGS by Clinical Cluster, hg19/gnomAD v2
# (legacy -- superseded by gnomAD v4)
# =============================================================================
# Stratifies carrier frequency and OR analyses by clinical clusters (C1-C3)
# plus TD and IDD controls. Also computes PGS odds ratios (top 25% vs
# bottom 25%) by cluster.
#
# Refactored from: script_zakaria/frequencies_or_clusters.py
# Input: carrier_annotations.tsv (1_carrier_annotation.py), cluster_assignments.csv (1_clustering),
#        PGS file, diag_listing.csv, individuals_metrics.tsv
# Output: cluster-level frequency/OR figures + PGS OR forest plot + tables
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import multipletests

_ORIGINAL_SCRIPTS_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
    "2_genetic_analysis", "scripts"))
sys.path.insert(0, _ORIGINAL_SCRIPTS_DIR)
from _config import (
    FIGURES_DIR, TABLES_DIR, DIAG_FILE, CLUSTERS_FILE, PGS_FILE,
    PALETTE_CLUSTERS, ORDER_CLUSTERS,
    LABELS_CONSTRAINED, LABELS_PGS,
    COLS_DELLOF_CONSTRAINED, PGS_TRAITS,
)

# Import reusable functions from the shared engine (stays in the main tree --
# also imported by the current v4 pipeline, see 2_genetic_analysis/README.md)
from importlib.util import spec_from_file_location, module_from_spec
_script_dir = os.path.dirname(os.path.abspath(__file__))
_spec = spec_from_file_location("script02", os.path.join(_ORIGINAL_SCRIPTS_DIR, "_carrier_freq_or_shared.py"))
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)
plot_carrier_frequencies = _mod.plot_carrier_frequencies
compute_and_plot_odds_ratios = _mod.compute_and_plot_odds_ratios

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)


def compute_pgs_odds_ratios(df, cluster_col, pgs_traits, palette, reference_group="NT"):
    """Compute PGS odds ratios (top 25% vs bottom 25%) for each cluster vs NT.

    Returns a DataFrame with OR, CI, and p-values for each cluster × trait.
    """
    results = []
    clusters = [c for c in df[cluster_col].dropna().unique() if c != reference_group]

    for cluster_name in clusters:
        df_sub = df[df[cluster_col].isin([reference_group, cluster_name])]

        for col in pgs_traits:
            if col not in df_sub.columns:
                continue
            top_thresh = df_sub[col].quantile(0.75)
            bot_thresh = df_sub[col].quantile(0.25)

            top_case = df_sub[
                (df_sub[cluster_col] == cluster_name) & (df_sub[col] >= top_thresh)
            ]["barcode"].nunique()
            tot_top = df_sub[df_sub[col] >= top_thresh]["barcode"].nunique()

            bot_case = df_sub[
                (df_sub[cluster_col] == cluster_name) & (df_sub[col] <= bot_thresh)
            ]["barcode"].nunique()
            tot_bot = df_sub[df_sub[col] <= bot_thresh]["barcode"].nunique()

            table = [
                [top_case, tot_top - top_case],
                [bot_case, tot_bot - bot_case],
            ]

            try:
                stat_table = Table2x2(table)
                OR = stat_table.oddsratio
                ci_low, ci_high = stat_table.oddsratio_confint()
                _, pvalue = fisher_exact(table)
            except Exception:
                OR, ci_low, ci_high, pvalue = np.nan, np.nan, np.nan, np.nan

            results.append({
                "cluster": cluster_name,
                "trait": col,
                "OR": OR,
                "CI_lower": ci_low,
                "CI_upper": ci_high,
                "fisher_pvalue": pvalue,
                "top_25_thresh": top_thresh,
                "bot_25_thresh": bot_thresh,
                "top_case": top_case,
                "tot_top": tot_top,
                "bot_case": bot_case,
                "tot_bot": tot_bot,
            })

    return pd.DataFrame(results)


def plot_pgs_forest(df_ors, palette, label_map, output_name):
    """Forest plot of PGS odds ratios by cluster."""
    cluster_order = ["C3", "C2", "C1"]
    df_plot = df_ors[~df_ors["cluster"].isin(["IDD"])].copy()
    df_plot["cluster_order"] = df_plot["cluster"].apply(
        lambda x: cluster_order.index(x) if x in cluster_order else len(cluster_order)
    )
    df_plot = df_plot.sort_values(
        by=["cluster_order", "OR"], ascending=[True, True]
    ).reset_index(drop=True)

    plt.figure(figsize=(12, 6), dpi=200)
    for idx, row in df_plot.iterrows():
        signif = row["fisher_pvalue"] < 0.05
        mfc = palette.get(row["cluster"], "black") if signif else "white"
        plt.errorbar(
            x=row["OR"], y=idx,
            xerr=[[row["OR"] - row["CI_lower"]], [row["CI_upper"] - row["OR"]]],
            fmt="o", mfc=mfc,
            mec=palette.get(row["cluster"], "black"),
            ecolor=palette.get(row["cluster"], "black"),
            markersize=5,
        )

    plt.axvline(1, color="grey", linestyle="--")
    y_labels = df_plot["trait"].map(label_map).fillna(df_plot["trait"])
    plt.yticks(range(len(df_plot)), y_labels, fontsize=10)
    plt.xlabel("Odds Ratio (top 25% vs bottom 25%)", fontsize=8)
    plt.xlim(0, 8)
    plt.tight_layout()
    plt.savefig(
        os.path.join(FIGURES_DIR, output_name),
        dpi=500, format="pdf", transparent=True, bbox_inches="tight",
    )
    plt.close()
    print(f"  Figure saved: {output_name}")


def main():
    # --- Load carrier annotations ---
    carrier_path = os.path.join(TABLES_DIR, "carrier_annotations.tsv")
    print(f"Loading carrier annotations from: {carrier_path}")
    df = pd.read_table(carrier_path)

    # Fix population labels
    df["PopulationS1"] = df["PopulationS1"].str.replace("ID", "IDD")
    df["Population1"] = df["Population1"].str.replace("ID", "IDD")
    df["Population_undiagnosed1"] = df["Population_undiagnosed1"].str.replace("ID", "IDD")
    # Rename TD → NT
    for col in ["PopulationS1", "Population1", "Population_undiagnosed1"]:
        df[col] = df[col].str.replace("TD", "NT")
    df = df[df["PopulationS1"] != "other"]
    df = df.drop_duplicates("ID")

    # --- Load cluster assignments + clinical data ---
    # Drop stale Cluster column from pre-processed carrier file (if present)
    if "Cluster" in df.columns:
        df = df.drop(columns=["Cluster"])

    print(f"Loading clusters from: {CLUSTERS_FILE}")
    clusters = pd.read_csv(CLUSTERS_FILE, low_memory=False)
    merge_cols = ["ID", "Cluster"]
    if "Relation_to_proposant" in clusters.columns:
        merge_cols.append("Relation_to_proposant")
    df = df.merge(clusters[merge_cols], how="left", on="ID")

    # If Relation_to_proposant wasn't in cluster file, load from clinical data
    if "Relation_to_proposant" not in df.columns:
        from _config import INDIVIDUALS_METRICS
        clinic = pd.read_table(INDIVIDUALS_METRICS)
        df = df.merge(clinic[["ID", "Relation_to_proposant"]], how="left", on="ID")

    # --- Build Cluster_with_NT column ---
    # C1/C2/C3 for Autism, plus NT and IDD participants
    def assign_cluster_nt(row):
        if row["PopulationS1"] == "NT" and row.get("Relation_to_proposant") == "participant":
            return "NT"
        elif row["PopulationS1"] == "IDD" and row.get("Relation_to_proposant") == "participant":
            return "IDD"
        elif pd.isna(row.get("Cluster")):
            return np.nan
        elif row["Cluster"] == "C1" and row["PopulationS1"] == "Autism":
            return "C1"
        elif row["Cluster"] == "C2" and row["PopulationS1"] == "Autism":
            return "C2"
        elif row["Cluster"] == "C3" and row["PopulationS1"] == "Autism":
            return "C3"
        return np.nan

    df["Cluster_with_NT"] = df.apply(assign_cluster_nt, axis=1)

    # --- Add diagnostic genetic flag ---
    if os.path.exists(DIAG_FILE):
        df_diag = pd.read_csv(DIAG_FILE, sep=";")
        df_diag = df_diag[df_diag["Individuals_GeneInclusionCriteria"].isna()]
        df["diag_genetic"] = df["barcode"].isin(df_diag["barcode"].tolist())
    else:
        df["diag_genetic"] = False

    # ====================================================================
    # Carrier frequency + OR by cluster
    # ====================================================================
    print("\n=== Cluster-level carrier frequencies (constrained) ===")
    plot_carrier_frequencies(
        df, "Cluster_with_NT", COLS_DELLOF_CONSTRAINED,
        PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_CONSTRAINED,
        "CLUSTER_TD_constraint_freq.pdf", ylim=(0, 1),
    )

    print("\n=== Cluster-level odds ratios (constrained) ===")
    compute_and_plot_odds_ratios(
        df, "Cluster_with_NT", COLS_DELLOF_CONSTRAINED,
        PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_CONSTRAINED,
        "CLUSTER_TD_constraint_or.pdf",
    )

    # ====================================================================
    # PGS analysis
    # ====================================================================
    print("\n=== PGS odds ratios by cluster ===")
    if os.path.exists(PGS_FILE):
        print(f"Loading PGS from: {PGS_FILE}")
        pgs_df = pd.read_table(PGS_FILE)
        pgs_cols = ["sample_id"] + [t for t in PGS_TRAITS if t in pgs_df.columns]
        pgs_df = pgs_df[pgs_cols].rename(columns={"sample_id": "barcode"})
        df = df.merge(pgs_df, how="left", on="barcode")

        df_ors = compute_pgs_odds_ratios(
            df, "Cluster_with_NT", PGS_TRAITS, PALETTE_CLUSTERS,
        )

        if len(df_ors) > 0:
            # FDR correction
            valid_pvals = df_ors["fisher_pvalue"].dropna().values
            if len(valid_pvals) > 0:
                reject, pvals_fdr, _, _ = multipletests(valid_pvals, method="fdr_bh")
                df_ors.loc[df_ors["fisher_pvalue"].notna(), "P_FDR"] = pvals_fdr
                df_ors.loc[df_ors["fisher_pvalue"].notna(), "Reject_FDR"] = reject

            # Save table
            df_ors.to_csv(
                os.path.join(TABLES_DIR, "pgs_odds_ratios_clusters.csv"), index=False
            )

            # Plot
            plot_pgs_forest(
                df_ors, PALETTE_CLUSTERS, LABELS_PGS,
                "CLUSTER_pgs_OR_top25_bottom25.pdf",
            )
        else:
            print("  No PGS OR results to plot.")
    else:
        print(f"Warning: PGS file not found ({PGS_FILE}), skipping PGS analysis")

    print("\n=== Done. Outputs in:", FIGURES_DIR, "===")


if __name__ == "__main__":
    main()
