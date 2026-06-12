#!/usr/bin/env python3
# =============================================================================
# 11 - Carrier Frequencies & Odds Ratios by Manual Cluster (hg19)
# =============================================================================
# Same as 02_gmm_clusters.py but uses Cluster_Manual instead of Cluster_GMM.
#
# Input:  carrier_annotations.tsv (2_genetic_analysis outputs)
#         df_multi_dataset_with_clusters.csv (Cluster_Manual column)
#         diag_listing.csv, individuals_metrics.tsv, PGS file
# Output: figures_manual/ and tables_manual/
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    CARRIER_ANNOTATIONS, CLUSTERS_GMM_FILE,
    INDIVIDUALS_METRICS, DIAG_FILE, PGS_FILE,
    FIGURES_MANUAL_DIR, TABLES_MANUAL_DIR,
    PALETTE_CLUSTERS, ORDER_CLUSTERS,
    COLS_DELLOF_CONSTRAINED, COLS_DELLOF_ALL,
    LABELS_CONSTRAINED, LABELS_ALL_GENES, LABELS_PGS,
    PGS_TRAITS,
)

from importlib.util import spec_from_file_location, module_from_spec
_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", ".."))
_spec = spec_from_file_location(
    "script02",
    os.path.join(_lib_dir, "LEAP-InovAND_resource", "2_genetic_analysis", "scripts", "02_carrier_freq_or.py"),
)
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)
plot_carrier_frequencies = _mod.plot_carrier_frequencies
compute_and_plot_odds_ratios = _mod.compute_and_plot_odds_ratios

os.makedirs(FIGURES_MANUAL_DIR, exist_ok=True)
os.makedirs(TABLES_MANUAL_DIR, exist_ok=True)


def compute_pgs_odds_ratios(df, cluster_col, pgs_traits, palette, reference_group="NT"):
    results = []
    clusters = [c for c in df[cluster_col].dropna().unique() if c != reference_group]
    for cluster_name in clusters:
        df_sub = df[df[cluster_col].isin([reference_group, cluster_name])]
        for col in pgs_traits:
            if col not in df_sub.columns:
                continue
            top_thresh = df_sub[col].quantile(0.75)
            bot_thresh = df_sub[col].quantile(0.25)
            top_case = df_sub[(df_sub[cluster_col] == cluster_name) & (df_sub[col] >= top_thresh)]["barcode"].nunique()
            tot_top = df_sub[df_sub[col] >= top_thresh]["barcode"].nunique()
            bot_case = df_sub[(df_sub[cluster_col] == cluster_name) & (df_sub[col] <= bot_thresh)]["barcode"].nunique()
            tot_bot = df_sub[df_sub[col] <= bot_thresh]["barcode"].nunique()
            table = [[top_case, tot_top - top_case], [bot_case, tot_bot - bot_case]]
            try:
                stat_table = Table2x2(table)
                OR = stat_table.oddsratio
                ci_low, ci_high = stat_table.oddsratio_confint()
                _, pvalue = fisher_exact(table)
            except Exception:
                OR, ci_low, ci_high, pvalue = np.nan, np.nan, np.nan, np.nan
            results.append({
                "cluster": cluster_name, "trait": col,
                "OR": OR, "CI_lower": ci_low, "CI_upper": ci_high,
                "fisher_pvalue": pvalue,
                "top_case": top_case, "tot_top": tot_top,
                "bot_case": bot_case, "tot_bot": tot_bot,
            })
    return pd.DataFrame(results)


def plot_pgs_forest(df_ors, palette, label_map, output_name, figures_dir):
    cluster_order = ["C3", "C2", "C1"]
    df_plot = df_ors[~df_ors["cluster"].isin(["IDD"])].copy()
    df_plot["cluster_order"] = df_plot["cluster"].apply(
        lambda x: cluster_order.index(x) if x in cluster_order else len(cluster_order)
    )
    df_plot = df_plot.sort_values(by=["cluster_order", "OR"], ascending=[True, True]).reset_index(drop=True)
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
    plt.savefig(os.path.join(figures_dir, output_name), dpi=500, format="pdf", transparent=True, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {output_name}")


def main():
    # --- Load carrier annotations ---
    print(f"Loading carrier annotations from: {CARRIER_ANNOTATIONS}")
    df = pd.read_table(CARRIER_ANNOTATIONS, low_memory=False)

    for col in ["PopulationS1", "Population1", "Population_undiagnosed1"]:
        df[col] = df[col].str.replace("ID", "IDD").str.replace("TD", "NT")
    df = df[df["PopulationS1"] != "other"].drop_duplicates("ID")

    # --- Load Manual cluster assignments ---
    print(f"Loading Manual clusters from: {CLUSTERS_GMM_FILE}")
    clusters = pd.read_csv(CLUSTERS_GMM_FILE, low_memory=False)
    df = df.merge(clusters[["ID", "Cluster_Manual"]], how="left", on="ID")

    clinic = pd.read_table(INDIVIDUALS_METRICS, low_memory=False)
    df = df.merge(clinic[["ID", "Relation_to_proposant"]], how="left", on="ID")

    # --- Build Cluster_with_NT (using Cluster_Manual) ---
    def assign_cluster_nt(row):
        if row["PopulationS1"] == "NT" and row.get("Relation_to_proposant") == "participant":
            return "NT"
        elif row["PopulationS1"] == "IDD" and row.get("Relation_to_proposant") == "participant":
            return "IDD"
        elif pd.isna(row.get("Cluster_Manual")):
            return np.nan
        elif row["Cluster_Manual"] == "C1" and row["PopulationS1"] == "Autism":
            return "C1"
        elif row["Cluster_Manual"] == "C2" and row["PopulationS1"] == "Autism":
            return "C2"
        elif row["Cluster_Manual"] == "C3" and row["PopulationS1"] == "Autism":
            return "C3"
        return np.nan

    df["Cluster_with_NT"] = df.apply(assign_cluster_nt, axis=1)
    print("Cluster_with_NT (Manual) distribution:")
    print(df["Cluster_with_NT"].value_counts(dropna=False).to_string())

    # --- Add diagnostic genetic flag ---
    if os.path.exists(DIAG_FILE):
        df_diag = pd.read_csv(DIAG_FILE, sep=";")
        df_diag = df_diag[df_diag["Individuals_GeneInclusionCriteria"].isna()]
        df["diag_genetic"] = df["barcode"].isin(df_diag["barcode"].tolist())
    else:
        df["diag_genetic"] = False

    # ====================================================================
    # Carrier frequency + OR by cluster (constrained)
    # ====================================================================
    print("\n=== Cluster-level carrier frequencies (constrained) ===")
    plot_carrier_frequencies(
        df, "Cluster_with_NT", COLS_DELLOF_CONSTRAINED,
        PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_CONSTRAINED,
        "manual_cluster_constraint_freq.pdf", ylim=(0, 1),
        figures_dir=FIGURES_MANUAL_DIR,
    )

    print("\n=== Cluster-level odds ratios (constrained) ===")
    compute_and_plot_odds_ratios(
        df, "Cluster_with_NT", COLS_DELLOF_CONSTRAINED,
        PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_CONSTRAINED,
        "manual_cluster_constraint_or.pdf",
        figures_dir=FIGURES_MANUAL_DIR, tables_dir=TABLES_MANUAL_DIR,
    )

    # ====================================================================
    # Carrier frequency + OR by cluster (all genes)
    # ====================================================================
    print("\n=== Cluster-level carrier frequencies (all genes) ===")
    cols_no_diag = [c for c in COLS_DELLOF_ALL if c != "diag_genetic"]
    labels_no_diag = {k: v for k, v in LABELS_ALL_GENES.items() if k != "diag_genetic"}
    plot_carrier_frequencies(
        df, "Cluster_with_NT", cols_no_diag,
        PALETTE_CLUSTERS, ORDER_CLUSTERS, labels_no_diag,
        "manual_cluster_allgenes_freq.pdf",
        figures_dir=FIGURES_MANUAL_DIR,
    )

    print("\n=== Cluster-level odds ratios (all genes) ===")
    compute_and_plot_odds_ratios(
        df, "Cluster_with_NT", cols_no_diag,
        PALETTE_CLUSTERS, ORDER_CLUSTERS, labels_no_diag,
        "manual_cluster_allgenes_or.pdf",
        figures_dir=FIGURES_MANUAL_DIR, tables_dir=TABLES_MANUAL_DIR,
    )

    # ====================================================================
    # PGS analysis
    # ====================================================================
    print("\n=== PGS odds ratios by Manual cluster ===")
    if os.path.exists(PGS_FILE):
        print(f"Loading PGS from: {PGS_FILE}")
        pgs_df = pd.read_table(PGS_FILE)
        pgs_cols = ["sample_id"] + [t for t in PGS_TRAITS if t in pgs_df.columns]
        pgs_df = pgs_df[pgs_cols].rename(columns={"sample_id": "barcode"})
        df_pgs = df.merge(pgs_df, how="left", on="barcode")

        df_ors = compute_pgs_odds_ratios(df_pgs, "Cluster_with_NT", PGS_TRAITS, PALETTE_CLUSTERS)

        if len(df_ors) > 0:
            valid_pvals = df_ors["fisher_pvalue"].dropna().values
            if len(valid_pvals) > 0:
                reject, pvals_fdr, _, _ = multipletests(valid_pvals, method="fdr_bh")
                df_ors.loc[df_ors["fisher_pvalue"].notna(), "P_FDR"] = pvals_fdr
                df_ors.loc[df_ors["fisher_pvalue"].notna(), "Reject_FDR"] = reject

            df_ors.to_csv(
                os.path.join(TABLES_MANUAL_DIR, "manual_pgs_or_clusters.csv"), index=False
            )
            plot_pgs_forest(df_ors, PALETTE_CLUSTERS, LABELS_PGS,
                            "manual_cluster_pgs_or.pdf", FIGURES_MANUAL_DIR)
        else:
            print("  No PGS OR results to plot.")
    else:
        print(f"Warning: PGS file not found ({PGS_FILE}), skipping")

    print(f"\n=== Done. Outputs in: {FIGURES_MANUAL_DIR} ===")


if __name__ == "__main__":
    main()
