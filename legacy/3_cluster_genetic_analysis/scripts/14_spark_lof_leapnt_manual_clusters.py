#!/usr/bin/env python3
# =============================================================================
# 14 - SPARK LoF Carrier Freq & OR with LEAP-InovAND NT (Manual Clusters)
# =============================================================================
# Supplementary analysis: SPARK autism clusters (Manual) tested against
# LEAP-InovAND NT subjects as reference (true neurotypicals, not siblings).
# Uses lof_ carrier columns common to both datasets (hg38).
#
# Input:  SPARK_lof_carrier_annotation.tsv, carrier_annotations_hg38.tsv,
#         df_multi_dataset_with_clusters.csv, individuals_metrics.tsv
# Output: figures_spark_leapnt/ and tables_spark_leapnt/
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "3_cluster_genetic_analysis", "scripts")))
from _config import (
    SPARK_LOF_CARRIER, CARRIER_ANNOTATIONS_HG38,
    CLUSTERS_GMM_FILE, INDIVIDUALS_METRICS, SPARK_PGS_FILE,
    FIGURES_SPARK_LEAPNT_MANUAL_DIR, TABLES_SPARK_LEAPNT_MANUAL_DIR,
    COLS_SPARK_LOF_CONSTRAINED, COLS_SPARK_LOF_ALL,
    LABELS_SPARK_LOF_CONSTRAINED, LABELS_SPARK_LOF_ALL,
    LABELS_PGS, PGS_TRAITS, PGS_FILE,
    PALETTE_CLUSTERS,
)

from importlib.util import spec_from_file_location, module_from_spec
_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", "..", ".."))
_spec = spec_from_file_location(
    "script02",
    os.path.join(_lib_dir, "LEAP-InovAND_resource", "2_genetic_analysis", "scripts", "_carrier_freq_or_shared.py"),
)
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)
plot_carrier_frequencies = _mod.plot_carrier_frequencies
compute_and_plot_odds_ratios = _mod.compute_and_plot_odds_ratios

os.makedirs(FIGURES_SPARK_LEAPNT_MANUAL_DIR, exist_ok=True)
os.makedirs(TABLES_SPARK_LEAPNT_MANUAL_DIR, exist_ok=True)

PALETTE_POP = {
    "Autism": "#8991FA",
    "Autism with IDD": "#324095",
    "Autism without IDD": "#5CAEE1",
    "NT": "#C1C2BC",
}
ORDER_POP = ["Autism with IDD", "Autism without IDD", "NT"]
ORDER_POP_OR = ["Autism with IDD", "Autism without IDD"]
ORDER_POP_SIMPLE = ["Autism", "NT"]
ORDER_POP_SIMPLE_OR = ["Autism"]

ORDER_CLUSTERS = ["NT", "C1", "C2", "C3"]


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
            top_case = df_sub[(df_sub[cluster_col] == cluster_name) & (df_sub[col] >= top_thresh)]["ID"].nunique()
            tot_top = df_sub[df_sub[col] >= top_thresh]["ID"].nunique()
            bot_case = df_sub[(df_sub[cluster_col] == cluster_name) & (df_sub[col] <= bot_thresh)]["ID"].nunique()
            tot_bot = df_sub[df_sub[col] <= bot_thresh]["ID"].nunique()
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
    df_plot = df_ors.copy()
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
    # --- Load SPARK carrier data ---
    print(f"Loading SPARK LoF carriers from: {SPARK_LOF_CARRIER}")
    df_spark = pd.read_table(SPARK_LOF_CARRIER, low_memory=False)
    df_spark["ID"] = df_spark["ID"].astype(str)
    print(f"  SPARK: {len(df_spark)} subjects")

    # --- Load LEAP-InovAND hg38 carrier data (NT only) ---
    print(f"Loading LEAP-InovAND hg38 carriers from: {CARRIER_ANNOTATIONS_HG38}")
    df_li = pd.read_table(CARRIER_ANNOTATIONS_HG38, low_memory=False)
    df_li["ID"] = df_li["ID"].astype(str)
    for col in ["PopulationS1", "Population1", "Population_undiagnosed1"]:
        df_li[col] = df_li[col].str.replace("ID", "IDD").str.replace("TD", "NT")

    # Get NT participants only
    clinic = pd.read_table(INDIVIDUALS_METRICS, low_memory=False)
    df_li = df_li.merge(clinic[["ID", "Relation_to_proposant"]], how="left", on="ID")
    df_nt = df_li[
        (df_li["PopulationS1"] == "NT") &
        (df_li["Relation_to_proposant"] == "participant")
    ].drop_duplicates("ID").copy()
    df_nt["population"] = "NT"
    print(f"  LEAP-InovAND NT: {len(df_nt)} subjects")

    # --- Load cluster assignments ---
    print(f"Loading cluster assignments from: {CLUSTERS_GMM_FILE}")
    clusters = pd.read_csv(CLUSTERS_GMM_FILE, low_memory=False)
    clusters["ID"] = clusters["ID"].astype(str)
    clusters_spark = clusters[clusters["Cohort"] == "SPARK"][["ID", "Cluster_Manual"]].copy()

    # --- Merge SPARK with clusters ---
    df_spark = df_spark.merge(clusters_spark, on="ID", how="left")
    df_spark["Cluster"] = df_spark["Cluster_Manual"]

    # Keep SPARK subjects with population
    df_spark_pop = df_spark[df_spark["population"].isin(["Autism with IDD", "Autism without IDD"])].copy()

    # --- Common lof_ carrier columns ---
    lof_carrier_cols = [c for c in COLS_SPARK_LOF_CONSTRAINED + COLS_SPARK_LOF_ALL
                        if c in df_spark.columns and c in df_nt.columns]
    lof_carrier_cols = list(dict.fromkeys(lof_carrier_cols))  # deduplicate

    # --- Combine: SPARK autism + LEAP-InovAND NT ---
    spark_subset = df_spark_pop[["ID", "population", "Cluster"] + lof_carrier_cols].copy()
    spark_subset["barcode"] = spark_subset["ID"]
    nt_subset = df_nt[["ID", "barcode", "population"] + lof_carrier_cols].copy()
    nt_subset["Cluster"] = np.nan

    df = pd.concat([spark_subset, nt_subset], ignore_index=True)
    print(f"\n  Combined dataset: {len(df)} subjects")
    print(f"    SPARK Autism with IDD: {(df['population'] == 'Autism with IDD').sum()}")
    print(f"    SPARK Autism without IDD: {(df['population'] == 'Autism without IDD').sum()}")
    print(f"    LEAP-InovAND NT: {(df['population'] == 'NT').sum()}")

    # Simple Autism grouping
    df["PopulationS1_simple"] = pd.Series(np.nan, index=df.index, dtype="object")
    df.loc[df["population"].isin(["Autism with IDD", "Autism without IDD"]), "PopulationS1_simple"] = "Autism"
    df.loc[df["population"] == "NT", "PopulationS1_simple"] = "NT"

    # Cluster_with_NT
    def assign_cluster_nt(row):
        if row["population"] == "NT":
            return "NT"
        elif pd.notna(row.get("Cluster")) and row["Cluster"] in ("C1", "C2", "C3"):
            return row["Cluster"]
        return np.nan

    df["Cluster_with_NT"] = df.apply(assign_cluster_nt, axis=1)
    print(f"\nCluster_with_NT distribution:")
    print(df["Cluster_with_NT"].value_counts(dropna=False).to_string())

    cols_constrained = [c for c in COLS_SPARK_LOF_CONSTRAINED if c in df.columns]
    cols_all = [c for c in COLS_SPARK_LOF_ALL if c in df.columns]

    # ====================================================================
    # 1. Population-level (Autism with/without IDD vs LEAP-InovAND NT)
    # ====================================================================
    print("\n=== Population-level frequencies (constrained) ===")
    plot_carrier_frequencies(
        df, "population", cols_constrained,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_leapnt_manual_lof_constraint_freq_population.pdf", ylim=(0, 0.8),
        figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )
    compute_and_plot_odds_ratios(
        df, "population", cols_constrained,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_leapnt_manual_lof_constraint_or_population.pdf",
        reference_group="NT",
        tables_dir=TABLES_SPARK_LEAPNT_MANUAL_DIR, figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )

    print("\n=== Population-level frequencies (all genes) ===")
    plot_carrier_frequencies(
        df, "population", cols_all,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_ALL,
        "SPARK_leapnt_manual_lof_allgenes_freq_population.pdf",
        figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )
    compute_and_plot_odds_ratios(
        df, "population", cols_all,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_ALL,
        "SPARK_leapnt_manual_lof_allgenes_or_population.pdf",
        reference_group="NT",
        tables_dir=TABLES_SPARK_LEAPNT_MANUAL_DIR, figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )

    # ====================================================================
    # 2. Autism grouped vs LEAP-InovAND NT
    # ====================================================================
    print("\n=== Autism vs LEAP-InovAND NT frequencies (constrained) ===")
    plot_carrier_frequencies(
        df, "PopulationS1_simple", cols_constrained,
        PALETTE_POP, ORDER_POP_SIMPLE, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_leapnt_manual_lof_constraint_freq_autism_vs_nt.pdf", ylim=(0, 0.8),
        figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_simple", cols_constrained,
        PALETTE_POP, ORDER_POP_SIMPLE_OR, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_leapnt_manual_lof_constraint_or_autism_vs_nt.pdf",
        reference_group="NT",
        tables_dir=TABLES_SPARK_LEAPNT_MANUAL_DIR, figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )

    print("\n=== Autism vs LEAP-InovAND NT frequencies (all genes) ===")
    plot_carrier_frequencies(
        df, "PopulationS1_simple", cols_all,
        PALETTE_POP, ORDER_POP_SIMPLE, LABELS_SPARK_LOF_ALL,
        "SPARK_leapnt_manual_lof_allgenes_freq_autism_vs_nt.pdf",
        figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_simple", cols_all,
        PALETTE_POP, ORDER_POP_SIMPLE_OR, LABELS_SPARK_LOF_ALL,
        "SPARK_leapnt_manual_lof_allgenes_or_autism_vs_nt.pdf",
        reference_group="NT",
        tables_dir=TABLES_SPARK_LEAPNT_MANUAL_DIR, figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
    )

    # ====================================================================
    # 3. Cluster-level (C1/C2/C3 vs LEAP-InovAND NT)
    # ====================================================================
    df_cluster = df[df["Cluster_with_NT"].notna()].copy()
    n_clustered = len(df_cluster)
    print(f"\n=== Cluster analysis: {n_clustered} subjects ===")

    if n_clustered > 0:
        print("\n=== Cluster-level frequencies (constrained) ===")
        plot_carrier_frequencies(
            df_cluster, "Cluster_with_NT", cols_constrained,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_CONSTRAINED,
            "SPARK_leapnt_manual_lof_cluster_constraint_freq.pdf", ylim=(0, 1),
            figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
        )
        compute_and_plot_odds_ratios(
            df_cluster, "Cluster_with_NT", cols_constrained,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_CONSTRAINED,
            "SPARK_leapnt_manual_lof_cluster_constraint_or.pdf",
            reference_group="NT",
            figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR, tables_dir=TABLES_SPARK_LEAPNT_MANUAL_DIR,
        )

        print("\n=== Cluster-level frequencies (all genes) ===")
        plot_carrier_frequencies(
            df_cluster, "Cluster_with_NT", cols_all,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_ALL,
            "SPARK_leapnt_manual_lof_cluster_allgenes_freq.pdf",
            figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR,
        )
        compute_and_plot_odds_ratios(
            df_cluster, "Cluster_with_NT", cols_all,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_ALL,
            "SPARK_leapnt_manual_lof_cluster_allgenes_or.pdf",
            reference_group="NT",
            figures_dir=FIGURES_SPARK_LEAPNT_MANUAL_DIR, tables_dir=TABLES_SPARK_LEAPNT_MANUAL_DIR,
        )

        # ================================================================
        # PGS analysis — SPARK PGS for clusters, LEAP-InovAND PGS for NT
        # ================================================================
        print("\n=== PGS odds ratios by Manual cluster ===")
        pgs_available = []
        df_pgs = df_cluster.copy()

        # SPARK PGS: merge on ID (spid == ID for SPARK subjects)
        if os.path.exists(SPARK_PGS_FILE):
            pgs_spark = pd.read_table(SPARK_PGS_FILE, low_memory=False)
            pgs_available = [t for t in PGS_TRAITS if t in pgs_spark.columns]
            pgs_spark = pgs_spark[["spid"] + pgs_available].rename(columns={"spid": "ID"})
            pgs_spark["ID"] = pgs_spark["ID"].astype(str)
            df_pgs = df_pgs.merge(pgs_spark, how="left", on="ID", suffixes=("", "_spark"))

        # LEAP-InovAND PGS: merge on barcode (sample_id == barcode for NT)
        if os.path.exists(PGS_FILE) and pgs_available:
            pgs_li = pd.read_table(PGS_FILE)
            pgs_li_cols = ["sample_id"] + [t for t in pgs_available if t in pgs_li.columns]
            pgs_li = pgs_li[pgs_li_cols].rename(columns={"sample_id": "barcode"})
            pgs_li["barcode"] = pgs_li["barcode"].astype(str)
            df_pgs = df_pgs.merge(pgs_li, how="left", on="barcode", suffixes=("", "_li"))
            for t in pgs_available:
                li_col = t + "_li" if t + "_li" in df_pgs.columns else None
                if li_col:
                    df_pgs[t] = df_pgs[t].fillna(df_pgs[li_col])

        if pgs_available:
            print(f"  PGS traits available: {len(pgs_available)}/{len(PGS_TRAITS)}")
            n_with_pgs = df_pgs[pgs_available[0]].notna().sum()
            print(f"  Subjects with PGS data: {n_with_pgs}/{len(df_pgs)}")

            df_ors = compute_pgs_odds_ratios(df_pgs, "Cluster_with_NT", pgs_available, PALETTE_CLUSTERS)
            if len(df_ors) > 0:
                valid_pvals = df_ors["fisher_pvalue"].dropna().values
                if len(valid_pvals) > 0:
                    reject, pvals_fdr, _, _ = multipletests(valid_pvals, method="fdr_bh")
                    df_ors.loc[df_ors["fisher_pvalue"].notna(), "P_FDR"] = pvals_fdr
                    df_ors.loc[df_ors["fisher_pvalue"].notna(), "Reject_FDR"] = reject
                df_ors.to_csv(
                    os.path.join(TABLES_SPARK_LEAPNT_MANUAL_DIR, "SPARK_leapnt_pgs_or_clusters.csv"), index=False
                )
                plot_pgs_forest(df_ors, PALETTE_CLUSTERS, LABELS_PGS,
                                "SPARK_leapnt_cluster_pgs_or.pdf", FIGURES_SPARK_LEAPNT_MANUAL_DIR)
        else:
            print("  PGS files not found, skipping")

    else:
        print("  No clustered subjects — skipping cluster analysis")

    print(f"\n=== Done. Outputs in: {FIGURES_SPARK_LEAPNT_MANUAL_DIR} ===")


if __name__ == "__main__":
    main()
