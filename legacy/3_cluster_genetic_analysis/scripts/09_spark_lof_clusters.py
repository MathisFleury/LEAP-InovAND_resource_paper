#!/usr/bin/env python3
# =============================================================================
# 09 - SPARK LoF Carrier Frequencies & Odds Ratios
# =============================================================================
# Carrier frequency and odds ratio analysis using SPARK LoF variants.
# Analyses by population (Autism with/without IDD, Sibs) and by GMM cluster.
#
# Input:  SPARK_lof_carrier_annotation.tsv
#         df_multi_dataset_with_clusters.csv (for cluster + Population1)
# Output: figures_spark/ and tables_spark/
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
    SPARK_LOF_CARRIER, CLUSTERS_GMM_FILE, SPARK_PGS_FILE,
    FIGURES_SPARK_DIR, TABLES_SPARK_DIR,
    COLS_SPARK_LOF_CONSTRAINED, COLS_SPARK_LOF_ALL,
    LABELS_SPARK_LOF_CONSTRAINED, LABELS_SPARK_LOF_ALL,
    LABELS_PGS, PGS_TRAITS,
)

# Import reusable plot functions from 2_genetic_analysis
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

os.makedirs(FIGURES_SPARK_DIR, exist_ok=True)
os.makedirs(TABLES_SPARK_DIR, exist_ok=True)

# Palettes and orders for SPARK populations (NT = unaffected siblings → "Sibs")
PALETTE_POP = {
    "Autism": "#8991FA",
    "Autism with IDD": "#324095",
    "Autism without IDD": "#5CAEE1",
    "Sibs": "#C1C2BC",
}
ORDER_POP = ["Autism with IDD", "Autism without IDD", "Sibs"]
ORDER_POP_OR = ["Autism with IDD", "Autism without IDD"]

ORDER_POP_SIMPLE = ["Autism", "Sibs"]
ORDER_POP_SIMPLE_OR = ["Autism"]

PALETTE_CLUSTERS = {
    "Sibs": "#C1C2BC",
    "C1": "#7A8B47",
    "C2": "#ff9fa0",
    "C3": "#e7ba52",
}
ORDER_CLUSTERS = ["Sibs", "C1", "C2", "C3"]


def compute_pgs_odds_ratios(df, cluster_col, pgs_traits, palette, reference_group="Sibs"):
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
    df = pd.read_table(SPARK_LOF_CARRIER, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    print(f"  Loaded: {len(df)} subjects")

    # --- Load cluster assignments ---
    print(f"Loading cluster assignments from: {CLUSTERS_GMM_FILE}")
    clusters = pd.read_csv(CLUSTERS_GMM_FILE, low_memory=False)
    clusters["ID"] = clusters["ID"].astype(str)
    clusters_spark = clusters[clusters["Cohort"] == "SPARK"].copy()
    print(f"  SPARK subjects in cluster DF: {len(clusters_spark)}")

    # Merge cluster + Population1
    df = df.merge(
        clusters_spark[["ID", "Cluster", "Population1"]],
        on="ID", how="left",
    )

    # --- Build population columns ---
    # SPARK "NT" are unaffected siblings, rename to "Sibs" for display
    df["population"] = df["population"].replace({"NT": "Sibs"})
    df["PopulationS1_simple"] = pd.Series(np.nan, index=df.index, dtype="object")
    df.loc[df["population"].isin(["Autism with IDD", "Autism without IDD"]), "PopulationS1_simple"] = "Autism"
    df.loc[df["population"] == "Sibs", "PopulationS1_simple"] = "Sibs"

    # Drop blank population (relatives without population label)
    df_pop = df[df["population"].isin(["Autism with IDD", "Autism without IDD", "Sibs"])].copy()
    print(f"  Subjects with known population: {len(df_pop)}")
    for pop in ORDER_POP:
        n = (df_pop["population"] == pop).sum()
        print(f"    {pop}: n={n}")

    # --- Build Cluster_with_Sibs ---
    def assign_cluster_sibs(row):
        if row["population"] == "Sibs":
            return "Sibs"
        elif pd.notna(row.get("Cluster")) and row["Cluster"] in ("C1", "C2", "C3"):
            return row["Cluster"]
        return np.nan

    df_pop["Cluster_with_Sibs"] = df_pop.apply(assign_cluster_sibs, axis=1)
    cluster_counts = df_pop["Cluster_with_Sibs"].value_counts(dropna=False)
    print(f"\nCluster_with_Sibs distribution:")
    print(cluster_counts.to_string())

    # Filter columns that exist
    cols_constrained = [c for c in COLS_SPARK_LOF_CONSTRAINED if c in df_pop.columns]
    cols_all = [c for c in COLS_SPARK_LOF_ALL if c in df_pop.columns]

    # ====================================================================
    # 1. Population-level analysis (Autism with/without IDD vs NT)
    # ====================================================================
    print("\n=== Population-level frequencies (constrained) ===")
    plot_carrier_frequencies(
        df_pop, "population", cols_constrained,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_freq_population.pdf", ylim=(0, 0.8),
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_pop, "population", cols_constrained,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_or_population.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    print("\n=== Population-level frequencies (all genes) ===")
    plot_carrier_frequencies(
        df_pop, "population", cols_all,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_freq_population.pdf",
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_pop, "population", cols_all,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_or_population.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    # ====================================================================
    # 2. Population-level (Autism grouped vs Sibs)
    # ====================================================================
    df_pop_simple = df_pop[df_pop["PopulationS1_simple"].notna()].copy()

    print("\n=== Autism vs Sibs frequencies (constrained) ===")
    plot_carrier_frequencies(
        df_pop_simple, "PopulationS1_simple", cols_constrained,
        PALETTE_POP, ORDER_POP_SIMPLE, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_freq_autism_vs_sibs.pdf", ylim=(0, 0.8),
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_pop_simple, "PopulationS1_simple", cols_constrained,
        PALETTE_POP, ORDER_POP_SIMPLE_OR, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_or_autism_vs_sibs.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    print("\n=== Autism vs Sibs frequencies (all genes) ===")
    plot_carrier_frequencies(
        df_pop_simple, "PopulationS1_simple", cols_all,
        PALETTE_POP, ORDER_POP_SIMPLE, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_freq_autism_vs_sibs.pdf",
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_pop_simple, "PopulationS1_simple", cols_all,
        PALETTE_POP, ORDER_POP_SIMPLE_OR, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_or_autism_vs_sibs.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    # ====================================================================
    # 3. EUR ancestry only
    # ====================================================================
    df_eur = df_pop[df_pop["EUR_ancestry"] == 1].copy()
    print(f"\n=== EUR ancestry: {len(df_eur)} subjects ===")

    print("\n=== EUR ancestry frequencies (constrained) ===")
    plot_carrier_frequencies(
        df_eur, "population", cols_constrained,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_freq_EUR.pdf", ylim=(0, 0.8),
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_eur, "population", cols_constrained,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_or_EUR.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    print("\n=== EUR ancestry frequencies (all genes) ===")
    plot_carrier_frequencies(
        df_eur, "population", cols_all,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_freq_EUR.pdf",
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_eur, "population", cols_all,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_or_EUR.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    # ====================================================================
    # 4. Non-EUR ancestry
    # ====================================================================
    df_noneur = df_pop[df_pop["EUR_ancestry"] != 1].copy()
    print(f"\n=== Non-EUR ancestry: {len(df_noneur)} subjects ===")

    print("\n=== Non-EUR ancestry frequencies (constrained) ===")
    plot_carrier_frequencies(
        df_noneur, "population", cols_constrained,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_freq_nonEUR.pdf",
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_noneur, "population", cols_constrained,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_CONSTRAINED,
        "SPARK_lof_constraint_or_nonEUR.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    print("\n=== Non-EUR ancestry frequencies (all genes) ===")
    plot_carrier_frequencies(
        df_noneur, "population", cols_all,
        PALETTE_POP, ORDER_POP, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_freq_nonEUR.pdf",
        figures_dir=FIGURES_SPARK_DIR,
    )
    compute_and_plot_odds_ratios(
        df_noneur, "population", cols_all,
        PALETTE_POP, ORDER_POP_OR, LABELS_SPARK_LOF_ALL,
        "SPARK_lof_allgenes_or_nonEUR.pdf",
        reference_group="Sibs",
        tables_dir=TABLES_SPARK_DIR, figures_dir=FIGURES_SPARK_DIR,
    )

    # ====================================================================
    # 5. Cluster-level analysis (C1/C2/C3 vs Sibs)
    # ====================================================================
    df_cluster = df_pop[df_pop["Cluster_with_Sibs"].notna()].copy()
    n_clustered = len(df_cluster)
    print(f"\n=== Cluster analysis: {n_clustered} subjects with cluster assignment ===")

    if n_clustered > 0:
        print("\n=== Cluster-level frequencies (constrained) ===")
        plot_carrier_frequencies(
            df_cluster, "Cluster_with_Sibs", cols_constrained,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_CONSTRAINED,
            "SPARK_lof_cluster_constraint_freq.pdf", ylim=(0, 1),
            figures_dir=FIGURES_SPARK_DIR,
        )
        compute_and_plot_odds_ratios(
            df_cluster, "Cluster_with_Sibs", cols_constrained,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_CONSTRAINED,
            "SPARK_lof_cluster_constraint_or.pdf",
            reference_group="Sibs",
            figures_dir=FIGURES_SPARK_DIR, tables_dir=TABLES_SPARK_DIR,
        )

        print("\n=== Cluster-level frequencies (all genes) ===")
        plot_carrier_frequencies(
            df_cluster, "Cluster_with_Sibs", cols_all,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_ALL,
            "SPARK_lof_cluster_allgenes_freq.pdf",
            figures_dir=FIGURES_SPARK_DIR,
        )
        compute_and_plot_odds_ratios(
            df_cluster, "Cluster_with_Sibs", cols_all,
            PALETTE_CLUSTERS, ORDER_CLUSTERS, LABELS_SPARK_LOF_ALL,
            "SPARK_lof_cluster_allgenes_or.pdf",
            reference_group="Sibs",
            figures_dir=FIGURES_SPARK_DIR, tables_dir=TABLES_SPARK_DIR,
        )
        # ================================================================
        # PGS analysis by cluster
        # ================================================================
        print("\n=== PGS odds ratios by cluster ===")
        if os.path.exists(SPARK_PGS_FILE):
            print(f"Loading SPARK PGS from: {SPARK_PGS_FILE}")
            pgs_df = pd.read_table(SPARK_PGS_FILE, low_memory=False)
            pgs_available = [t for t in PGS_TRAITS if t in pgs_df.columns]
            pgs_cols = ["spid"] + pgs_available
            pgs_df = pgs_df[pgs_cols].rename(columns={"spid": "ID"})
            pgs_df["ID"] = pgs_df["ID"].astype(str)
            df_pgs = df_cluster.merge(pgs_df, how="left", on="ID")
            print(f"  PGS traits available: {len(pgs_available)}/{len(PGS_TRAITS)}")

            df_ors = compute_pgs_odds_ratios(df_pgs, "Cluster_with_Sibs", pgs_available, PALETTE_CLUSTERS)

            if len(df_ors) > 0:
                valid_pvals = df_ors["fisher_pvalue"].dropna().values
                if len(valid_pvals) > 0:
                    reject, pvals_fdr, _, _ = multipletests(valid_pvals, method="fdr_bh")
                    df_ors.loc[df_ors["fisher_pvalue"].notna(), "P_FDR"] = pvals_fdr
                    df_ors.loc[df_ors["fisher_pvalue"].notna(), "Reject_FDR"] = reject

                df_ors.to_csv(
                    os.path.join(TABLES_SPARK_DIR, "SPARK_lof_pgs_or_clusters.csv"), index=False
                )
                plot_pgs_forest(df_ors, PALETTE_CLUSTERS, LABELS_PGS,
                                "SPARK_lof_cluster_pgs_or.pdf", FIGURES_SPARK_DIR)
            else:
                print("  No PGS OR results to plot.")
        else:
            print(f"Warning: SPARK PGS file not found ({SPARK_PGS_FILE}), skipping")

    else:
        print("  No clustered subjects — skipping cluster analysis")

    print(f"\n=== Done. Outputs in: {FIGURES_SPARK_DIR} ===")


if __name__ == "__main__":
    main()
