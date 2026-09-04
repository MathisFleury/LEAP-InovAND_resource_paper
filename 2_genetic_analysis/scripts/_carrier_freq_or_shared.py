#!/usr/bin/env python3
# =============================================================================
# 3 - Carrier Frequencies and Odds Ratios (Population-Level)
# =============================================================================
# Computes carrier frequencies and odds ratios (vs TD) for rare variants
# across multiple gene lists, ancestries, cohorts, and variant types.
#
# Refactored from: script_zakaria/frequencies_or.py (~3,100 lines → parameterized)
# Input: carrier_annotations.tsv (from script 1), diag_listing.csv
# Output: frequency barplots (PDF), OR forest plots (PDF), OR tables (CSV)
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    FIGURES_DIR, TABLES_DIR, DIAG_FILE, DUP_FILE,
    PALETTE_FREQ, ORDER_FREQ, ORDER_OR,
    PALETTE_POPULATION1, ORDER_POPULATION1,
    LABELS_CONSTRAINED, LABELS_ALL_GENES,
    LABELS_DELLOFMISS_CONSTRAINED, LABELS_DELLOFMISS_ALL,
    COLS_DELLOF_CONSTRAINED, COLS_DELLOF_ALL,
    COLS_DELLOFMISS_CONSTRAINED, COLS_DELLOFMISS_ALL,
    COLS_DUP_CONSTRAINED, LABELS_DUP_CONSTRAINED,
    COLS_DUP_ALL, LABELS_DUP_ALL,
    COMBO_PANEL_KEYS, COMBO_PANEL_LABELS,
    COMBO_SOURCES_DELLOF, COMBO_SOURCES_DELLOFDUP, COMBO_SOURCES_DELLOFDUPMISS,
)

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# ============================================================================
# Reusable analysis functions
# ============================================================================

def plot_carrier_frequencies(df, group_col, columns_to_plot, palette, order,
                             label_map, output_name, title="", ylim=(0, 0.8),
                             figures_dir=None):
    """Barplot of carrier frequencies by group. Saves PDF to figures_dir."""
    if not columns_to_plot:
        print(f"  No columns to plot for {output_name}")
        return
    if figures_dir is None:
        figures_dir = FIGURES_DIR
    df_sub = df[df[group_col].notna()].copy()
    df_sub["Group"] = df_sub[group_col]
    df_melted = pd.melt(
        df_sub,
        id_vars=["Group"],
        value_vars=columns_to_plot,
        var_name="GENELIST",
        value_name="Value",
    )

    g = sns.catplot(
        data=df_melted, x="GENELIST", y="Value", hue="Group",
        kind="bar", height=4, aspect=10 / 4,
        errorbar=("ci", 95), hue_order=order, palette=palette,
    )
    g.set_axis_labels("", "Proportion of carriers")
    g.fig.suptitle(title, y=0.95)

    # Add n= to legend labels
    if g._legend is not None:
        group_counts = df_sub.groupby("Group").size().to_dict()
        for text in g._legend.texts:
            label = text.get_text()
            n = group_counts.get(label, "?")
            text.set_text(f"{label} (n={n})")

    g.set_xticklabels(
        [label_map.get(t.get_text(), t.get_text()) for t in g.ax.get_xticklabels()],
        rotation=0, ha="center",
    )
    g.set(ylim=ylim)
    plt.tight_layout()
    plt.savefig(
        os.path.join(figures_dir, output_name),
        dpi=500, format="pdf", transparent=True, bbox_inches="tight",
    )
    plt.close()
    print(f"  Figure saved: {output_name}")


def make_combo_carrier_cols(df, prefix, sources, panel_keys, panel_labels):
    """Build combined carrier columns by OR-ing the per-panel source columns.
    Returns (cols, labels) suitable for plot_carrier_frequencies / compute_and_plot_odds_ratios.
    Missing source columns are skipped (their absence is treated as no contribution).
    """
    cols = []
    labels = {}
    for key in panel_keys:
        missing = [c for c in sources[key] if c not in df.columns]
        if missing:
            print(f"  Warning: combo {prefix}/{key} missing sources: {missing}")
        srcs = [c for c in sources[key] if c in df.columns]
        out_col = f"combo_{prefix}_{key}"
        if not srcs:
            df[out_col] = False
        else:
            df[out_col] = df[srcs].fillna(False).astype(bool).any(axis=1)
        cols.append(out_col)
        labels[out_col] = panel_labels[key]
    return cols, labels


def compute_and_plot_odds_ratios(df, group_col, columns_to_plot, palette, order,
                                  label_map, output_name, reference_group="NT",
                                  title="", fid_col="FID", save_table=True,
                                  figures_dir=None, tables_dir=None):
    """Compute carrier counts, Fisher exact OR vs reference, FDR correction.
    Saves forest plot (PDF) and optionally OR table (CSV).
    """
    if figures_dir is None:
        figures_dir = FIGURES_DIR
    if tables_dir is None:
        tables_dir = TABLES_DIR
    # --- Deduplicate by family (keep carrier preferentially) ---
    df_filtered = df[df[group_col].notna()].copy()
    # if fid_col in df_filtered.columns:
    #     df_filtered["_has_gene_true"] = df_filtered[columns_to_plot].any(axis=1)
    #     df_filtered = (
    #         df_filtered
    #         .sort_values(by=[group_col, "_has_gene_true"], ascending=[True, False])
    #         .drop_duplicates(subset=[group_col, fid_col], keep="first")
    #         .drop(columns=["_has_gene_true"])
    #     )

    # --- Compute counts ---
    population_names = [p for p in order if p in df_filtered[group_col].unique()]
    # Ensure reference group is always included for OR computation
    if reference_group not in population_names and reference_group in df_filtered[group_col].unique():
        population_names.append(reference_group)
    data_rows = []
    for gene in columns_to_plot:
        for pop in population_names:
            df_pop = df_filtered[df_filtered[group_col] == pop]
            n_carriers = int(df_pop[gene].sum())
            total = len(df_pop)
            data_rows.append({
                "GENELIST": gene, "group": pop,
                "n_carriers": n_carriers,
                "n_non_carriers": total - n_carriers,
                "Frequency": n_carriers / total if total > 0 else 0,
            })
    df_counts = pd.DataFrame(data_rows)

    # --- Fisher exact tests vs reference ---
    odds_ratio_results = []
    for gene in columns_to_plot:
        gene_data = df_counts[df_counts["GENELIST"] == gene]
        ref_row = gene_data[gene_data["group"] == reference_group]
        if ref_row.empty:
            continue
        ref_carriers = ref_row["n_carriers"].iloc[0]
        ref_non = ref_row["n_non_carriers"].iloc[0]

        for _, row in gene_data.iterrows():
            pop = row["group"]
            if pop == reference_group:
                odds_ratio_results.append({
                    "GENELIST": gene, "group": pop,
                    "Odds_Ratio": 1.0, "CI_Lower": 1.0, "CI_Upper": 1.0,
                    "P_Value": np.nan,
                })
                continue

            table = np.array([
                [row["n_carriers"], row["n_non_carriers"]],
                [ref_carriers, ref_non],
            ])
            try:
                stat_table = Table2x2(table)
                or_val = stat_table.oddsratio
                ci_lower, ci_upper = stat_table.oddsratio_confint(alpha=0.05)
            except (ValueError, ZeroDivisionError):
                or_val = ci_lower = ci_upper = np.nan
            try:
                _, p_value = fisher_exact(table)
            except ValueError:
                p_value = np.nan

            odds_ratio_results.append({
                "GENELIST": gene, "group": pop,
                "Odds_Ratio": or_val, "CI_Lower": ci_lower,
                "CI_Upper": ci_upper, "P_Value": p_value,
            })

    df_or = pd.DataFrame(odds_ratio_results)

    # Merge carrier counts into OR table
    if len(df_or) > 0 and len(df_counts) > 0:
        df_or = df_or.merge(
            df_counts[["GENELIST", "group", "n_carriers", "n_non_carriers"]],
            on=["GENELIST", "group"], how="left",
        )

    if len(df_or) == 0:
        print(f"  No OR results computed for {output_name}")
        return df_or

    # FDR correction
    df_or_plot = df_or.dropna(subset=["Odds_Ratio", "CI_Lower", "CI_Upper", "P_Value"])
    df_or_plot = df_or_plot[df_or_plot["group"] != reference_group].copy()

    if len(df_or_plot) > 0:
        valid_pvals = df_or_plot["P_Value"].dropna().values
        if len(valid_pvals) > 0:
            reject, pvals_fdr, _, _ = multipletests(valid_pvals, method="fdr_bh")
            df_or_plot.loc[df_or_plot["P_Value"].notna(), "P_FDR"] = pvals_fdr
            df_or_plot.loc[df_or_plot["P_Value"].notna(), "Reject_FDR"] = reject

    # --- Save OR table ---
    if save_table:
        table_name = output_name.replace(".pdf", ".csv")
        df_or.to_csv(os.path.join(tables_dir, table_name), index=False)

    # --- Forest plot ---
    if len(df_or_plot) == 0:
        print(f"  No OR data to plot for {output_name}")
        return df_or

    group_order_map = {g: i for i, g in enumerate(order)}
    outcome_order_map = {o: i for i, o in enumerate(columns_to_plot)}
    offset = 0.17

    plt.figure(figsize=(10, 3))
    for _, row in df_or_plot.iterrows():
        x_pos = outcome_order_map.get(row["GENELIST"], 0) + \
                (group_order_map.get(row["group"], 0) - 1) * offset
        is_ns = row["CI_Lower"] <= 1 <= row["CI_Upper"]
        plt.errorbar(
            x_pos, row["Odds_Ratio"],
            yerr=[[row["Odds_Ratio"] - row["CI_Lower"]],
                  [row["CI_Upper"] - row["Odds_Ratio"]]],
            fmt="o", color=palette.get(row["group"], "black"),
            capsize=5, linewidth=1.5,
            markerfacecolor="none" if is_ns else palette.get(row["group"], "black"),
            markeredgewidth=2 if is_ns else 1,
        )

    plt.xticks(ticks=np.arange(len(columns_to_plot)), labels=columns_to_plot, rotation=0)
    group_counts = df_filtered.groupby(group_col).size().to_dict()
    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                    markerfacecolor=palette.get(g, "black"), markersize=8,
                    label=f"{g} (n={group_counts.get(g, '?')})")
        for g in order if g != reference_group
    ]
    plt.legend(handles=legend_handles, title="Group")
    plt.axhline(y=1, linestyle="--", color="gray")
    plt.yscale("log")
    plt.ylim(0.01, 100)
    plt.ylabel("Odds Ratio (OR)")
    plt.xlabel("")
    plt.title(title)

    ax = plt.gca()
    ax.set_xticklabels(
        [label_map.get(t.get_text(), t.get_text()) for t in ax.get_xticklabels()],
        rotation=0, ha="center",
    )
    ax.spines["bottom"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(
        os.path.join(figures_dir, output_name),
        dpi=500, format="pdf", transparent=True, bbox_inches="tight",
    )
    plt.close()
    print(f"  Figure saved: {output_name}")
    return df_or


# ============================================================================
# Main
# ============================================================================

def main():
    # --- Load data ---
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

    # Add diagnostic genetic flag
    if os.path.exists(DIAG_FILE):
        df_diag = pd.read_csv(DIAG_FILE, sep=";")
        df_diag = df_diag[df_diag["Individuals_GeneInclusionCriteria"].isna()]
        df["diag_genetic"] = df["barcode"].isin(df_diag["barcode"].tolist())
    else:
        print(f"Warning: {DIAG_FILE} not found, skipping diag_genetic flag")
        df["diag_genetic"] = False

    # Build population frequency column (PopulationS1_freq)
    def assign_pop_freq(row):
        if row["PopulationS1"] == "Autism":
            return "Autism"
        elif row["PopulationS1"] == "NT":
            return "NT"
        elif row["PopulationS1"] == "IDD":
            return "IDD"
        elif row["Population_undiagnosed1"] == "Undiagnosed siblings":
            return "Undiagnosed siblings"
        elif row["Population_undiagnosed1"] == "Undiagnosed parents":
            return "Undiagnosed parents"
        return np.nan

    df["PopulationS1_freq"] = df.apply(assign_pop_freq, axis=1)

    # Build Population1 with undiagnosed breakdown
    def assign_pop1(row):
        if row["PopulationS1_freq"] == "Undiagnosed siblings":
            return "Undiagnosed siblings"
        elif row["PopulationS1_freq"] == "Undiagnosed parents":
            return "Undiagnosed parents"
        elif row["Population1"] == "IDD":
            return "IDD"
        elif row["Population1"] == "NT":
            return "NT"
        elif row["Population1"] == "Autism to exclude":
            return np.nan
        elif row["Population1"] == "Autism with IDD":
            return "Autism with IDD"
        elif row["Population1"] == "Autism without IDD":
            return "Autism without IDD"
        return np.nan

    df["population1_new"] = df.apply(assign_pop1, axis=1)

    order_population1_full = [
        "Autism without IDD", "Autism with IDD", "IDD",
        "Undiagnosed siblings", "Undiagnosed parents", "NT",
    ]

    # --- Label maps for constrained (without diag_genetic prefix) ---
    labels_constrained_no_diag = {
        k: v for k, v in LABELS_CONSTRAINED.items() if k != "diag_genetic"
    }

    # ====================================================================
    # Run all analysis blocks
    # ====================================================================

    # --- 1. PAN ancestry, constrained, Population1 (5 groups) ---
    print("\n=== PAN ancestry, constrained, Population1 ===")
    plot_carrier_frequencies(
        df, "population1_new", COLS_DELLOF_CONSTRAINED,
        PALETTE_FREQ, order_population1_full, LABELS_CONSTRAINED,
        "PAN_constraint_freq_population1.pdf",
    )
    compute_and_plot_odds_ratios(
        df, "population1_new", COLS_DELLOF_CONSTRAINED,
        PALETTE_FREQ, order_population1_full, LABELS_CONSTRAINED,
        "PAN_constraint_or_population1.pdf",
    )

    # --- 2. PAN ancestry, constrained, PopulationS1 (Autism grouped) ---
    print("\n=== PAN ancestry, constrained, PopulationS1 ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", COLS_DELLOF_CONSTRAINED,
        PALETTE_FREQ, ORDER_FREQ, LABELS_CONSTRAINED,
        "PAN_constraint_freq.pdf", ylim=(0, 0.75),
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", COLS_DELLOF_CONSTRAINED,
        PALETTE_FREQ, ORDER_FREQ, LABELS_CONSTRAINED,
        "PAN_constraint_or.pdf",
    )

    # --- 3. PAN ancestry, all genes, PopulationS1 ---
    print("\n=== PAN ancestry, all genes, PopulationS1 ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES,
        "all_ancestries_freq.pdf",
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES,
        "all_ancestries_or.pdf",
        title="LoF + DEL - Multi-ancestry - all constraint",
    )

    # --- 4. EUR ancestry, constrained ---
    print("\n=== EUR ancestry, constrained ===")
    df_eur = df[df["EUR_ancestry"] == True]
    cols_eur_constrained = [c for c in COLS_DELLOF_CONSTRAINED if c != "diag_genetic"]
    plot_carrier_frequencies(
        df_eur, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "EUR_constraint_freq.pdf", ylim=(0, 1),
    )
    compute_and_plot_odds_ratios(
        df_eur, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "EUR_constraint_or.pdf",
    )

    # --- 5. EUR ancestry, all genes ---
    print("\n=== EUR ancestry, all genes ===")
    plot_carrier_frequencies(
        df_eur, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES,
        "EUR_freq.pdf",
    )
    compute_and_plot_odds_ratios(
        df_eur, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES,
        "EUR_or.pdf",
    )

    # --- 6. Non-EUR ancestry, constrained ---
    print("\n=== Non-EUR ancestry, constrained ===")
    df_noneur = df[df["EUR_ancestry"] != True]
    plot_carrier_frequencies(
        df_noneur, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "nonEUR_constraint_freq.pdf",
    )
    compute_and_plot_odds_ratios(
        df_noneur, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "nonEUR_constraint_or.pdf",
    )

    # --- 7. Non-EUR ancestry, all genes ---
    print("\n=== Non-EUR ancestry, all genes ===")
    plot_carrier_frequencies(
        df_noneur, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES,
        "nonEUR_freq.pdf",
    )
    compute_and_plot_odds_ratios(
        df_noneur, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES,
        "nonEUR_or.pdf",
    )

    # --- 8. InovAND cohort, constrained ---
    print("\n=== InovAND cohort, constrained ===")
    df_inovand = df[df["cohort"] == "INOVAND"]
    plot_carrier_frequencies(
        df_inovand, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_freq_inovand.pdf", ylim=(0, 0.7),
    )
    compute_and_plot_odds_ratios(
        df_inovand, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_or_inovand.pdf",
    )

    # --- 9. InovAND cohort, all genes ---
    print("\n=== InovAND cohort, all genes ===")
    plot_carrier_frequencies(
        df_inovand, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES,
        "PAN_freq_inovand.pdf",
    )
    compute_and_plot_odds_ratios(
        df_inovand, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES,
        "PAN_or_inovand.pdf",
    )

    # --- 10. LEAP cohort, constrained ---
    print("\n=== LEAP cohort, constrained ===")
    df_leap = df[df["cohort"] == "LEAP"]
    plot_carrier_frequencies(
        df_leap, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_freq_leap.pdf",
    )
    compute_and_plot_odds_ratios(
        df_leap, "PopulationS1_freq", cols_eur_constrained,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_or_leap.pdf",
    )

    # --- 11. LEAP cohort, all genes ---
    print("\n=== LEAP cohort, all genes ===")
    plot_carrier_frequencies(
        df_leap, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES,
        "PAN_freq_leap.pdf",
    )
    compute_and_plot_odds_ratios(
        df_leap, "PopulationS1_freq", COLS_DELLOF_ALL,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES,
        "PAN_or_leap.pdf",
    )

    # --- 12. DEL+LoF+Miss, constrained ---
    print("\n=== DEL+LoF+Miss, constrained ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", COLS_DELLOFMISS_CONSTRAINED,
        PALETTE_FREQ, ORDER_FREQ, LABELS_DELLOFMISS_CONSTRAINED,
        "PAN_freq_dellofmiss_constraint.pdf", ylim=(0, 1),
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", COLS_DELLOFMISS_CONSTRAINED,
        PALETTE_FREQ, ORDER_FREQ, LABELS_DELLOFMISS_CONSTRAINED,
        "PAN_or_dellofmiss_constraint.pdf",
    )

    # --- 13. DEL+LoF+Miss, all genes ---
    print("\n=== DEL+LoF+Miss, all genes ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", COLS_DELLOFMISS_ALL,
        PALETTE_FREQ, ORDER_FREQ, LABELS_DELLOFMISS_ALL,
        "PAN_freq_dellofmiss_allconstraint.pdf",
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", COLS_DELLOFMISS_ALL,
        PALETTE_FREQ, ORDER_OR, LABELS_DELLOFMISS_ALL,
        "PAN_or_dellofmiss_allconstraint.pdf",
    )

    # --- 14. Duplications, constrained, PopulationS1 ---
    print("\n=== Duplications, constrained, PopulationS1 ===")
    dup_constraint_cols = [c for c in COLS_DUP_CONSTRAINED if c != "dup_any_constraint_carrier"]
    dup_all_cols = [c for c in COLS_DUP_ALL if c != "dup_any_constraint_carrier"]
    dup_real_cols = dup_constraint_cols + dup_all_cols
    df_dup = pd.read_table(DUP_FILE)[["ID"] + dup_real_cols]
    df_with_dup = df.merge(df_dup, on="ID", how="left")
    missing_dup = df_with_dup[dup_real_cols[0]].isna().sum()
    if missing_dup:
        print(f"  Warning: {missing_dup} individuals missing dup annotation; treating as non-carrier")
    for c in dup_real_cols:
        df_with_dup[c] = df_with_dup[c].fillna(False).astype(bool)
    # TODO: drop this synthetic OR once DUP_FILE provides an overall "any
    # constrained gene" carrier flag — read that column directly instead.
    df_with_dup["dup_any_constraint_carrier"] = df_with_dup[dup_constraint_cols].any(axis=1)

    plot_carrier_frequencies(
        df_with_dup, "PopulationS1_freq", COLS_DUP_CONSTRAINED,
        PALETTE_FREQ, ORDER_FREQ, LABELS_DUP_CONSTRAINED,
        "PAN_dup_constraint_freq.pdf",
    )
    compute_and_plot_odds_ratios(
        df_with_dup, "PopulationS1_freq", COLS_DUP_CONSTRAINED,
        PALETTE_FREQ, ORDER_OR, LABELS_DUP_CONSTRAINED,
        "PAN_dup_constraint_or.pdf",
    )

    # --- 14b. Duplications, all genes (no constraint filter), PopulationS1 ---
    print("\n=== Duplications, all genes, PopulationS1 ===")
    plot_carrier_frequencies(
        df_with_dup, "PopulationS1_freq", COLS_DUP_ALL,
        PALETTE_FREQ, ORDER_FREQ, LABELS_DUP_ALL,
        "PAN_dup_allconstraint_freq.pdf",
    )
    compute_and_plot_odds_ratios(
        df_with_dup, "PopulationS1_freq", COLS_DUP_ALL,
        PALETTE_FREQ, ORDER_OR, LABELS_DUP_ALL,
        "PAN_dup_allconstraint_or.pdf",
    )

    # --- 15. Combined-variant panels (DEL+LoF / +DUP / +DUP+MISS) ---
    print("\n=== Combined-variant panels (DEL+LoF | +DUP | +DUP+MISS) ===")
    for prefix, sources, freq_name, or_name in [
        ("dellof",        COMBO_SOURCES_DELLOF,
         "PAN_dellof_constraint_freq.pdf",        "PAN_dellof_constraint_or.pdf"),
        ("dellofdup",     COMBO_SOURCES_DELLOFDUP,
         "PAN_dellofdup_constraint_freq.pdf",     "PAN_dellofdup_constraint_or.pdf"),
        ("dellofdupmiss", COMBO_SOURCES_DELLOFDUPMISS,
         "PAN_dellofdupmiss_constraint_freq.pdf", "PAN_dellofdupmiss_constraint_or.pdf"),
    ]:
        cols, labels = make_combo_carrier_cols(
            df_with_dup, prefix, sources, COMBO_PANEL_KEYS, COMBO_PANEL_LABELS,
        )
        plot_carrier_frequencies(
            df_with_dup, "PopulationS1_freq", cols,
            PALETTE_FREQ, ORDER_FREQ, labels, freq_name,
        )
        compute_and_plot_odds_ratios(
            df_with_dup, "PopulationS1_freq", cols,
            PALETTE_FREQ, ORDER_OR, labels, or_name,
        )

    print("\n=== Done. All figures and tables in:", FIGURES_DIR, "===")


if __name__ == "__main__":
    main()
