#!/usr/bin/env python3
# =============================================================================
# 04 - Carrier Frequencies & Odds Ratios (hg38 / GRCh38)
# =============================================================================
# Same analyses as script 02 but using the GRCh38-based carrier annotation
# file, which includes additional gene lists: metabolic and SynGO+ChromEpiTF.
#
# Input: hg38 carrier file (external drive), diag_listing.csv
# Output: figures_hg38/ and tables_hg38/ in outputs/
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    FIGURES_DIR_HG38, TABLES_DIR_HG38, CARRIER_HG38, DIAG_FILE, TABLES_DIR, DUP_FILE,
    PALETTE_FREQ, ORDER_FREQ, ORDER_OR,
    PALETTE_POPULATION1, ORDER_POPULATION1,
    LABELS_CONSTRAINED_EXT, LABELS_ALL_GENES_EXT,
    COLS_DELLOF_CONSTRAINED_EXT, COLS_DELLOF_ALL_EXT,
    COLS_MISS_CONSTRAINED, COLS_MISS_ALL,
    LABELS_MISS_CONSTRAINED, LABELS_MISS_ALL,
    COLS_DUP_CONSTRAINED, LABELS_DUP_CONSTRAINED,
    COMBO_PANEL_KEYS, COMBO_PANEL_LABELS,
    COMBO_SOURCES_DELLOF, COMBO_SOURCES_DELLOFDUP, COMBO_SOURCES_DELLOFDUPMISS_HG38,
)

# Import reusable functions from script 02
from importlib.util import spec_from_file_location, module_from_spec
_script_dir = os.path.dirname(os.path.abspath(__file__))
_spec = spec_from_file_location("script02", os.path.join(_script_dir, "02_carrier_freq_or.py"))
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)
plot_carrier_frequencies = _mod.plot_carrier_frequencies
compute_and_plot_odds_ratios = _mod.compute_and_plot_odds_ratios
make_combo_carrier_cols = _mod.make_combo_carrier_cols

os.makedirs(FIGURES_DIR_HG38, exist_ok=True)
os.makedirs(TABLES_DIR_HG38, exist_ok=True)


def main():
    # --- Load hg38 carrier data ---
    # Prefer generated file from script 01b; fall back to pre-annotated file
    generated_path = os.path.join(TABLES_DIR_HG38, "carrier_annotations_hg38.tsv")
    if os.path.exists(generated_path):
        print(f"Loading generated hg38 carrier file: {generated_path}")
        df = pd.read_table(generated_path, low_memory=False)
    elif os.path.exists(CARRIER_HG38):
        print(f"Loading pre-annotated hg38 carrier file: {CARRIER_HG38}")
        df = pd.read_table(CARRIER_HG38, low_memory=False)
    else:
        print(f"Error: No hg38 carrier file found.")
        print(f"  Tried: {generated_path}")
        print(f"  Tried: {CARRIER_HG38}")
        return

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
    if "diag_genetic" not in df.columns:
        if os.path.exists(DIAG_FILE):
            df_diag = pd.read_csv(DIAG_FILE, sep=";")
            df_diag = df_diag[df_diag["Individuals_GeneInclusionCriteria"].isna()]
            df["diag_genetic"] = df["barcode"].isin(df_diag["barcode"].tolist())
        else:
            print(f"Warning: {DIAG_FILE} not found, skipping diag_genetic flag")
            df["diag_genetic"] = False

    # Build population frequency column
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

    # Filter columns that actually exist in the hg38 file
    def filter_existing(cols):
        return [c for c in cols if c in df.columns]

    cols_constrained = filter_existing(COLS_DELLOF_CONSTRAINED_EXT)
    cols_all = filter_existing(COLS_DELLOF_ALL_EXT)
    cols_miss_constrained = filter_existing(COLS_MISS_CONSTRAINED)
    cols_miss_all = filter_existing(COLS_MISS_ALL)

    # Constrained labels without diag_genetic (for subgroup analyses)
    labels_constrained_no_diag = {
        k: v for k, v in LABELS_CONSTRAINED_EXT.items() if k != "diag_genetic"
    }
    cols_constrained_no_diag = [c for c in cols_constrained if c != "diag_genetic"]

    # ====================================================================
    # Run all analysis blocks (figures → figures_hg38/)
    # ====================================================================

    # --- 1. PAN ancestry, constrained, Population1 ---
    print("\n=== PAN ancestry, constrained, Population1 ===")
    plot_carrier_frequencies(
        df, "population1_new", cols_constrained,
        PALETTE_FREQ, order_population1_full, LABELS_CONSTRAINED_EXT,
        "PAN_constraint_freq_population1.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df, "population1_new", cols_constrained,
        PALETTE_FREQ, order_population1_full, LABELS_CONSTRAINED_EXT,
        "PAN_constraint_or_population1.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 2. PAN ancestry, constrained, PopulationS1 ---
    print("\n=== PAN ancestry, constrained, PopulationS1 ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", cols_constrained,
        PALETTE_FREQ, ORDER_FREQ, LABELS_CONSTRAINED_EXT,
        "PAN_constraint_freq.pdf", ylim=(0, 0.75),
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", cols_constrained,
        PALETTE_FREQ, ORDER_FREQ, LABELS_CONSTRAINED_EXT,
        "PAN_constraint_or.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 3. PAN ancestry, all genes, PopulationS1 ---
    print("\n=== PAN ancestry, all genes, PopulationS1 ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES_EXT,
        "all_ancestries_freq.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES_EXT,
        "all_ancestries_or.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 4. EUR ancestry, constrained ---
    print("\n=== EUR ancestry, constrained ===")
    df_eur = df[df["EUR_ancestry"] == True]
    plot_carrier_frequencies(
        df_eur, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "EUR_constraint_freq.pdf", ylim=(0, 1),
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_eur, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "EUR_constraint_or.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 5. EUR ancestry, all genes ---
    print("\n=== EUR ancestry, all genes ===")
    plot_carrier_frequencies(
        df_eur, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES_EXT,
        "EUR_freq.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_eur, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES_EXT,
        "EUR_or.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 6. Non-EUR ancestry, constrained ---
    print("\n=== Non-EUR ancestry, constrained ===")
    df_noneur = df[df["EUR_ancestry"] != True]
    plot_carrier_frequencies(
        df_noneur, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "nonEUR_constraint_freq.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_noneur, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "nonEUR_constraint_or.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 7. Non-EUR ancestry, all genes ---
    print("\n=== Non-EUR ancestry, all genes ===")
    plot_carrier_frequencies(
        df_noneur, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES_EXT,
        "nonEUR_freq.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_noneur, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES_EXT,
        "nonEUR_or.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 8. InovAND cohort, constrained ---
    print("\n=== InovAND cohort, constrained ===")
    df_inovand = df[df["cohort"] == "INOVAND"]
    plot_carrier_frequencies(
        df_inovand, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_freq_inovand.pdf", ylim=(0, 0.7),
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_inovand, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_or_inovand.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 9. InovAND cohort, all genes ---
    print("\n=== InovAND cohort, all genes ===")
    plot_carrier_frequencies(
        df_inovand, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES_EXT,
        "PAN_freq_inovand.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_inovand, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES_EXT,
        "PAN_or_inovand.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 10. LEAP cohort, constrained ---
    print("\n=== LEAP cohort, constrained ===")
    df_leap = df[df["cohort"] == "LEAP"]
    plot_carrier_frequencies(
        df_leap, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_freq_leap.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_leap, "PopulationS1_freq", cols_constrained_no_diag,
        PALETTE_FREQ, ORDER_FREQ, labels_constrained_no_diag,
        "PAN_constraint_or_leap.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 11. LEAP cohort, all genes ---
    print("\n=== LEAP cohort, all genes ===")
    plot_carrier_frequencies(
        df_leap, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_FREQ, LABELS_ALL_GENES_EXT,
        "PAN_freq_leap.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_leap, "PopulationS1_freq", cols_all,
        PALETTE_FREQ, ORDER_OR, LABELS_ALL_GENES_EXT,
        "PAN_or_leap.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 12. DEL+LoF+Miss, constrained ---
    print("\n=== DEL+LoF+Miss, constrained ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", cols_miss_constrained,
        PALETTE_FREQ, ORDER_FREQ, LABELS_MISS_CONSTRAINED,
        "PAN_freq_dellofmiss_constraint.pdf", ylim=(0, 1),
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", cols_miss_constrained,
        PALETTE_FREQ, ORDER_FREQ, LABELS_MISS_CONSTRAINED,
        "PAN_or_dellofmiss_constraint.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 13. DEL+LoF+Miss, all genes ---
    print("\n=== DEL+LoF+Miss, all genes ===")
    plot_carrier_frequencies(
        df, "PopulationS1_freq", cols_miss_all,
        PALETTE_FREQ, ORDER_FREQ, LABELS_MISS_ALL,
        "PAN_freq_dellofmiss_allconstraint.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df, "PopulationS1_freq", cols_miss_all,
        PALETTE_FREQ, ORDER_OR, LABELS_MISS_ALL,
        "PAN_or_dellofmiss_allconstraint.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 14. Duplications, constrained, PopulationS1 ---
    print("\n=== Duplications, constrained, PopulationS1 ===")
    dup_real_cols = [c for c in COLS_DUP_CONSTRAINED if c != "dup_any_constraint_carrier"]
    df_dup = pd.read_table(DUP_FILE)[["ID"] + dup_real_cols]
    df_with_dup = df.merge(df_dup, on="ID", how="left")
    missing_dup = df_with_dup[dup_real_cols[0]].isna().sum()
    if missing_dup:
        print(f"  Warning: {missing_dup} individuals missing dup annotation; treating as non-carrier")
    for c in dup_real_cols:
        df_with_dup[c] = df_with_dup[c].fillna(False).astype(bool)
    # TODO: drop this synthetic OR once DUP_FILE provides an overall "any
    # constrained gene" carrier flag — read that column directly instead.
    df_with_dup["dup_any_constraint_carrier"] = df_with_dup[dup_real_cols].any(axis=1)

    plot_carrier_frequencies(
        df_with_dup, "PopulationS1_freq", COLS_DUP_CONSTRAINED,
        PALETTE_FREQ, ORDER_FREQ, LABELS_DUP_CONSTRAINED,
        "PAN_dup_constraint_freq.pdf",
        figures_dir=FIGURES_DIR_HG38,
    )
    compute_and_plot_odds_ratios(
        df_with_dup, "PopulationS1_freq", COLS_DUP_CONSTRAINED,
        PALETTE_FREQ, ORDER_OR, LABELS_DUP_CONSTRAINED,
        "PAN_dup_constraint_or.pdf",
        tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
    )

    # --- 15. Combined-variant panels (DEL+LoF / +DUP / +DUP+MISS) ---
    print("\n=== Combined-variant panels (DEL+LoF | +DUP | +DUP+MISS) ===")
    for prefix, sources, freq_name, or_name in [
        ("dellof",        COMBO_SOURCES_DELLOF,
         "PAN_dellof_constraint_freq.pdf",        "PAN_dellof_constraint_or.pdf"),
        ("dellofdup",     COMBO_SOURCES_DELLOFDUP,
         "PAN_dellofdup_constraint_freq.pdf",     "PAN_dellofdup_constraint_or.pdf"),
        ("dellofdupmiss", COMBO_SOURCES_DELLOFDUPMISS_HG38,
         "PAN_dellofdupmiss_constraint_freq.pdf", "PAN_dellofdupmiss_constraint_or.pdf"),
    ]:
        cols, labels = make_combo_carrier_cols(
            df_with_dup, prefix, sources, COMBO_PANEL_KEYS, COMBO_PANEL_LABELS,
        )
        plot_carrier_frequencies(
            df_with_dup, "PopulationS1_freq", cols,
            PALETTE_FREQ, ORDER_FREQ, labels, freq_name,
            figures_dir=FIGURES_DIR_HG38,
        )
        compute_and_plot_odds_ratios(
            df_with_dup, "PopulationS1_freq", cols,
            PALETTE_FREQ, ORDER_OR, labels, or_name,
            tables_dir=TABLES_DIR_HG38, figures_dir=FIGURES_DIR_HG38,
        )

    print("\n=== Done. All hg38 figures and tables in:", FIGURES_DIR_HG38, "===")


if __name__ == "__main__":
    main()
