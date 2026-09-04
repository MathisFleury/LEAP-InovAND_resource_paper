#!/usr/bin/env python3
# =============================================================================
# 07 - LOEUF score vs Cortical Thickness + Subcortical Volumes (hg19)
# =============================================================================
# Regresses cortical thickness, surface area, and subcortical volumes against
# -log10(LOEUF best score) per gene list using Pearson r (= standardized β
# for a single z-scored predictor).
# FDR correction (Benjamini-Hochberg) per subset.
#
# Input:  carrier_annotations.tsv (2_genetic_analysis/outputs/tables)
#         MRI_ANAT_INOVAND_LEAP_COMBAT.tsv
# Output: figures_genetics/, tables_genetics/
# =============================================================================

import os
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", "..", ".."))

CARRIER_HG19 = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "legacy", "2_genetic_analysis", "outputs", "tables",
    "carrier_annotations.tsv",
)
MRI_ANAT_FILE = os.path.join(
    _lib_dir, "imaging2genet", "0_input", "dataframes",
    "MRI_ANAT_INOVAND_LEAP_COMBAT.tsv",
)
OUTPUT_FIGURES = os.path.join(_script_dir, "..", "outputs", "figures_genetics")
OUTPUT_TABLES  = os.path.join(_script_dir, "..", "outputs", "tables_genetics")
os.makedirs(OUTPUT_FIGURES, exist_ok=True)
os.makedirs(OUTPUT_TABLES,  exist_ok=True)

LOEUF_COLS = {
    "constrained": [
        "dellof_constraint_genes_best_score",
        "dellof_hcndddomv6_xlinked_boyz_constraint_genes_best_score",
        "dellof_sparksfari1_constraint_genes_best_score",
        "dellof_eagles_ds_constraint_genes_best_score",
        "dellof_syngo_constraint_genes_best_score",
        "dellof_chromepitf_constraint_genes_best_score",
        "dellof_syngo_chromepitf_constraint_genes_best_score",
    ],
    "all": [
        "dellof_proteinconding_genes_best_score",
        "dellof_hcndddomv6_xlinked_boyz_genes_best_score",
        "dellof_sparksfari1_genes_best_score",
        "dellof_eagles_ds_genes_best_score",
        "dellof_syngo_genes_best_score",
        "dellof_chromepitf_genes_best_score",
        "dellof_syngo_chromepitf_genes_best_score",
    ],
}

LOEUF_LABELS = {
    "dellof_constraint_genes_best_score":                          "All (constrained)",
    "dellof_hcndddomv6_xlinked_boyz_constraint_genes_best_score": "HCNDD DOM+XL\n(constrained)",
    "dellof_sparksfari1_constraint_genes_best_score":              "SPARK SFARI\n(constrained)",
    "dellof_eagles_ds_constraint_genes_best_score":                "EAGLE\n(constrained)",
    "dellof_syngo_constraint_genes_best_score":                    "SynGO\n(constrained)",
    "dellof_chromepitf_constraint_genes_best_score":               "ChromEpiTF\n(constrained)",
    "dellof_syngo_chromepitf_constraint_genes_best_score":         "SynGO+ChromEpiTF\n(constrained)",
    "dellof_proteinconding_genes_best_score":                      "Protein coding",
    "dellof_hcndddomv6_xlinked_boyz_genes_best_score":            "HCNDD DOM+XL",
    "dellof_sparksfari1_genes_best_score":                         "SPARK SFARI",
    "dellof_eagles_ds_genes_best_score":                           "EAGLE",
    "dellof_syngo_genes_best_score":                               "SynGO",
    "dellof_chromepitf_genes_best_score":                          "ChromEpiTF",
    "dellof_syngo_chromepitf_genes_best_score":                    "SynGO+ChromEpiTF",
}

SUBCORTICAL_COLS = [
    "Left-Hippocampus", "Right-Hippocampus",
    "Left-Amygdala",    "Right-Amygdala",
    "Left-Caudate",     "Right-Caudate",
    "Left-Putamen",     "Right-Putamen",
    "Left-Pallidum",    "Right-Pallidum",
    "Left-Thalamus",    "Right-Thalamus",
    "Left-Accumbens-area", "Right-Accumbens-area",
]

SUBCORTICAL_LABELS = {
    c: c.replace("Left-", "L ").replace("Right-", "R ").replace("-area", " area")
    for c in SUBCORTICAL_COLS
}


def load_data():
    print(f"Loading carrier annotations: {CARRIER_HG19}")
    df = pd.read_table(CARRIER_HG19, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    for col in ["PopulationS1", "Population1", "Population_undiagnosed1"]:
        if col in df.columns:
            df[col] = df[col].str.replace("ID", "IDD").str.replace("TD", "NT")
    df = df[df["PopulationS1"] != "other"].drop_duplicates("ID")

    print(f"Loading MRI data: {MRI_ANAT_FILE}")
    df_m = pd.read_csv(MRI_ANAT_FILE, sep="\t", low_memory=False)
    df_m["ID"] = df_m["ID"].astype(str)
    cortical_cols = [c for c in df_m.columns
                     if (c.endswith("_thickness") or c.endswith("_area"))
                     and (c.startswith("lh_") or c.startswith("rh_"))
                     and "MeanThickness" not in c and "WhiteSurfArea" not in c]
    subco_cols = [c for c in SUBCORTICAL_COLS if c in df_m.columns]

    df = df.merge(df_m[["ID"] + cortical_cols + subco_cols], on="ID", how="inner")
    print(f"Merged: {len(df)} individuals")
    return df, cortical_cols, subco_cols


def get_mri_label(col):
    m = re.match(r"^(lh|rh)_(.+)_(thickness|area)$", col)
    if m:
        side = "L" if m.group(1) == "lh" else "R"
        return f"{side} {m.group(2).replace('_',' ')} ({m.group(3)[:5]})"
    return SUBCORTICAL_LABELS.get(col, col)


def compute_stats(df, loeuf_cols, mri_cols):
    rows = []
    for lc in loeuf_cols:
        if lc not in df.columns:
            print(f"    Column missing: {lc} — skipping")
            continue
        raw = df[lc].dropna()
        raw = raw[raw > 0]
        x = -np.log10(raw)

        for mc in mri_cols:
            if mc not in df.columns:
                continue
            mri = df[mc].dropna()
            common = x.index.intersection(mri.index)
            if len(common) < 10:
                continue

            x_z = scipy_stats.zscore(x.loc[common].values)
            y_z = scipy_stats.zscore(mri.loc[common].values)

            try:
                r, p = pearsonr(x_z, y_z)
            except Exception:
                r, p = np.nan, np.nan
            if np.isnan(r):
                continue

            rows.append({
                "loeuf_col":         lc,
                "loeuf_label":       LOEUF_LABELS.get(lc, lc),
                "mri_col":           mc,
                "mri_label":         get_mri_label(mc),
                "n":                 len(common),
                "standardized_beta": r,
                "p_value":           p,
            })

    df_s = pd.DataFrame(rows)
    if len(df_s) > 0:
        valid = df_s["p_value"].notna()
        _, p_fdr, _, _ = multipletests(df_s.loc[valid, "p_value"].values, method="fdr_bh")
        df_s.loc[valid, "p_fdr"] = p_fdr
    return df_s


def plot_heatmap(df_s, mri_cols_subset, output_name, title):
    df_plot = df_s[df_s["mri_col"].isin(mri_cols_subset)]
    if df_plot.empty:
        return
    loeuf_labels = list(dict.fromkeys(df_plot["loeuf_label"].tolist()))
    mri_labels   = list(dict.fromkeys(df_plot["mri_label"].tolist()))

    beta_mat = np.full((len(loeuf_labels), len(mri_labels)), np.nan)
    pfdr_mat = np.full((len(loeuf_labels), len(mri_labels)), np.nan)

    for _, row in df_plot.iterrows():
        i = loeuf_labels.index(row["loeuf_label"])
        j = mri_labels.index(row["mri_label"])
        beta_mat[i, j] = row["standardized_beta"]
        pfdr_mat[i, j] = row.get("p_fdr", np.nan)

    nrows = len(loeuf_labels)
    ncols = len(mri_labels)
    fig, ax = plt.subplots(figsize=(max(6, ncols * 0.5 + 2), nrows * 0.7 + 2))
    sns.heatmap(
        beta_mat, mask=np.isnan(beta_mat),
        annot=False, cmap="RdBu_r", center=0,
        linewidths=0.3,
        cbar_kws={"shrink": 0.6, "label": "Standardized β (-log10 LOEUF ~ MRI)"},
        xticklabels=mri_labels, yticklabels=loeuf_labels, ax=ax,
    )
    for i in range(nrows):
        for j in range(ncols):
            if np.isnan(pfdr_mat[i, j]):
                continue
            p = pfdr_mat[i, j]
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            if star:
                ax.text(j + 0.5, i + 0.7, star, ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")

    ax.set_title(title, fontsize=11, fontweight="bold", pad=12)
    plt.xticks(rotation=60, ha="right", fontsize=6)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()
    path = os.path.join(OUTPUT_FIGURES, output_name)
    plt.savefig(path, dpi=500, format="pdf", transparent=True, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {output_name}")


def main():
    df, cortical_cols, subco_cols = load_data()

    for subset_name, cols in LOEUF_COLS.items():
        cols_exist = [c for c in cols if c in df.columns]
        print(f"\n=== hg19 {subset_name} ({len(cols_exist)} gene lists) ===")
        df_s = compute_stats(df, cols_exist, cortical_cols + subco_cols)
        if df_s.empty:
            print("  No results — skipping")
            continue

        thick_cols = [c for c in cortical_cols if c.endswith("_thickness")]
        plot_heatmap(df_s, thick_cols,
                     f"loeuf_thickness_heatmap_{subset_name}.pdf",
                     f"LOEUF vs Cortical Thickness (hg19, {subset_name})")

        area_cols = [c for c in cortical_cols if c.endswith("_area")]
        plot_heatmap(df_s, area_cols,
                     f"loeuf_area_heatmap_{subset_name}.pdf",
                     f"LOEUF vs Surface Area (hg19, {subset_name})")

        plot_heatmap(df_s, subco_cols,
                     f"loeuf_subcortical_heatmap_{subset_name}.pdf",
                     f"LOEUF vs Subcortical Volume (hg19, {subset_name})")

        df_s.to_csv(
            os.path.join(OUTPUT_TABLES, f"loeuf_mri_{subset_name}.csv"), index=False
        )
        print(f"  Table saved: loeuf_mri_{subset_name}.csv")

    print(f"\n=== Done. Outputs in: {OUTPUT_FIGURES} ===")


if __name__ == "__main__":
    main()
