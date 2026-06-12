#!/usr/bin/env python3
# =============================================================================
# 10 - Standardized Beta Coefficients in FDR-Significant Regions
# =============================================================================
# For brain regions that are FDR-significant (p < 0.05) under either
# SynGO-constrained or ChromEpiTF-constrained LOEUF analysis, compares the
# standardized OLS beta coefficients (β from MRI_z ~ β · LOEUF_z, equivalent
# to Pearson r for a single predictor) across three gene lists:
#   - ChromEpiTF constrained   (dellof_chromepitf_constraint_genes_best_score)
#   - SynGO constrained        (dellof_syngo_constraint_genes_best_score)
#   - Protein coding           (dellof_proteinconding_genes_best_score)
#
# Outputs:
#   - figures_genetics/beta_coeff_forest_plot.pdf
#   - tables_genetics/beta_coeff_sign_regions.csv
# =============================================================================

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

from _mri_loader import load_mri

_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", ".."))

CLUSTERS_COMPLETE_FILE = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "df_clusters_complete_kmeans.csv",
)
CORR_STATS_FILE = os.path.join(
    _script_dir, "..", "outputs", "tables_genetics",
    "genetics_mri_correlation_statistics.csv",
)
OUTPUT_FIGURES = os.path.join(_script_dir, "..", "outputs", "figures_genetics")
OUTPUT_TABLES  = os.path.join(_script_dir, "..", "outputs", "tables_genetics")

GENE_LISTS = {
    "ChromEpiTF\n(constrained)": "dellof_chromepitf_constraint_genes_best_score",
    "SynGO\n(constrained)":      "dellof_syngo_constraint_genes_best_score",
    "Protein coding":             "dellof_proteinconding_genes_best_score",
}

PALETTE = {
    "ChromEpiTF\n(constrained)": "#8991FA",
    "SynGO\n(constrained)":      "#9AD5D3",
    "Protein coding":             "#F5A623",
}


def load_data():
    print(f"Loading LOEUF data: {CLUSTERS_COMPLETE_FILE}")
    df_g = pd.read_csv(CLUSTERS_COMPLETE_FILE, low_memory=False)
    df_g["ID"] = df_g["ID"].astype(str)
    drop_mri = [c for c in df_g.columns if c.startswith("lh_") or c.startswith("rh_")]
    df_g = df_g.drop(columns=drop_mri)

    print(f"Loading MRI data...")
    df_m = load_mri()

    cortical_cols = [c for c in df_m.columns
                     if (c.endswith("_thickness") or c.endswith("_area"))
                     and (c.startswith("lh_") or c.startswith("rh_"))
                     and "MeanThickness" not in c and "WhiteSurfArea" not in c]

    df = df_g.merge(df_m[["ID"] + cortical_cols], on="ID", how="inner")
    print(f"Merged: {len(df)} individuals, {len(cortical_cols)} cortical regions")
    return df, cortical_cols


def get_significant_mri_cols():
    """Return MRI columns FDR-significant (p<0.05) under SynGO or ChromEpiTF constrained."""
    print(f"\nLoading correlation statistics: {CORR_STATS_FILE}")
    df_s = pd.read_csv(CORR_STATS_FILE)

    sig_features = [
        "dellof_syngo_constraint_genes_best_score",
        "dellof_chromepitf_constraint_genes_best_score",
    ]
    sig = df_s[
        df_s["genetic_feature"].isin(sig_features) &
        (df_s["p_fdr"] < 0.05) &
        (df_s["mri_type"].isin(["thickness", "area"]))
    ]
    # Reconstruct mri_col from hemisphere + region + mri_type
    sig = sig.copy()
    sig["mri_col"] = (
        sig["hemisphere"].str.replace("left", "lh").str.replace("right", "rh")
        + "_" + sig["region"]
        + "_" + sig["mri_type"]
    )
    mri_cols = sorted(sig["mri_col"].unique())
    print(f"  FDR-significant regions (SynGO or CHROM constrained): {len(mri_cols)}")
    return mri_cols


def compute_standardized_beta(df, loeuf_col, mri_col):
    """
    Compute standardized OLS β from MRI_z ~ β · LOEUF_z.
    For a single predictor, β_std = Pearson r.
    Returns (beta, p_value).
    """
    raw = df[loeuf_col].dropna()
    raw = raw[raw > 0]
    x = -np.log10(raw)

    mri = df[mri_col].dropna()
    common = x.index.intersection(mri.index)
    if len(common) < 10:
        return np.nan, np.nan

    x_vals = x.loc[common].values
    y_vals = mri.loc[common].values

    # Standardize
    x_z = scipy_stats.zscore(x_vals)
    y_z = scipy_stats.zscore(y_vals)

    try:
        r, p = pearsonr(x_z, y_z)   # β_std = r for single predictor
        return r, p
    except Exception:
        return np.nan, np.nan


def plot_forest(df_beta):
    """
    Forest / point plot comparing standardized beta across gene lists
    for FDR-significant regions.

    One panel per gene list (3 panels in a row).
    Each panel:
      - Y axis: brain region (sorted by mean beta across gene lists)
      - X axis: standardized beta
      - Point: observed beta, colored by sign (positive=blue, negative=red)
      - Vertical dashed line at 0
    """
    # Sort regions by mean beta across gene lists (descending)
    all_mri_cols = sorted(df_beta["mri_col"].unique())
    mean_beta = {}
    for mc in all_mri_cols:
        vals = df_beta[df_beta["mri_col"] == mc]["beta"].values
        mean_beta[mc] = np.nanmean(vals)
    sorted_cols = sorted(all_mri_cols, key=lambda c: mean_beta.get(c, 0), reverse=True)
    rl_map = df_beta.set_index("mri_col")["region_label"].to_dict()
    region_labels = [rl_map.get(c, c) for c in sorted_cols]
    n_regions = len(sorted_cols)

    gene_list_order = list(GENE_LISTS.keys())
    ncols = 3
    figsize = (ncols * 4, max(4, n_regions * 0.4 + 1.5))
    fig, axes = plt.subplots(1, ncols, figsize=figsize, sharey=True)

    for col_idx, gl in enumerate(gene_list_order):
        ax = axes[col_idx]
        ax.axvline(0, color="grey", lw=0.8, ls="--", zorder=0)

        gl_rows = df_beta[df_beta["gene_list"] == gl].set_index("mri_col")

        for i, mc in enumerate(sorted_cols):
            if mc not in gl_rows.index:
                continue
            row = gl_rows.loc[mc]
            beta  = row["beta"]
            p_fdr = row.get("p_fdr", 1.0)

            color = "#4C72B0" if beta > 0 else "#C44E52"
            alpha = 1.0 if p_fdr < 0.05 else 0.35

            ax.scatter([beta], [i], color=color, s=25,
                       alpha=alpha, zorder=3)

        ax.set_yticks(range(n_regions))
        if col_idx == 0:
            ax.set_yticklabels(region_labels, fontsize=6)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("Standardized \u03b2", fontsize=8)
        ax.set_title(gl.replace("\n", " "), fontsize=9, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#4C72B0",
                   markersize=5, alpha=1.0, label="Significant (p\u22060.05)"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
                   markersize=5, alpha=0.35, label="Not significant"),
        ]
        ax.legend(handles=legend_elements, fontsize=6, loc="lower right",
                  frameon=False)

    fig.suptitle(
        "Standardized \u03b2 in FDR-Significant Regions\n"
        "(ChromEpiTF constrained or SynGO constrained, FDR p < 0.05)",
        fontsize=10, fontweight="bold",
    )
    plt.tight_layout()
    path = os.path.join(OUTPUT_FIGURES, "beta_coeff_forest_plot.pdf")
    plt.savefig(path, dpi=300, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"  Saved: beta_coeff_forest_plot.pdf")


def region_label(mri_col):
    m = re.match(r"^(lh|rh)_(.+)_(thickness|area)$", mri_col)
    if m:
        side = "L" if m.group(1) == "lh" else "R"
        return f"{side} {m.group(2).replace('_', ' ')} ({m.group(3)[:5]})"
    return mri_col


def main():
    df, cortical_cols = load_data()

    sig_mri_cols = get_significant_mri_cols()
    # Restrict to columns that exist in dataset
    sig_mri_cols = [c for c in sig_mri_cols if c in df.columns]
    print(f"  Available in dataset: {len(sig_mri_cols)}")

    rows_beta = []

    print(f"\nComputing standardized \u03b2 "
          f"({len(GENE_LISTS)} gene lists \u00d7 {len(sig_mri_cols)} regions)...")

    for gl_label, loeuf_col in GENE_LISTS.items():
        if loeuf_col not in df.columns:
            print(f"  Column missing: {loeuf_col} \u2014 skipping")
            continue
        for mri_col in sig_mri_cols:
            beta, p_val = compute_standardized_beta(df, loeuf_col, mri_col)
            if np.isnan(beta):
                continue

            rows_beta.append({
                "gene_list":    gl_label,
                "loeuf_col":    loeuf_col,
                "mri_col":      mri_col,
                "region_label": region_label(mri_col),
                "beta":         beta,
                "p_value":      p_val,
            })

        print(f"  {gl_label.replace(chr(10), ' ')}: {sum(1 for r in rows_beta if r['gene_list'] == gl_label)} regions")

    df_beta = pd.DataFrame(rows_beta)
    if df_beta.empty:
        print("No results \u2014 check input files.")
        return

    # FDR per gene list
    for gl in df_beta["gene_list"].unique():
        idx = df_beta["gene_list"] == gl
        valid = idx & df_beta["p_value"].notna()
        if valid.sum() > 0:
            _, pfdr, _, _ = multipletests(df_beta.loc[valid, "p_value"].values, method="fdr_bh")
            df_beta.loc[valid, "p_fdr"] = pfdr

    # Save CSV
    os.makedirs(OUTPUT_TABLES, exist_ok=True)
    os.makedirs(OUTPUT_FIGURES, exist_ok=True)
    df_beta.to_csv(os.path.join(OUTPUT_TABLES, "beta_coeff_sign_regions.csv"), index=False)
    print(f"\n  Table saved: beta_coeff_sign_regions.csv")

    # --- Forest plot ---
    plot_forest(df_beta)

    print(f"\n=== Done. Outputs in: {OUTPUT_FIGURES} ===")


if __name__ == "__main__":
    main()
