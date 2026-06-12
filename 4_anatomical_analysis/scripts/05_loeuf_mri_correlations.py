#!/usr/bin/env python3
# =============================================================================
# 05 - LOEUF vs MRI Correlation Statistics (hg19 legacy)
# =============================================================================
# Computes Pearson correlations between -log10(LOEUF best score) and ALL
# cortical thickness / surface area / subcortical volume regions.
# Applies FDR correction.
#
# Outputs:
#   - tables_genetics/genetics_mri_correlation_statistics.csv  (feeds R ggseg script)
#
# Input:  df_clusters_complete_kmeans.csv   (LOEUF best-score columns + Autism subset)
#         MRI_ANAT_INOVAND_LEAP_COMBAT.tsv
# =============================================================================

import os
import sys
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", ".."))

CLUSTERS_COMPLETE_FILE = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "df_clusters_complete_kmeans.csv",
)
MRI_ANAT_FILE = os.path.join(
    _lib_dir, "imaging2genet", "0_input", "dataframes",
    "MRI_ANAT_INOVAND_LEAP_COMBAT.tsv",
)
OUTPUT_TABLES = os.path.join(_script_dir, "..", "outputs", "tables_genetics")
OUTPUT_FIGURES = os.path.join(_script_dir, "..", "outputs", "figures_genetics")
os.makedirs(OUTPUT_TABLES, exist_ok=True)
os.makedirs(OUTPUT_FIGURES, exist_ok=True)

# LOEUF features to analyse (dellof + miss, constrained and all, SynGO/CHROM focus)
LOEUF_FEATURES = {
    "dellof_syngo_constraint_genes_best_score":      "dellof_constrained",
    "dellof_chromepitf_constraint_genes_best_score": "dellof_constrained",
    "miss_syngo_constraint_genes_best_score":        "miss_constrained",
    "miss_chromepitf_constraint_genes_best_score":   "miss_constrained",
    "dellof_syngo_genes_best_score":                 "dellof_all",
    "dellof_chromepitf_genes_best_score":            "dellof_all",
    "miss_syngo_genes_best_score":                   "miss_all",
    "miss_chromepitf_genes_best_score":              "miss_all",
    # Additional gene lists
    "dellof_hcndddomv6_constraint_genes_best_score": "dellof_constrained",
    "dellof_sparksfari1_constraint_genes_best_score":"dellof_constrained",
    "dellof_constraint_genes_best_score":            "dellof_constrained",
    "dellof_proteinconding_genes_best_score":        "dellof_all",
}

# Subcortical structures of interest
SUBCORTICAL_COLS = [
    "Left-Hippocampus", "Right-Hippocampus",
    "Left-Amygdala",    "Right-Amygdala",
    "Left-Caudate",     "Right-Caudate",
    "Left-Putamen",     "Right-Putamen",
    "Left-Pallidum",    "Right-Pallidum",
    "Left-Thalamus",    "Right-Thalamus",
    "Left-Accumbens-area", "Right-Accumbens-area",
]


def parse_mri_col(col):
    """Parse hemisphere, region, mri_type from a FreeSurfer column name."""
    # Cortical: lh_superiortemporal_thickness / rh_bankssts_area
    m = re.match(r"^(lh|rh)_(.+)_(thickness|area)$", col)
    if m:
        hemi = "left" if m.group(1) == "lh" else "right"
        return hemi, m.group(2), m.group(3)
    # Subcortical: Left-Hippocampus / Right-Amygdala
    m2 = re.match(r"^(Left|Right)-(.+)$", col)
    if m2:
        hemi = "left" if m2.group(1) == "Left" else "right"
        region = m2.group(2).lower().replace("-", "_")
        return hemi, region, "volume"
    return None, None, None


def load_data():
    print(f"Loading LOEUF data: {CLUSTERS_COMPLETE_FILE}")
    df_g = pd.read_csv(CLUSTERS_COMPLETE_FILE, low_memory=False)
    df_g["ID"] = df_g["ID"].astype(str)
    # Drop any MRI columns already present (use COMBAT-harmonised file instead)
    drop_mri = [c for c in df_g.columns if c.startswith("lh_") or c.startswith("rh_")]
    df_g = df_g.drop(columns=drop_mri)

    print(f"Loading MRI data: {MRI_ANAT_FILE}")
    df_m = pd.read_csv(MRI_ANAT_FILE, sep="\t", low_memory=False)
    df_m["ID"] = df_m["ID"].astype(str)

    # Select all thickness, area, and subcortical columns that exist
    cortical_cols = [c for c in df_m.columns
                     if (c.endswith("_thickness") or c.endswith("_area"))
                     and (c.startswith("lh_") or c.startswith("rh_"))
                     and "MeanThickness" not in c and "WhiteSurfArea" not in c]
    subco_cols = [c for c in SUBCORTICAL_COLS if c in df_m.columns]
    mri_cols = cortical_cols + subco_cols

    df = df_g.merge(df_m[["ID"] + mri_cols], on="ID", how="inner")
    print(f"Merged: {len(df)} individuals, {len(mri_cols)} MRI regions")
    return df, mri_cols


def compute_correlations(df, mri_cols):
    rows = []

    available = {col: vtype for col, vtype in LOEUF_FEATURES.items() if col in df.columns}
    print(f"LOEUF features available: {len(available)}")

    for loeuf_col, variant_type in available.items():
        raw = df[loeuf_col].dropna()
        raw = raw[raw > 0]
        neg_log = -np.log10(raw)

        for mri_col in mri_cols:
            if mri_col not in df.columns:
                continue
            mri_vals = df[mri_col].dropna()
            common = neg_log.index.intersection(mri_vals.index)
            if len(common) < 10:
                continue

            x = neg_log.loc[common].values
            y = mri_vals.loc[common].values

            try:
                r, p = pearsonr(x, y)
            except Exception:
                continue

            hemi, region, mri_type = parse_mri_col(mri_col)
            if hemi is None:
                continue

            rows.append({
                "genetic_type":    "LOEUF",
                "genetic_feature": loeuf_col,
                "variant_type":    variant_type,
                "mri_col":         mri_col,
                "mri_type":        mri_type,
                "hemisphere":      hemi,
                "region":          region,
                "correlation":     r,
                "p_value":         p,
                "sample_size":     len(common),
            })

        print(f"  {loeuf_col}: {len([r for r in rows if r['genetic_feature'] == loeuf_col])} regions computed")

    df_corr = pd.DataFrame(rows)
    if len(df_corr) > 0:
        # FDR per variant_type × mri_type combination
        for (vtype, mtype), grp in df_corr.groupby(["variant_type", "mri_type"]):
            valid = grp["p_value"].notna()
            if valid.sum() == 0:
                continue
            _, p_fdr, _, _ = multipletests(grp.loc[valid, "p_value"].values, method="fdr_bh")
            df_corr.loc[grp[valid].index, "p_fdr"] = p_fdr

    return df_corr


def main():
    df, mri_cols = load_data()

    print(f"\nComputing correlations...")
    df_corr = compute_correlations(df, mri_cols)

    print(f"\nTotal correlations: {len(df_corr)}")
    if "p_fdr" in df_corr.columns:
        sig_fdr = (df_corr["p_fdr"] < 0.05).sum()
        print(f"Significant (FDR < 0.05): {sig_fdr}")

    # Save main statistics file (feeds the R ggseg script)
    stats_path = os.path.join(OUTPUT_TABLES, "genetics_mri_correlation_statistics.csv")
    df_corr.drop(columns=["mri_col"], errors="ignore").to_csv(stats_path, index=False)
    print(f"\n  Table saved: genetics_mri_correlation_statistics.csv")

    print(f"\n=== Done. Outputs in: {OUTPUT_TABLES} ===")


if __name__ == "__main__":
    main()
