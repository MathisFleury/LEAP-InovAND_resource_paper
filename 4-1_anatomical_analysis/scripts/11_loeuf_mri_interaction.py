#!/usr/bin/env python3
# =============================================================================
# 11 - LOEUF × MRI Interaction: Scatter, Forest, IQ-SRS plots (hg19 + hg38)
# =============================================================================
# Produces scatter plots (LOEUF vs cortical thickness in highlighted regions),
# forest plots (beta across all CT regions), and IQ vs SRS-2 scatter with
# LOEUF sizing and carrier outlines.
#
# Gene lists: dellof_ + miss_ (constrained and all), SynGO + ChromEpiTF.
# Correlation statistics for brain maps are read from script 05's CSV (hg19)
# or computed here (hg38, where miss_ is not applicable).
#
# Inputs (hg19):  df_clusters_complete_kmeans.csv, freesurfer_zscore_demo_only.tsv (via load_mri())
#                 tables_genetics/genetics_mri_correlation_statistics.csv
# Inputs (hg38):  carrier_annotations_hg38.tsv, df_clusters_complete_kmeans.csv,
#                 freesurfer_zscore_demo_only.tsv (via load_mri())
#
# Outputs:
#   figures_genetics/loeuf_mri_scatter_hg19.pdf
#   figures_genetics/loeuf_mri_forest_hg19.pdf
#   figures_genetics/loeuf_iq_srs_hg19.pdf
#   tables_genetics/loeuf_mri_interaction_statistics_hg19.csv
#   figures_genetics_hg38/loeuf_mri_scatter_hg38.pdf
#   figures_genetics_hg38/loeuf_mri_forest_hg38.pdf
#   tables_genetics_hg38/loeuf_mri_interaction_statistics_hg38.csv
# =============================================================================

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

from _mri_loader import load_mri

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir    = os.path.normpath(os.path.join(_script_dir, "..", "..", ".."))

CLUSTERS_FILE  = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "df_clusters_complete_kmeans.csv")
CARRIER_HG38   = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "2_genetic_analysis", "outputs",
    "tables_hg38", "carrier_annotations_hg38.tsv")

OUT_FIG_HG19   = os.path.join(_script_dir, "..", "outputs", "figures_genetics")
OUT_TAB_HG19   = os.path.join(_script_dir, "..", "outputs", "tables_genetics")
OUT_FIG_HG38   = os.path.join(_script_dir, "..", "outputs", "figures_genetics_hg38")
OUT_TAB_HG38   = os.path.join(_script_dir, "..", "outputs", "tables_genetics_hg38")
for d in [OUT_FIG_HG19, OUT_TAB_HG19, OUT_FIG_HG38, OUT_TAB_HG38]:
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# Gene lists
# ---------------------------------------------------------------------------
# hg19 constrained: dellof + miss for SynGO and ChromEpiTF
GENE_LISTS_HG19_CONSTRAINED = [
    ("dellof_proteinconding_genes_best_score",         "Protein coding",              "gray"),
    ("dellof_chromepitf_constraint_genes_best_score",  "CHROM delloF (constrained)",  "#8E44AD"),
    ("miss_chromepitf_constraint_genes_best_score",    "CHROM Miss (constrained)",    "#5B2C6F"),
    ("dellof_syngo_constraint_genes_best_score",       "SynGO delloF (constrained)",  "#4A90E2"),
    ("miss_syngo_constraint_genes_best_score",         "SynGO Miss (constrained)",    "#1A5276"),
]

# hg19 all-carriers: dellof + miss for SynGO and ChromEpiTF
GENE_LISTS_HG19_ALL = [
    ("dellof_proteinconding_genes_best_score",  "Protein coding", "gray"),
    ("dellof_chromepitf_genes_best_score",      "CHROM delloF",   "#8E44AD"),
    ("miss_chromepitf_genes_best_score",        "CHROM Miss",     "#5B2C6F"),
    ("dellof_syngo_genes_best_score",           "SynGO delloF",   "#4A90E2"),
    ("miss_syngo_genes_best_score",             "SynGO Miss",     "#1A5276"),
]

# hg38 constrained: combined dellof (lof + del, most constrained gene wins)
GENE_LISTS_HG38_CONSTRAINED = [
    ("dellof_syngo_constraint_best_score_hg38",    "SynGO delloF (constrained)",  "#4A90E2"),
    ("dellof_chromepitf_constraint_best_score_hg38", "CHROM delloF (constrained)", "#8E44AD"),
]

# hg38 all-carriers: combined dellof unconstrained
GENE_LISTS_HG38_ALL = [
    ("dellof_syngo_best_score_hg38",       "SynGO delloF",  "#4A90E2"),
    ("dellof_chromepitf_best_score_hg38",  "CHROM delloF",  "#8E44AD"),
]

# Regions to highlight in scatter plots
HIGHLIGHT_REGIONS = [
    ("lh_superiortemporal_thickness",  "Left STG\n(thickness)"),
    ("lh_lateraloccipital_thickness",  "Left Lateral Occipital\n(thickness)"),
]

CLUSTER_COLORS = {
    "C1": "#7A8B47",
    "C2": "#ff9fa0",
    "C3": "#e7ba52",
    "NT": "#C1C2BC",
}

# LOEUF sizing tiers for IQ-SRS plot
LOEUF_COLUMN_IQ  = "dellof_syngo_chromepitf_constraint_genes_best_score"
GENE_NAME_COL_IQ = "dellof_syngo_chromepitf_constraint_genes_best_gene"
SIZE_BINS  = [(0.40, np.inf, 1, 60),   (0.25, 0.40, 7, 200),
              (0.10, 0.25, 10, 400),   (0.00, 0.10, 100, 800)]
SIZE_MAP   = {1: 60, 7: 200, 10: 400, 100: 800, 350: 900}

SYNGO_CARRIER_COL = "dellof_syngo_contraint_carrier"
CHROM_CARRIER_COL = "dellof_chromepitf_contraint_carrier"


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------
def neg_log10(series):
    s = series.dropna()
    return -np.log10(s[s > 0])


def loeuf_size_col(score_col, df):
    s = df[score_col].copy().fillna(np.inf).clip(upper=350)
    result = pd.Series(1, index=df.index)
    for lo, hi, code, _ in SIZE_BINS:
        mask = (s > lo) & (s <= hi)
        result[mask] = code
    result[s >= 350] = 350
    return result


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_hg19():
    print(f"Loading hg19 data from df_clusters_complete...")
    df = pd.read_csv(CLUSTERS_FILE, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df["_pop"] = df.get("PopulationS1", pd.Series("", index=df.index))
    df.loc[df["_pop"] == "TD", "Cluster"] = "NT"
    df = df[df["_pop"] != "other"].drop_duplicates("ID")

    df_m = load_mri()
    ct_cols = [c for c in df_m.columns
               if c.endswith("_thickness") and (c.startswith("lh_") or c.startswith("rh_"))
               and "MeanThickness" not in c]
    df = df.merge(df_m[["ID"] + ct_cols], on="ID", how="inner", suffixes=("", "_mri"))
    for c in ct_cols:
        if c + "_mri" in df.columns:
            df[c] = df[c + "_mri"]
            df.drop(columns=[c + "_mri"], inplace=True)

    print(f"  hg19 merged: {len(df)} individuals, {len(ct_cols)} CT regions")
    return df, ct_cols


def load_hg38():
    print(f"Loading hg38 data...")
    df38 = pd.read_table(CARRIER_HG38, low_memory=False)
    df38["ID"] = df38["ID"].astype(str)
    for col in ["PopulationS1", "Population1"]:
        if col in df38.columns:
            df38[col] = df38[col].str.replace("ID", "IDD").str.replace("TD", "NT")
    df38 = df38[df38["PopulationS1"] != "other"].drop_duplicates("ID")

    df19 = pd.read_csv(CLUSTERS_FILE, low_memory=False, usecols=lambda c: c in
                       ["ID", "Cluster", "IQ", "SRS_tscore",
                        "dellof_syngo_contraint_carrier", "dellof_chromepitf_contraint_carrier",
                        "dellof_syngo_chromepitf_constraint_genes_best_score",
                        "dellof_syngo_chromepitf_constraint_genes_best_gene"])
    df19["ID"] = df19["ID"].astype(str)
    df38 = df38.merge(df19, on="ID", how="left")
    df38.loc[df38["PopulationS1"] == "TD", "Cluster"] = "NT"

    # Combine lof + del into dellof equivalents (min LOEUF = most constrained gene)
    _COMBINE = [
        ("dellof_syngo_constraint_best_score_hg38",
         "lof_syngo_contraint_best_score",       "del_syngo_constraint_genes_best_score"),
        ("dellof_chromepitf_constraint_best_score_hg38",
         "lof_chromepitf_contraint_best_score",  "del_chromepitf_constraint_genes_best_score"),
        ("dellof_syngo_best_score_hg38",
         "lof_syngo_best_score",                 "del_syngo_genes_best_score"),
        ("dellof_chromepitf_best_score_hg38",
         "lof_chromepitf_best_score",            "del_chromepitf_genes_best_score"),
    ]
    for new_col, lof_col, del_col in _COMBINE:
        lof_s = df38[lof_col] if lof_col in df38.columns else pd.Series(np.nan, index=df38.index)
        del_s = df38[del_col] if del_col in df38.columns else pd.Series(np.nan, index=df38.index)
        df38[new_col] = np.fmin(lof_s.values, del_s.values)  # min ignoring NaN

    df_m = load_mri()
    ct_cols = [c for c in df_m.columns
               if (c.endswith("_thickness") or c.endswith("_area"))
               and (c.startswith("lh_") or c.startswith("rh_"))
               and "MeanThickness" not in c and "WhiteSurfArea" not in c]
    df38 = df38.merge(df_m[["ID"] + ct_cols], on="ID", how="inner")

    print(f"  hg38 merged: {len(df38)} individuals, {len(ct_cols)} cortical regions")
    return df38, ct_cols


# ---------------------------------------------------------------------------
# Correlation analysis
# ---------------------------------------------------------------------------
def run_correlation_analysis(df, gene_lists, ct_cols, label="hg19"):
    print(f"\nRunning LOEUF-MRI correlation analysis ({label})...")
    rows = []
    for col, label_gl, _ in gene_lists:
        if col not in df.columns:
            print(f"  Skipping {col}: column not found")
            continue
        x_full = neg_log10(df[col])
        n_carriers = len(x_full)
        for mri_col in ct_cols:
            if mri_col not in df.columns:
                continue
            mri = df[mri_col].dropna()
            common = x_full.index.intersection(mri.index)
            if len(common) < 10:
                continue
            x = x_full.loc[common].values
            y = mri.loc[common].values
            xz = scipy_stats.zscore(x)
            yz = scipy_stats.zscore(y)
            r, p = pearsonr(xz, yz)
            p_perm = np.nan
            m = re.match(r"^(lh|rh)_(.+)_(thickness|area)$", mri_col)
            hemisphere = ("left" if m.group(1) == "lh" else "right") if m else ""
            region     = m.group(2) if m else mri_col
            mri_type   = m.group(3) if m else ""
            rows.append({
                "genetic_type":    "LOEUF",
                "genetic_feature": col,
                "gene_list_label": label_gl,
                "mri_col":         mri_col,
                "mri_type":        mri_type,
                "hemisphere":      hemisphere,
                "region":          region,
                "sample_size":     len(common),
                "n_total_carriers": n_carriers,
                "correlation":     r,
                "p_value":         p,
                "p_permutation":   p_perm,
            })
        print(f"  {label_gl}: {n_carriers} carriers")

    df_stats = pd.DataFrame(rows)
    if df_stats.empty:
        return df_stats

    for col in df_stats["genetic_feature"].unique():
        idx = df_stats["genetic_feature"] == col
        valid = idx & df_stats["p_value"].notna()
        if valid.sum() > 0:
            _, pfdr, _, _ = multipletests(df_stats.loc[valid, "p_value"].values,
                                          method="fdr_bh")
            df_stats.loc[valid, "p_fdr"] = pfdr
    return df_stats


# ---------------------------------------------------------------------------
# Scatter plots for highlighted regions
# ---------------------------------------------------------------------------
def _scatter_region(ax, df, loeuf_col, loeuf_label, color, mri_col, r, p_val):
    x_raw = neg_log10(df[loeuf_col])
    mri   = df[mri_col].dropna()
    common = x_raw.index.intersection(mri.index)
    if len(common) < 5:
        return
    xs = x_raw.loc[common].values
    ys = mri.loc[common].values
    slope, intercept, *_ = scipy_stats.linregress(xs, ys)
    x_range = np.linspace(xs.min(), xs.max(), 100)
    ax.plot(x_range, slope * x_range + intercept,
            color=color, linewidth=2, zorder=4,
            label=f"{loeuf_label.replace(chr(10), ' ')} "
                  f"r={r:.3f}, p={p_val:.3f} (n={len(common)})")


def plot_scatter_panels(df, gene_lists, df_stats, out_dir, suffix="hg19"):
    n_regions = len(HIGHLIGHT_REGIONS)
    fig, axes = plt.subplots(1, n_regions, figsize=(6 * n_regions, 6), sharey=False)
    if n_regions == 1:
        axes = [axes]

    protein_col = gene_lists[0][0]

    for ax_idx, (mri_col, region_display) in enumerate(HIGHLIGHT_REGIONS):
        ax = axes[ax_idx]

        if protein_col in df.columns and mri_col in df.columns:
            x_bg = neg_log10(df[protein_col])
            y_bg = df[mri_col].dropna()
            common_bg = x_bg.index.intersection(y_bg.index)
            if len(common_bg) > 0:
                bg_df = pd.DataFrame({
                    "x": x_bg.loc[common_bg].values,
                    "y": y_bg.loc[common_bg].values,
                    "Cluster": df.loc[common_bg, "Cluster"].fillna("NT").values,
                })
                for cl, grp in bg_df.groupby("Cluster"):
                    ax.scatter(grp["x"], grp["y"],
                               color=CLUSTER_COLORS.get(cl, "#999999"),
                               alpha=0.6, s=25, linewidths=0, label=cl, zorder=2)

        for loeuf_col, loeuf_label, color in gene_lists:
            if loeuf_col not in df.columns or mri_col not in df.columns:
                continue
            row = df_stats[
                (df_stats["genetic_feature"] == loeuf_col) &
                (df_stats["mri_col"] == mri_col)
            ]
            if row.empty:
                continue
            r     = row["correlation"].values[0]
            p_val = row["p_value"].values[0]
            _scatter_region(ax, df, loeuf_col, loeuf_label, color, mri_col, r, p_val)

        ax.set_xlabel("-log10(LOEUF)", fontsize=10)
        ax.set_ylabel("Cortical thickness (mm)", fontsize=10)
        ax.set_title(region_display, fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=7, loc="upper left", frameon=False)

    plt.suptitle(f"LOEUF x Cortical Thickness ({suffix})",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(out_dir, f"loeuf_mri_scatter_{suffix}.pdf")
    plt.savefig(path, dpi=300, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"  Saved: loeuf_mri_scatter_{suffix}.pdf")


# ---------------------------------------------------------------------------
# Forest plot across all CT regions
# ---------------------------------------------------------------------------
def plot_forest(df_stats, gene_lists, out_dir, suffix="hg19"):
    ct_stats = df_stats[df_stats["mri_type"] == "thickness"].copy()
    if ct_stats.empty:
        return

    all_cols = sorted(ct_stats["mri_col"].unique())
    mean_r = {c: ct_stats[ct_stats["mri_col"] == c]["correlation"].mean() for c in all_cols}
    sorted_cols = sorted(all_cols, key=lambda c: mean_r.get(c, 0), reverse=True)
    region_labels = []
    for c in sorted_cols:
        m = re.match(r"^(lh|rh)_(.+)_thickness$", c)
        if m:
            side = "L" if m.group(1) == "lh" else "R"
            region_labels.append(f"{side} {m.group(2).replace('_', ' ')}")
        else:
            region_labels.append(c)

    n_gl  = len(gene_lists)
    n_reg = len(sorted_cols)
    fig, axes = plt.subplots(1, n_gl, figsize=(n_gl * 4, max(5, n_reg * 0.35 + 1.5)),
                             sharey=True)
    if n_gl == 1:
        axes = [axes]

    for ax_idx, (col, gl_label, color) in enumerate(gene_lists):
        ax = axes[ax_idx]
        ax.axvline(0, color="grey", lw=0.8, ls="--", zorder=0)
        sub = ct_stats[ct_stats["genetic_feature"] == col].set_index("mri_col")
        for i, mc in enumerate(sorted_cols):
            if mc not in sub.index:
                continue
            row   = sub.loc[mc]
            r     = row["correlation"]
            p_p   = row.get("p_fdr", 1.0)
            alpha = 1.0 if p_p < 0.05 else 0.3
            pt_color = "#4C72B0" if r > 0 else "#C44E52"
            ax.scatter([r], [i], color=pt_color, s=20, alpha=alpha, zorder=3)
            ax.plot([0, r], [i, i], color=pt_color, lw=0.6, alpha=alpha, zorder=2)

        ax.set_yticks(range(n_reg))
        if ax_idx == 0:
            ax.set_yticklabels(region_labels, fontsize=5.5)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("Pearson r", fontsize=8)
        ax.set_title(gl_label.replace("\n", " "), fontsize=9, fontweight="bold",
                     color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(f"LOEUF x Cortical Thickness — Pearson r across regions ({suffix})",
                 fontsize=10, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(out_dir, f"loeuf_mri_forest_{suffix}.pdf")
    plt.savefig(path, dpi=300, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"  Saved: loeuf_mri_forest_{suffix}.pdf")


# ---------------------------------------------------------------------------
# IQ vs SRS-2 scatter (hg19 only)
# ---------------------------------------------------------------------------
def plot_iq_srs(df, out_dir):
    required = ["IQ", "SRS_tscore", "Cluster", LOEUF_COLUMN_IQ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"  IQ-SRS plot: skipping (missing columns: {missing})")
        return

    df_plot = df.dropna(subset=["IQ", "SRS_tscore"]).copy()
    df_plot["_size_code"] = loeuf_size_col(LOEUF_COLUMN_IQ, df_plot)
    df_plot["_was_nan"]   = df[LOEUF_COLUMN_IQ].loc[df_plot.index].isna()

    fig, ax = plt.subplots(figsize=(9, 7))

    for cl, grp in df_plot.groupby("Cluster"):
        valid = grp[~grp["_was_nan"]]
        nan   = grp[grp["_was_nan"]]
        if len(valid) > 0:
            sizes = valid["_size_code"].map(SIZE_MAP)
            ax.scatter(valid["SRS_tscore"], valid["IQ"],
                       c=CLUSTER_COLORS.get(cl, "#999999"),
                       s=sizes, alpha=0.7, linewidths=0, zorder=2, label=cl)
        if len(nan) > 0:
            ax.scatter(nan["SRS_tscore"], nan["IQ"],
                       c=CLUSTER_COLORS.get(cl, "#999999"),
                       s=40, alpha=0.7, linewidths=0, marker="D", zorder=2)

    syngo_cols = SYNGO_CARRIER_COL in df_plot.columns
    chrom_cols = CHROM_CARRIER_COL in df_plot.columns

    if syngo_cols and chrom_cols:
        syngo = df_plot[df_plot[SYNGO_CARRIER_COL] == 1]
        chrom = df_plot[df_plot[CHROM_CARRIER_COL] == 1]
        both  = df_plot[(df_plot[SYNGO_CARRIER_COL] == 1) &
                        (df_plot[CHROM_CARRIER_COL] == 1)]

        def outline(sub, color, label):
            valid = sub[~sub["_was_nan"]]
            nan_s = sub[sub["_was_nan"]]
            if len(valid) > 0:
                sizes = valid["_size_code"].map(SIZE_MAP)
                ax.scatter(valid["SRS_tscore"], valid["IQ"],
                           facecolors="none", edgecolors=color,
                           s=sizes, linewidths=1.5, zorder=5, label=label)
            if len(nan_s) > 0:
                ax.scatter(nan_s["SRS_tscore"], nan_s["IQ"],
                           facecolors="none", edgecolors=color,
                           s=50, linewidths=1.5, marker="D", zorder=5)

        syngo_only = syngo[~syngo.index.isin(both.index)]
        chrom_only = chrom[~chrom.index.isin(both.index)]
        outline(syngo_only, "#4A90E2", "SynGO carrier")
        outline(chrom_only, "#8E44AD", "CHROM carrier")
        outline(both,       "#27ae60", "SynGO + CHROM")

    if GENE_NAME_COL_IQ in df_plot.columns:
        annotate = df_plot[
            df_plot["_size_code"].isin([10, 100]) &
            df_plot[GENE_NAME_COL_IQ].notna()
        ]
        for _, row in annotate.iterrows():
            ax.text(row["SRS_tscore"] + 1, row["IQ"] + 1,
                    str(row[GENE_NAME_COL_IQ]),
                    fontsize=7, style="italic", color="black",
                    ha="left", va="bottom", zorder=6)

    for srs_thresh in [60, 75]:
        ax.axvline(srs_thresh, color="black", linestyle="--", linewidth=0.8)
    for iq_thresh in [70, 130]:
        ax.axhline(iq_thresh, color="black", linestyle="--", linewidth=0.8)

    ax.set_xlabel("SRS-2 Total T-score", fontsize=11)
    ax.set_ylabel("Measured IQ", fontsize=11)
    ax.set_title("IQ vs SRS-2 with LOEUF constraint and carrier status",
                 fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles, labels_l = ax.get_legend_handles_labels()
    unique = dict(zip(labels_l, handles))
    ax.legend(unique.values(), unique.keys(),
              title="Cluster / Carrier", fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)

    plt.tight_layout()
    path = os.path.join(out_dir, "loeuf_iq_srs_hg19.pdf")
    plt.savefig(path, dpi=300, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"  Saved: loeuf_iq_srs_hg19.pdf")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    # ------------------------------------------------------------------
    # hg19
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("hg19 analysis")
    print("=" * 60)
    df19, ct_cols19 = load_hg19()

    # Constrained gene lists (dellof + miss)
    df_stats_constrained = run_correlation_analysis(
        df19, GENE_LISTS_HG19_CONSTRAINED, ct_cols19,
        label="hg19 constrained (dellof + miss)")
    df_stats_constrained.to_csv(
        os.path.join(OUT_TAB_HG19, "loeuf_mri_interaction_statistics_hg19.csv"),
        index=False)
    print(f"  Statistics saved: loeuf_mri_interaction_statistics_hg19.csv")

    # All-carriers gene lists for scatter plots
    df_stats_all = run_correlation_analysis(
        df19, GENE_LISTS_HG19_ALL, ct_cols19,
        label="hg19 all carriers (dellof + miss)")

    plot_scatter_panels(df19, GENE_LISTS_HG19_ALL, df_stats_all, OUT_FIG_HG19, suffix="hg19")
    plot_forest(df_stats_constrained, GENE_LISTS_HG19_CONSTRAINED, OUT_FIG_HG19, suffix="hg19")
    plot_iq_srs(df19, OUT_FIG_HG19)

    # ------------------------------------------------------------------
    # hg38
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("hg38 analysis")
    print("=" * 60)
    try:
        df38, ct_cols38 = load_hg38()

        df_stats38_c = run_correlation_analysis(
            df38, GENE_LISTS_HG38_CONSTRAINED, ct_cols38, label="hg38 constrained")
        df_stats38_a = run_correlation_analysis(
            df38, GENE_LISTS_HG38_ALL, ct_cols38, label="hg38 all carriers")

        df_stats38 = pd.concat([df_stats38_c, df_stats38_a], ignore_index=True)
        df_stats38.to_csv(
            os.path.join(OUT_TAB_HG38, "loeuf_mri_interaction_statistics_hg38.csv"),
            index=False)
        print(f"  Statistics saved: loeuf_mri_interaction_statistics_hg38.csv")

        plot_scatter_panels(df38, GENE_LISTS_HG38_CONSTRAINED, df_stats38_c,
                            OUT_FIG_HG38, suffix="hg38")
        plot_forest(df_stats38_c, GENE_LISTS_HG38_CONSTRAINED, OUT_FIG_HG38,
                    suffix="hg38")
    except Exception as e:
        print(f"  hg38 failed: {e}")

    print("\n=== Outputs complete ===")
    print(f"  hg19 figures : {OUT_FIG_HG19}")
    print(f"  hg38 figures : {OUT_FIG_HG38}")


if __name__ == "__main__":
    main()
