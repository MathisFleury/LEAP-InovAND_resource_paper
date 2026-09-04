#!/usr/bin/env python3
# =============================================================================
# 04 - IQ x SRS scatter enriched with clusters, genetics and adaptive behaviour
# =============================================================================
# Four-panel scatter (patchworklib) of SRS-2 t-score (x) vs full-scale IQ (y):
#   - coloured by cluster or by Population1
#   - marker size = LOEUF of the best constrained DEL/LoF gene (binned on the
#     gnomAD LOEUF percentile cut-offs; largest = most constrained), or the
#     VABS-II ABC adaptive score
#   - HCNDD dom/rec-X carriers get a coloured outline; the most constrained
#     genes are annotated in italics (labels repelled off the circles)
# Plus a single extra panel where the carrier outlines mark SynGO / ChromEpiTF
# / both gene sets. Diamonds mark individuals missing the size variable; dashed
# lines are the clinical cut-offs (SRS 60/75, IQ 70/130).
#
# Clinical + clusters: the curated k-means solution written by
#   1_clustering/scripts/04_run_clustering.py (in-project output).
# Genetics: gnomAD v4 carrier matrix (CARRIER_V4) — LOEUF-scored DEL/LoF calls.
# Output:   figures/figure_final_IQ_SRS.pdf, figures/figure_IQ_SRS_geneset.pdf
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import patchworklib as pw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    CLUSTERS_CURATED_FILE, CARRIER_V4, FIGURES_DIR, TABLES_DIR,
    PALETTE_CLUSTERS, PALETTE_POPULATION1,
)

# HCNDD dominant / rec-X constrained carriers -> outline colour
CARRIER_TYPES = [("dellof_hcndddomv6_contraint_carrier", "HCNDD-DOM-X-REC", "green")]

# Gene-set panel. Columns gs_syngo/gs_chrom/gs_both are derived in load_data()
# so the three groups are mutually exclusive.
GENESET_CARRIERS = [
    ("gs_syngo", "SynGO", "#4a90e2"),
    ("gs_chrom", "ChromEpiTF", "#8e44ad"),
    ("gs_both", "SynGO + ChromEpiTF", "#008000"),
]

# curated table uses the old TD/ID spelling; align to the project's NT/IDD
POP_RELABEL = {
    "TD": "NT", "ID": "IDD",
    "Autism with ID": "Autism with IDD", "Autism without ID": "Autism without IDD",
}

# LOEUF bands (band 5 = most constrained = biggest circle). Stricter than the
# raw gnomAD percentiles so only genuinely constrained genes get a big circle.
LOEUF_EDGES = [0.10, 0.20, 0.36, 0.45]  # ascending; band = 5,4,3,2 else 1
LOEUF_SIZE_MAP = {1: 60, 2: 120, 3: 220, 4: 450, 5: 850}
LOEUF_LABELS = {
    5: "LOEUF < 0.10",
    4: "0.10 ≤ LOEUF < 0.20",
    3: "0.20 ≤ LOEUF < 0.36",
    2: "0.36 ≤ LOEUF < 0.45",
    1: "LOEUF ≥ 0.45",
}
LOEUF_ANNOTATE_BELOW = 0.36  # label constrained carrier genes below this LOEUF

# VABS-II ABC adaptive score bands (bigger circle = more impaired)
VIN_SIZE_MAP = {40: 900, 60: 450, 80: 225, 100: 100}
VIN_LABELS = {40: "< 50", 60: "50–79", 80: "80–99", 100: "≥ 100"}


def load_data():
    """Build the plotting frame: curated k-means clusters + v4 genetics.

    Clinical scores and the k-means Cluster label come from the in-project
    curated clustering output; genetics from the gnomAD v4 carrier matrix.
    Matching a row against the v4 matrix (inner join) is what restricts the
    figure to the WGS-genotyped cohort."""
    df = pd.read_csv(CLUSTERS_CURATED_FILE, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df = df.drop_duplicates("ID")
    df["VINELANDt"] = df["vabsabcabc_standard"]  # in-project table has no VINELANDt col
    df["Population1"] = df["Population1"].replace(POP_RELABEL)
    # show NT participants in their own colour instead of a k-means C-label
    df.loc[df["Population1"] == "NT", "Cluster"] = "NT"

    gen = pd.read_table(CARRIER_V4, low_memory=False)
    gen["ID"] = gen["ID"].astype(str)
    gene_cols = ["dellof_hcndddomv6_constraint_genes_best_gene",
                 "dellof_hcndddomv6_constraint_genes_best_score",
                 "dellof_constraint_genes_best_score",
                 "dellof_constraint_genes_best_gene",
                 # gene-set-specific best genes so a ringed carrier is labelled
                 # with the gene of its ring (not a more-constrained other gene)
                 "dellof_syngo_constraint_genes_best_gene",
                 "dellof_syngo_constraint_genes_best_score",
                 "dellof_chromepitf_constraint_genes_best_gene",
                 "dellof_chromepitf_constraint_genes_best_score",
                 "dellof_syngo_chromepitf_constraint_genes_best_gene",
                 "dellof_syngo_chromepitf_constraint_genes_best_score"]
    carrier_cols = [c[0] for c in CARRIER_TYPES] + [
        "dellof_syngo_contraint_carrier", "dellof_chromepitf_contraint_carrier"]
    # inner join = keep only WGS-genotyped individuals (present in the v4 matrix)
    df = df.merge(gen[["ID"] + gene_cols + carrier_cols], on="ID", how="inner")

    df = df.dropna(subset=["IQ", "SRS_tscore", "Cluster"])
    print(f"  Individuals with complete data: {len(df)}")

    # mutually exclusive gene-set carrier groups
    sy = df["dellof_syngo_contraint_carrier"] == 1
    ch = df["dellof_chromepitf_contraint_carrier"] == 1
    df["gs_syngo"] = (sy & ~ch).astype(int)
    df["gs_chrom"] = (ch & ~sy).astype(int)
    df["gs_both"] = (sy & ch).astype(int)

    # default = overall most-constrained gene (used for non-carrier points)
    df["gene"] = df["dellof_constraint_genes_best_gene"]
    df["LOEUF_raw"] = df["dellof_constraint_genes_best_score"].replace([np.inf, -np.inf], np.nan)
    # gene-set carriers: label by the gene of their ring so the two agree
    # (a SynGO ring shows its SynGO gene, e.g. ITSN1 — not a non-SynGO gene).
    for mask, gcol, scol in [
        (df["gs_syngo"] == 1, "dellof_syngo_constraint_genes_best_gene",
         "dellof_syngo_constraint_genes_best_score"),
        (df["gs_chrom"] == 1, "dellof_chromepitf_constraint_genes_best_gene",
         "dellof_chromepitf_constraint_genes_best_score"),
        (df["gs_both"] == 1, "dellof_syngo_chromepitf_constraint_genes_best_gene",
         "dellof_syngo_chromepitf_constraint_genes_best_score"),
    ]:
        df.loc[mask, "gene"] = df.loc[mask, gcol]
        df.loc[mask, "LOEUF_raw"] = df.loc[mask, scol].replace([np.inf, -np.inf], np.nan)
    df["was_nan_LOEUF"] = df["LOEUF_raw"].isna()
    lo = df["LOEUF_raw"].fillna(2.0)  # missing -> smallest band
    df["LOEUF"] = np.select(
        [lo < LOEUF_EDGES[0], lo < LOEUF_EDGES[1], lo < LOEUF_EDGES[2], lo < LOEUF_EDGES[3]],
        [5, 4, 3, 2], default=1,
    )

    df["was_nan_VINELAND"] = df["VINELANDt"].isna()
    vin = df["VINELANDt"].fillna(1)
    df["VINELANDt"] = np.select(
        [vin < 50, vin < 80, vin < 100, vin >= 100], [40, 60, 80, 100], default=1,
    )
    return df


def _legend(ax, df, hue, palette, size_map, size_labels, size_title, carriers):
    """Build an explicit legend: hue groups, size bands (grey circles), carriers."""
    handles = [Line2D([], [], ls="", label=hue, alpha=0)]  # section header
    for g in palette:
        if (df[hue] == g).any():
            handles.append(Line2D([], [], marker="o", ls="", markerfacecolor=palette[g],
                                  markeredgecolor="none", markersize=8, label=str(g)))
    handles.append(Line2D([], [], ls="", label=size_title, alpha=0))
    for key in sorted(size_map):
        handles.append(Line2D([], [], marker="o", ls="", markerfacecolor="0.6",
                              markeredgecolor="none", markersize=size_map[key] ** 0.5,
                              label=size_labels[key]))
    handles.append(Line2D([], [], ls="", label="Carrier", alpha=0))
    for _, label, color in carriers:
        handles.append(Line2D([], [], marker="o", ls="", markerfacecolor="none",
                              markeredgecolor=color, markeredgewidth=1.5, markersize=9, label=label))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.05, 1.0),
              fontsize=7.5, labelspacing=1.1, borderpad=0.8, frameon=False)


def _overlay(ax, df, size_col, size_map, nan_col, carriers, annotate):
    """Carrier outlines + gene annotations + clinical cut-off lines on one axis."""
    for col, label, color in carriers:
        sub = df[df[col] == 1]
        valid = sub[~sub[nan_col]]
        ax.scatter(valid["SRS_tscore"], valid["IQ"], facecolors="none", edgecolors=color,
                   s=valid[size_col].map(size_map), linewidths=1.5, marker="o", label=label)
        miss = sub[sub[nan_col]]
        ax.scatter(miss["SRS_tscore"], miss["IQ"], facecolors="none", edgecolors=color,
                   s=min(size_map.values()), linewidths=1.2, marker="D")

    if annotate:
        # label only this panel's outlined carriers, most-constrained band,
        # deduplicated by gene — labelling every constrained gene is unreadable
        is_carrier = np.logical_or.reduce([df[c[0]] == 1 for c in carriers])
        # label every constrained carrier point (no per-gene de-duplication)
        genes = (df[is_carrier & (df["LOEUF_raw"] < LOEUF_ANNOTATE_BELOW)]
                 .dropna(subset=["gene"]))
        texts, targets = [], []
        for _, r in genes.iterrows():
            texts.append(ax.text(r["SRS_tscore"], r["IQ"], r["gene"], style="italic",
                                 fontsize=7, color="black", zorder=6))
            # matplotlib scatter s -> marker DIAMETER is sqrt(s) pts, so radius = sqrt(s)/2
            targets.append((r["SRS_tscore"], r["IQ"], 0.5 * size_map[r["LOEUF"]] ** 0.5))
        if texts:
            try:
                from adjustText import adjust_text
                pts = df.dropna(subset=["SRS_tscore", "IQ"])
                # move labels first (no arrows), then draw connectors that stop
                # at each circle's edge via a per-target shrinkB = its radius
                adjust_text(texts, x=pts["SRS_tscore"].tolist(), y=pts["IQ"].tolist(), ax=ax,
                            expand=(1.6, 2.0), force_text=(0.8, 1.2), force_static=(0.4, 0.6))
                # connector under the text (zorder<6, text drawn on top); shrinkB =
                # circle radius + margin so it stops outside the circle
                for t, (cx, cy, r_pts) in zip(texts, targets):
                    ax.annotate("", xy=(cx, cy), xytext=t.get_position(), zorder=4,
                                arrowprops=dict(arrowstyle="-", color="0.35", lw=0.7,
                                                shrinkA=6, shrinkB=r_pts + 2))
            except Exception as e:  # keep the figure even if label repel fails
                print(f"  WARNING: adjust_text failed ({e}); labels left un-repelled.")

    for xv in (60, 75):
        ax.axvline(xv, color="black", ls="--", lw=1)
    for yv in (70, 130):
        ax.axhline(yv, color="black", ls="--", lw=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("SRS-2 total T-score")
    ax.set_ylabel("Measured IQ")
    ax.set_title("")


def _panel(df, hue, palette, size_col, missing_size, carriers=CARRIER_TYPES):
    """One pw.Brick scatter: round markers for known size, diamonds for missing."""
    unmapped = sorted(set(df[hue].dropna().unique()) - set(palette))
    if unmapped:
        n_drop = int(df[hue].isin(unmapped).sum())
        print(f"  NOTE: dropping {n_drop} rows with {hue} values not in the "
              f"palette (no colour defined): {unmapped}")
        df = df[~df[hue].isin(unmapped)]

    is_loeuf = size_col == "LOEUF"
    size_map = LOEUF_SIZE_MAP if is_loeuf else VIN_SIZE_MAP
    size_labels = LOEUF_LABELS if is_loeuf else VIN_LABELS
    size_title = "LOEUF" if is_loeuf else "VABS-II ABC"
    nan_col = "was_nan_LOEUF" if is_loeuf else "was_nan_VINELAND"

    ax = pw.Brick(figsize=(10, 7))
    sns.scatterplot(data=df[~df[nan_col]], x="SRS_tscore", y="IQ", hue=hue, size=size_col,
                    sizes=size_map, palette=palette, ax=ax, alpha=0.7, edgecolor="none",
                    marker="o", legend=False)
    miss = df[df[nan_col]]
    sns.scatterplot(data=miss, x="SRS_tscore", y="IQ", hue=hue, palette=palette, s=missing_size,
                    ax=ax, alpha=0.7, edgecolor="none", marker="D", legend=False)
    _overlay(ax, df, size_col, size_map, nan_col, carriers, annotate=is_loeuf)
    _legend(ax, df, hue, palette, size_map, size_labels, size_title, carriers)
    return ax


def export_r_frame(df):
    """Dump the plotting frame for the R (ggrepel) version of the figure.
    Keeps the tested merge/band logic here so R only has to plot."""
    cols = ["ID", "SRS_tscore", "IQ", "Cluster", "Population1",
            "LOEUF", "was_nan_LOEUF", "LOEUF_raw", "gene",
            "VINELANDt", "was_nan_VINELAND",
            "dellof_hcndddomv6_contraint_carrier", "gs_syngo", "gs_chrom", "gs_both",
            # HCNDD-specific best gene/score (for the HCNDD presentation variant)
            "dellof_hcndddomv6_constraint_genes_best_gene",
            "dellof_hcndddomv6_constraint_genes_best_score"]
    out = os.path.join(TABLES_DIR, "iq_srs_plot_frame.csv")
    os.makedirs(TABLES_DIR, exist_ok=True)
    df[cols].to_csv(out, index=False)
    print(f"  Saved: {out}")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    df = load_data()
    export_r_frame(df)

    b1 = _panel(df, "Cluster", PALETTE_CLUSTERS, "LOEUF", 40)
    b2 = _panel(df, "Population1", PALETTE_POPULATION1, "LOEUF", 40)
    b3 = _panel(df, "Cluster", PALETTE_CLUSTERS, "VINELANDt", 30)
    b4 = _panel(df, "Population1", PALETTE_POPULATION1, "VINELANDt", 30)

    final = (b3 | b4) / (b1 | b2)
    out = os.path.join(FIGURES_DIR, "figure_final_IQ_SRS.pdf")
    final.savefig(out, dpi=300)
    print(f"  Saved: {out}")

    # Extra single panel: LOEUF-sized circles, SynGO / ChromEpiTF / both outlines
    gs = _panel(df, "Cluster", PALETTE_CLUSTERS, "LOEUF", 40, carriers=GENESET_CARRIERS)
    out_gs = os.path.join(FIGURES_DIR, "figure_IQ_SRS_geneset.pdf")
    gs.savefig(out_gs, dpi=300)
    print(f"  Saved: {out_gs}")


if __name__ == "__main__":
    main()
