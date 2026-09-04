#!/usr/bin/env python3
# =============================================================================
# 10 - Cluster scatter panels (IQ x SRS), coloured 4 ways
# =============================================================================
# Four small sklearn-KMeans-style scatter plots of the clustering space
# (Measured IQ x SRS-2 t-score), the same points coloured by:
#   (a) k-means Cluster (+ centroids)   (b) Autism / IDD / NT
#   (c) ADHD vs no ADHD                 (d) Male / Female
#
# Curated is the priority regime → reads the curated k-means table, writes to
# outputs/curated/figures/. Saves the 2x2 panel + each panel as its own PDF.
#
# Usage:  python3 10_cluster_scatter_panels.py [frozen]
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_dir = os.path.normpath(os.path.join(_script_dir, "..", ".."))
sys.path.insert(0, _project_dir)
import config  # noqa: E402

FROZEN = "frozen" in sys.argv[1:]
TABLE = config.CLUSTERS_KMEANS_FROZEN if FROZEN else config.CLUSTERS_CURATED
FIG_DIR = os.path.join(_script_dir, "..", "outputs",
                       *(() if FROZEN else ("curated",)), "figures")
FIG_DIR = os.path.normpath(FIG_DIR)
os.makedirs(FIG_DIR, exist_ok=True)

# --- Colours -----------------------------------------------------------------
GREY = "#d9d9d9"
COL_CLUSTER = {"C1": "#7A8B47", "C2": "#ff9fa0", "C3": "#e7ba52"}
COL_POP = {"Autism without IDD": "#5CAEE1", "Autism with IDD": "#324095",
           "IDD": "#D8A4CB", "NT": "#C1C2BC", "Other": GREY}
COL_ADHD = {"ADHD": "#ff7f0e", "No ADHD": "#C1C2BC", "Unknown": GREY}  # orange (ref: figures_papers)
COL_SEX = {"Male": "#5B9E6F", "Female": "#8E44AD"}  # muted green / purple — clear of the ADHD orange

# --- Data (clustered cohort only) --------------------------------------------
df = pd.read_csv(TABLE, low_memory=False)
df = df[df["Cluster"].notna() & df["IQ"].notna() & df["SRS_tscore"].notna()].copy()

# Derive the categorical labels for each colouring.
df["_cluster"] = df["Cluster"]
_POP_KEEP = ("Autism without IDD", "Autism with IDD", "IDD", "NT")
# fold the lone 'Autism to exclude' subject (IQ 94, no IDD) into Autism w/o IDD
df["_pop"] = (df["Population1"].replace({"Autism to exclude": "Autism without IDD"})
              .map(lambda v: v if v in _POP_KEEP else "Other"))
df["_adhd"] = df["pheno_ADHD"].map({2.0: "ADHD", 1.0: "No ADHD"}).fillna("Unknown")
# Sex: 1 = Male, 0 = Female (autism cohort is male-skewed, ~2.8:1 here).
# ponytail: coding assumed from the M/F ratio — flip the map if the codebook says otherwise.
df["_sex"] = df["Sex"].map({1.0: "Male", 0.0: "Female"})

PANELS = [
    ("k-means clusters", "_cluster", COL_CLUSTER, ["C1", "C2", "C3"], True),
    ("NT / Autism w/o IDD / IDD", "_pop", COL_POP,
     ["Autism without IDD", "Autism with IDD", "IDD", "NT"], False),
    ("ADHD", "_adhd", COL_ADHD, ["ADHD", "No ADHD"], False),
    ("Sex", "_sex", COL_SEX, ["Male", "Female"], False),
]

X, Y = "SRS_tscore", "IQ"
xlim = (df[X].min() - 4, df[X].max() + 4)
ylim = (df[Y].min() - 6, df[Y].max() + 6)

# Clinical reference thresholds: IQ 70 / 130, SRS-2 60 / 75.
IQ_LINES = [70, 130]
SRS_LINES = [60, 75]


def draw_panel(ax, col, palette, order, centroids):
    present = [g for g in order if (df[col] == g).any()]
    for g in present:
        sub = df[df[col] == g]
        ax.scatter(sub[X], sub[Y], s=9, color=palette[g], alpha=0.65,
                   edgecolor="none", zorder=2, label=f"{g} (n={len(sub)})")
    if centroids:  # k-means centroids in (SRS, IQ) space
        for g in present:
            sub = df[df[col] == g]
            ax.scatter(sub[X].mean(), sub[Y].mean(), marker="o", s=90,
                       color=palette[g], edgecolor="black", linewidth=1.5, zorder=5)
    # Clinical reference lines (labelled via ticks at the threshold values).
    for xv in SRS_LINES:
        ax.axvline(xv, color="#555555", ls=(0, (4, 3)), lw=0.9, zorder=4)
    for yv in IQ_LINES:
        ax.axhline(yv, color="#555555", ls=(0, (4, 3)), lw=0.9, zorder=4)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xticks(SRS_LINES); ax.set_yticks(IQ_LINES)
    ax.tick_params(labelsize=7)
    ax.set_aspect(abs((xlim[1] - xlim[0]) / (ylim[1] - ylim[0])))
    for s in ax.spines.values():
        s.set_edgecolor("#999999")
    ax.legend(loc="upper right", fontsize=6.5, frameon=False,
              handletextpad=0.2, borderaxespad=0.2, labelspacing=0.25)


# --- Combined 2x2 ------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(7.6, 7.8))
for ax, (title, col, pal, order, cen) in zip(axes.ravel(), PANELS):
    draw_panel(ax, col, pal, order, cen)
    ax.set_title(title, fontsize=12, fontweight="bold")
fig.supxlabel("SRS-2 t-score", fontsize=11)
fig.supylabel("Measured IQ", fontsize=11)
fig.tight_layout()
combined = os.path.join(FIG_DIR, "cluster_scatter_panels.pdf")
fig.savefig(combined, dpi=300, bbox_inches="tight")
plt.close(fig)

# --- Combined single column (cluster -> population -> ADHD -> Sex) -----------
figc, axesc = plt.subplots(len(PANELS), 1, figsize=(3.4, 3.2 * len(PANELS)))
for ax, (title, col, pal, order, cen) in zip(axesc, PANELS):
    draw_panel(ax, col, pal, order, cen)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("SRS-2 t-score", fontsize=8)   # per-panel axis labels
    ax.set_ylabel("Measured IQ", fontsize=8)
figc.tight_layout()
figc.savefig(os.path.join(FIG_DIR, "cluster_scatter_panels_1col.pdf"),
             dpi=300, bbox_inches="tight")
plt.close(figc)

# --- Each panel as its own little figure -------------------------------------
slug = {"k-means clusters": "cluster", "NT / Autism w/o IDD / IDD": "population",
        "ADHD": "adhd", "Sex": "sex"}
for title, col, pal, order, cen in PANELS:
    f, ax = plt.subplots(figsize=(3.0, 3.0))
    draw_panel(ax, col, pal, order, cen)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("SRS-2 t-score", fontsize=8)
    ax.set_ylabel("Measured IQ", fontsize=8)
    out = os.path.join(FIG_DIR, f"cluster_scatter_{slug[title]}.pdf")
    f.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(f)

print(f"{'[frozen]' if FROZEN else '[curated]'} N={len(df)}  ->  {FIG_DIR}")
print("  cluster_scatter_panels.pdf + cluster_scatter_{cluster,population,adhd,sex}.pdf")
