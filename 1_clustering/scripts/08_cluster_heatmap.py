#!/usr/bin/env python3
# =============================================================================
# 08 - Cluster heatmap (Figure 4a style, dendrogram removed)
# =============================================================================
# Recreates the paper's Figure-4 panel-a heatmap but WITHOUT the dendrogram:
# rows are ordered by cluster label (k-means primary) instead of by a Ward
# tree, since the primary method is no longer hierarchical.
#
# Layout (left→right): phenotype annotation strips (Autism / IDD / ADHD,
# grey = absent, black = unknown), z-scored "Measured IQ values" and
# "SRS-2 t-score" heatmap columns (RdBu_r, ±2), cluster strip on the right.
#
# Usage:  python3 08_cluster_heatmap.py [kmeans|ward|gmm]   (default kmeans)
# Output: outputs/figures/cluster_heatmap[ _ward | _gmm ].pdf
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize

_script_dir = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.normpath(os.path.join(_script_dir, "..", "outputs", "tables"))
FIGURES_DIR = os.path.normpath(os.path.join(_script_dir, "..", "outputs", "figures"))

method = sys.argv[1] if len(sys.argv) > 1 else "kmeans"
suffix = "" if method == "kmeans" else f"_{method}"
TABLE = os.path.join(TABLES_DIR, f"individuals_metrics_with_clusters{suffix}.csv")
OUT = os.path.join(FIGURES_DIR, f"cluster_heatmap{suffix}.pdf")

# --- Colours ---
PALETTE_CLUSTERS = {"C1": "#7A8B47", "C2": "#ff9fa0", "C3": "#e7ba52"}
PHENO_COLORS = {"Autism": "#6a5acd", "IDD": "#ff69b4", "ADHD": "#c71585"}
ABSENT, UNKNOWN = "#d9d9d9", "#000000"

# --- Data ---
df = pd.read_csv(TABLE, low_memory=False)
df = df[df["Cluster"].notna()].copy()

# z-score the two heatmap features across the whole cohort
z_iq = (df["IQ"] - df["IQ"].mean()) / df["IQ"].std()
z_srs = (df["SRS_tscore"] - df["SRS_tscore"].mean()) / df["SRS_tscore"].std()
df["z_IQ"], df["z_SRS"] = z_iq, z_srs

# Cluster order top→bottom by descending mean IQ (C1/C2 high, C3 low),
# ponytail: within-cluster sort by IQ gives the smooth gradient the
# dendrogram used to provide — no tree needed.
cluster_order = df.groupby("Cluster")["IQ"].mean().sort_values(ascending=False).index.tolist()
df["_ck"] = df["Cluster"].map({c: i for i, c in enumerate(cluster_order)})
df = df.sort_values(["_ck", "IQ"], ascending=[True, False]).reset_index(drop=True)

n = len(df)

# --- Phenotype annotation matrix (RGB) ---
pheno_cols = [("Autism", "pheno_ASD"), ("IDD", "pheno_ID"), ("ADHD", "pheno_ADHD")]
from matplotlib.colors import to_rgb

anno = np.ones((n, len(pheno_cols), 3))
for j, (name, col) in enumerate(pheno_cols):
    for i, v in enumerate(df[col].values):
        if v == -9 or pd.isna(v):
            anno[i, j] = to_rgb(UNKNOWN)
        elif v == 2:                       # phenotype present
            anno[i, j] = to_rgb(PHENO_COLORS[name])
        else:                              # absent
            anno[i, j] = to_rgb(ABSENT)

# --- Figure ---
fig = plt.figure(figsize=(6, 9))
# width ratios: 3 pheno strips | 2 z-score cols | cluster strip
gs = fig.add_gridspec(1, 3, width_ratios=[3, 2.2, 0.4], wspace=0.08)

# Panel 1: phenotype strips
ax0 = fig.add_subplot(gs[0])
ax0.imshow(anno, aspect="auto", interpolation="nearest")
ax0.set_xticks(range(len(pheno_cols)))
ax0.set_xticklabels([p[0] for p in pheno_cols], rotation=90, fontsize=9)
ax0.set_yticks([])
for s in ax0.spines.values():
    s.set_visible(False)

# Panel 2: z-score heatmap (IQ, SRS)
ax1 = fig.add_subplot(gs[1])
Z = df[["z_IQ", "z_SRS"]].values
norm = Normalize(vmin=-2, vmax=2)
im = ax1.imshow(Z, aspect="auto", cmap="RdBu_r", norm=norm, interpolation="nearest")
ax1.set_xticks([0, 1])
ax1.set_xticklabels(["Measured IQ values", "SRS-2 t-score"], rotation=90, fontsize=9)
ax1.set_yticks([])
for s in ax1.spines.values():
    s.set_visible(False)

# Panel 3: cluster strip + labels
ax2 = fig.add_subplot(gs[2])
cvals = df["Cluster"].map({c: i for i, c in enumerate(cluster_order)}).values
cmap = ListedColormap([PALETTE_CLUSTERS[c] for c in cluster_order])
ax2.imshow(cvals.reshape(-1, 1), aspect="auto", cmap=cmap, interpolation="nearest")
ax2.set_xticks([])
ax2.set_yticks([])
# cluster name at the vertical centre of each block
start = 0
for c in cluster_order:
    size = int((df["Cluster"] == c).sum())
    ax2.text(0.7, start + size / 2, c, va="center", ha="left",
             fontsize=12, fontweight="bold")
    start += size
for s in ax2.spines.values():
    s.set_visible(False)

# --- Legends ---
from matplotlib.patches import Patch
legend_items = [Patch(facecolor=PHENO_COLORS["Autism"], label="Phenotype Autism"),
                Patch(facecolor=PHENO_COLORS["IDD"], label="Phenotype IDD"),
                Patch(facecolor=PHENO_COLORS["ADHD"], label="Phenotype ADHD"),
                Patch(facecolor=UNKNOWN, label="Unknown")]
ax0.legend(handles=legend_items, loc="lower left", bbox_to_anchor=(0, 1.02),
           fontsize=8, frameon=False, ncol=2)

# colorbar
cbar = fig.colorbar(im, ax=ax2, fraction=0.6, pad=1.4, ticks=[-2, -1, 0, 1, 2])
cbar.set_label("z-score", fontsize=9)
cbar.ax.tick_params(labelsize=8)

plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"N={n}  cluster order (top→bottom): {cluster_order}")
print(f"saved: {OUT}")

# ponytail self-check: annotation matrix is fully classified, no white gaps
assert not np.any(np.all(anno == 1.0, axis=2)), "unclassified phenotype cell"
