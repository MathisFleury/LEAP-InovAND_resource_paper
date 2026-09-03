#!/usr/bin/env python3
# =============================================================================
# 09 - Cluster concept map (IQ x SRS, k-means)
# =============================================================================
# Conceptual, minimalist view of the k-means solution: each cluster is a soft
# density "territory" in IQ x SRS space. No points, no marginal distributions,
# no numeric clutter — just the three regions, their centroids, and directional
# axes. Cluster palette is the project-wide PALETTE_CLUSTERS.
#
# Usage:  python3 09_cluster_jointplot.py [kmeans|ward|gmm]   (default kmeans)
# Output: outputs/figures/cluster_conceptmap[ _ward | _gmm ].pdf
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, to_rgb
import seaborn as sns

_script_dir = os.path.dirname(os.path.abspath(__file__))
# Moved to legacy/1_clustering/: TABLES_DIR still points at the real
# 1_clustering/outputs/tables/ (04_run_clustering.py's shared output --
# a genuine read dependency, not just historical co-location). FIGURES_DIR
# is self-contained here since nothing else reads this script's output.
_ORIG_SECTION = os.path.normpath(os.path.join(_script_dir, "..", "..", "1_clustering"))
TABLES_DIR = os.path.join(_ORIG_SECTION, "outputs", "tables")
FIGURES_DIR = os.path.normpath(os.path.join(_script_dir, "outputs", "figures"))

method = sys.argv[1] if len(sys.argv) > 1 else "kmeans"
suffix = "" if method == "kmeans" else f"_{method}"
TABLE = os.path.join(TABLES_DIR, f"individuals_metrics_with_clusters{suffix}.csv")
OUT = os.path.join(FIGURES_DIR, f"cluster_conceptmap{suffix}.pdf")

PALETTE_CLUSTERS = {"C1": "#7A8B47", "C2": "#ff9fa0", "C3": "#e7ba52"}
INK, MUTED = "#2b2b2b", "#8a8a8a"

df = pd.read_csv(TABLE, low_memory=False)
df = df[df["Cluster"].notna() & df["IQ"].notna() & df["SRS_tscore"].notna()].copy()
order = sorted(PALETTE_CLUSTERS)

plt.rcParams.update({"font.family": "sans-serif", "text.color": INK})

fig, ax = plt.subplots(figsize=(7.5, 7.5))

for c in order:
    sub = df[df["Cluster"] == c]
    col = PALETTE_CLUSTERS[c]
    x, y = sub["SRS_tscore"], sub["IQ"]
    # soft organic blob: white->colour ramp, no outline
    cmap = LinearSegmentedColormap.from_list(c, ["white", col])
    sns.kdeplot(x=x, y=y, ax=ax, fill=True, cmap=cmap, levels=8,
                thresh=0.10, bw_adjust=1.4, alpha=0.55, zorder=1)
    # centroid dot + big label
    cx, cy = x.mean(), y.mean()
    ax.scatter(cx, cy, s=90, color=col, edgecolor="white", linewidth=2, zorder=3)
    ax.annotate(c, (cx, cy), fontsize=22, fontweight="bold", color=col,
                ha="center", va="center", xytext=(0, 30), textcoords="offset points",
                path_effects=[pe.withStroke(linewidth=4, foreground="white")], zorder=4)

# --- minimalist directional axes (arrows, no ticks) ---
ax.set_xlabel(""); ax.set_ylabel("")          # drop seaborn's auto labels
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)
xmin, xmax = df["SRS_tscore"].min() - 5, df["SRS_tscore"].max() + 5
ymin, ymax = df["IQ"].min() - 8, df["IQ"].max() + 8
ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
arrow = dict(arrowstyle="-|>", color=MUTED, lw=1.6)
ax.annotate("", xy=(xmax, ymin), xytext=(xmin, ymin), arrowprops=arrow)
ax.annotate("", xy=(xmin, ymax), xytext=(xmin, ymin), arrowprops=arrow)
ax.text(xmax, ymin - (ymax - ymin) * 0.05, "higher SRS-2 t-score  →",
        ha="right", va="top", fontsize=13, color=MUTED)
ax.text(xmin - (xmax - xmin) * 0.035, ymax, "higher IQ  →", rotation=90,
        ha="center", va="top", fontsize=13, color=MUTED)

plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"N={len(df)}  sizes: {df['Cluster'].value_counts().sort_index().to_dict()}")
print(f"saved: {OUT}")
