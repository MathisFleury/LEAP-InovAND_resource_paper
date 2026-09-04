#!/usr/bin/env python3
# =============================================================================
# 6 - Unified cluster-stability figure (k-means / Ward / GMM)
# =============================================================================
# Merges the three per-method stability figures into one supplementary figure:
#   row A = k-means (primary), row B = Ward, row C = GMM.
# Each row has the same four panels (per-individual modal probability,
# subsampling, noise, sample-size curve), read from the CSVs 5_cluster_
# stability writes.
#
# Curated is the priority regime → default reads the flat outputs/tables/.
# Usage:  python3 6_stability_figure_merged.py [frozen]
# Output: outputs/figures/cluster_stability_merged.pdf (frozen: self-
#         contained under legacy/1_clustering/outputs/, same as
#         5_cluster_stability.py's frozen wrapper).
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

_script_dir = os.path.dirname(os.path.abspath(__file__))
FROZEN = "frozen" in sys.argv[1:]
OUT_BASE = (os.path.join(_script_dir, "..", "..", "legacy", "1_clustering", "outputs") if FROZEN
            else os.path.join(_script_dir, "..", "outputs"))
OUT_BASE = os.path.normpath(OUT_BASE)
TABLES_DIR = os.path.join(OUT_BASE, "tables")
FIG_DIR = os.path.join(OUT_BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

METHODS = [("kmeans", "K-means (primary)"), ("ward", "Ward"), ("gmm", "GMM")]
COLOR = "#324095"


def _read(method, kind):
    return pd.read_csv(os.path.join(TABLES_DIR, f"stability_{kind}_{method}.csv"))


def _errbar(ax, x, d, dashed=None):
    yerr = np.vstack([d["ARI_mean"] - d["ARI_p2.5"], d["ARI_p97.5"] - d["ARI_mean"]])
    ax.errorbar(x, d["ARI_mean"], yerr=yerr, marker="o", ms=5, lw=1.6,
                color=COLOR, capsize=3)
    if dashed is not None:
        ax.axhline(dashed, ls="--", color="grey", lw=0.8)
    ax.set_ylim(0, 1.02)


COL_TITLES = ["Per-individual stability", "Subsampling", "Noise", "Sample-size curve"]

fig, axes = plt.subplots(len(METHODS), 4, figsize=(18, 11))

for r, (m, label) in enumerate(METHODS):
    # (1) per-individual modal probability histogram
    ax = axes[r][0]
    pi = _read(m, "per_individual")
    ax.hist(pi["modal_assignment_prob"], bins=np.linspace(0, 1, 31),
            color=COLOR, alpha=0.85, edgecolor="black", linewidth=0.4)
    ax.axvline(0.90, ls="--", color="grey", lw=0.8)
    ax.set_xlabel("Modal cluster assignment probability")
    ax.set_ylabel("Individuals")

    # (2) subsampling — x = subsample size
    ax = axes[r][1]
    sub = _read(m, "subsampling")
    _errbar(ax, sub["n_subset"], sub, dashed=0.75)
    ax.set_xlabel("Subsample size")
    ax.set_ylabel("ARI vs reference")

    # (3) noise — x = noise SD
    ax = axes[r][2]
    noi = _read(m, "noise")
    _errbar(ax, noi["noise_sd"], noi, dashed=0.75)
    ax.set_xlabel("Noise SD (std-space units)")
    ax.set_ylabel("ARI vs reference")

    # (4) sample-size curve — x = N
    ax = axes[r][3]
    ss = _read(m, "sample_size")
    _errbar(ax, ss["n"], ss, dashed=0.75)
    ax.set_xlabel("N (training subset)")
    ax.set_ylabel("ARI vs reference")

    # column titles on the top row only
    if r == 0:
        for c, t in enumerate(COL_TITLES):
            axes[r][c].set_title(t, fontsize=12)

    # big row label (A / B / C) + method name at the left of each row
    letter = chr(ord("A") + r)
    axes[r][0].annotate(letter, xy=(-0.28, 1.02), xycoords="axes fraction",
                        fontsize=26, fontweight="bold", va="bottom", ha="right")
    axes[r][0].annotate(label, xy=(-0.28, 0.5), xycoords="axes fraction",
                        fontsize=14, fontweight="bold", rotation=90,
                        va="center", ha="right")

fig.tight_layout(rect=(0.04, 0, 1, 1))
out = os.path.join(FIG_DIR, "cluster_stability_merged.pdf")
fig.savefig(out, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"{'[frozen]' if FROZEN else '[curated]'} saved: {out}")
