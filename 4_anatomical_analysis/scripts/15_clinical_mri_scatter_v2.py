#!/usr/bin/env python3
"""
Scatter / correlation plots for selected MRI x clinical associations (autism).
Reuses the revised-MRI + curated-clinical autism dataset from
12_clinical_mri_correlations_v2.build_dataset().

Panels:
  1. L superior-temporal thickness  vs  SRS-2 t-score
  2. L superior-temporal thickness  vs  full-scale IQ
  3. L cerebellum (cortex) volume    vs  full-scale IQ

MRI values are the revised z-scored (ComBat + age/sex/eTIV/euler-regressed) units.
Output: ../outputs/figures_clinical_mri/clinical_mri_scatter.pdf (+ .png)
"""

import os
import sys
from pathlib import Path
from importlib import import_module

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent))
_m12 = import_module("12_clinical_mri_correlations_v2")
import _config  # noqa: E402

OUT_DIR = Path(_config.OUTPUT_BASE) / "figures_clinical_mri"

# (x column, y column, x label, y label, title, point colour)
PANELS = [
    ("lh_superiortemporal_thickness", "SRS_tscore",
     "L superior-temporal thickness (z)", "SRS-2 t-score", "STG thickness vs SRS-2", "#4A90E2"),
    ("lh_superiortemporal_thickness", "total_IQ",
     "L superior-temporal thickness (z)", "Full-scale IQ", "STG thickness vs FSIQ", "#7A8B47"),
    ("Left-Cerebellum-Cortex", "total_IQ",
     "L cerebellum volume (z)", "Full-scale IQ", "L cerebellum volume vs FSIQ", "#8E44AD"),
]


def main():
    df, _, _ = _m12.build_dataset()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    for ax, (xc, yc, xl, yl, ttl, col) in zip(axes, PANELS):
        d = df[[xc, yc]].apply(lambda s: __import__("pandas").to_numeric(s, errors="coerce")).dropna()
        r, p = pearsonr(d[xc], d[yc])
        sns.regplot(x=xc, y=yc, data=d, ax=ax, color=col,
                    scatter_kws=dict(s=16, alpha=0.55, edgecolor="none"),
                    line_kws=dict(color="black", lw=1.6))
        ax.set_xlabel(xl, fontsize=10); ax.set_ylabel(yl, fontsize=10)
        ax.set_title(ttl, fontsize=11, fontweight="bold")
        pstr = f"p = {p:.3f}" if p >= 1e-3 else f"p = {p:.1e}"
        ax.text(0.04, 0.96, f"r = {r:.2f}\n{pstr}\nn = {len(d)}", transform=ax.transAxes,
                va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.85))
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("MRI × clinical associations — autism (revised anat, curated clinical)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = OUT_DIR / "clinical_mri_scatter.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
