#!/usr/bin/env python3
"""
Per-ROI beta-coefficient figures (LOEUF x MRI, hg38 regperm).

For each ROI in ROIS, plot the OLS beta of (ROI metric ~ -log10 LOEUF best-score)
for the three gene lists (Protein Coding / CHROM / SynGO), with +/-1 SE bars and
permutation-test significance stars. Styled after
eeg_mri-pipeline/.../genetics_mri_analysis/beta_coefficients_iq_loeuf_final.pdf.

Reads the tables written by 8_loeuf_mri_regression_permutation.py
  (outputs/tables_genetics_hg38_regperm/loeuf_reg_*.csv)
and writes to outputs/figures_genetics_hg38_regperm/.

Run:  python3 10_beta_coefficients_roi.py
"""

import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _config  # noqa: E402

TABLES_DIR = Path(_config.OUTPUT_BASE) / "tables_genetics_hg38_regperm"
FIG_DIR = Path(_config.OUTPUT_BASE) / "figures_genetics_hg38_regperm"

# ROIs to feature (significant in the permutation test). Edit freely.
ROIS = [
    ("rh_pericalcarine_thickness",            "R pericalcarine thickness"),
    ("lh_caudalanteriorcingulate_thickness", "L caudal anterior cingulate thickness"),
]

# x order + colours (match the reference figure).
ORDER = ["Protein coding", "CHROM", "SynGO"]
XLABEL = {"Protein coding": "All Genes", "CHROM": "CHROM", "SynGO": "SYNGO"}
COLORS = {"Protein coding": "gray", "CHROM": "#8E44AD", "SynGO": "#4A90E2"}


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def main():
    files = glob.glob(str(TABLES_DIR / "loeuf_reg_*.csv"))
    if not files:
        sys.exit(f"No regperm tables in {TABLES_DIR} — run 08 first.")
    alldf = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    for mri_col, nice in ROIS:
        sub = alldf[alldf["mri_col"] == mri_col]
        if sub.empty:
            print(f"skip {mri_col}: not in tables"); continue
        rows = {r["gene_list_label"]: r for _, r in sub.iterrows()}
        present = [g for g in ORDER if g in rows]

        fig, ax = plt.subplots(figsize=(4.2, 4.6))
        ax.axhline(0, ls="--", color="grey", lw=1)
        for i, g in enumerate(present):
            r = rows[g]
            beta = r.get("beta_robust", r["beta"])  # robust (bootstrap-median) beta
            se, p = r["se"], r["p_permutation"]
            sig = p < 0.05
            ax.errorbar(i, beta, yerr=se, fmt="o", color=COLORS[g],  # +/-1 SE (paper)
                        capsize=5, capthick=2, markersize=12 if sig else 10,
                        markeredgecolor="black" if sig else COLORS[g],
                        markeredgewidth=1.5 if sig else 0.5,
                        elinewidth=2.5 if sig else 2, alpha=0.85, zorder=3)
            if sig:  # p-value on top of the significant point(s)
                ax.text(i, beta + se + 0.03, f"p = {p:.3f}", ha="center",
                        va="bottom", fontsize=10, fontweight="bold", color="black")

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax + 0.10 * (ymax - ymin))   # headroom for sig stars
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels([XLABEL[g] for g in present])
        ax.set_xlim(-0.6, len(present) - 0.4)
        ax.set_xlabel("Gene List")
        ax.set_ylabel(f"Beta Coefficient\n({nice} vs -log LOEUF)", fontsize=10)
        ax.set_title(f"Beta Coefficients\n{nice} vs -log LOEUF", fontsize=12, fontweight="bold")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

        # legend
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker="o", color="w", label=XLABEL[g],
                          markerfacecolor=COLORS[g], markersize=10,
                          markeredgecolor="black" if rows[g]["p_permutation"] < 0.05 else COLORS[g])
                   for g in present]
        ax.legend(handles=handles, title="Gene List", loc="center left",
                  bbox_to_anchor=(1.02, 0.5), frameon=True, fontsize=9)

        fig.tight_layout()
        slug = mri_col.replace("_thickness", "").replace("_", "-")
        out = FIG_DIR / f"beta_coefficients_{slug}_loeuf.pdf"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {out}")
        for g in present:
            r = rows[g]
            print(f"    {g:14s} beta={r.get('beta_robust', r['beta']):+.3f} se={r['se']:.3f} "
                  f"p_perm={r['p_permutation']:.4f} {stars(r['p_permutation'])}")


if __name__ == "__main__":
    main()
