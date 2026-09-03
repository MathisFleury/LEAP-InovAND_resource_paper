#!/usr/bin/env python3
"""
Build the composite figure:
  a) thickness brain maps (SynGO, CHROM)      -> left column
  b) two significant-permutation boxplots      -> right column
       row1: SynGO  brain  <-> lh caudal ant. cingulate (SynGO  perm p<0.05)
       row2: CHROM  brain  <-> rh transversetemporal    (CHROM  perm p<0.05)
  c) IQ x SRS cluster/LOEUF panel (reframed 1.35x wide) -> centered below

Boxplots reuse 08's regperm tables; x-ticks relabeled (All Genes->Protein Coding),
main title + "Gene List" x-title dropped, legend dropped (x-ticks carry the labels).

Prereqs (run first):
  Rscript 18_thickness_maps_for_composite.R
  Rscript ../../10_clinical_analysis/scripts/04e_plot_iq_srs_cluster_reframe.R

Run:  python3 19_build_composite.py
"""
import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import fitz  # PyMuPDF

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _config  # noqa: E402

TABLES = Path(_config.OUTPUT_BASE) / "tables_genetics_hg38_regperm"
COMP   = Path(_config.OUTPUT_BASE) / "figures_genetics_hg38_regperm" / "composite"
COMP.mkdir(parents=True, exist_ok=True)

PANEL_C = (Path(__file__).resolve().parents[3] / "10_clinical_analysis" / "outputs"
           / "figures" / "figure_IQ_SRS_cluster_loeuf_reframe.pdf")

CM = 28.3465  # pt per cm

# ROI -> nice label. Order mirrors panel a (ChromEpiTF/CHROM left, SynGO right).
ROIS = [
    ("rh_transversetemporal_thickness",      "R transverse temporal thickness"),
    ("lh_caudalanteriorcingulate_thickness", "L caudal anterior cingulate thickness"),
]
ORDER  = ["Protein coding", "CHROM", "SynGO"]
XLABEL = {"Protein coding": "Protein Coding", "CHROM": "CHROM", "SynGO": "SYNGO"}
COLORS = {"Protein coding": "gray", "CHROM": "#8E44AD", "SynGO": "#4A90E2"}


def make_boxplot(alldf, mri_col, nice):
    sub = alldf[alldf["mri_col"] == mri_col]
    rows = {r["gene_list_label"]: r for _, r in sub.iterrows()}
    present = [g for g in ORDER if g in rows]

    fig, ax = plt.subplots(figsize=(3.4, 3.6))
    ax.axhline(0, ls="--", color="grey", lw=1)
    for i, g in enumerate(present):
        r = rows[g]
        beta = r.get("beta_robust", r["beta"])
        se, p = r["se"], r["p_permutation"]
        sig = p < 0.05
        ax.errorbar(i, beta, yerr=se, fmt="o", color=COLORS[g],
                    capsize=5, capthick=2, markersize=12 if sig else 10,
                    markeredgecolor="black" if sig else COLORS[g],
                    markeredgewidth=1.5 if sig else 0.5,
                    elinewidth=2.5 if sig else 2, alpha=0.85, zorder=3)
        if sig:
            ax.text(i, beta + se + 0.03, f"p = {p:.3f}", ha="center",
                    va="bottom", fontsize=10, fontweight="bold", color="black")

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 0.10 * (ymax - ymin))
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([XLABEL[g] for g in present])
    ax.set_xlim(-0.6, len(present) - 0.4)
    ax.set_ylabel(f"Beta Coefficient\n({nice} vs -log LOEUF)", fontsize=9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    out = COMP / f"box_{mri_col}.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def place(page, src_path, x, y_top, height):
    """Place src PDF's first page with top-left at (x, y_top), scaled to `height`,
    aspect preserved. Returns (width, height) occupied."""
    d = fitz.open(src_path)
    p0 = d[0]
    sw, sh = p0.rect.width, p0.rect.height
    w = height * sw / sh
    page.show_pdf_page(fitz.Rect(x, y_top, x + w, y_top + height), d, 0)
    d.close()
    return w, height


def label(page, txt, x, y):
    page.insert_text((x, y), txt, fontsize=16, fontname="Helvetica-Bold")


def main():
    files = glob.glob(str(TABLES / "loeuf_reg_*.csv"))
    alldf = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    boxes = {mc: make_boxplot(alldf, mc, nice) for mc, nice in ROIS}

    panel_a = COMP / "panel_a_thickness.pdf"     # ChromEpiTF | SynGO, shared bar
    for pth in [panel_a, PANEL_C]:
        if not Path(pth).exists():
            sys.exit(f"Missing prerequisite: {pth}\nRun the two Rscripts first (see docstring).")

    M = 0.7 * CM          # page margin
    GAP_AB = 0.6 * CM     # gap between panel a and panel b
    GAP_BX = 0.3 * CM     # gap between the two boxplots
    GAP_ROW = 1.0 * CM    # gap between top band and panel c
    LBL = 0.6 * CM        # label band above each panel

    H_A = 5.5 * CM        # panel a (maps) height
    H_BOX = 6.5 * CM      # boxplot height
    H_C = 18 * CM         # panel c height — dominant

    def wof(path, h):
        d = fitz.open(path); r = d[0].rect; d.close()
        return h * r.width / r.height

    # --- top band geometry: [ a ]  [ box1 box2 ] ---
    a_w = wof(panel_a, H_A)
    bx_w = [wof(boxes[mc], H_BOX) for mc, _ in ROIS]
    b_w = sum(bx_w) + GAP_BX
    band_h = max(H_A, H_BOX)
    top_w = a_w + GAP_AB + b_w
    c_w = wof(PANEL_C, H_C)

    content_w = max(top_w, c_w)
    page_w = content_w + 2 * M
    page_h = M + LBL + band_h + GAP_ROW + LBL + H_C + M

    doc = fitz.open()
    page = doc.new_page(width=page_w, height=page_h)

    x0, y0 = M, M + LBL
    # panel a (vertically centered in the band)
    place(page, panel_a, x0, y0 + (band_h - H_A) / 2, H_A)
    label(page, "a", x0, y0 - 4)
    # panel b: two boxplots side by side
    xb = x0 + a_w + GAP_AB
    label(page, "b", xb, y0 - 4)
    for (mc, _), w in zip(ROIS, bx_w):
        place(page, boxes[mc], xb, y0 + (band_h - H_BOX) / 2, H_BOX)
        xb += w + GAP_BX

    # panel c centered, big
    y_c = M + LBL + band_h + GAP_ROW + LBL
    x_c = M + (content_w - c_w) / 2
    place(page, PANEL_C, x_c, y_c, H_C)
    label(page, "c", M, y_c - 4)

    out = COMP / "figure_composite.pdf"
    doc.save(out)
    print(f"Saved: {out}  ({page_w/CM:.1f} x {page_h/CM:.1f} cm)")


if __name__ == "__main__":
    main()
