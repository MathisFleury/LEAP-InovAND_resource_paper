#!/usr/bin/env python3.11
"""
Composite QC figure — Framewise Displacement, Euler number, ComBat AUC-ROC.

2x2 grid of pre-rendered panel PDFs (each already produced by its own
analysis script — this script only lays them out):
  a. Framewise displacement (Power et al.) by phenotype group
       <- 6_functional_analysis/concat/scripts/12_motion_confound_check.py
  b. FreeSurfer Euler number by phenotype group (structural QC companion to a)
       <- 4_anatomical_analysis/scripts/17_euler_confound_check.py
  c. ComBat AUC-ROC (fMRI): SVM scanner-classification AUC before vs after
     ComBat, per scanner
       <- 6_functional_analysis/concat/scripts/13_sequence_combat_svm.py (by-machine)
  d. ComBat AUC-ROC (structural MRI): same check on the anatomical pipeline
     (structural companion to c)
       <- orphaned output at 4_anatomical_analysis/outputs/qc_sequence/
          anat_sequence_svm_auc_bymachine.pdf (generating script not in repo)

Re-run the source scripts first if the underlying data changed; this script
just assembles their current PDF outputs.
"""
from pathlib import Path

import fitz  # PyMuPDF

ROOT = Path(__file__).parent.parent
ROWS = [
    [("a", ROOT / "6_functional_analysis/concat/outputs/qc_motion/fd_violin.pdf"),
     ("b", ROOT / "4_anatomical_analysis/outputs/figures/qc_euler/euler_violin.pdf")],
    [("c", ROOT / "6_functional_analysis/concat/outputs/qc_sequence/sequence_svm_auc_bymachine.pdf"),
     ("d", ROOT / "4_anatomical_analysis/outputs/qc_sequence/anat_sequence_svm_auc_bymachine.pdf")],
]
OUT = Path(__file__).parent / "Figure_QC_FD_Euler_ComBatROC.pdf"

CM = 28.3465  # pt per cm
MARGIN = 0.7 * CM
LABEL_BAND = 0.7 * CM
GAP_COL = 0.5 * CM
GAP_ROW = 0.6 * CM
CONTENT_W = 18 * CM


def main():
    missing = [p for row in ROWS for _, p in row if not p.exists()]
    if missing:
        raise SystemExit("Missing panel PDF(s), run the source script(s) first:\n"
                          + "\n".join(f"  {p}" for p in missing))

    row_heights = []
    for row in ROWS:
        n = len(row)
        col_w = (CONTENT_W - GAP_COL * (n - 1)) / n
        heights = []
        for _, p in row:
            d = fitz.open(p)
            r = d[0].rect
            heights.append(col_w * r.height / r.width)
            d.close()
        row_heights.append((col_w, heights))

    page_w = CONTENT_W + 2 * MARGIN
    page_h = (MARGIN + sum(LABEL_BAND + max(hs) for _, hs in row_heights)
              + GAP_ROW * (len(ROWS) - 1) + MARGIN)

    doc = fitz.open()
    page = doc.new_page(width=page_w, height=page_h)

    y = MARGIN
    for row, (col_w, heights) in zip(ROWS, row_heights):
        row_h = max(heights)
        x = MARGIN
        for (letter, path), h in zip(row, heights):
            page.insert_text((x, y + LABEL_BAND - 4), letter,
                              fontsize=16, fontname="Helvetica-Bold")
            src = fitz.open(path)
            y_panel = y + LABEL_BAND + (row_h - h) / 2  # vertically centre in the row band
            page.show_pdf_page(fitz.Rect(x, y_panel, x + col_w, y_panel + h), src, 0)
            src.close()
            x += col_w + GAP_COL
        y += LABEL_BAND + row_h + GAP_ROW

    doc.save(OUT)
    print(f"Saved: {OUT}  ({page_w / CM:.1f} x {page_h / CM:.1f} cm)")


if __name__ == "__main__":
    main()
