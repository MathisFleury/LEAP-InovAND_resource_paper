#!/usr/bin/env python3.11
"""
Standalone figure — diagnosis x age / x sex interaction sanity check
(reviewer-requested robustness panel), full sample.

Side by side:
  a. Diagnosis x age interaction, combined ggseg brain maps (4 metrics)
       <- 11_age_sex_analysis/scripts/08_combined_interaction_brain_figure.R
  b. Diagnosis x sex interaction, combined ggseg brain maps (4 metrics)
       <- same script

Re-run 08_combined_interaction_brain_figure.R first if the underlying
interaction stats changed; this script just assembles its current outputs.
"""
from pathlib import Path

import fitz  # PyMuPDF

ROOT = Path(__file__).parent.parent
FIG_DIR = ROOT / "11_age_sex_analysis/outputs/figures"
PANELS = [
    ("a", FIG_DIR / "combined_interaction_brain_maps_age.pdf"),
    ("b", FIG_DIR / "combined_interaction_brain_maps_sex.pdf"),
]
OUT = Path(__file__).parent / "Figure_AgeSex_Interaction.pdf"

CM = 28.3465  # pt per cm
MARGIN = 0.7 * CM
LABEL_BAND = 0.7 * CM
GAP = 0.5 * CM
PANEL_H = 22 * CM


def main():
    missing = [p for _, p in PANELS if not p.exists()]
    if missing:
        raise SystemExit("Missing panel PDF(s), run 08_combined_interaction_brain_figure.R first:\n"
                          + "\n".join(f"  {p}" for p in missing))

    widths = []
    for _, p in PANELS:
        d = fitz.open(p)
        r = d[0].rect
        widths.append(PANEL_H * r.width / r.height)
        d.close()

    page_w = sum(widths) + GAP * (len(PANELS) - 1) + 2 * MARGIN
    page_h = MARGIN + LABEL_BAND + PANEL_H + MARGIN

    doc = fitz.open()
    page = doc.new_page(width=page_w, height=page_h)

    x = MARGIN
    y = MARGIN
    for (letter, path), w in zip(PANELS, widths):
        page.insert_text((x, y + LABEL_BAND - 4), letter,
                          fontsize=16, fontname="Helvetica-Bold")
        src = fitz.open(path)
        page.show_pdf_page(fitz.Rect(x, y + LABEL_BAND, x + w, y + LABEL_BAND + PANEL_H), src, 0)
        src.close()
        x += w + GAP

    doc.save(OUT)
    print(f"Saved: {OUT}  ({page_w / CM:.1f} x {page_h / CM:.1f} cm)")


if __name__ == "__main__":
    main()
