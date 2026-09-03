#!/usr/bin/env python3
"""
Subcortical Volume Brain Visualization — Autism vs NT (REVISION, R3)

v2 of ../../scripts/04_plot_subcortical_yabplot.py.  Renders Cohen's d
on the aseg atlas using a fixed [-0.4, 0.4] color scale to match the
cortical brain maps emitted by 02_plot_anatomical_mri_brain_v2.R and the
per-cluster v2 plots in 5_cluster_anatomical_analysis/curated/.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import pyvista as pv

pv.global_theme.notebook = True

from yabplot import plot_subcortical, get_atlas_regions

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
INPUT_DIR = _SECTION_DIR / "outputs" / "figures" / "r_input_files"
OUTPUT_DIR = _SECTION_DIR / "outputs" / "figures"
INPUT_FILE = INPUT_DIR / "t_stat_anat_aseg_volume_mri_autism_vs_control.csv"

# Cohen's d colorbar limit is computed from the data (= max |d|, with a
# 0.1 floor) so no region is clipped to NA.
D_FLOOR = 0.1

LABEL_MAP = {
    "Left-Thalamus":     "Left_Thalamus",
    "Right-Thalamus":    "Right_Thalamus",
    "Left-Caudate":      "Left_Caudate",
    "Right-Caudate":     "Right_Caudate",
    "Left-Putamen":      "Left_Putamen",
    "Right-Putamen":     "Right_Putamen",
    "Left-Pallidum":     "Left_Pallidum",
    "Right-Pallidum":    "Right_Pallidum",
    "Left-Hippocampus":  "Left_Hippocampus",
    "Right-Hippocampus": "Right_Hippocampus",
    "Left-Amygdala":     "Left_Amygdala",
    "Right-Amygdala":    "Right_Amygdala",
}

ALL_REGIONS = get_atlas_regions("aseg", "subcortical")


def load_data(filepath):
    return pd.read_csv(filepath)


def make_yabplot_data(df, fdr_only=False):
    data = {}
    for _, row in df.iterrows():
        yab_label = LABEL_MAP.get(row["label"])
        if not yab_label:
            continue
        if fdr_only and row.get("p_fdr", 1.0) >= 0.05:
            data[yab_label] = float("nan")
        else:
            data[yab_label] = float(row["cohens_d"])
    for r in ALL_REGIONS:
        data.setdefault(r, float("nan"))
    return data


def plot_and_save(data, vmax, out_path, title=""):
    plotter = plot_subcortical(
        data=data, atlas="aseg", cmap="RdBu_r",
        vminmax=[-vmax, vmax], display_type="object",
        figsize=(1200, 500),
    )
    if title:
        plotter.add_title(title, font_size=12)
    plotter.screenshot(str(out_path))
    plotter.close()
    print(f"  Saved: {out_path.name}")


def main():
    print("=" * 60)
    print("Subcortical yabplot — Autism vs NT  (Cohen's d, REVISION)")
    print("=" * 60)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT_FILE.exists():
        print(f"ERROR: input file not found ({INPUT_FILE}). "
              "Run 01_anatomical_mri_autism_nt_v2.py first.")
        return 1
    df = load_data(INPUT_FILE)
    d_max = max(float(np.nanmax(np.abs(df["cohens_d"].values))), D_FLOOR)
    print(f"  d_max (data-driven): {d_max:.3f}")

    data_all = make_yabplot_data(df, fdr_only=False)
    plot_and_save(
        data_all, d_max,
        OUTPUT_DIR / "autism_vs_nt_volume_aseg_yabplot.png",
        title=f"Autism vs NT - Subcortical Volume (Cohen's d, max={d_max:.2f})",
    )
    data_fdr = make_yabplot_data(df, fdr_only=True)
    plot_and_save(
        data_fdr, d_max,
        OUTPUT_DIR / "autism_vs_nt_volume_aseg_yabplot_fdr.png",
        title=f"Autism vs NT - Subcortical Volume (FDR<0.05, Cohen's d, max={d_max:.2f})",
    )
    print(f"\nDone. Outputs in {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
