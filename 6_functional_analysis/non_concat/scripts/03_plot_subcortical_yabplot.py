#!/usr/bin/env python3
"""
Subcortical Functional Connectivity — Yabplot Visualization
Autism vs NT

Uses the custom 4S156 subcortical atlas (CIT168Subcortical, SubcorticalHCP,
Cerebellum, ThalamusHCP) built from the BIDS NIfTI.

For each region present in the custom atlas, counts the number of significant
FDR-corrected edges involving that region.

Run _resources/build_custom_subcortical_atlas.py once before this script.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import pyvista as pv

pv.global_theme.notebook = True

from yabplot import plot_subcortical

# =============================================================================
# PATHS
# =============================================================================
import os

_SCRIPT_DIR  = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
# Honor AUTISM_TD_OUTPUT_DIR (same env var as 01/02) so legacy/revised
# subcortical PNGs don't overwrite each other.
_DEFAULT_DIR = _SECTION_DIR / 'outputs' / 'figures' / 'revised'
INPUT_DIR    = Path(os.environ.get('AUTISM_TD_OUTPUT_DIR', _DEFAULT_DIR))
OUTPUT_DIR   = INPUT_DIR

CUSTOM_ATLAS_DIR = str(
    Path(__file__).parent.parent.parent.parent /
    '_resources' / 'custom_subcortical_atlas_4S156'
)

HYPER_FILE = INPUT_DIR / 'autism_vs_td_hyperconnectivity_full_data.csv'
HYPO_FILE  = INPUT_DIR / 'autism_vs_td_hypoconnectivity_full_data.csv'


def get_custom_atlas_regions(atlas_dir: str) -> set:
    """Return set of region names from VTK files in the custom atlas directory."""
    return {p.stem for p in Path(atlas_dir).glob('*.vtk')}


def count_edges_per_region(df: pd.DataFrame, atlas_regions: set) -> dict:
    """
    Count how many edges involve each custom atlas region
    (source or target endpoint matches a region in the atlas).
    """
    counts = {r: 0 for r in atlas_regions}
    for _, row in df.iterrows():
        for endpoint in (row['source'], row['target']):
            if endpoint in counts:
                counts[endpoint] += 1
    # 0 → NaN so those regions appear grey
    return {k: (float(v) if v > 0 else float('nan')) for k, v in counts.items()}


def plot_and_save(data: dict, n_max: float, out_path: Path, title: str, colormap: str):
    plotter = plot_subcortical(
        data=data,
        custom_atlas_path=CUSTOM_ATLAS_DIR,
        cmap=colormap,
        vminmax=[1, n_max],
        display_type='object',
        figsize=(1200, 500),
    )
    if title:
        plotter.add_title(title, font_size=12)
    plotter.screenshot(str(out_path))
    plotter.close()
    print(f"  Saved: {out_path.name}")


def main():
    print("=" * 60)
    print("SUBCORTICAL CONNECTIVITY YABPLOT — Autism vs NT")
    print("=" * 60)
    print(f"Custom atlas: {CUSTOM_ATLAS_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    atlas_regions = get_custom_atlas_regions(CUSTOM_ATLAS_DIR)
    if not atlas_regions:
        print(f"ERROR: No VTK meshes found in {CUSTOM_ATLAS_DIR}")
        print("Run _resources/build_custom_subcortical_atlas.py first.")
        return 1
    print(f"Custom atlas regions available: {len(atlas_regions)}\n")

    for label, filepath, colormap in [
        ('Hyperconnectivity', HYPER_FILE, 'Reds'),
        ('Hypoconnectivity',  HYPO_FILE,  'Blues'),
    ]:
        print(f"{label}:")
        if not filepath.exists():
            print(f"  WARNING: {filepath.name} not found, skipping.")
            continue

        df = pd.read_csv(filepath)
        print(f"  Loaded {len(df)} edges")

        counts = count_edges_per_region(df, atlas_regions)
        valid  = {k: v for k, v in counts.items() if not np.isnan(v)}
        if valid:
            n_max = max(valid.values())
            print(f"  Regions with edges: {len(valid)}, max: {int(n_max)}")
            for region, n in sorted(valid.items(), key=lambda x: -x[1]):
                print(f"    {region}: {int(n)}")
            title = f'Autism vs NT — {label} (No. of edges, FDR < 0.05)'
        else:
            n_max = 1.0
            print(f"  No subcortical edges found for {label.lower()} "
                  f"— writing empty-atlas PNG for parity with other variants.")
            title = f'Autism vs NT — {label} (no subcortical edges, FDR < 0.05)'

        suffix = 'hyper' if 'Hyper' in label else 'hypo'
        plot_and_save(
            counts, n_max,
            OUTPUT_DIR / f'subcortical_connectivity_autism_vs_nt_{suffix}_yabplot.png',
            title=title,
            colormap=colormap,
        )

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
