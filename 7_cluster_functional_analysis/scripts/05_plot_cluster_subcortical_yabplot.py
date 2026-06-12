#!/usr/bin/env python3
"""
Subcortical Functional Connectivity — Yabplot Visualization
Per autism cluster vs NT

Uses the custom 4S156 subcortical atlas (CIT168Subcortical, SubcorticalHCP,
Cerebellum, ThalamusHCP) built from the BIDS NIfTI.

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
_SCRIPT_DIR  = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
INPUT_DIR    = _SECTION_DIR / 'outputs' / 'figures'
OUTPUT_DIR   = _SECTION_DIR / 'outputs' / 'figures'

CUSTOM_ATLAS_DIR = str(
    Path(__file__).parent.parent.parent /
    '_resources' / 'custom_subcortical_atlas_4S156'
)

CLUSTERS = ['C1', 'C2', 'C3']


def get_custom_atlas_regions(atlas_dir: str) -> set:
    return {p.stem for p in Path(atlas_dir).glob('*.vtk')}


def count_edges_per_region(df: pd.DataFrame, atlas_regions: set) -> dict:
    """
    Parse source/target from 'region' column (format: con_source/target),
    count edges involving each custom atlas region.
    """
    counts = {r: 0 for r in atlas_regions}
    for region in df['region']:
        clean = region.removeprefix('con_')
        parts = clean.split('/', 1)
        source = parts[0]
        target = parts[1] if len(parts) > 1 else ''
        for endpoint in (source, target):
            if endpoint in counts:
                counts[endpoint] += 1
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
    print(f"    Saved: {out_path.name}")


def main():
    print("=" * 60)
    print("SUBCORTICAL CONNECTIVITY YABPLOT — Clusters vs NT")
    print("=" * 60)
    print(f"Custom atlas: {CUSTOM_ATLAS_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    atlas_regions = get_custom_atlas_regions(CUSTOM_ATLAS_DIR)
    if not atlas_regions:
        print(f"ERROR: No VTK meshes found in {CUSTOM_ATLAS_DIR}")
        print("Run _resources/build_custom_subcortical_atlas.py first.")
        return 1
    print(f"Custom atlas regions available: {len(atlas_regions)}\n")

    for cluster in CLUSTERS:
        print(f"--- Cluster {cluster} ---")
        csv_path = INPUT_DIR / f'cluster_{cluster}_connectivity_autism_vs_td.csv'
        if not csv_path.exists():
            print(f"  WARNING: {csv_path.name} not found, skipping.")
            continue

        df = pd.read_csv(csv_path)
        df_sig = df[df['p_fdr'] < 0.05]
        # Hyper / hypo by Cohen's d direction (Reviewer 3); falls back to t_stat
        # if cohens_d is absent in legacy outputs.
        dir_col = 'cohens_d' if 'cohens_d' in df_sig.columns else 't_stat'
        df_hyper = df_sig[df_sig[dir_col] > 0]
        df_hypo  = df_sig[df_sig[dir_col] < 0]
        print(f"  FDR-sig: {len(df_sig)} total, {len(df_hyper)} hyper, {len(df_hypo)} hypo")

        for direction, df_dir, colormap, label in [
            ('hyper', df_hyper, 'Reds',  'Hyperconnectivity'),
            ('hypo',  df_hypo,  'Blues', 'Hypoconnectivity'),
        ]:
            if df_dir.empty:
                print(f"  No {direction} edges — skipping.")
                continue

            counts = count_edges_per_region(df_dir, atlas_regions)
            valid  = {k: v for k, v in counts.items() if not np.isnan(v)}
            if not valid:
                print(f"  No subcortical {direction} edges found.")
                continue

            n_max = max(valid.values())
            print(f"  {label}: {len(valid)} regions, max {int(n_max)} edges")
            for region, n in sorted(valid.items(), key=lambda x: -x[1]):
                print(f"    {region}: {int(n)}")

            plot_and_save(
                counts, n_max,
                OUTPUT_DIR / f'cluster_{cluster}_subcortical_connectivity_{direction}_yabplot.png',
                title=f'Cluster {cluster} vs NT — {label} (No. of edges, FDR < 0.05)',
                colormap=colormap,
            )

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
