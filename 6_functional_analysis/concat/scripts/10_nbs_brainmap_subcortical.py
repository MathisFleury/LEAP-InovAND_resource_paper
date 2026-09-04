#!/usr/bin/env python3
"""
Whole-group Autism-vs-NT NBS subcortical insets (figure_6b style, section 6-2).

Adapted from 7_cluster_functional_analysis/scripts/08_figure6b_nbs_subcortical.py
but for the single whole-group contrast (Autism vs NT) instead of per cluster.
Reads nbs_edges.csv (region, direction, t_stat) written by 08_nbs_brainmap_edges.py.

For each direction (hypo/hyper) counts NBS-component edges incident to each
subcortical region and renders a yabplot inset. Merged thalamus (LH/RH_Thalamus
built by unioning the 7+7 nucleus meshes), cerebellum dropped, 2 views (right|left).

Input : outputs/<AUTISM_TD_OUTPUT_DIR>/autism_vs_td_<dir>connectivity_full_data.csv (see AUTISM_TD_OUTPUT_DIR)
Output: <AUTISM_TD_OUTPUT_DIR>/autism_vs_nt_nbs_subcortical_<dir>.png
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import pyvista as pv

pv.global_theme.notebook = True
from yabplot import plot_subcortical

_SCRIPT_DIR = Path(__file__).parent
_SECTION = _SCRIPT_DIR.parent
# Read the NBS largest-component edges from the nbs_double dir (04's full_data,
# source,target) — the 6-2-canonical whole-group component, consistent with the
# double network matrix (06) and the combined figure (11). Set via IO_DIR /
# AUTISM_TD_OUTPUT_DIR (defaults to the same convention as 04/06).
NBS_DIR = Path(os.environ.get("AUTISM_TD_OUTPUT_DIR",
               _SECTION / "outputs" / "figures" / os.environ.get("NBS_SUBDIR", "nbs_wholegroup")))
_FULL_ATLAS_DIR = _SECTION.parent.parent / "_resources" / "custom_subcortical_atlas_4S156"

DIRS = [("hypo", "Blues", "Hypoconnectivity"),
        ("hyper", "Reds", "Hyperconnectivity")]
VIEWS = ["right_lateral", "left_lateral"]   # two views only: right | left
EXCLUDE_PREFIX = "Cerebellar_"              # cerebellum dropped (atlas too coarse)
_THAL_NUCLEI = ["Pulvinar", "Anterior", "Medio_Dorsal", "Ventral_Latero_Dorsal",
                "Central_Lateral-Lateral_Posterior-Medial_Pulvinar",
                "Ventral_Anterior", "Ventral_Latero_Ventral"]


def make_filtered_atlas_dir() -> str:
    """Drop cerebellar meshes; replace 14 thalamic-nucleus meshes with 2 merged
    LH_Thalamus / RH_Thalamus (union of the 7 per side) to match the v0.11
    merged-thalamus connectivity."""
    dst = _SECTION.parent.parent / "_resources" / "custom_subcortical_atlas_4S156_merged"
    dst.mkdir(parents=True, exist_ok=True)
    for f in dst.glob("*.vtk"):
        f.unlink()
    thal_names = {f"{h}-{n}" for h in ("LH", "RH") for n in _THAL_NUCLEI}
    kept = 0
    for f in _FULL_ATLAS_DIR.glob("*.vtk"):
        if f.name.startswith(EXCLUDE_PREFIX) or f.stem in thal_names:
            continue
        (dst / f.name).symlink_to(f)
        kept += 1
    for hemi, out in (("LH", "LH_Thalamus"), ("RH", "RH_Thalamus")):
        parts = [_FULL_ATLAS_DIR / f"{hemi}-{n}.vtk" for n in _THAL_NUCLEI]
        parts = [p for p in parts if p.exists()]
        if not parts:
            continue
        merged = pv.MultiBlock([pv.read(str(p)) for p in parts]).combine()
        merged.save(str(dst / f"{out}.vtk"))
        kept += 1
    print(f"filtered atlas (no cerebellum, merged thalamus): {kept} meshes -> {dst.name}")
    return str(dst)


def atlas_regions(atlas_dir: str) -> set:
    return {p.stem for p in Path(atlas_dir).glob("*.vtk")}


def count_edges_per_region(src_series, tgt_series, atlas_set) -> dict:
    """Count NBS-component edges incident to each atlas region (source+target)."""
    counts = {r: 0 for r in atlas_set}
    for s, t in zip(src_series, tgt_series):
        for endpoint in (str(s), str(t)):
            if endpoint in counts:
                counts[endpoint] += 1
    return {k: (float(v) if v > 0 else float("nan")) for k, v in counts.items()}


def main() -> int:
    if not any(_FULL_ATLAS_DIR.glob("*.vtk")):
        print(f"ERROR: no VTK meshes in {_FULL_ATLAS_DIR}; run "
              "_resources/build_custom_subcortical_atlas.py first (needs mount).")
        return 1

    atlas_dir = make_filtered_atlas_dir()
    regs = atlas_regions(atlas_dir)
    print(f"Custom atlas regions (no cerebellum): {len(regs)}")

    # full_data source,target per direction (04's output), from the nbs_double dir
    per_cell, gmax = {}, 0
    for key, _, _ in DIRS:
        fp = NBS_DIR / f"autism_vs_td_{key}connectivity_full_data.csv"
        if not fp.exists():
            print(f"  {key}: {fp.name} not found — skipping"); continue
        d = pd.read_csv(fp)
        counts = count_edges_per_region(d["source"], d["target"], regs)
        per_cell[key] = counts
        vals = [v for v in counts.values() if not np.isnan(v)]
        if vals:
            gmax = max(gmax, max(vals))
    n_max = max(gmax, 2.0)
    print(f"shared subcortical edge-count max = {int(n_max)}")

    for key, counts in per_cell.items():
        cmap = dict((k, c) for k, c, _ in DIRS)[key]
        vals = [v for v in counts.values() if not np.isnan(v)]
        if not vals:
            print(f"  {key}: no subcortical edges — skipping")
            continue
        out = NBS_DIR / f"autism_vs_nt_nbs_subcortical_{key}.png"
        plotter = plot_subcortical(
            data=counts, custom_atlas_path=atlas_dir,
            cmap=cmap, vminmax=[1, n_max], display_type="object",
            views=VIEWS, layout=(1, 2), figsize=(1000, 460),
        )
        plotter.screenshot(str(out))
        plotter.close()
        print(f"  Saved: {out.name}  ({len(vals)} regions, max {int(max(vals))})")

    print(f"\nAll subcortical insets in {NBS_DIR}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
