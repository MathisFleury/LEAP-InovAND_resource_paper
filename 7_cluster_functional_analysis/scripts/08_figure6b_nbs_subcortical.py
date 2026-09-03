#!/usr/bin/env python3
"""
Figure 6b — subcortical yabplot insets from the per-cluster NBS components.

Subcortical analogue / companion of 07_figure6b_nbs_cortical.R. For each autism
cluster (C1/C2/C3) and each direction (hypo = cluster<NT, hyper = cluster>NT),
counts NBS-component edges incident to each custom-atlas subcortical region
(CIT168 striatum/basal ganglia, HCP hippocampus/amygdala, 10 cerebellar
parcels, HCP thalamic nuclei) and renders a yabplot subcortical figure.

Input : outputs/figures/nbs/cluster_<C>_nbs_edges.csv  (region, direction, t_stat)
        (written by 06_cluster_nbs.py)
Output: outputs/figures/nbs/cluster_<C>_nbs_subcortical_<dir>.png

Cerebellum: EXCLUDED here. The 4S156 atlas carves the whole cerebellum into
only 10 coarse parcels (vs 100 cortical), too messy for this panel — the
cerebellar meshes are dropped so yabplot neither colours nor draws them.

Run _resources/build_custom_subcortical_atlas.py once first (needs the mount).
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
NBS_DIR = _SECTION / "outputs" / "figures" / os.environ.get("NBS_SUBDIR", "nbs")
_FULL_ATLAS_DIR = _SECTION.parent / "_resources" / "custom_subcortical_atlas_4S156"

CLUSTERS = ["C1", "C2", "C3"]
DIRS = [("hypo", "Blues", "Hypoconnectivity"),
        ("hyper", "Reds", "Hyperconnectivity")]
VIEWS = ["right_lateral", "left_lateral"]   # two views only: right | left
EXCLUDE_PREFIX = "Cerebellar_"              # cerebellum dropped (atlas too coarse)
# 14 thalamic nuclei meshes → merged into LH_Thalamus / RH_Thalamus (the v0.11
# connectivity uses the 2 merged thalamic ROIs; the dseg has no voxels for them,
# so we union the nucleus meshes to draw them).
_THAL_NUCLEI = ["Pulvinar", "Anterior", "Medio_Dorsal", "Ventral_Latero_Dorsal",
                "Central_Lateral-Lateral_Posterior-Medial_Pulvinar",
                "Ventral_Anterior", "Ventral_Latero_Ventral"]


def make_filtered_atlas_dir() -> str:
    """Build a sibling atlas dir that (a) drops cerebellar meshes and (b) replaces
    the 14 thalamic-nucleus meshes with 2 merged LH_Thalamus / RH_Thalamus meshes
    (union of the 7 per side), matching the merged-thalamus v0.11 connectivity."""
    import pyvista as _pv
    dst = _SECTION.parent / "_resources" / "custom_subcortical_atlas_4S156_merged"
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
    # merged thalamus meshes (union of the per-hemisphere nucleus meshes)
    for hemi, out in (("LH", "LH_Thalamus"), ("RH", "RH_Thalamus")):
        parts = [_FULL_ATLAS_DIR / f"{hemi}-{n}.vtk" for n in _THAL_NUCLEI]
        parts = [p for p in parts if p.exists()]
        if not parts:
            continue
        merged = _pv.MultiBlock([_pv.read(str(p)) for p in parts]).combine()
        merged.save(str(dst / f"{out}.vtk"))
        kept += 1
    print(f"filtered atlas (no cerebellum, merged thalamus): {kept} meshes -> {dst.name}")
    return str(dst)


def atlas_regions(atlas_dir: str) -> set:
    return {p.stem for p in Path(atlas_dir).glob("*.vtk")}


def count_edges_per_region(regions_series, atlas_set) -> dict:
    counts = {r: 0 for r in atlas_set}
    for region in regions_series:
        clean = region.removeprefix("con_")
        parts = clean.split("/", 1)
        for endpoint in (parts[0], parts[1] if len(parts) > 1 else ""):
            if endpoint in counts:
                counts[endpoint] += 1
    # NaN for zero so those regions render as the grey context mesh
    return {k: (float(v) if v > 0 else float("nan")) for k, v in counts.items()}


def main() -> int:
    if not any(_FULL_ATLAS_DIR.glob("*.vtk")):
        print(f"ERROR: no VTK meshes in {_FULL_ATLAS_DIR}; run "
              "_resources/build_custom_subcortical_atlas.py first (needs mount).")
        return 1
    atlas_dir = make_filtered_atlas_dir()   # cerebellum-free
    regs = atlas_regions(atlas_dir)
    print(f"Custom atlas regions (no cerebellum): {len(regs)}")

    # Shared max across every cluster x direction so all insets share a scale.
    per_cell = {}
    gmax = 0
    for C in CLUSTERS:
        fp = NBS_DIR / f"cluster_{C}_nbs_edges.csv"
        if not fp.exists():
            print(f"  missing {fp.name}; skipping {C}")
            continue
        d = pd.read_csv(fp)
        for key, _, _ in DIRS:
            # regs excludes cerebellum, so cerebellar edges are dropped here too
            counts = count_edges_per_region(d[d.direction == key].region, regs)
            per_cell[(C, key)] = counts
            vals = [v for v in counts.values() if not np.isnan(v)]
            if vals:
                gmax = max(gmax, max(vals))
    n_max = max(gmax, 2.0)
    print(f"shared subcortical edge-count max = {int(n_max)}")

    for (C, key), counts in per_cell.items():
        cmap = dict((k, c) for k, c, _ in DIRS)[key]
        label = dict((k, lb) for k, _, lb in DIRS)[key]
        vals = [v for v in counts.values() if not np.isnan(v)]
        if not vals:
            print(f"  {C} {key}: no subcortical edges — skipping")
            continue
        out = NBS_DIR / f"cluster_{C}_nbs_subcortical_{key}.png"
        # No overlaid title (the composite labels each inset). Two views only.
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
