#!/usr/bin/env python3
"""
Generate the cluster-anatomical MRI outputs from the curated clustering.

Runs run_cluster_anatomical_analysis.py once per method in CLUSTER_METHODS
(default: kmeans, the only method kept in this tree — set
CLUSTER_METHODS=kmeans,ward,gmm to reproduce the full sensitivity sweep,
whose ward/gmm results belong in ../legacy/), each time pointing the
pipeline at the curated cluster-label file via the CLUSTER_FILE env var
(read by 1_generate_cluster_mri_inputs.py). No file swapping.

kmeans with no OUT_SUFFIX lands directly in outputs/ (the canonical
result, same name the pipeline already stages into). Any other
method/suffix lands in outputs_curated_<method><suffix>/ instead, so it
never overwrites outputs/. kmeans always runs last so its result isn't
wiped by another method's staging cleanup.

Note: this duplicates the curated-method sensitivity check that
1_clustering/scripts/run_anatomical_with_other_methods.py also runs
(file-swap based) for this same pipeline. Kept for its extra
CLUSTER_METHODS/OUT_SUFFIX capabilities (method subsetting, no-Euler
sensitivity variant) that the generic driver doesn't have.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

PY = sys.executable
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent          # 5_cluster_anatomical_analysis/
PIPELINE = _SCRIPT_DIR / "run_cluster_anatomical_analysis.py"

CURATED_DIR = _SECTION_DIR.parent / "1_clustering" / "outputs" / "tables"
_ALL_METHOD_FILES = {
    "kmeans": CURATED_DIR / "individuals_metrics_with_clusters_curated.csv",
    "ward":   CURATED_DIR / "individuals_metrics_with_clusters_curated_ward.csv",
    "gmm":    CURATED_DIR / "individuals_metrics_with_clusters_curated_gmm.csv",
}
# CLUSTER_METHODS restricts which methods to run. OUT_SUFFIX is appended to
# outputs_curated_<method> so a variant (e.g. the no-Euler MRI table via
# MRI_NDD_FILE) lands in outputs_curated_<method>_noeuler without clobbering
# the main results. MRI_NDD_FILE flows through os.environ.
_WANT = [m.strip() for m in os.environ.get("CLUSTER_METHODS", "kmeans").split(",") if m.strip()]
OUT_SUFFIX = os.environ.get("OUT_SUFFIX", "")
# kmeans last: its target is the same dir every method stages into, so an
# earlier position would get wiped by the next method's pre-run cleanup.
_ORDER = sorted((m for m in _WANT if m in _ALL_METHOD_FILES), key=lambda m: m == "kmeans")


def target_dir(method):
    if method == "kmeans" and not OUT_SUFFIX:
        return _SECTION_DIR / "outputs"
    return _SECTION_DIR / f"outputs_curated_{method}{OUT_SUFFIX}"


def main() -> bool:
    out_dir = _SECTION_DIR / "outputs"
    ok_all = True
    for method in _ORDER:
        cluster_file = _ALL_METHOD_FILES[method]
        print(f"\n{'=' * 60}\nCURATED {method.upper()}\n{'=' * 60}")
        if not cluster_file.exists():
            print(f"  SKIP — curated file missing: {cluster_file}")
            ok_all = False
            continue

        target = target_dir(method)
        if out_dir.exists() and out_dir != target:
            shutil.rmtree(out_dir)

        env = {**os.environ, "CLUSTER_FILE": str(cluster_file)}
        res = subprocess.run([PY, str(PIPELINE)], env=env)
        # Judge success on the core per-cluster output, not the exit code
        # (optional plotting steps can fail without invalidating the stats).
        core = (out_dir / "tables" / "r_input_files"
                / "t_stat_cluster_C1_thickness_dk_mri_cluster_vs_td.csv")
        if not core.exists():
            print(f"  FAILED (exit {res.returncode}, no core output) — skipping {method}")
            ok_all = False
            continue
        if res.returncode != 0:
            print(f"  NOTE: exit {res.returncode} — optional steps failed; core stats complete.")

        if target != out_dir:
            if target.exists():
                shutil.rmtree(target)
            out_dir.rename(target)
        print(f"  curated results -> {target.name}/")

    print("\nDone.")
    return ok_all


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
