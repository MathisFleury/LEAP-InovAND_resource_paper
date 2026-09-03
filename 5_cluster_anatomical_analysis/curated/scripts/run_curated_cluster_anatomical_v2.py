#!/usr/bin/env python3
"""
Generate the v2 cluster-anatomical MRI outputs from the CURATED clustering.

Runs run_cluster_anatomical_analysis_v2.py once per clustering method
(k-means / Ward / GMM), each time pointing the pipeline at the curated
cluster-label file via the CLUSTER_FILE env var (read by
01_generate_cluster_mri_inputs_v2.py, which prefers the curated
`population_group` phenotype column). No file swapping.

The paper outputs/ are preserved: stashed aside during the run and restored
afterwards. Curated results land in dedicated folders:

  curated/outputs_curated_kmeans/   curated/outputs_curated_ward/
  curated/outputs_curated_gmm/      curated/outputs/  <- paper, untouched
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

PY = sys.executable
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent          # 5_cluster_anatomical_analysis/curated/
PIPELINE = _SCRIPT_DIR / "run_cluster_anatomical_analysis_v2.py"

CURATED_DIR = Path(
    "/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/"
    "1_clustering/outputs/curated/tables"
)
_ALL_METHOD_FILES = {
    "kmeans": CURATED_DIR / "individuals_metrics_with_clusters_curated.csv",
    "ward":   CURATED_DIR / "individuals_metrics_with_clusters_curated_ward.csv",
    "gmm":    CURATED_DIR / "individuals_metrics_with_clusters_curated_gmm.csv",
}
# CLUSTER_METHODS restricts which methods to run (e.g. "kmeans" for a sensitivity
# variant). OUT_SUFFIX is appended to outputs_curated_<method> so a variant (e.g.
# the no-Euler MRI table via MRI_NDD_FILE) lands in outputs_curated_<method>_noeuler
# without clobbering the main results. MRI_NDD_FILE flows through os.environ.
_WANT = [m.strip() for m in os.environ.get("CLUSTER_METHODS", "kmeans,ward,gmm").split(",") if m.strip()]
METHOD_FILES = {m: _ALL_METHOD_FILES[m] for m in _WANT if m in _ALL_METHOD_FILES}
OUT_SUFFIX = os.environ.get("OUT_SUFFIX", "")


def main() -> bool:
    out_dir = _SECTION_DIR / "outputs"
    stash = _SECTION_DIR / ".outputs_paper_stash"

    # Stash the paper outputs/ so the curated runs never erase them.
    had_paper = out_dir.exists()
    if had_paper:
        if stash.exists():
            shutil.rmtree(stash)
        out_dir.rename(stash)

    ok_all = True
    try:
        for method, cluster_file in METHOD_FILES.items():
            print(f"\n{'=' * 60}\nCURATED {method.upper()}\n{'=' * 60}")
            if not cluster_file.exists():
                print(f"  SKIP — curated file missing: {cluster_file}")
                ok_all = False
                continue

            if out_dir.exists():
                shutil.rmtree(out_dir)

            env = {**os.environ, "CLUSTER_FILE": str(cluster_file)}
            res = subprocess.run([PY, str(PIPELINE)], env=env)
            # Judge success on the core per-cluster output, not the exit code
            # (optional plotting steps can fail without invalidating the stats).
            core = (out_dir / "figures" / "r_input_files"
                    / "t_stat_cluster_C1_thickness_dk_mri_cluster_vs_td.csv")
            if not core.exists():
                print(f"  FAILED (exit {res.returncode}, no core output) — skipping {method}")
                ok_all = False
                continue
            if res.returncode != 0:
                print(f"  NOTE: exit {res.returncode} — optional steps failed; core stats complete.")

            target = _SECTION_DIR / f"outputs_curated_{method}{OUT_SUFFIX}"
            if target.exists():
                shutil.rmtree(target)
            out_dir.rename(target)
            print(f"  curated results -> {target.name}/")
    finally:
        # Restore the paper outputs/ regardless of outcome.
        if out_dir.exists():
            shutil.rmtree(out_dir)
        if had_paper:
            stash.rename(out_dir)
            print("\npaper results restored -> outputs/")

    print("\nDone.")
    return ok_all


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
