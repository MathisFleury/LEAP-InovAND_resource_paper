#!/usr/bin/env python3
"""
Run the gnomad v4 cluster genetic analysis (15_v4_cluster_clusters.py, all 5
clustering strategies) once per CURATED clustering solution -- legacy.

Retired: 15_v4_cluster_clusters.py's k-means branch is now covered by the
current tree's 1_carrier_freq_or_clusters.py, whose curated-method
comparison is already handled by 1_clustering/scripts/
run_genetics_with_other_methods.py (file-swap based). This script is kept
for reproducing the pre-consolidation legacy runs across all 5 strategies.

Genetic analogue of
  5_cluster_anatomical_analysis/curated/scripts/run_curated_cluster_anatomical_v2.py
: re-runs 15_v4_cluster_clusters.py for each curated re-clustering method
(k-means / Ward / GMM), pointing it at that method's [ID, Cluster] label file
via env vars. Each method's constrained carrier-freq + OR figures/tables land in
  legacy/3_cluster_genetic_analysis/outputs_curated_<method>/figures_<method>_v4/
  legacy/3_cluster_genetic_analysis/outputs_curated_<method>/tables_<method>_v4/
Uses the gnomad v4 carrier matrix (CARRIER_V4) and never touches the frozen
paper-label outputs produced by `15_v4_cluster_clusters.py` standalone.
"""

import os
import subprocess
import sys
from pathlib import Path

PY = os.environ.get("LEAP_INOVAND_PYTHON", "python3.11")
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent                       # legacy/3_cluster_genetic_analysis/
_ROOT = _SECTION_DIR.parent.parent                      # LEAP-InovAND_resource/
ANALYSIS = _SCRIPT_DIR / "15_v4_cluster_clusters.py"

# Curated re-clustering label files (1_clustering/scripts/4_run_clustering.py)
CURATED_DIR = _ROOT / "1_clustering" / "outputs" / "tables"
METHOD_FILES = {
    "kmeans": CURATED_DIR / "individuals_metrics_with_clusters_curated.csv",
    "ward":   CURATED_DIR / "individuals_metrics_with_clusters_curated_ward.csv",
    "gmm":    CURATED_DIR / "individuals_metrics_with_clusters_curated_gmm.csv",
}


def main() -> bool:
    ok_all = True
    for method, cluster_file in METHOD_FILES.items():
        print(f"\n{'=' * 60}\nCURATED {method.upper()} (v4 genetic)\n{'=' * 60}")
        if not cluster_file.exists():
            print(f"  SKIP — curated file missing: {cluster_file}")
            ok_all = False
            continue

        env = {
            **os.environ,
            "GENETIC_CLUSTER_FILE":  str(cluster_file),
            "GENETIC_CLUSTER_LABEL": method,
            "GENETIC_CLUSTER_COL":   "Cluster",
            "GENETIC_OUTPUT_BASE":   str(_SECTION_DIR / f"outputs_curated_{method}"),
        }
        res = subprocess.run([PY, str(ANALYSIS)], env=env)
        if res.returncode != 0:
            print(f"  FAILED (exit {res.returncode})")
            ok_all = False

    print("\nDone." if ok_all else "\nDone (with failures).")
    return ok_all


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
