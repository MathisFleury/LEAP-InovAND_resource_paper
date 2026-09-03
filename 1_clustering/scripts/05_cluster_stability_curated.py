#!/usr/bin/env python3
"""
Cluster stability analysis on the CURATED clustering.

Runs 05_cluster_stability.py on the curated clustering's feature set
(cluster_assignments_curated.csv → ID, IQ, SRS_tscore) for all three
methods, writing to 1_clustering/outputs/curated/. The frozen-paper run
(05_cluster_stability.py with no env) is unaffected.
"""

import os
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
STABILITY = _SCRIPT_DIR / "05_cluster_stability.py"
FEATURES_CSV = (_SCRIPT_DIR / ".." / "outputs" / "curated" / "tables"
                / "cluster_assignments_curated.csv").resolve()

if not FEATURES_CSV.exists():
    sys.exit(f"Curated features not found: {FEATURES_CSV}\n"
             "Run 04_run_clustering_curated.py first.")

env = {**os.environ,
       "STABILITY_FEATURES_CSV": str(FEATURES_CSV),
       "STABILITY_OUTPUT_SUBDIR": "curated"}
sys.exit(subprocess.run([sys.executable, str(STABILITY)], env=env).returncode)
