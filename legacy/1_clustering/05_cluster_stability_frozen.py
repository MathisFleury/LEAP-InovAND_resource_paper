#!/usr/bin/env python3
# =============================================================================
# 05 - Cluster stability, FROZEN / paper-reproduction
# =============================================================================
# Thin invoker for 1_clustering/scripts/05_cluster_stability.py's shared
# engine, in its original frozen mode: derives the same ready-CSV shape
# (ID, IQ, SRS_tscore) the engine expects from the DEPRECATED
# individuals_metrics.tsv (see root CLAUDE.md), then subprocess-calls the
# engine. Curated is the default/only mode shown in 1_clustering/scripts/
# itself; this wrapper (plus the ID-exclusion list, which was frozen-data-
# specific) is what keeps the deprecated-data reference out of that main tree.
#
# Run: cd legacy/1_clustering && python3.11 05_cluster_stability_frozen.py
# Output: 1_clustering/outputs/ (same place the frozen run always wrote to).
# =============================================================================
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_THIS_DIR = Path(__file__).parent.resolve()
_REPO_ROOT = _THIS_DIR.parent.parent
ENGINE = _REPO_ROOT / "1_clustering" / "scripts" / "05_cluster_stability.py"

DATA_PATH = os.environ.get(
    "LEAP_INOVAND_DATA", str(_REPO_ROOT.parent / "imaging2genet" / "0_input" / "dataframes")
)
INDIVIDUALS_METRICS = os.path.join(DATA_PATH, "individuals_metrics.tsv")

# Data-quality exclusions specific to the frozen dataset (never applied to
# curated, which has its own QC upstream).
EXCLUDED_IDS = {
    "429385763020",
    "C0733-011-137-001",
    "C0733-011-155-001",
}


def build_ready_csv() -> Path:
    """Derive the (ID, IQ, SRS_tscore) ready-CSV shape the engine expects,
    from the raw frozen TSV -- same logic the engine used to run inline."""
    df = pd.read_csv(INDIVIDUALS_METRICS, sep="\t", low_memory=False)
    df = df.drop_duplicates(subset=["ID"]).replace({999: np.nan, 998: np.nan})
    df = df[df["Relation_to_proposant"] == "participant"]
    df["IQ"] = df["total_IQ"].fillna(df["performance_IQ"])
    df_clust = df[["IQ", "SRS_tscore", "ID"]].dropna().reset_index(drop=True)
    df_clust = df_clust[~df_clust["ID"].astype(str).isin(EXCLUDED_IDS)].reset_index(drop=True)

    out_dir = _THIS_DIR / "outputs" / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "frozen_features_for_stability.csv"
    df_clust.to_csv(out, index=False)
    print(f"Built frozen ready CSV: {out}  (N={len(df_clust)})")
    return out


def main():
    ready_csv = build_ready_csv()
    env = {**os.environ, "STABILITY_FEATURES_CSV": str(ready_csv), "STABILITY_OUTPUT_SUBDIR": ""}
    sys.exit(subprocess.run([sys.executable, str(ENGINE)], env=env).returncode)


if __name__ == "__main__":
    main()
