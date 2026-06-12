#!/usr/bin/env python3
"""
Run the cluster-anatomical pipelines (5_cluster_anatomical_analysis legacy
and /revision/ v2) with K-means, Ward, and GMM cluster labels (curated data).

Optional flag: --autism-only — uses the autism-only clusterings produced by
running 04_run_clustering_curated.py with --autism-only (i.e. cluster file
basenames become df_clusters_complete_curated_autism{,_ward,_gmm}.csv).

For each method ∈ METHODS:
  1. Swap the sibling-repo and local cluster files into the live paths.
  2. Run legacy + v2 pipelines (writes to outputs/).
  3. Rename outputs/ → outputs_curated[_autism]_<method>/.
  4. Filter outputs to keep only **combined** figures.

Final state in each anatomical dir (when both modes have been run):
  outputs_paper/                        — paper-data k-means
  outputs_curated_{kmeans,ward,gmm}/    — curated all-pop clusters
  outputs_curated_autism_{kmeans,ward,gmm}/  — curated autism-only clusters
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--autism-only", action="store_true",
                     help="Use autism-only cluster files (df_clusters_complete_curated_autism*.csv).")
_parser.add_argument("--nt-c1", action="store_true",
                     help="Pass --nt-c1 to the anatomical pipelines: restrict the NT "
                          "control pool to NT subjects with Cluster == 'C1' (NT-like).")
_args = _parser.parse_args()
AUTISM_ONLY = _args.autism_only
NT_C1 = _args.nt_c1
# CLUSTER_TAG identifies which cluster file to load (depends on clustering mode).
# OUT_TAG identifies the output directory (also reflects the --nt-c1 analysis variant).
CLUSTER_TAG = "curated_autism" if AUTISM_ONLY else "curated"
OUT_TAG = CLUSTER_TAG + ("_nt_c1" if NT_C1 else "")
TAG = CLUSTER_TAG  # kept for backwards compat with the rest of the script
print(f"=== Mode: {'AUTISM-ONLY' if AUTISM_ONLY else 'ALL POPULATIONS'}"
      f"{', NT=C1' if NT_C1 else ''} (tag={OUT_TAG}) ===")

PY = "/usr/local/bin/python3.11"
ROOT = Path("/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource")
SIBLING = Path("/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataset_paper/dataframes")
LIVE = SIBLING / "df_clusters_complete_kmeans.csv"

# v2 anatomical reads a LOCAL file rather than the sibling — must also swap.
LOCAL_INDIV = ROOT / "1_clustering" / "outputs" / "tables" / "individuals_metrics_with_clusters.csv"
CURATED_DIR = ROOT / "1_clustering" / "outputs" / TAG / "tables"

PIPELINES = [
    #(ROOT / "5_cluster_anatomical_analysis", "scripts/run_cluster_anatomical_analysis.py"),
    (ROOT / "5_cluster_anatomical_analysis" / "revision", "scripts/run_cluster_anatomical_analysis_v2.py"),
]

METHODS = ["kmeans", "ward", "gmm"]


def swap_sibling_to(method):
    """Back up live sibling + local files and overwrite with curated <method>."""
    method_suffix = "" if method == "kmeans" else f"_{method}"
    file_tag = TAG  # "curated" or "curated_autism"
    sibling_src = SIBLING / f"df_clusters_complete_{file_tag}{method_suffix}.csv"
    sibling_backup = SIBLING / "df_clusters_complete_kmeans_paperdata.csv"
    if not sibling_src.exists():
        raise FileNotFoundError(f"Curated sibling file missing: {sibling_src}")
    if not sibling_backup.exists():
        shutil.copy2(LIVE, sibling_backup)
        print(f"  backed up: {LIVE.name} → {sibling_backup.name}")
    shutil.copy2(sibling_src, LIVE)
    print(f"  installed: {sibling_src.name} → {LIVE.name}")

    # Also swap the LOCAL file the v2 anatomical reads.
    local_src = CURATED_DIR / f"individuals_metrics_with_clusters_{file_tag}{method_suffix}.csv"
    local_backup = LOCAL_INDIV.with_name(LOCAL_INDIV.stem + "_paperdata" + LOCAL_INDIV.suffix)
    if not local_src.exists():
        raise FileNotFoundError(f"Curated local file missing: {local_src}")
    if not local_backup.exists():
        shutil.copy2(LOCAL_INDIV, local_backup)
        print(f"  backed up: {LOCAL_INDIV.name} → {local_backup.name}")
    shutil.copy2(local_src, LOCAL_INDIV)
    print(f"  installed: {local_src.name} → {LOCAL_INDIV.name}")


def restore_sibling():
    sibling_backup = SIBLING / "df_clusters_complete_kmeans_paperdata.csv"
    if sibling_backup.exists():
        shutil.copy2(sibling_backup, LIVE)
        sibling_backup.unlink()
        print(f"  restored: {LIVE.name}")
    local_backup = LOCAL_INDIV.with_name(LOCAL_INDIV.stem + "_paperdata" + LOCAL_INDIV.suffix)
    if local_backup.exists():
        shutil.copy2(local_backup, LOCAL_INDIV)
        local_backup.unlink()
        print(f"  restored: {LOCAL_INDIV.name}")


def run_one_pipeline(section_dir: Path, script_rel: str) -> bool:
    script = section_dir / script_rel
    print(f"\n  → {script.relative_to(ROOT)}")
    cmd = [PY, str(script)]
    if NT_C1:
        cmd.append("--nt-c1")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"    FAILED (exit {res.returncode}). stderr tail:\n    " + res.stderr.strip()[-400:])
        return False
    tail = "\n    ".join(res.stdout.strip().splitlines()[-2:])
    print(f"    OK\n    {tail}")
    return True


def filter_combined_only(out_dir: Path):
    """Keep combined_*.pdf, *_cohens_d_forest.{pdf,csv}, r_input_files/, and
    cluster_counts_mri.csv. Delete per-cluster cluster_C*_*.{pdf,png} files."""
    figs = out_dir / "figures"
    if not figs.exists():
        print(f"    (no figures/ under {out_dir.name})")
        return
    kept, removed = 0, 0
    for p in figs.rglob("*"):
        if not p.is_file():
            continue
        if p.parent.name == "r_input_files":
            kept += 1
            continue
        name = p.name
        if (name.startswith("combined_clusters_")
                or name.endswith("_cohens_d_forest.pdf")
                or name.endswith("_cohens_d_forest.csv")
                or name == "cluster_counts_mri.csv"):
            kept += 1
        else:
            p.unlink()
            removed += 1
    print(f"    figures filtered: kept {kept}, removed {removed}")


def main():
    print("=== Running anatomical pipelines for Ward and GMM (curated) ===")
    for method in METHODS:
        print(f"\n=========================================================")
        print(f"  METHOD: {method.upper()}")
        print(f"=========================================================")
        try:
            swap_sibling_to(method)
            for section_dir, script_rel in PIPELINES:
                out_dir = section_dir / "outputs"
                # Ensure outputs/ doesn't exist before the run (rename old if it does)
                if out_dir.exists():
                    shutil.rmtree(out_dir)
                ok = run_one_pipeline(section_dir, script_rel)
                if not ok:
                    print(f"    pipeline failed — leaving outputs/ in place for inspection")
                    continue
                target = section_dir / f"outputs_{OUT_TAG}_{method}"
                if target.exists():
                    shutil.rmtree(target)
                out_dir.rename(target)
                print(f"    moved outputs/ → {target.name}/")
                filter_combined_only(target)
        finally:
            restore_sibling()

    print("\n=========================================================")
    print("  DONE")
    print("=========================================================")
    for section_dir, _ in PIPELINES:
        print(f"\n  {section_dir.relative_to(ROOT)}/")
        for d in sorted(section_dir.glob("outputs_*")):
            n_pdfs = sum(1 for _ in d.rglob("*.pdf"))
            n_pngs = sum(1 for _ in d.rglob("*.png"))
            n_csv  = sum(1 for _ in d.rglob("*.csv"))
            print(f"    {d.name:<30}  {n_pdfs:>3} pdfs · {n_pngs:>3} pngs · {n_csv:>3} csvs")


if __name__ == "__main__":
    main()
