#!/usr/bin/env python3
"""
Run every cluster-aware downstream pipeline using the curated cluster
assignments (LEAP + INOVAND + INFOR, see 04_run_clustering_curated.py).

Strategy (preserves paper-data outputs):
  1. For each section, rename outputs/ → outputs_paper/ (idempotent).
  2. Swap the cluster files used by downstream readers:
        - 1_clustering/outputs/tables/cluster_assignments.csv      (Group B)
        - 1_clustering/outputs/tables/individuals_metrics_with_clusters.csv
        - eeg_mri-pipeline/.../df_clusters_complete_kmeans.csv    (Group A)
     Each is backed up to *_paperdata.csv and overwritten by its curated
     counterpart for the duration of the runs.
  3. Run every downstream pipeline / script.
  4. Restore each swapped file from its *_paperdata.csv backup.

After completion:
  - outputs_paper/   = original paper-data results (from before).
  - outputs/         = curated-data results (new).

The cluster file backups (*_paperdata.csv) are removed after restoration.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

PY = "/usr/local/bin/python3.11"
RSCRIPT = "Rscript"
ROOT = Path("/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource")
SIBLING = Path("/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataset_paper/dataframes")

# Files we will swap (existing → backup, curated → existing)
SWAPS = [
    {
        "live": ROOT / "1_clustering" / "outputs" / "tables" / "cluster_assignments.csv",
        "curated": ROOT / "1_clustering" / "outputs" / "curated" / "tables" / "cluster_assignments_curated.csv",
    },
    {
        "live": ROOT / "1_clustering" / "outputs" / "tables" / "individuals_metrics_with_clusters.csv",
        "curated": ROOT / "1_clustering" / "outputs" / "curated" / "tables" / "individuals_metrics_with_clusters_curated.csv",
    },
    {
        "live": SIBLING / "df_clusters_complete_kmeans.csv",
        "curated": SIBLING / "df_clusters_complete_curated.csv",
    },
]

# Section output dirs that should be moved to outputs_paper/ before re-running
SECTION_OUTPUT_DIRS = [
    ROOT / "2_genetic_analysis" / "outputs",
    ROOT / "3_cluster_genetic_analysis" / "outputs",
    ROOT / "5_cluster_anatomical_analysis" / "outputs",
    ROOT / "5_cluster_anatomical_analysis" / "curated" / "outputs",
    ROOT / "7_cluster_functional_analysis" / "outputs",
    ROOT / "9_cluster_eeg_analysis" / "outputs",
    ROOT / "10_clinical_analysis" / "outputs",
]

# Downstream commands to run, in order
COMMANDS = [
    # 2_genetic_analysis (cluster-aware only)
    [PY, str(ROOT / "2_genetic_analysis" / "scripts" / "03_carrier_freq_or_clusters.py")],
    # 3_cluster_genetic_analysis (no orchestrator — run all 14)
    *[[PY, str(p)] for p in sorted((ROOT / "3_cluster_genetic_analysis" / "scripts").glob("0?_*.py"))
      if p.name.startswith(("01_", "02_", "03_", "04_", "05_", "06_", "07_", "08_", "09_"))],
    *[[PY, str(p)] for p in sorted((ROOT / "3_cluster_genetic_analysis" / "scripts").glob("1?_*.py"))],
    # 5_cluster_anatomical (orchestrators)
    [PY, str(ROOT / "5_cluster_anatomical_analysis" / "scripts" / "run_cluster_anatomical_analysis.py")],
    [PY, str(ROOT / "5_cluster_anatomical_analysis" / "curated" / "scripts" / "run_cluster_anatomical_analysis_v2.py")],
    # 7_cluster_functional
    [PY, str(ROOT / "7_cluster_functional_analysis" / "scripts" / "run_cluster_functional_analysis.py")],
    # 9_cluster_eeg
    [PY, str(ROOT / "9_cluster_eeg_analysis" / "scripts" / "run_cluster_eeg_analysis.py")],
    # 10_clinical
    [PY, str(ROOT / "10_clinical_analysis" / "scripts" / "01_plot_psychomotor_milestones.py")],
    [PY, str(ROOT / "10_clinical_analysis" / "scripts" / "02_plot_verbal_status.py")],
    # NOTE: this used to also run a top-level revision/scripts/{05_synGO_chromEpiTF_per_cluster.py,
    # 07_geneset_provenance.py, 09_cluster_vs_cohort_sMRI.py} chain. That path never
    # existed in this repo (checked: no such folder or filenames anywhere in the
    # tree) -- the analysis it referred to was apparently never migrated in from
    # wherever this script was adapted from. Removed rather than left dangling;
    # if that gene-set/provenance/cluster-vs-cohort analysis still matters, it
    # needs to be written (or ported) fresh, not pointed at a path that never existed.
]


def move_outputs_to_paper(section_dir: Path):
    paper_dir = section_dir.parent / "outputs_paper"
    if section_dir.exists() and not paper_dir.exists():
        section_dir.rename(paper_dir)
        print(f"  moved {section_dir.relative_to(ROOT)} → outputs_paper/")
    elif paper_dir.exists():
        print(f"  outputs_paper/ already exists for {section_dir.parent.name} — leaving in place")
        if section_dir.exists():
            # Move whatever is in outputs/ to an outputs_orphan/ so we start clean
            orphan = section_dir.parent / "outputs_orphan"
            if orphan.exists():
                shutil.rmtree(orphan)
            section_dir.rename(orphan)
            print(f"  also moved current outputs/ → outputs_orphan/ (will be removed)")
            shutil.rmtree(orphan, ignore_errors=True)


def swap_files():
    print("\n=== Swapping cluster files ===")
    for s in SWAPS:
        live, cur = Path(s["live"]), Path(s["curated"])
        if not cur.exists():
            print(f"  SKIP (curated missing): {cur}")
            continue
        if not live.exists():
            print(f"  SKIP (live missing):    {live}")
            continue
        backup = live.with_name(live.stem + "_paperdata" + live.suffix)
        if not backup.exists():
            shutil.copy2(live, backup)
            print(f"  backed up: {live.name} → {backup.name}")
        shutil.copy2(cur, live)
        print(f"  installed: {cur.name} → {live.name}")


def restore_files():
    print("\n=== Restoring paper cluster files ===")
    for s in SWAPS:
        live = Path(s["live"])
        backup = live.with_name(live.stem + "_paperdata" + live.suffix)
        if backup.exists():
            shutil.copy2(backup, live)
            backup.unlink()
            print(f"  restored: {backup.name} → {live.name}  (backup removed)")
        else:
            print(f"  SKIP (no backup): {live.name}")


def run(cmd):
    print(f"\n--- {' '.join(cmd[-2:])} ---")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  FAILED (exit {res.returncode})")
        print("  STDERR tail:", res.stderr.strip()[-500:])
        return False
    # Print last few stdout lines
    tail = "\n".join(res.stdout.strip().splitlines()[-3:])
    print(f"  OK\n  {tail}")
    return True


def main():
    print("== Move paper outputs to outputs_paper/ ==")
    for d in SECTION_OUTPUT_DIRS:
        move_outputs_to_paper(d)

    swap_files()

    print("\n== Running downstream pipelines with curated cluster labels ==")
    results = []
    try:
        for cmd in COMMANDS:
            results.append((cmd[-1], run(cmd)))
    finally:
        restore_files()

    print("\n== Summary ==")
    failed = [c for c, ok in results if not ok]
    for c, ok in results:
        print(f"  [{'OK' if ok else 'FAIL'}] {c}")
    print(f"\n{len(results) - len(failed)}/{len(results)} succeeded.")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
