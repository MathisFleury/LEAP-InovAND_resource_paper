#!/usr/bin/env python3
"""
Run the cluster-genetic pipelines (3_cluster_genetic_analysis, 14 scripts;
plus the cluster-aware 2_genetic_analysis/03_carrier_freq_or_clusters.py)
for each clustering method ∈ {kmeans, ward, gmm} from the curated data.

Optional --autism-only flag uses the autism-only clusterings produced by
04_run_clustering_curated.py --autism-only.

For each method:
  1. Swap `1_clustering/outputs/tables/cluster_assignments.csv` to the
     curated version for that method (the local file is what these scripts
     read).  Also swaps the sibling-repo df_clusters_complete_kmeans.csv
     (used by 2_genetic_analysis/03 only).
  2. Run the 15 cluster-aware genetic scripts.
  3. Move outputs:
        3_cluster_genetic_analysis/outputs → outputs_<tag>_<method>/
        2_genetic_analysis/outputs/figures → outputs_<tag>_<method>/figures
        (the upstream carrier_annotations.tsv stays in 2_genetic/outputs/tables/)
  4. Restore live files on exit (try/finally).

Final state after running both modes (--autism-only and default):
  3_cluster_genetic_analysis/
    outputs_paper/                              # original paper-data k-means
    outputs_curated_{kmeans,ward,gmm}/          # curated all-pop clusters
    outputs_curated_autism_{kmeans,ward,gmm}/   # curated autism-only clusters
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--autism-only", action="store_true",
                     help="Use autism-only cluster files.")
_args = _parser.parse_args()
AUTISM_ONLY = _args.autism_only
TAG = "curated_autism" if AUTISM_ONLY else "curated"
print(f"=== Mode: {'AUTISM-ONLY' if AUTISM_ONLY else 'ALL POPULATIONS'} (tag={TAG}) ===")

PY = "/usr/local/bin/python3.11"
ROOT = Path("/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource")
SIBLING = Path("/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataset_paper/dataframes")

LIVE_SIBLING = SIBLING / "df_clusters_complete_kmeans.csv"
LIVE_LOCAL_ASSIGN = ROOT / "1_clustering" / "outputs" / "tables" / "cluster_assignments.csv"
LIVE_LOCAL_INDIV  = ROOT / "1_clustering" / "outputs" / "tables" / "individuals_metrics_with_clusters.csv"
CURATED_DIR = ROOT / "1_clustering" / "outputs" / TAG / "tables"

METHODS = ["kmeans", "ward", "gmm"]

# Scripts that actually read our cluster_assignments.csv (and thus produce
# different outputs under different clustering methods). The other scripts
# in 3_cluster_genetic_analysis use external clusterings (GMM from
# large_clustering_datasets, manual cutoffs, SPARK-LOF) and produce
# identical outputs regardless of which method we ran — they're not
# re-executed here.
GENETIC_SCRIPTS = [
    ROOT / "2_genetic_analysis/scripts/03_carrier_freq_or_clusters.py",
    ROOT / "3_cluster_genetic_analysis/scripts/01_legacy_clusters.py",
    ROOT / "3_cluster_genetic_analysis/scripts/03_legacy_clusters_nt_c1.py",
    ROOT / "3_cluster_genetic_analysis/scripts/05_hg38_legacy_clusters.py",
    ROOT / "3_cluster_genetic_analysis/scripts/07_hg38_legacy_clusters_nt_c1.py",
]

# Subdir renames applied to each method's output dir after the scripts run.
LEGACY_SUFFIXES = ["", "_nt_c1", "_hg38", "_hg38_nt_c1"]


def rename_legacy_subdirs(out_dir: Path, method: str):
    """Rename figures_legacy{,_nt_c1,_hg38,_hg38_nt_c1} → figures_method_<method>{,_nt_c1,...}.
    Same for tables_*. Makes it obvious which clustering method produced the files."""
    if not out_dir.exists():
        return
    for kind in ("figures", "tables"):
        for suf in LEGACY_SUFFIXES:
            src = out_dir / f"{kind}_legacy{suf}"
            dst = out_dir / f"{kind}_method_{method}{suf}"
            if src.exists():
                if dst.exists():
                    shutil.rmtree(dst)
                src.rename(dst)
                print(f"    renamed: {src.name}/ → {dst.name}/")

# Output dirs to relocate to outputs_<tag>_<method>/ after each method run.
SECTION_OUTPUTS = [
    ROOT / "3_cluster_genetic_analysis" / "outputs",
    # 2_genetic_analysis/outputs is mostly upstream (carrier annotations);
    # we won't move the whole dir — instead, after each method run we copy
    # only the *_cluster_* files into a per-method subdir below.
]


def _backup_path(p: Path) -> Path:
    return p.with_name(p.stem + "_paperdata" + p.suffix)


def swap_to(method: str):
    method_suffix = "" if method == "kmeans" else f"_{method}"
    file_tag = TAG

    swaps = [
        (LIVE_SIBLING,
         SIBLING / f"df_clusters_complete_{file_tag}{method_suffix}.csv"),
        (LIVE_LOCAL_ASSIGN,
         CURATED_DIR / f"cluster_assignments_{file_tag}{method_suffix}.csv"),
        (LIVE_LOCAL_INDIV,
         CURATED_DIR / f"individuals_metrics_with_clusters_{file_tag}{method_suffix}.csv"),
    ]
    for live, src in swaps:
        if not src.exists():
            raise FileNotFoundError(f"Source missing: {src}")
        bak = _backup_path(live)
        if not bak.exists():
            shutil.copy2(live, bak)
            print(f"  backed up: {live.name}")
        shutil.copy2(src, live)
        print(f"  installed: {src.name} → {live.name}")


def restore_all():
    for live in (LIVE_SIBLING, LIVE_LOCAL_ASSIGN, LIVE_LOCAL_INDIV):
        bak = _backup_path(live)
        if bak.exists():
            shutil.copy2(bak, live)
            bak.unlink()
            print(f"  restored: {live.name}")


def run_script(script: Path) -> bool:
    print(f"  → {script.name}")
    res = subprocess.run([PY, str(script)], capture_output=True, text=True)
    if res.returncode != 0:
        print(f"    FAILED. stderr tail: {res.stderr.strip()[-300:]}")
        return False
    tail = res.stdout.strip().splitlines()[-1:] if res.stdout.strip() else []
    print(f"    OK  {tail[0] if tail else ''}")
    return True


def collect_2genetic_cluster_outputs(target_dir: Path, method: str):
    """Move CLUSTER_* files from 2_genetic_analysis/outputs/{figures,tables}
    into target_dir/{figures,tables}_method_<method>/ — same subdir as the
    rest of the per-method cluster artifacts (PGS forest, carrier OR plots,
    etc.). Removes the source files so they don't bleed into the next
    method's run."""
    src = ROOT / "2_genetic_analysis" / "outputs"
    if not src.exists():
        return
    moved = 0
    for fig in (src / "figures").glob("CLUSTER_*"):
        dest = target_dir / f"figures_method_{method}"
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fig, dest / fig.name)
        fig.unlink()
        moved += 1
    for tab in (src / "tables").glob("CLUSTER_*"):
        dest = target_dir / f"tables_method_{method}"
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tab, dest / tab.name)
        tab.unlink()
        moved += 1
    if moved:
        print(f"    captured {moved} 2_genetic CLUSTER_* files → figures/tables_method_{method}/")


def main():
    results = []
    try:
        for method in METHODS:
            print(f"\n=========================================================")
            print(f"  METHOD: {method.upper()}  ({TAG})")
            print(f"=========================================================")
            swap_to(method)
            ok_count = 0
            for script in GENETIC_SCRIPTS:
                ok = run_script(script)
                ok_count += int(ok)
                results.append((method, script.name, ok))
            print(f"\n  {ok_count}/{len(GENETIC_SCRIPTS)} scripts succeeded")

            # Move 3_cluster_genetic_analysis/outputs → outputs_<tag>_<method>/
            section3 = ROOT / "3_cluster_genetic_analysis"
            src = section3 / "outputs"
            tgt = section3 / f"outputs_{TAG}_{method}"
            if tgt.exists():
                shutil.rmtree(tgt)
            if src.exists():
                src.rename(tgt)
                print(f"  moved 3_cluster_genetic/outputs → {tgt.name}/")
            # Rename figures_legacy* → figures_method_<method>* (and tables_*)
            rename_legacy_subdirs(tgt, method)
            collect_2genetic_cluster_outputs(tgt, method)
    finally:
        print("\n--- Restoring live files ---")
        restore_all()

    print("\n=========================================================")
    print("  DONE")
    print("=========================================================")
    failed = [r for r in results if not r[2]]
    print(f"\n  {len(results) - len(failed)}/{len(results)} scripts succeeded across all methods.")
    if failed:
        for m, n, _ in failed:
            print(f"    [{m}] {n}")
    section3 = ROOT / "3_cluster_genetic_analysis"
    for d in sorted(section3.glob("outputs_*")):
        n_pdfs = sum(1 for _ in d.rglob("*.pdf"))
        n_csv  = sum(1 for _ in d.rglob("*.csv"))
        print(f"  {d.name:<40}  {n_pdfs:>4} pdfs · {n_csv:>4} csvs")


if __name__ == "__main__":
    sys.exit(0 if main() is None else 1)
