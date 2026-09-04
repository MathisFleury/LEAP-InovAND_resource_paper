#!/usr/bin/env python3
"""
Cluster Anatomical MRI Analysis Pipeline — curated (v2)

Runs the v2 per-cluster anatomical pipeline rebuilt on the QC + ComBat +
age/sex/eTIV-regressed FreeSurfer table, with LEAP_W1 / INOVAND_T1 wave
selection (see 01_generate_cluster_mri_inputs_v2.py for details).

Steps:
  1. Per-cluster MRI inputs — CT / SA / subcortical contrasts vs pooled
     NT + pairwise cluster comparisons (t / Cohen's d / FDR).
  2. Cortical brain maps — ggseg (R), fill = Cohen's d, FDR outlines.
  3. Subcortical brain map — yabplot, Cohen's d.
  4. Figure 6a composite grid (ggseg t-stat maps, 3 modalities x 3 clusters).
  5. Global structural measures by cluster (violin: thickness/eTIV/area).
  6. Publication gt tables (per-cluster MRI stats).
"""

import argparse
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--nt-c1", action="store_true",
    help="Pass --nt-c1 to step 1 (restrict NT pool to NT subjects with Cluster == 'C1').",
)
ARGS, _ = _parser.parse_known_args()


def run_script(script_path: Path, description: str) -> bool:
    print(f"\n{'=' * 60}\nRunning: {description}\n{'=' * 60}")
    try:
        if script_path.suffix == ".py":
            cmd = [sys.executable, str(script_path)]
            if ARGS.nt_c1 and script_path.name.startswith("01_"):
                cmd.append("--nt-c1")
        elif script_path.suffix == ".R":
            cmd = ["Rscript", str(script_path)]
        else:
            raise ValueError(f"Unsupported script type: {script_path.suffix}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stdout)
        if result.stderr:
            print("Warnings/Errors:", result.stderr[:500])
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {description}:")
        print(f"  STDOUT: {e.stdout[:500]}")
        print(f"  STDERR: {e.stderr[:500]}")
        return False


def main() -> bool:
    print("Cluster Anatomical MRI Analysis Pipeline — curated (v2)")
    print("=" * 60)

    steps = [
        (_SCRIPT_DIR / "01_generate_cluster_mri_inputs_v2.py",
         "Per-cluster MRI inputs — Cohen's d / FDR", True),
        (_SCRIPT_DIR / "02_plot_cluster_mri_brain_v2.R",
         "Cortical brain maps (ggseg, Cohen's d)", False),
        (_SCRIPT_DIR / "04_combined_cluster_brain_figure_v2.R",
         "Combined single-page cluster brain figure (per cluster, FDR)", False),
        (_SCRIPT_DIR / "03_plot_subcortical_cluster_yabplot_v2.py",
         "Subcortical brain map (yabplot, Cohen's d)", False),
        (_SCRIPT_DIR / "05_figure6a_anatomical_v2.R",
         "Figure 6a composite grid (3 modalities x 3 clusters)", False),
        (_SCRIPT_DIR / "06_cluster_global_measures_violin_v2.py",
         "Global structural measures by cluster (violin)", False),
        (_SCRIPT_DIR / "07_gt_tables_cluster_v2.R",
         "Publication gt tables (per-cluster MRI stats)", False),
    ]

    results = {}
    for script_path, description, required in steps:
        if not script_path.exists():
            msg = f"Script not found: {script_path}"
            if required:
                print(f"ERROR: {msg}")
                return False
            print(f"WARNING: {msg} (skipping)")
            continue
        success = run_script(script_path, description)
        results[description] = success
        if not success and required:
            print(f"Required step failed: {description}")
            return False

    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    for step, success in results.items():
        status = "SUCCESS" if success else "FAILED"
        print(f"  [{status}] {step}")

    return all(results.values())


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
