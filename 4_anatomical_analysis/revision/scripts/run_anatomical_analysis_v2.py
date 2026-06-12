#!/usr/bin/env python3
"""
Anatomical MRI Analysis Pipeline — REVISION (v2)

Runs the v2 anatomical pipeline rebuilt on the QC + ComBat +
age/sex/eTIV-regressed FreeSurfer table, with LEAP_W1 / INOVAND_T1
wave selection (see 01_anatomical_mri_autism_nt_v2.py for details).

Steps:
  1. Autism vs NT contrasts — t / Cohen's d / FDR per region.
  2. Cortical brain maps — ggseg (R), fill = Cohen's d, FDR outlines.
  3. Subcortical brain map — yabplot, Cohen's d.
  4. Euler number by population — surface-reconstruction QC summary.
  5. LOEUF (hg38) × MRI Pearson correlations — all-carriers only, per-panel FDR.
  6. LOEUF (hg38) brain maps — ggseg cortical + aseg subcortical; raw p<0.01 outlines.
"""

import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent


def run_script(script_path: Path, description: str) -> bool:
    print(f"\n{'=' * 60}\nRunning: {description}\n{'=' * 60}")
    try:
        if script_path.suffix == ".py":
            cmd = [sys.executable, str(script_path)]
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
    print("Anatomical MRI Analysis Pipeline — REVISION (v2)")
    print("=" * 60)

    steps = [
        (_SCRIPT_DIR / "01_anatomical_mri_autism_nt_v2.py",
         "Autism vs NT — t / Cohen's d / FDR", True),
        (_SCRIPT_DIR / "02_plot_anatomical_mri_brain_v2.R",
         "Cortical brain maps (ggseg, Cohen's d)", False),
        (_SCRIPT_DIR / "03_plot_subcortical_yabplot_v2.py",
         "Subcortical brain map (yabplot, Cohen's d)", False),
        (_SCRIPT_DIR / "04_euler_by_population.py",
         "Euler number by population (QC summary)", False),
        (_SCRIPT_DIR / "05_loeuf_mri_correlations_hg38_v2.py",
         "LOEUF (hg38) × MRI Pearson correlations (per-panel FDR)", False),
        (_SCRIPT_DIR / "06_loeuf_mri_brain_maps_hg38_v2.R",
         "LOEUF (hg38) brain maps (cortical + subcortical, p<0.01)", False),
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
