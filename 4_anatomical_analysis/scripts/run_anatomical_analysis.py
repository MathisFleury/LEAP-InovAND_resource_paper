#!/usr/bin/env python3
"""
Anatomical MRI Analysis Pipeline

Runs the complete structural MRI analysis:
1. Anatomical MRI Autism vs TD comparison (t-statistics, Cohen's d)
2. Structural MRI-IQ-SRS-Vineland correlations
3. Brain visualizations using R (ggseg)
"""

import os
import sys
import subprocess
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent

def run_script(script_path, description):
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")
    try:
        if script_path.suffix == '.py':
            result = subprocess.run([sys.executable, str(script_path)],
                                    capture_output=True, text=True, check=True)
        elif script_path.suffix == '.R':
            result = subprocess.run(['Rscript', str(script_path)],
                                    capture_output=True, text=True, check=True)
        else:
            raise ValueError(f"Unsupported script type: {script_path.suffix}")
        print(result.stdout)
        if result.stderr:
            print("Warnings/Errors:", result.stderr[:500])
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {description}:")
        print(f"  STDOUT: {e.stdout[:500]}")
        print(f"  STDERR: {e.stderr[:500]}")
        return False


def main():
    print("Anatomical MRI Analysis Pipeline")
    print("=" * 60)

    steps = [
        (_SCRIPT_DIR / '01_anatomical_mri_autism_td_analysis.py',
         'Anatomical MRI Autism vs TD Analysis', True),
        (_SCRIPT_DIR / '03_plot_anatomical_mri_brain_visualizations.R',
         'Brain Visualizations (ggseg)', False),
        (_SCRIPT_DIR / '02_structural_mri_iq_srs_analysis.py',
         'Structural MRI-IQ-SRS Correlations', False),
    ]

    results = {}
    for script_path, description, required in steps:
        if not script_path.exists():
            msg = f"Script not found: {script_path}"
            if required:
                print(f"ERROR: {msg}")
                return False
            else:
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
    success = main()
    sys.exit(0 if success else 1)
