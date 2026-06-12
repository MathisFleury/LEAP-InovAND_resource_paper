#!/usr/bin/env python3
"""
fMRI Functional Connectivity Analysis Pipeline

Runs the complete fMRI Autism vs TD analysis:
1. Connectivity analysis (t-tests, FDR, hyper/hypo)
2. Schaefer atlas brain visualizations (R)
"""

import os
import sys
import subprocess
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent


def run_script(script_path, description, python_executable=None):
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")
    try:
        if script_path.suffix == '.py':
            cmd = [python_executable or sys.executable, str(script_path)]
        elif script_path.suffix == '.R':
            cmd = ['Rscript', str(script_path)]
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


def main():
    print("fMRI Functional Connectivity Analysis Pipeline")
    print("=" * 60)

    steps = [
        (_SCRIPT_DIR / '01_autism_td_connectivity_analysis.py', 'Autism vs TD Connectivity Analysis', True, None),
        (_SCRIPT_DIR / '02_plot_autism_td_schaefer.R', 'Schaefer Atlas Brain Visualization', False, None),
        (_SCRIPT_DIR / '03_plot_subcortical_yabplot.py', 'Subcortical Connectivity Yabplot', False, '/usr/local/bin/python3.11'),
    ]

    results = {}
    for script_path, description, required, py_exec in steps:
        if not script_path.exists():
            msg = f"Script not found: {script_path}"
            if required:
                print(f"ERROR: {msg}")
                return False
            else:
                print(f"WARNING: {msg} (skipping)")
                continue
        success = run_script(script_path, description, py_exec)
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
    start = time.time()
    success = main()
    print(f"\nPipeline duration: {time.time() - start:.2f} seconds")
    sys.exit(0 if success else 1)
