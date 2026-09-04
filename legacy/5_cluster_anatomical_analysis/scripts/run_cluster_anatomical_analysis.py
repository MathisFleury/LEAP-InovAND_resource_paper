#!/usr/bin/env python3
"""
Cluster Anatomical Analysis Pipeline

Runs the complete cluster-based structural MRI analysis:
  1. Generate cluster MRI input CSVs (t-statistics, Cohen's d)
  2. Cortical brain visualizations using R (ggseg)
  3. Subcortical brain visualization (yabplot)
  4. Combined single-page brain figure (per cluster, FDR)
"""

import argparse
import sys
import subprocess
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--nt-c1", action="store_true",
    help="Pass --nt-c1 to step 1 (restrict NT pool to NT subjects with Cluster == 'C1').",
)
ARGS, _ = _parser.parse_known_args()


def run_script(script_path, description, python_executable=None):
    """
    Run a Python or R script as a subprocess.

    Parameters
    ----------
    script_path : Path
        Absolute path to the script to execute.
    description : str
        Human-readable label printed in the header.
    python_executable : str or None
        Explicit Python interpreter path for .py scripts.
        Defaults to sys.executable when None.

    Returns
    -------
    bool
        True on success, False on failure.
    """
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")

    script_path = Path(script_path)
    if not script_path.exists():
        print(f"ERROR: Script not found: {script_path}")
        return False

    try:
        if script_path.suffix == '.py':
            interpreter = python_executable if python_executable else sys.executable
            cmd = [interpreter, str(script_path)]
            if ARGS.nt_c1 and script_path.name.startswith('01_'):
                cmd.append('--nt-c1')
        elif script_path.suffix == '.R':
            cmd = ['Rscript', str(script_path)]
        else:
            print(f"ERROR: Unsupported script type: {script_path.suffix}")
            return False

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stdout)
        if result.stderr:
            print("Warnings/Errors:")
            print(result.stderr[:1000])
        return True

    except subprocess.CalledProcessError as e:
        print(f"ERROR running {description}:")
        print(f"  STDOUT: {e.stdout[:500]}")
        print(f"  STDERR: {e.stderr[:500]}")
        return False


def main():
    print("Cluster Anatomical Analysis Pipeline")
    print("=" * 60)

    # Each entry: (script_path, description, required, python_executable)
    # required=True means a failure in this step halts the pipeline.
    # python_executable=None falls back to sys.executable for .py scripts.
    steps = [
        (
            _SCRIPT_DIR / '01_generate_cluster_mri_inputs.py',
            'MRI Cluster Input Generation',
            True,
            '/usr/local/bin/python3.11',
        ),
        (
            _SCRIPT_DIR / '02_plot_cluster_mri_brain.R',
            'MRI Cluster Brain Visualizations',
            False,
            None,
        ),
        (
            _SCRIPT_DIR / '03_plot_subcortical_cluster_yabplot.py',
            'Subcortical Cluster Yabplot',
            False,
            '/usr/local/bin/python3.11',
        ),
        (
            _SCRIPT_DIR / '04_combined_cluster_brain_figure_frozen.R',
            'Combined single-page cluster brain figure (per cluster, FDR)',
            False,
            None,
        ),
    ]

    results = {}
    for script_path, description, required, python_executable in steps:
        if not Path(script_path).exists():
            msg = f"Script not found: {script_path}"
            if required:
                print(f"ERROR: {msg}")
                return False
            else:
                print(f"WARNING: {msg} (skipping)")
                continue

        # For .R scripts, python_executable is unused; run_script handles dispatch
        success = run_script(script_path, description, python_executable)
        results[description] = success

        if not success and required:
            print(f"\nRequired step failed: {description}")
            print("Pipeline aborted.")
            return False

    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    for step_name, success in results.items():
        status = "SUCCESS" if success else "FAILED"
        print(f"  [{status}] {step_name}")

    overall = all(results.values())
    if overall:
        print("\nAll steps completed successfully.")
    else:
        print("\nOne or more steps failed.")
    return overall


if __name__ == '__main__':
    start = time.time()
    success = main()
    print(f"\nPipeline duration: {time.time() - start:.2f} seconds")
    sys.exit(0 if success else 1)
