#!/usr/bin/env python3
"""
Cluster fMRI Functional Analysis Pipeline

Orchestrates the full cluster-level fMRI connectivity analysis:
  1. Generate per-cluster connectivity inputs (Python)
  2. Plot Schaefer atlas visualisations for each cluster (R)
"""

import sys
import subprocess
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent


def run_script(script_path: Path, description: str, python_executable: str | None = None) -> bool:
    """
    Execute a Python (.py) or R (.R) script as a subprocess.

    Parameters
    ----------
    script_path        : Path to the script file.
    description        : Human-readable step name (used in printed headers).
    python_executable  : Path to the Python interpreter to use for .py files.
                         Falls back to the currently running interpreter when None.

    Returns
    -------
    True on success, False on failure.
    """
    print(f"\n{'=' * 60}")
    print(f"Running: {description}")
    print(f"{'=' * 60}")

    if not script_path.exists():
        print(f"ERROR: Script not found: {script_path}")
        return False

    if script_path.suffix == '.py':
        executable = python_executable or sys.executable
        cmd = [executable, str(script_path)]
    elif script_path.suffix == '.R':
        cmd = ['Rscript', str(script_path)]
    else:
        print(f"ERROR: Unsupported script type: {script_path.suffix}")
        return False

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stdout)
        if result.stderr:
            # R sends normal messages to stderr; truncate to avoid noise
            print("Warnings/Info:\n" + result.stderr[:1000])
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR running '{description}':")
        print(f"  STDOUT: {e.stdout[:800]}")
        print(f"  STDERR: {e.stderr[:800]}")
        return False


def main() -> bool:
    print("=" * 60)
    print("CLUSTER fMRI FUNCTIONAL ANALYSIS PIPELINE")
    print("=" * 60)

    # Each entry: (script_path, description, is_required, python_executable)
    steps = [
        (
            _SCRIPT_DIR / '01_generate_cluster_fmri_inputs.py',
            'fMRI Cluster Input Generation',
            True,
            '/usr/local/bin/python3.11',
        ),
        (
            _SCRIPT_DIR / '02_plot_cluster_fmri_schaefer.R',
            'fMRI Cluster Schaefer Visualization',
            False,
            None,
        ),
        (
            _SCRIPT_DIR / '03_cluster_hyper_hypo_analysis.py',
            'Hyper/Hypo Edge Count Analysis',
            False,
            '/usr/local/bin/python3.11',
        ),
        (
            _SCRIPT_DIR / '04_plot_cluster_network_matrix.R',
            'Network Connectivity Matrix',
            False,
            None,
        ),
        (
            _SCRIPT_DIR / '05_plot_cluster_subcortical_yabplot.py',
            'Subcortical Connectivity Yabplot',
            False,
            '/usr/local/bin/python3.11',
        ),
        (
            _SCRIPT_DIR / '06_cluster_nbs.py',
            'Per-cluster NBS (Figure 6b data)',
            False,
            '/usr/local/bin/python3.11',
        ),
        (
            # Must run before 07 -- the R composite embeds its PNG insets.
            _SCRIPT_DIR / '08_figure6b_nbs_subcortical.py',
            'Figure 6b subcortical yabplot insets',
            False,
            '/usr/local/bin/python3.11',
        ),
        (
            _SCRIPT_DIR / '07_figure6b_nbs_cortical.R',
            'Figure 6b full composite (cortical + subcortical insets)',
            False,
            None,
        ),
        (
            _SCRIPT_DIR / '09_nbs_prominent_networks.py',
            'Most prominent networks/regions per cluster x direction',
            False,
            '/usr/local/bin/python3.11',
        ),
    ]

    results: dict[str, bool] = {}

    for script_path, description, required, py_exec in steps:
        if not script_path.exists():
            msg = f"Script not found: {script_path}"
            if required:
                print(f"ERROR: {msg}")
                return False
            else:
                print(f"WARNING: {msg} (skipping optional step)")
                continue

        success = run_script(script_path, description, py_exec)
        results[description] = success

        if not success and required:
            print(f"\nRequired step failed: '{description}' — aborting pipeline.")
            return False

    # Summary
    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    for step, success in results.items():
        status = "SUCCESS" if success else "FAILED "
        print(f"  [{status}] {step}")

    return all(results.values())


if __name__ == '__main__':
    start = time.time()
    success = main()
    print(f"\nPipeline duration: {time.time() - start:.2f} seconds")
    sys.exit(0 if success else 1)
