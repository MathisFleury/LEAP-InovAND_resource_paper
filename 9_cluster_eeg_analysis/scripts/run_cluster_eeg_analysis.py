#!/usr/bin/env python3
"""
EEG Cluster Analysis Pipeline

Orchestrates the complete EEG cluster analysis:
1. Alpha peak — Autism clusters vs TD
2. Power bands — Autism clusters vs TD
"""

import sys
import subprocess
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent


def run_script(script_path: Path, description: str,
               python_executable: str = None) -> bool:
    """Run a Python script as a subprocess and report success/failure.

    Parameters
    ----------
    script_path        : Absolute path to the .py script.
    description        : Human-readable name shown in the console header.
    python_executable  : Python interpreter to use (defaults to sys.executable).

    Returns
    -------
    True on success, False on failure.
    """
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Script : {script_path}")
    print(f"{'='*60}")

    if not script_path.exists():
        print(f"ERROR: script not found: {script_path}")
        return False

    executable = python_executable or sys.executable

    try:
        result = subprocess.run(
            [executable, str(script_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("Warnings/Info:", result.stderr[:1000])
        return True
    except subprocess.CalledProcessError as exc:
        print(f"ERROR running '{description}':")
        if exc.stdout:
            print(f"  STDOUT: {exc.stdout[:1000]}")
        if exc.stderr:
            print(f"  STDERR: {exc.stderr[:1000]}")
        return False
    except FileNotFoundError:
        print(f"ERROR: Python executable not found: {executable}")
        return False


def main() -> bool:
    """Run all cluster EEG analysis steps sequentially.

    Returns True if all required steps succeeded, False otherwise.
    """
    print("EEG Cluster Analysis Pipeline")
    print("=" * 60)

    steps = [
        (
            _SCRIPT_DIR / '01_eeg_alpha_peak_clusters.py',
            'EEG Alpha Peak Cluster Analysis',
            True,
            '/usr/local/bin/python3.11',
        ),
        (
            _SCRIPT_DIR / '02_eeg_power_bands_clusters.py',
            'EEG Power Bands Cluster Analysis',
            True,
            '/usr/local/bin/python3.11',
        ),
    ]

    results = {}
    for script_path, description, required, py_exec in steps:
        success = run_script(script_path, description, py_exec)
        results[description] = success
        if not success and required:
            print(f"\nRequired step failed: '{description}' — aborting pipeline.")
            break

    # Summary
    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    for step_name, step_success in results.items():
        status = "SUCCESS" if step_success else "FAILED"
        print(f"  [{status}] {step_name}")

    overall = all(results.values()) if results else False
    print(f"\nOverall: {'SUCCESS' if overall else 'FAILED'}")
    return overall


if __name__ == '__main__':
    start = time.time()
    success = main()
    print(f"\nPipeline duration: {time.time() - start:.2f} seconds")
    sys.exit(0 if success else 1)
