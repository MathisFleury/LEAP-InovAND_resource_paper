#!/usr/bin/env python3
"""
EEG Autism vs NT Analysis Pipeline (curated clinical, k-means).

0. Rebuild corrected alpha peak from raw $IMG5 features (preprocessing)
1. Alpha peak — Autism vs NT
2. Multi-band power — Autism vs NT
"""

import sys
import subprocess
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_PREP = _SCRIPT_DIR.parent / "preprocessing"
PY = "/usr/local/bin/python3.11"


def run_script(script_path, description, python_executable=PY):
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")
    try:
        cmd = [python_executable or sys.executable, str(script_path)]
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
    print("EEG Autism vs TD Analysis Pipeline")
    print("=" * 60)

    steps = [
        (_PREP / 'build_alpha_peak_corrected.py', 'Preprocessing: rebuild corrected alpha peak (raw)', True, PY),
        (_SCRIPT_DIR / '01_eeg_alpha_peak_autism_td.py', 'EEG Alpha Peak Analysis', True, PY),
        (_SCRIPT_DIR / '02_eeg_power_bands_autism_td.py', 'EEG Multi-Band Power Analysis', True, PY),
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
