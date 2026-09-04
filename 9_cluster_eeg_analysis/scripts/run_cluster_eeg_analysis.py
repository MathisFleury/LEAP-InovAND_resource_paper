#!/usr/bin/env python3
"""
EEG Cluster Analysis Pipeline (curated clinical, k-means clusters vs NT).

0. Rebuild corrected alpha peak from raw $IMG5 (shared 8_eeg preprocessing)
1. Alpha peak — Autism k-means clusters vs NT
2. Power bands — Autism k-means clusters vs NT
"""

import sys
import subprocess
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_PREP = _SCRIPT_DIR.parent.parent / "8_eeg_analysis" / "preprocessing"  # shared
PY = "/usr/local/bin/python3.11"


def run_script(script_path: Path, description: str) -> bool:
    print(f"\n{'='*60}\nRunning: {description}\n{'='*60}")
    if not script_path.exists():
        print(f"ERROR: script not found: {script_path}")
        return False
    try:
        result = subprocess.run([PY, str(script_path)], capture_output=True, text=True, check=True)
        print(result.stdout)
        if result.stderr:
            print("Warnings/Info:", result.stderr[:1000])
        return True
    except subprocess.CalledProcessError as exc:
        print(f"ERROR running '{description}':\n  STDOUT: {exc.stdout[-1000:]}\n  STDERR: {exc.stderr[-1000:]}")
        return False


def main() -> bool:
    print("EEG Cluster Analysis Pipeline (k-means, curated)\n" + "=" * 60)
    steps = [
        (_PREP / 'build_alpha_peak_corrected.py', 'Preprocessing: rebuild corrected alpha peak (raw)'),
        (_SCRIPT_DIR / '01_eeg_alpha_peak_clusters.py', 'EEG Alpha Peak — clusters vs NT'),
        (_SCRIPT_DIR / '02_eeg_power_bands_clusters.py', 'EEG Power Bands — clusters vs NT'),
    ]
    results = {}
    for script_path, description in steps:
        results[description] = run_script(script_path, description)
        if not results[description]:
            print(f"\nRequired step failed: '{description}' — aborting.")
            break

    print("\n" + "=" * 60 + "\nPIPELINE SUMMARY\n" + "=" * 60)
    for name, ok in results.items():
        print(f"  [{'SUCCESS' if ok else 'FAILED'}] {name}")
    return all(results.values()) and len(results) == len(steps)


if __name__ == '__main__':
    start = time.time()
    success = main()
    print(f"\nPipeline duration: {time.time() - start:.2f} seconds")
    sys.exit(0 if success else 1)
