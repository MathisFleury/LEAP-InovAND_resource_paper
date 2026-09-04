#!/usr/bin/env python3
"""
fMRI Functional Connectivity Analysis Pipeline (Autism vs NT)

Runs the full analysis suite on a harmonised connectivity table:
  01  edge-level t-tests + BH-FDR (hyper/hypo)            [python]
  02  network-block label-permutation test (FDR)         [python]
  03  NBS — component-level FWER (Zalesky)                [python]
  04  NBS largest-component edge export                   [python]
  05  Schaefer-atlas brain visualisation                 [R]
  06  double network-matrix figure (hyper/hypo)          [R]
  07  subcortical connectivity yabplot                   [python]

Targeting a specific dataset/variant: set the env vars the scripts honour
(inherited by every subprocess here), e.g. to run on a 6-2 v0.11 variant:

  AUTISM_TD_FMRI_CSV=.../6_functional_analysis/concat/preprocessing/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv \
  AUTISM_TD_OUTPUT_DIR=.../6_functional_analysis/concat/outputs/nogsr_concat \
  python3.11 run_functional_analysis.py

With no env vars set it runs on this section's default preprocessing output.
"""

import os
import sys
import subprocess
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
PY = "/usr/local/bin/python3.11"   # project Python (see CLAUDE.md)


def run_script(script_path, description, python_executable=None, env_override=None):
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
        env = {**os.environ, **env_override} if env_override else None
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, env=env)
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

    if os.environ.get("AUTISM_TD_FMRI_CSV"):
        print(f"  input : {os.environ['AUTISM_TD_FMRI_CSV']}")
    if os.environ.get("AUTISM_TD_OUTPUT_DIR"):
        print(f"  output: {os.environ['AUTISM_TD_OUTPUT_DIR']}")

    # NBS double-matrix renders into a subdir of the run's output dir: 04 exports
    # the largest-component edges there, then 06 is re-run against it.
    base_out = os.environ.get("AUTISM_TD_OUTPUT_DIR",
                              str(_SCRIPT_DIR.parent / "outputs" / "figures" / "revised"))
    nbs_double_dir = str(Path(base_out) / "nbs_double")

    # (script, description, required, python_exec, env_override)  — 01 must run
    # first (it writes the *_for_r.csv / *_full_data.csv the R + yabplot steps
    # consume). 04 must run before the second 06 (NBS double matrix).
    steps = [
        (_SCRIPT_DIR / '01_autism_td_connectivity_analysis.py', 'Edge t-tests + FDR (hyper/hypo)', True, PY, None),
        (_SCRIPT_DIR / '02_network_block_permutation.py', 'Network-block permutation (FDR)', False, PY, None),
        (_SCRIPT_DIR / '03_nbs_test.py', 'NBS component-level FWER', False, PY, None),
        (_SCRIPT_DIR / '05_plot_autism_td_schaefer.R', 'Schaefer atlas brain visualisation', False, None, None),
        (_SCRIPT_DIR / '06_plot_double_network_matrix.R', 'Double network-matrix figure', False, None, None),
        (_SCRIPT_DIR / '04_nbs_component_edges.py', 'NBS largest-component edge export', False, PY, None),
        (_SCRIPT_DIR / '06_plot_double_network_matrix.R', 'NBS double network-matrix figure', False, None,
         {"AUTISM_TD_OUTPUT_DIR": nbs_double_dir,
          "AUTISM_TD_TABLES_DIR": nbs_double_dir, "AUTISM_TD_FIGURES_DIR": nbs_double_dir}),
        (_SCRIPT_DIR / '07_plot_subcortical_yabplot.py', 'Subcortical connectivity yabplot', False, PY, None),
    ]

    results = {}
    for script_path, description, required, py_exec, env_override in steps:
        if not script_path.exists():
            msg = f"Script not found: {script_path}"
            if required:
                print(f"ERROR: {msg}")
                return False
            else:
                print(f"WARNING: {msg} (skipping)")
                continue
        success = run_script(script_path, description, py_exec, env_override)
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
