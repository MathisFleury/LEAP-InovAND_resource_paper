#!/usr/bin/env python3
"""
fMRI Functional Connectivity Analysis Pipeline (Autism vs NT)

Runs the full analysis suite on a harmonised connectivity table:
  01  edge-level t-tests + BH-FDR (hyper/hypo)            [python]
  06  network-block label-permutation test (FDR)         [python]
  07  NBS — component-level FWER (Zalesky)                [python]
  02  Schaefer-atlas brain visualisation                 [R]
  05  double network-matrix figure (hyper/hypo)          [R]
  03  subcortical connectivity yabplot                   [python]

Then two more blocks: sensitivity analyses (01b/01c cutoff + covariate
variants, and their comparison summary/figures against the primary above),
and a legacy-vs-revised preprocessing QA comparison (12_plot_legacy_vs_revised_schaefer.R).

Targeting a specific dataset/variant: set the env vars the scripts honour
(inherited by every subprocess here), e.g. to run on a 6-2 v0.11 variant:

  AUTISM_TD_FMRI_CSV=.../6-2_functional_analysis/preprocessing/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv \
  AUTISM_TD_OUTPUT_DIR=.../6-2_functional_analysis/outputs/nogsr_concat \
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

    # NBS double-matrix renders into a subdir of the run's output dir: 08 exports
    # the largest-component edges there, then 05 is re-run against it.
    base_out = os.environ.get("AUTISM_TD_OUTPUT_DIR",
                              str(_SCRIPT_DIR.parent / "outputs" / "figures" / "revised"))
    nbs_double_dir = str(Path(base_out) / "nbs_double")

    # (script, description, required, python_exec, env_override)  — 01 must run
    # first (it writes the *_for_r.csv / *_full_data.csv the R + yabplot steps
    # consume). 08 must run before the second 05 (NBS double matrix).
    steps = [
        (_SCRIPT_DIR / '01_autism_td_connectivity_analysis.py', 'Edge t-tests + FDR (hyper/hypo)', True, PY, None),
        (_SCRIPT_DIR / '06_network_block_permutation.py', 'Network-block permutation (FDR)', False, PY, None),
        (_SCRIPT_DIR / '07_nbs_test.py', 'NBS component-level FWER', False, PY, None),
        (_SCRIPT_DIR / '02_plot_autism_td_schaefer.R', 'Schaefer atlas brain visualisation', False, None, None),
        (_SCRIPT_DIR / '05_plot_double_network_matrix.R', 'Double network-matrix figure', False, None, None),
        (_SCRIPT_DIR / '08_nbs_component_edges.py', 'NBS largest-component edge export', False, PY, None),
        (_SCRIPT_DIR / '05_plot_double_network_matrix.R', 'NBS double network-matrix figure', False, None,
         {"AUTISM_TD_OUTPUT_DIR": nbs_double_dir}),
        (_SCRIPT_DIR / '03_plot_subcortical_yabplot.py', 'Subcortical connectivity yabplot', False, PY, None),
    ]

    # Sensitivity analyses (6-min/mean_fd cutoff, covariate-adjusted) -- run
    # after the primary chain above since 02/03 compare against its output.
    sensitivity_steps = [
        (_SCRIPT_DIR / '01b_autism_td_sensitivity_6min.py', 'Sensitivity: >=6 min + mean_fd<1mm cutoff', False, PY, None),
        (_SCRIPT_DIR / '01c_autism_td_covariate_adjusted.py', 'Sensitivity: OLS group + mean_fd + minutes covariates', False, PY, None),
        (_SCRIPT_DIR / '09_sensitivity_comparison_summary.py', 'Sensitivity: primary vs 6-min vs covariate-adjusted summary', False, PY, None),
        (_SCRIPT_DIR / '10_plot_sensitivity_comparison.py', 'Sensitivity: n-significant / sample-size bar charts', False, PY, None),
        (_SCRIPT_DIR / '11_plot_autism_td_schaefer_sensitivity.R', 'Sensitivity: Schaefer atlas maps (6-min + covariate-adjusted)', False, None, None),
    ]

    # Legacy-vs-revised preprocessing comparison (QA diagnostic, not a
    # sensitivity analysis) -- only meaningful if outputs/figures/legacy/
    # (preprocessing_cohort.py's old output) still exists locally.
    legacy_compare_steps = [
        (_SCRIPT_DIR / '12_plot_legacy_vs_revised_schaefer.R', 'QA: legacy vs revised preprocessing Schaefer comparison', False, None, None),
    ]

    results = {}
    for block_name, block_steps in [("primary", steps), ("sensitivity", sensitivity_steps),
                                     ("legacy-comparison", legacy_compare_steps)]:
        if block_name != "primary":
            print(f"\n{'#' * 60}\n# {block_name.upper()} BLOCK\n{'#' * 60}")
        if not _run_block(block_steps, results, run_script):
            print(f"Required step failed in {block_name} block -- aborting.")
            return False

    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    for step, success in results.items():
        status = "SUCCESS" if success else "FAILED"
        print(f"  [{status}] {step}")

    return all(results.values())


def _run_block(steps, results, run_script):
    """Run one ordered list of steps, recording each result into `results`.
    Returns False (abort the whole pipeline) iff a required step is missing
    or fails; otherwise True even if optional steps failed."""
    for script_path, description, required, py_exec, env_override in steps:
        if not script_path.exists():
            msg = f"Script not found: {script_path}"
            if required:
                print(f"ERROR: {msg}")
                results[description] = False
                return False
            print(f"WARNING: {msg} (skipping)")
            continue
        success = run_script(script_path, description, py_exec, env_override)
        results[description] = success
        if not success and required:
            print(f"Required step failed: {description}")
            return False
    return True


if __name__ == "__main__":
    start = time.time()
    success = main()
    print(f"\nPipeline duration: {time.time() - start:.2f} seconds")
    sys.exit(0 if success else 1)
