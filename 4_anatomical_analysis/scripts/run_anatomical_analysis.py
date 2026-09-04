#!/usr/bin/env python3
"""
Anatomical MRI Analysis Pipeline (current, priority regime)

Runs the v2 anatomical pipeline rebuilt on the QC + ComBat +
age/sex/eTIV-regressed FreeSurfer table, with LEAP_W1 / INOVAND_T1
wave selection (see 1_anatomical_mri_autism_nt.py for details).

Steps:
  1. Autism vs NT contrasts — t / Cohen's d / FDR per region.
  2. Cortical brain maps — ggseg (R), fill = Cohen's d, FDR outlines.
  3. Subcortical brain map — yabplot, Cohen's d.
  4. Euler number by population — surface-reconstruction QC summary.
  5. LOEUF (hg38) × MRI Pearson correlations — all-carriers only, per-panel FDR.
  6. LOEUF (hg38) brain maps — ggseg cortical + aseg subcortical; raw p<0.01 outlines.
  7. LOEUF (hg38) OLS beta + permutation, per-ROI beta figures, gt tables.
  8. Clinical x MRI correlations (curated clinical) + brain maps + scatter panels.
  9. QC/robustness: per-site effect size, age-imbalance sensitivity, Euler
     confound check (whole-cohort + all-features whole-brain figure).
  10. Figure 7 composite (thickness maps + IQxSRS panel + beta boxplots) and
      the supplementary brain-map grid.

Step 19 (build_composite) ALSO requires
../../10_clinical_analysis/scripts/04e_plot_iq_srs_cluster_reframe.R to have
been run first (a different section) -- it checks for that prerequisite
itself and exits with a clear message if missing, rather than failing silently.
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
    print("Anatomical MRI Analysis Pipeline (current, priority regime)")
    print("=" * 60)

    steps = [
        (_SCRIPT_DIR / "1_anatomical_mri_autism_nt.py",
         "Autism vs NT — t / Cohen's d / FDR", True),
        (_SCRIPT_DIR / "2_plot_anatomical_mri_brain.R",
         "Cortical brain maps (ggseg, Cohen's d)", False),
        (_SCRIPT_DIR / "7_combined_brain_figure.R",
         "Combined single-page brain figure (4 features, FDR, ROI labels)", False),
        (_SCRIPT_DIR / "3_plot_subcortical_yabplot.py",
         "Subcortical brain map (yabplot, Cohen's d)", False),
        (_SCRIPT_DIR / "4_euler_by_population.py",
         "Euler number by population (QC summary)", False),
        (_SCRIPT_DIR / "5_loeuf_mri_correlations_hg38.py",
         "LOEUF (v4) × MRI Pearson correlations (per pathway×feature FDR)", False),
        (_SCRIPT_DIR / "6_loeuf_mri_brain_maps_hg38.R",
         "LOEUF (v4) correlation brain maps (per-pathway combined figure)", False),
        (_SCRIPT_DIR / "8_loeuf_mri_regression_permutation.py",
         "LOEUF (v4) × MRI OLS β + permutation (gene-list-specific null)", False),
        (_SCRIPT_DIR / "9_loeuf_combined_brain_figure.R",
         "LOEUF (v4) regression β brain maps (per-pathway, perm p<0.05)", False),
        (_SCRIPT_DIR / "10_beta_coefficients_roi.py",
         "Per-ROI beta-coefficient figures (LOEUF x MRI, hg38 regperm)", False),
        (_SCRIPT_DIR / "11_age_imbalance_sensitivity.py",
         "Age-imbalance sensitivity (Autism vs NT, ComBat file)", False),
        (_SCRIPT_DIR / "12_clinical_mri_correlations.py",
         "Clinical x MRI Pearson correlations (curated clinical)", False),
        (_SCRIPT_DIR / "13_clinical_mri_brain_maps.R",
         "Clinical x MRI brain maps (ggseg)", False),
        (_SCRIPT_DIR / "14_site_effect_size_stg.py",
         "Per-site effect size, L superior temporal thickness", False),
        (_SCRIPT_DIR / "15_clinical_mri_scatter.py",
         "Clinical x MRI scatter panels", False),
        (_SCRIPT_DIR / "16_gt_tables_anat.R",
         "Publication gt tables (Autism vs NT, curated)", False),
        (_SCRIPT_DIR / "17_euler_confound_check.py",
         "Euler-number (image-quality) confound check", False),
        (_SCRIPT_DIR / "18_euler_phenotype_all_features.py",
         "Euler x phenotype whole-brain figure (all 232 features)", False),
        (_SCRIPT_DIR / "20_thickness_maps_for_composite.R",
         "Thickness brain maps for the Figure 7 composite", False),
        (_SCRIPT_DIR / "19_build_composite.py",
         "Figure 7 composite (also needs 10_clinical_analysis's 04e R script "
         "run first -- see module docstring)", False),
        (_SCRIPT_DIR / "21_supp_brain_grid.R",
         "Supplementary brain-map grid (4 features x 3 gene lists)", False),
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
