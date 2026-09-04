# 4. Anatomical MRI Analysis

Structural MRI (FreeSurfer-derived) for autism vs NT, IQ/SRS correlates,
LOEUF carrier effects, and brain visualisations, on the QC+ComBat+age/sex/
eTIV-regressed FreeSurfer table used in the paper.

This tree runs the **current, priority regime**. The earlier ComBat-only
(no QC/regression) pipeline is superseded and lives in
`../legacy/4_anatomical_analysis/` (local only, not part of this public
release) with its own, fully self-contained `run_all.sh`.

## Scripts

Numbered in run order (see `run_anatomical_analysis.py`'s own `steps` list
for the exact sequence — a couple of steps run out of numeric order, e.g.
`7` runs right after `2`, and `20` before `19`; the numbering otherwise
follows the analysis narrative).

| Script | Description |
|--------|-------------|
| `_config.py` | Shared configuration: MRI table, cluster roster, curated-clinical toggle |
| `curated_clinical.py` | Loads the curated per-cohort clinical TSVs (population_group/age/sex) |
| `run_anatomical_analysis.py` | Orchestrator — runs all 21 steps below in dependency order |
| `1_anatomical_mri_autism_nt_v2.py` | Autism vs NT contrasts — t / Cohen's d / FDR per region |
| `2_plot_anatomical_mri_brain_v2.R` | Cortical brain maps (ggseg, Cohen's d) |
| `_combined_brain_figure_shared.R` | Shared single-page brain-figure engine (Cohen's d or t-stat, parameterised via `ANAT_*` env vars). Used by `7` here *and* by the legacy frozen pipeline (`legacy/4_anatomical_analysis/scripts/13_combined_brain_figure_frozen.R`), so it stays in this tree rather than being duplicated. |
| `3_plot_subcortical_yabplot_v2.py` | Subcortical brain map (yabplot, Cohen's d) |
| `4_euler_by_population.py` | Euler number by population (surface-reconstruction QC) |
| `5_loeuf_mri_correlations_hg38_v2.py` | LOEUF (gnomAD v4) x MRI Pearson correlations, per-panel FDR |
| `6_loeuf_mri_brain_maps_hg38_v2.R` | LOEUF correlation brain maps (per-pathway) |
| `7_combined_brain_figure_v2.R` | Wrapper around `_combined_brain_figure_shared.R` — combined single-page brain figure (4 features, FDR, ROI labels) |
| `8_loeuf_mri_regression_permutation_v2.py` | LOEUF x MRI OLS beta + permutation (gene-list-specific null) |
| `9_loeuf_combined_brain_figure_v2.R` | LOEUF regression-beta brain maps (per-pathway, perm p<0.05) |
| `10_beta_coefficients_roi_v2.py` | Per-ROI beta-coefficient figures |
| `11_age_imbalance_sensitivity.py` | Age-imbalance sensitivity (Autism vs NT) |
| `12_clinical_mri_correlations_v2.py` | Clinical x MRI Pearson correlations (curated clinical) |
| `13_clinical_mri_brain_maps_v2.R` | Clinical x MRI brain maps (ggseg) |
| `14_site_effect_size_stg.py` | Per-site effect size, L superior temporal thickness |
| `15_clinical_mri_scatter_v2.py` | Clinical x MRI scatter panels |
| `16_gt_tables_anat_v2.R` | Publication gt tables (Autism vs NT) |
| `17_euler_confound_check.py` | Euler-number (image-quality) confound check |
| `18_euler_phenotype_all_features.py` | Euler x phenotype whole-brain figure (all 232 features) |
| `19_build_composite.py` | Figure 7 composite (thickness maps + IQxSRS panel + beta boxplots). Also needs `10_clinical_analysis/scripts/04e_plot_iq_srs_cluster_reframe.R` run first — checks for that prerequisite itself. |
| `20_thickness_maps_for_composite.R` | Thickness brain maps feeding the Figure 7 composite |
| `21_supp_brain_grid.R` | Supplementary brain-map grid (4 features x 3 gene lists) |

## Run

```bash
# From 4_anatomical_analysis/
./run_all.sh

# Legacy paper-reproduction pipeline (local only, not part of this public release)
cd ../legacy/4_anatomical_analysis && ./run_all.sh
```

## Data

Set `LEAP_INOVAND_DATA` and mount `$IMG5` as described in the root README.
A few section-specific overrides:

- `CLUSTER_FILE` — override the autism-vs-NT cluster/phenotype roster (defaults to the frozen paper roster; curated runs point this at a curated cluster file)
- `CURATED_CLINICAL` — `1` (default) loads curated per-cohort clinical TSVs via `curated_clinical.py`; `0` reproduces the frozen paper roster instead
- `LOEUF_OUT_DIR` — override the output base directory name (default `outputs`)
