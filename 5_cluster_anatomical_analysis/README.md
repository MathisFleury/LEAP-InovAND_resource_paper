# 5. Cluster Anatomical MRI Analysis

Structural MRI contrasts by clinical cluster (t-stats, Cohen's d), ggseg
and subcortical visualisations, the Figure 6a composite grid, and
publication gt tables — on the current curated k-means clustering.

This tree runs the **current, priority regime** — curated k-means only.
Every other clustering method (Ward, GMM, NT-restricted, autism-only) run
during review, plus the earlier frozen per-cluster pipeline, are superseded
and live in `../legacy/5_cluster_anatomical_analysis/` (local only, not
part of this public release) with their own, fully self-contained
orchestrator.

## Scripts

| Script | Description |
|--------|-------------|
| `run_cluster_anatomical_analysis.py` | Orchestrator — runs steps 1-7 below in dependency order (`4` runs before `3`) |
| `1_generate_cluster_mri_inputs.py` | Per-cluster MRI inputs — CT/SA/subcortical contrasts vs pooled NT + pairwise cluster comparisons (t / Cohen's d / FDR) |
| `2_plot_cluster_mri_brain.R` | Cortical brain maps (ggseg, Cohen's d) |
| `_combined_cluster_brain_figure_shared.R` | Shared single-page cluster-brain-figure engine, parameterised via `ANAT_*` env vars. Used by `4` here *and* by the legacy frozen pipeline (`legacy/5_cluster_anatomical_analysis/scripts/04_combined_cluster_brain_figure_frozen.R`), so it stays in this tree rather than being duplicated. |
| `3_plot_subcortical_cluster_yabplot.py` | Subcortical brain map (yabplot, Cohen's d) |
| `4_combined_cluster_brain_figure.R` | Wrapper around `_combined_cluster_brain_figure_shared.R` — combined single-page cluster brain figure (per cluster, FDR) |
| `5_figure6a_anatomical.R` | Figure 6a composite grid (ggseg t-stat maps, 3 modalities x 3 clusters) |
| `6_cluster_global_measures_violin.py` | Global structural measures by cluster (violin: thickness/eTIV/area) |
| `7_gt_tables_cluster.R` | Publication gt tables (per-cluster MRI stats) |
| `run_curated_cluster_anatomical.py` | Reruns the pipeline once per curated clustering method, env-var based (`CLUSTER_METHODS=kmeans`, the default and only one kept in this tree — set `CLUSTER_METHODS=kmeans,ward,gmm` to reproduce the full sensitivity sweep, whose ward/gmm results live in `../legacy/`). All 7 steps stage into `outputs/` via `CLUSTER_OUT_DIR` (default `"outputs"`); for kmeans that's also the final resting name (no rename needed), other methods get renamed to `outputs_curated_<method>/`. Overlaps with `1_clustering/scripts/run_anatomical_with_other_methods.py` (file-swap based) for the same core sensitivity check, but adds `OUT_SUFFIX` (e.g. a no-Euler MRI-table variant) that the generic driver doesn't have. |

## Outputs

`outputs/figures/` and `outputs/tables/` are kept separate: `tables/` holds
all CSV data (including the `r_input_files/` intermediate stats every
plotting step reads) and gt-table renders (`table_cluster_*.{html,pdf,png}`
— table renders, not figures, even though some formats are rasterised);
`figures/` holds only the actual plots (brain maps, yabplots, combined
figures, Figure 6a, the violin plot).

## Run

```bash
cd 5_cluster_anatomical_analysis/scripts
python3.11 run_cluster_anatomical_analysis.py
# Optional: restrict NT pool to NT subjects with Cluster == 'C1'
python3.11 run_cluster_anatomical_analysis.py --nt-c1

# Or, equivalently, via the curated wrapper (also supports ward/gmm sensitivity runs):
python3.11 run_curated_cluster_anatomical.py

# Legacy frozen + other-method pipelines (local only, not part of this public release)
cd ../../legacy/5_cluster_anatomical_analysis/scripts
python3.11 run_cluster_anatomical_analysis.py
```

## Data

Reads the curated cluster labels from `1_clustering/outputs/tables/` and
the QC+ComBat+regressed FreeSurfer table (see root `README.md` for the
`$IMG5` mount and `LEAP_INOVAND_DATA` setup).
