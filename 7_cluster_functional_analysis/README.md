# 7. Cluster Functional MRI Analysis

Functional connectivity, autism *clusters* (C1/C2/C3) vs pooled NT — the
per-cluster analogue of `6_functional_analysis`'s whole-group analysis (see
its `METHODS_functional_6_vs_7.md` for how the two relate).

## Primary regime

Reads the **same connectivity table as section 6's primary regime**:
`6_functional_analysis/concat/preprocessing/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv`.

```bash
cd 7_cluster_functional_analysis && ./run_all.sh
```

Runs the per-cluster Schaefer / hyper-hypo / network-matrix / yabplot chain
(01–05), the 6-minute-cutoff and covariate-adjusted sensitivity checks
(01b/01c), and the Figure 6b per-cluster NBS chain (06–09).

Until this fix, `01_generate_cluster_fmri_inputs.py` had **no** connectivity
override at all — hardcoded to `6_functional_analysis/non_concat` (the
legacy XCP-D 0.8 pipeline), so steps 01–05's figures were built from the
wrong data. Only step 06 (NBS) supported a `FMRI_CONN_FILE` override, and
only via a manual export never captured in any committed script. Both are
now fixed: `01` defaults to the concat primary regime (still overridable via
`FMRI_CONN_FILE`, same env var `06_cluster_nbs.py` / `_sensitivity_utils`
honour), and `run_all.sh` pins it explicitly for a real run.

## Outputs

- `outputs/figures/` — plots only (Schaefer, network-matrix, yabplot,
  `nbs_nogsr_concat/` — Figure 6b composite + subcortical insets).
- `outputs/tables/` — all CSVs, including `r_input_files/` and
  `nbs_nogsr_concat/` (NBS edge lists, summary, prominent networks/regions).
- `outputs/sensitivity/` — an older NBS sweep across denoising variants
  (`nbs`, `nbs_gsr`, `nbs_gsr_concat`, `nbs_nogsr`; mirrors section 6's
  gsr/nogsr/concat sweep) — exploratory, not individually reported.
  Not to be confused with `outputs/figures/sensitivity/`, the 01b/01c
  6-minute-cutoff / covariate-adjusted checks' own output folder — same
  word, two unrelated sensitivity analyses at different folder depths.

An earlier, frozen paper-vintage run lives in
`../legacy/7_cluster_functional_analysis/outputs_paper/` (local only, not
part of this public release).
