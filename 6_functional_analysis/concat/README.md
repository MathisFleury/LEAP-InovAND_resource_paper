# 6. Functional MRI — concat (current, priority regime)

Functional connectivity, autism vs NT, on XCP-D v0.11. Self-contained — no
dependency on `../non_concat/` (the legacy XCP-D 0.8 pipeline, kept for the
non-concatenated-runs analysis; see `../README.md`).

## Primary regime

**`nogsr_concat`**: no global-signal-regression, runs concatenated per
subject (≥6 min usable data). This is what the manuscript reports (n=393
autism / 327 NT — the "reported whole-cohort N") and what
`11_age_sex_analysis` imports directly.

```bash
cd 6_functional_analysis/concat && ./run_all.sh
```

Runs preprocessing on `nogsr_concat` (build connectivity → ComBat → regress
→ z-score), the core edge/network-block/NBS analysis suite, the Figure 6b/3
brain-map chain, and QC (motion confound, sequence/ComBat SVM, FD x
connectivity, network proportions).

## Variant glossary

`run_gsr_all_modes.sh` (run manually — a long sweep, not part of
`run_all.sh`) reruns the same suite across 9 sensitivity variants. Variant
names are `<denoising>[_<modifier>]`:

| Denoising | Meaning |
|---|---|
| `nogsr` | 24P/8mm, no global-signal regression (default base) |
| `gsr` | 36P/6mm, with global-signal regression |
| `nogsr6` | 24P/6mm |

| Modifier suffix | Meaning |
|---|---|
| `_concat` | runs concatenated per subject before connectivity (the primary regime uses this) |
| `_mixed` | mixed run-selection strategy (see `step_01_build_connectivity_mixed.py`) |
| `_regfirst` | regress-then-ComBat order, instead of ComBat-then-regress |
| `_INOVAND` / `_LEAP` | restricted to one cohort (own ComBat batch structure) |
| `_anatQCpass` | filtered to subjects passing anatomical-MRI QC (`qc/anat_qc_crossref.py`) |

e.g. `gsr_concat_INOVAND` = 36P/6mm GSR, run-concatenated, INOVAND only.

## Outputs

- `outputs/nogsr_concat/` — the primary regime, actively rebuilt by `run_all.sh`.
- `outputs/figures/` — Figure 6b/3 brain-map chain, one `nbs_wholegroup_<variant>/`
  subfolder per variant it's been run for, plus cross-variant comparison PDFs.
- `outputs/qc_motion/`, `outputs/qc_sequence/` — QC, rebuilt by `run_all.sh` steps 12/13.
- `outputs/qc_connectivity_flatness/` — a one-off manual flat-timeseries
  investigation (no script regenerates it; see its own README inside).
- `outputs/sensitivity/` — the 9-variant GSR sweep from `run_gsr_all_modes.sh`
  (run once, 2026-07-08, not rebuilt since) plus a handful of earlier
  one-off nogsr-family checks (2026-06-29/30). **Only `gsr_concat` is
  individually reported**: its brain-map output
  (`outputs/tables/nbs_wholegroup_gsr_concat/`) feeds the manuscript's GSR
  sensitivity supplementary figure via `_resources/build_gsr_sensitivity_figure.R`.
  Everything else in `sensitivity/` was explored during revision and isn't
  cited in the manuscript. Note `gsr_concat` was last rebuilt 2026-07-16,
  before `nogsr_concat`'s latest (2026-09-03) rebuild — worth re-running if
  the two need to be compared on identical upstream data.

## Run

```bash
cd 6_functional_analysis/concat && ./run_all.sh          # primary regime (nogsr_concat)
./run_gsr_all_modes.sh                              # optional: full 9-variant sweep
```
