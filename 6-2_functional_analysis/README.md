# 6-2. Functional MRI (current pipeline)

Functional connectivity, autism vs NT, on XCP-D v0.11. Self-contained — no
dependency on `6_functional_analysis/` (the legacy XCP-D 0.8 pipeline).

## Primary regime

**`nogsr_concat`**: no global-signal-regression, runs concatenated per
subject (≥6 min usable data). This is what the manuscript reports (n=393
autism / 327 NT — the "reported whole-cohort N") and what
`11_age_sex_analysis` imports directly.

```bash
cd 6-2_functional_analysis && ./run_all.sh
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

## Run

```bash
cd 6-2_functional_analysis && ./run_all.sh          # primary regime (nogsr_concat)
./run_gsr_all_modes.sh                              # optional: full 9-variant sweep
```
