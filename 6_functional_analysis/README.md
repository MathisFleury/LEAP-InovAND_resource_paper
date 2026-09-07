# 6. Functional MRI Analysis

Functional connectivity, whole-group autism vs NT. Two parallel analyses,
differing in preprocessing (XCP-D version) and in whether runs are
concatenated per subject before computing connectivity:

| | `concat/` | `non_concat/` |
|---|---|---|
| XCP-D version | v0.11 | v0.8 (legacy) |
| Runs | concatenated per subject (≥6 min usable data) | single best run per subject |
| Status | **current (primary)** | kept as a secondary/historical analysis |
| Self-contained | yes, no dependency on `non_concat/` | yes |

See `concat/README.md` for the primary pipeline (also the `nogsr_concat`
sensitivity-variant glossary — GSR, mixed run-selection, regress-order,
per-cohort, anat-QC — via `concat/run_gsr_all_modes.sh`).

## Downstream consumers

`11_age_sex_analysis` and `7_cluster_functional_analysis` read `concat/`'s
`nogsr_concat` output directly, but still import shared statistical helpers
(`_sensitivity_utils.py`) from `non_concat/scripts/` for the sensitivity
(6-minute-cutoff / motion-covariate) re-runs — a pre-existing cross-folder
coupling, not introduced by this reorganisation.
