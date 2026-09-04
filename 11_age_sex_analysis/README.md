# 11. Age x sex interaction analysis

Diagnosis x age and diagnosis x sex interaction models (anatomical +
functional), age-binned sanity checks, and age-composition sensitivity
analyses — the manuscript's "The Effects of Age and Sex" section.

See `age_sex_writeup.md` for the full Methods/Results write-up. Run order
(`run_all.sh`) follows the scripts' own dependency chain: `01` is a hard
dependency of `05`, `06`, `09`, `10` (imported directly, not just run
first), and reads `4_anatomical_analysis/`'s config and script 01
— run section 4's curated pipeline before this one. `10` already invokes
`11_age_bin_brain_grid.R` itself; `run_all.sh` does not call it separately.

## Run

```bash
cd 11_age_sex_analysis && ./run_all.sh
```

## Outputs

`outputs/tables/` — all CSVs, including `r_input_files/` and the per-bin
`r_input_bin_<bin>/` / sensitivity `r_input_<age_5_22,age_matched,full>/`
ggseg-input folders, and `nbs/` (edge lists, summaries). `outputs/figures/`
— all plots, including `nbs/` (double-matrix + combined figures) and
`composite/` (the age-binned brain-map grids). A PDF's matching source-data
CSV (e.g. `anat_interaction_*.csv` next to `anat_interaction_*.pdf`) stays
alongside it in `figures/` rather than moving to `tables/`.
